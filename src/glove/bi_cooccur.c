//  Tool to calculate word-word cooccurrence statistics for sentence aligned
//  corpus of 2 languages.
//
//  Copyright (c) 2014 The Board of Trustees of
//  The Leland Stanford Junior University. All Rights Reserved.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
//
//  For more information, bug reports, fixes, contact:
//    Jilei Wang (wangjileiRUC@gmail.com)


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// #include "bisearch_mapping_table.h"
#include "hash_mapping_table.h"

#define TSIZE 1048576
#define SEED 1159241
#define HASHFN bitwisehash

/* in count, most sentence are 3 to 60 words, the longest sentence in my corpus has 341 words*/
#define MAX_SENTENCE_LENGTH 500

/************************************************************
  Structures, Public Variables, and Helping Functions
************************************************************/

/** Original structures */

static const int MAX_STRING_LENGTH = 1000;
typedef double real;

typedef struct cooccur_rec {
    int word1;
    int word2;
    real val;
} CREC;

typedef struct cooccur_rec_id {
    int word1;
    int word2;
    real val;
    int id;
} CRECID;

typedef struct hashrec {
    char *word;
    long long id;
    struct hashrec *next;
} HASHREC;

//-----------------------------------------------------------

/** sentence of 2 languages*/

HASHREC **sentence_1, **sentence_2;
int sentence_length_1 = 0, sentence_length_2 = 0;

//-----------------------------------------------------------

/* variables from get_bi_cooccurrence */

char overflow_filename[200];
int fidcounter = 1, ind = 0;
FILE *foverflow;
CREC *cr;
real *bigram_table;
long long *lookup;

//-----------------------------------------------------------

/* CJ-Glo parameters */

char *mapping_table_file_path = "data/kanji_sc_perline.txt";
int cjglo = 0;
real sentence_rate = 0.0;
real cjglo_rate = 1.0;
int cj_matched_count = 0;

/* Original parameters */

int verbose = 3; // 0, 1, 2, or 3 (debug)
long long max_product; // Cutoff for product of word frequency ranks below which cooccurrence counts will be stored in a compressed full array
long long overflow_length; // Number of cooccurrence records whose product exceeds max_product to store in memory before writing to disk
int window_size = 15; // default context window size
int symmetric = 1; // 0: asymmetric, 1: symmetric
real memory_limit = 3; // soft limit, in gigabytes, used to estimate optimal array sizes
char *vocab_file, *file_head, *corpus_file_1, *corpus_file_2;


//-----------------------------------------------------------

/* Efficient string comparison */
int scmp( char *s1, char *s2 ) {
    while (*s1 != '\0' && *s1 == *s2) {s1++; s2++;}
    return(*s1 - *s2);
}

/* Move-to-front hashing and hash function from Hugh Williams, http://www.seg.rmit.edu.au/code/zwh-ipl/ */

/* Simple bitwise hash function */
unsigned int bitwisehash(char *word, int tsize, unsigned int seed) {
    char c;
    unsigned int h;
    h = seed;
    for (; (c =* word) != '\0'; word++) h ^= ((h << 5) + c + (h >> 2));
    return((unsigned int)((h&0x7fffffff) % tsize));
}

/* Create hash table, initialise pointers to NULL */
HASHREC ** inithashtable() {
    int	i;
    HASHREC **ht;
    ht = (HASHREC **) malloc( sizeof(HASHREC *) * TSIZE );
    for (i = 0; i < TSIZE; i++) ht[i] = (HASHREC *) NULL;
    return(ht);
}

/* Search hash table for given string, return record if found, else NULL */
HASHREC *hashsearch(HASHREC **ht, char *w) {
    HASHREC	*htmp, *hprv;
    unsigned int hval = HASHFN(w, TSIZE, SEED);
    for (hprv = NULL, htmp=ht[hval]; htmp != NULL && scmp(htmp->word, w) != 0; hprv = htmp, htmp = htmp->next);
    if ( htmp != NULL && hprv!=NULL ) { // move to front on access
        hprv->next = htmp->next;
        htmp->next = ht[hval];
        ht[hval] = htmp;
    }
    return(htmp);
}

/* Insert string in hash table, check for duplicates which should be absent */
void hashinsert(HASHREC **ht, char *w, long long id) {
    HASHREC	*htmp, *hprv;
    unsigned int hval = HASHFN(w, TSIZE, SEED);
    for (hprv = NULL, htmp = ht[hval]; htmp != NULL && scmp(htmp->word, w) != 0; hprv = htmp, htmp = htmp->next);
    if (htmp == NULL) {
        htmp = (HASHREC *) malloc(sizeof(HASHREC));
        htmp->word = (char *) malloc(strlen(w) + 1);
        strcpy(htmp->word, w);
        htmp->id = id;
        htmp->next = NULL;
        if (hprv == NULL) ht[hval] = htmp;
        else hprv->next = htmp;
    }
    else fprintf(stderr, "Error, duplicate entry located: %s.\n",htmp->word);
    return;
}

/* Read word from input stream */
int get_word(char *word, FILE *fin) {
    int i = 0, ch;
    while (!feof(fin)) {
        ch = fgetc(fin);
        if (ch == 13) continue;
        if ((ch == ' ') || (ch == '\t') || (ch == '\n')) {
            if (i > 0) {
                if (ch == '\n') ungetc(ch, fin);
                break;
            }
            if (ch == '\n') return 1;
            else continue;
        }
        word[i++] = ch;
        if (i >= MAX_STRING_LENGTH - 1) i--;   // truncate words that exceed max length
    }
    word[i] = 0;
    return 0;
}

/* Write sorted chunk of cooccurrence records to file, accumulating duplicate entries */
int write_chunk(CREC *cr, long long length, FILE *fout) {
    if (length == 0) return 0;

    long long a = 0;
    CREC old = cr[a];
    
    for (a = 1; a < length; a++) {
        if (cr[a].word1 == old.word1 && cr[a].word2 == old.word2) {
            old.val += cr[a].val;
            continue;
        }
        fwrite(&old, sizeof(CREC), 1, fout);
        old = cr[a];
    }
    fwrite(&old, sizeof(CREC), 1, fout);
    return 0;
}

/* Check if two cooccurrence records are for the same two words, used for qsort */
int compare_crec(const void *a, const void *b) {
    int c;
    if ( (c = ((CREC *) a)->word1 - ((CREC *) b)->word1) != 0) return c;
    else return (((CREC *) a)->word2 - ((CREC *) b)->word2);
    
}

/* Check if two cooccurrence records are for the same two words */
int compare_crecid(CRECID a, CRECID b) {
    int c;
    if ( (c = a.word1 - b.word1) != 0) return c;
    else return a.word2 - b.word2;
}

/* Swap two entries of priority queue */
void swap_entry(CRECID *pq, int i, int j) {
    CRECID temp = pq[i];
    pq[i] = pq[j];
    pq[j] = temp;
}

/* Insert entry into priority queue */
void insert(CRECID *pq, CRECID new, int size) {
    int j = size - 1, p;
    pq[j] = new;
    while ( (p=(j-1)/2) >= 0 ) {
        if (compare_crecid(pq[p],pq[j]) > 0) {swap_entry(pq,p,j); j = p;}
        else break;
    }
}

/* Delete entry from priority queue */
void delete(CRECID *pq, int size) {
    int j, p = 0;
    pq[p] = pq[size - 1];
    while ( (j = 2*p+1) < size - 1 ) {
        if (j == size - 2) {
            if (compare_crecid(pq[p],pq[j]) > 0) swap_entry(pq,p,j);
            return;
        }
        else {
            if (compare_crecid(pq[j], pq[j+1]) < 0) {
                if (compare_crecid(pq[p],pq[j]) > 0) {swap_entry(pq,p,j); p = j;}
                else return;
            }
            else {
                if (compare_crecid(pq[p],pq[j+1]) > 0) {swap_entry(pq,p,j+1); p = j + 1;}
                else return;
            }
        }
    }
}

//-----------------------------------------------------------

/* Write top node of priority queue to file, accumulating duplicate entries */
int merge_write(CRECID new, CRECID *old, FILE *fout) {
    if (new.word1 == old->word1 && new.word2 == old->word2) {
        old->val += new.val;
        return 0; // Indicates duplicate entry
    }
    fwrite(old, sizeof(CREC), 1, fout);
    *old = new;
    return 1; // Actually wrote to file
}

/* Merge [num] sorted files of cooccurrence records */
int merge_files(int num) {
    int i, size;
    long long counter = 0;
    CRECID *pq, new, old;
    char filename[200];
    FILE **fid, *fout;
    fid = malloc(sizeof(FILE) * num);
    pq = malloc(sizeof(CRECID) * num);
    fout = stdout;
    if (verbose > 1) fprintf(stderr, "Merging cooccurrence files: processed 0 lines.");
    
    /* Open all files and add first entry of each to priority queue */
    for (i = 0; i < num; i++) {
        sprintf(filename,"%s_%04d.bin",file_head,i);
        fid[i] = fopen(filename,"rb");
        if (fid[i] == NULL) {fprintf(stderr, "Unable to open file %s.\n",filename); return 1;}
        fread(&new, sizeof(CREC), 1, fid[i]);
        new.id = i;
        insert(pq,new,i+1);
    }
    
    /* Pop top node, save it in old to see if the next entry is a duplicate */
    size = num;
    old = pq[0];
    i = pq[0].id;
    delete(pq, size);
    fread(&new, sizeof(CREC), 1, fid[i]);
    if (feof(fid[i])) size--;
    else {
        new.id = i;
        insert(pq, new, size);
    }
    
    /* Repeatedly pop top node and fill priority queue until files have reached EOF */
    while (size > 0) {
        counter += merge_write(pq[0], &old, fout); // Only count the lines written to file, not duplicates
        if ((counter%100000) == 0) if (verbose > 1) fprintf(stderr,"\033[39G%lld lines.",counter);
        i = pq[0].id;
        delete(pq, size);
        fread(&new, sizeof(CREC), 1, fid[i]);
        if (feof(fid[i])) size--;
        else {
            new.id = i;
            insert(pq, new, size);
        }
    }
    fwrite(&old, sizeof(CREC), 1, fout);
    fprintf(stderr,"\033[0GMerging cooccurrence files: processed %lld lines.\n",++counter);
    for (i=0;i<num;i++) {
        sprintf(filename,"%s_%04d.bin",file_head,i);
        remove(filename);
    }
    fprintf(stderr,"\n");
    return 0;
}

/************************************************************
  CJ-GLO: Bilingual Cooccurence Matrix Count Function
************************************************************/

int check_overflow() {
    if (ind >= overflow_length - window_size) { // If overflow buffer is (almost) full, sort it and write it to temporary file
        qsort(cr, ind, sizeof(CREC), compare_crec);
        write_chunk(cr,ind,foverflow);
        fclose(foverflow);
        fidcounter++;
        sprintf(overflow_filename,"%s_%04d.bin",file_head,fidcounter);
        foverflow = fopen(overflow_filename,"w");
        ind = 0;
    }
    return 0;
}

int add_to_table(int w1, int w2, real increment) {
    int i;
    if ( w1 < max_product/w2 ) { 
        // Product is small enough to store in a full array
        // Weight by inverse of distance between words
        bigram_table[lookup[w1-1] + w2 - 2] += increment; 
    }
    else { 
        // Product is too big, data is likely to be sparse. 
        // Store these entries in a temporary buffer to be sorted, 
        // merged (accumulated), and written to file when it gets full.
        cr[ind].word1 = w1;
        cr[ind].word2 = w2;
        cr[ind].val = increment;
        ind++; // Keep track of how full temporary buffer is
        check_overflow();
    }
    return 0;
}

/**
 * Calculate the words in the dual sentence's influence to a certain word w1 
 */
int calculate_dual_sentence(HASHREC **sentence, int sentence_length, int w1) {
    int i, w2;
    if (sentence_rate <= 0.000001) return 0;
    for (i = 0; i < sentence_length; i++) {
        w2 = sentence[i]->id;
        add_to_table(w1, w2, sentence_rate);
    }
    return 0;
}

/** 
 * Calculate the words in a certain window's influence to a certain word w1,
 * The window and target word can be either in one sentence or in dual sentences.
 */
int calculate_window(HASHREC **sentence, int sentence_length, int w1, int center_pos, real rate) {
    int begin, end, k, w2, i;
    real increment;
    begin = (center_pos - window_size > 0) ? center_pos - window_size : 0;
    end  = (center_pos + window_size < sentence_length) ? center_pos + window_size : sentence_length - 1;
    if (verbose > 2) fprintf(stderr, "center_pos: %d, begin: %d, end: %d\n", center_pos, begin, end);
    // add the left and right window in sentence to cooccur table
    for (k = begin; k <= end; k++) {
        if (k == center_pos) {
            if (symmetric > 0) continue;
            else break;
        }
        w2 = sentence[k]->id; // Context word (frequency rank)
        increment = rate * 1.0/fabs((real)(center_pos - k));
        add_to_table(w1, w2, increment);
    }
    return 0;
}

/**
 * Find the corresponding words in the correspongding sentence, 
 * which include common characters with target word
 * lang_id : 0 - k2sc, 1 - sc2k
 */
int calculate_cjglo(HASHREC **sentence, int sentence_length, HASHREC *target_word, int lang_id) {
    int w1, i;
    char *word1, *word2;
    w1 = target_word->id;
    // do not include the prefix
    word1 = target_word->word + 1;
    for (i = 0; i < sentence_length; i++) {
        // do not include the prefix
        word2 = sentence[i]->word + 1;
        if (has_common_character(word1, word2, lang_id)) {
            calculate_window(sentence, sentence_length, w1, i, cjglo_rate);
        }
    }
    return 0;
}

/**
 * read a sentence from corpus file, save corresponding HASHREC into array, 
 * return sentence length 
 */
int read_sentence(HASHREC **sentence, FILE *fin, HASHREC **vocab_hash, char prefix) {
    int length = 0, flag;
    char str_prefix[MAX_STRING_LENGTH + 2], *str;
    str_prefix[0] = prefix;
    str = str_prefix + 1;
    while (1) {
        flag = get_word(str, fin);
        /* EOF or new line */
        if (feof(fin)) return -1;
        if (flag == 1) return length;
        sentence[length] = hashsearch(vocab_hash, str_prefix);
        if (sentence[length] == NULL) continue; // Skip out-of-vocabulary words
        length += 1;
    }
}

/**
 * Collect word-word cooccurrence counts from sentence aligned bilingual corpus 
 */
int get_bi_cooccurrence() {
    int flag, x, y, sentence_counter = 0;;
    int pos_1, pos_2;
    long long a, i, j = 0, k, id, counter = 0, vocab_size, w1, w2, w3, *history;
    char format[20], str[MAX_STRING_LENGTH + 1];
    FILE *fid, *fin_1, *fin_2;
    real r;
    HASHREC *htmp, **vocab_hash = inithashtable();
    
    cr = malloc(sizeof(CREC) * (overflow_length + 1));
    history = malloc(sizeof(long long) * window_size);

    /* inti sentence */
    sentence_1 = (HASHREC **) malloc( sizeof(HASHREC *) * MAX_SENTENCE_LENGTH );
    sentence_2 = (HASHREC **) malloc( sizeof(HASHREC *) * MAX_SENTENCE_LENGTH );
    for (i = 0; i < MAX_SENTENCE_LENGTH; i++) {
        sentence_1[i] = (HASHREC *) NULL;
        sentence_2[i] = (HASHREC *) NULL;
    }
    
    /* print parameters */
    fprintf(stderr, "COUNTING COOCCURRENCES\n");
    if (verbose > 0) {
        fprintf(stderr, "window size: %d\n", window_size);
        if (symmetric == 0) fprintf(stderr, "context: asymmetric\n");
        else fprintf(stderr, "context: symmetric\n");
    }
    if (verbose > 1) fprintf(stderr, "max product: %lld\n", max_product);
    if (verbose > 1) fprintf(stderr, "overflow length: %lld\n", overflow_length);
    sprintf(format,"%%%ds %%lld", MAX_STRING_LENGTH); // Format to read from vocab file, which has (irrelevant) frequency data
    if (verbose > 1) fprintf(stderr, "Reading vocab from file \"%s\"...", vocab_file);
    
    /* read vocab from file, save the rank */
    fid = fopen(vocab_file,"r");
    if (fid == NULL) {fprintf(stderr,"Unable to open vocab file %s.\n",vocab_file); return 1;}
    while (fscanf(fid, format, str, &id) != EOF) hashinsert(vocab_hash, str, ++j); // Here id is not used: inserting vocab words into hash table with their frequency rank, j
    fclose(fid);
    vocab_size = j;

    j = 0;
    if (verbose > 1) fprintf(stderr, "loaded %lld words.\nBuilding lookup table...", vocab_size);
    
    /* Build auxiliary lookup table used to index into bigram_table */
    lookup = (long long *)calloc( vocab_size + 1, sizeof(long long) );
    if (lookup == NULL) {
        fprintf(stderr, "Couldn't allocate memory!");
        return 1;
    }
    lookup[0] = 1;
    for (a = 1; a <= vocab_size; a++) {
        if ((lookup[a] = max_product / a) < vocab_size) lookup[a] += lookup[a-1];
        else lookup[a] = lookup[a-1] + vocab_size;
    }
    if (verbose > 1) fprintf(stderr, "table contains %lld elements.\n",lookup[a-1]);
    
    /* Allocate memory for full array which will store all cooccurrence counts for words whose product of frequency ranks is less than max_product */
    bigram_table = (real *)calloc( lookup[a-1] , sizeof(real) );
    if (bigram_table == NULL) {
        fprintf(stderr, "Couldn't allocate memory!");
        return 1;
    }
    fprintf(stderr, "Allocate memory for bigram_table.\n");

    /* Open 2 coupus files */
    fin_1 = fopen(corpus_file_1, "r");
    fin_2 = fopen(corpus_file_2, "r");
    sprintf(format,"%%%ds",MAX_STRING_LENGTH);
    sprintf(overflow_filename,"%s_%04d.bin",file_head, fidcounter);
    
    /* overflow file */
    foverflow = fopen(overflow_filename,"w");
    
    /* if open CJ-CLO modle, read the Chinese-Japanses common character mapping table */
    if (cjglo > 0) {
        load_mapping_tables();
    }

    /* For each token in input stream, calculate a weighted cooccurrence sum within window_size */
    while (1) {
        /* Read each sentence pairs */
        sentence_length_1 = read_sentence(sentence_1, fin_1, vocab_hash, '1');
        sentence_length_2 = read_sentence(sentence_2, fin_2, vocab_hash, '2');
        sentence_counter++;
        if (sentence_counter % 7000 == 0) fprintf(stderr, "read %d sentence pairs\n", sentence_counter);

        /* EOF */
        if (sentence_length_1 < 0 || sentence_length_2 < 0) break;
        /* no words in hastable in a sentence */
        if (sentence_length_1 == 0 || sentence_length_2 == 0) continue;

        // for every word in sentence 1
        // if cjglo, simple 2 kanji, lang_id = 1
        for (pos_1 = 0; pos_1 < sentence_length_1; pos_1++) {
            w1 = sentence_1[pos_1]->id;
            calculate_window(sentence_1, sentence_length_1, w1, pos_1, 1);
            calculate_dual_sentence(sentence_2, sentence_length_2, w1);
            if (cjglo > 0) calculate_cjglo(sentence_2, sentence_length_2, sentence_1[pos_1], 1);
        }

        // for every word in sentence 2
        // if cjglo, kanji 2 simple chinese, lang_id = 0
        for (pos_2 = 0; pos_2 < sentence_length_2; pos_2++) {
            w1 = sentence_2[pos_2]->id;
            calculate_window(sentence_2, sentence_length_2, w1, pos_2, 1);
            calculate_dual_sentence(sentence_1, sentence_length_1, w1);
            if (cjglo > 0) calculate_cjglo(sentence_1, sentence_length_1, sentence_2[pos_2], 0);
        }
    }
    
    /* Write out temp buffer for the final time (it may not be full) */
    if (verbose > 1) fprintf(stderr,"\033[0GProcessed %lld tokens.\n",counter);
    qsort(cr, ind, sizeof(CREC), compare_crec);
    write_chunk(cr,ind,foverflow);
    sprintf(overflow_filename,"%s_0000.bin",file_head);
    
    /* Write out full bigram_table, skipping zeros */
    if (verbose > 1) fprintf(stderr, "Writing cooccurrences to disk");
    fid = fopen(overflow_filename,"w");
    j = 1e6;
    for (x = 1; x <= vocab_size; x++) {
        if ( (long long) (0.75*log(vocab_size / x)) < j) {j = (long long) (0.75*log(vocab_size / x)); if (verbose > 1) fprintf(stderr,".");} // log's to make it look (sort of) pretty
        for (y = 1; y <= (lookup[x] - lookup[x-1]); y++) {
            if ((r = bigram_table[lookup[x-1] - 2 + y]) != 0) {
                fwrite(&x, sizeof(int), 1, fid);
                fwrite(&y, sizeof(int), 1, fid);
                fwrite(&r, sizeof(real), 1, fid);
            }
        }
    }
    
    if (verbose > 1) fprintf(stderr,"%d files in total.\n",fidcounter + 1);
    fclose(fid);
    fclose(foverflow);
    free(cr);
    free(lookup);
    free(bigram_table);
    free(vocab_hash);
    return merge_files(fidcounter + 1); // Merge the sorted temporary files
}

/************************************************************
                        Main Function
************************************************************/

int find_arg(char *str, int argc, char **argv) {
    int i;
    for (i = 1; i < argc; i++) {
        if (!scmp(str, argv[i])) {
            if (i == argc - 1) {
                printf("No argument given for %s\n", str);
                exit(1);
            }
            return i;
        }
    }
    return -1;
}

int main(int argc, char **argv) {
    int i;
    real rlimit, n = 1e5;
    vocab_file = malloc(sizeof(char) * MAX_STRING_LENGTH);
    file_head = malloc(sizeof(char) * MAX_STRING_LENGTH);
    corpus_file_1 = malloc(sizeof(char) * MAX_STRING_LENGTH);
    corpus_file_2 = malloc(sizeof(char) * MAX_STRING_LENGTH);
    
    if (argc == 1) {
        printf("Tool to calculate word-word cooccurrence statistics\n");
        printf("Author: Jeffrey Pennington (jpennin@stanford.edu)\n\n");
        printf("Usage options:\n");
        printf("\t-verbose <int>\n");
        printf("\t\tSet verbosity: 0, 1, or 2 (default)\n");
        printf("\t-cjglo <int>\n");
        printf("\t\tSet cjglo model: 0 (default) close model, 1 open model\n");
        printf("\t-cjglo-rate <float>\n");
        printf("\t\tSet cjglo rate as a factor of how cjglo word pairs contribute to the cooccurrence matrix\n");
        printf("\t-sentence-rate <float>\n");
        printf("\t\tSet sentence rate as a factor of the whole corresponding sentence as a context of the target word\n");
        printf("\t-symmetric <int>\n");
        printf("\t\tIf <int> = 0, only use left context; if <int> = 1 (default), use left and right\n");
        printf("\t-window-size <int>\n");
        printf("\t\tNumber of context words to the left (and to the right, if symmetric = 1); default 15\n");
        printf("\t-vocab-file <file>\n");
        printf("\t\tFile containing vocabulary (truncated unigram counts, produced by 'vocab_count'); default vocab.txt\n");
        printf("\t-corpus-file-1 <file>\n");
        printf("\t\tFile that each line is a sentence, each sentence contains vocabularys seperated by space.\n");
        printf("\t-corpus-file-2 <file>\n");
        printf("\t\tFile that each line is a sentence, each sentence contains vocabularys seperated by space.\n");
        printf("\t-memory <float>\n");
        printf("\t\tSoft limit for memory consumption, in GB -- based on simple heuristic, so not extremely accurate; default 4.0\n");
        printf("\t-max-product <int>\n");
        printf("\t\tLimit the size of dense cooccurrence array by specifying the max product <int> of the frequency counts of the two cooccurring words.\n\t\tThis value overrides that which is automatically produced by '-memory'. Typically only needs adjustment for use with very large corpora.\n");
        printf("\t-overflow-length <int>\n");
        printf("\t\tLimit to length <int> the sparse overflow array, which buffers cooccurrence data that does not fit in the dense array, before writing to disk. \n\t\tThis value overrides that which is automatically produced by '-memory'. Typically only needs adjustment for use with very large corpora.\n");
        printf("\t-overflow-file <file>\n");
        printf("\t\tFilename, excluding extension, for temporary files; default overflow\n");

        printf("\nExample usage:\n");
        printf("./cooccur -verbose 2 -symmetric 0 -window-size 10 -vocab-file vocab.txt -memory 8.0 -overflow-file tempoverflow -corpus-file-1 corpus.zh.txt -corpus-file-2 corpus.ja.txt > cooccurrences.bin\n\n");
        return 0;
    }

    if ((i = find_arg((char *)"-corpus-file-1", argc, argv)) > 0) strcpy(corpus_file_1, argv[i + 1]);
    //else return 1;
    if ((i = find_arg((char *)"-corpus-file-2", argc, argv)) > 0) strcpy(corpus_file_2, argv[i + 1]);
    //else return 1;
    if ((i = find_arg((char *)"-cjglo", argc, argv)) > 0) cjglo = atoi(argv[i + 1]);
    if ((i = find_arg((char *)"-cjglo-rate", argc, argv)) > 0) cjglo_rate = atof(argv[i + 1]);
    if ((i = find_arg((char *)"-sentence-rate", argc, argv)) > 0) sentence_rate = atof(argv[i + 1]);
    if ((i = find_arg((char *)"-verbose", argc, argv)) > 0) verbose = atoi(argv[i + 1]);
    if ((i = find_arg((char *)"-symmetric", argc, argv)) > 0) symmetric = atoi(argv[i + 1]);
    if ((i = find_arg((char *)"-window-size", argc, argv)) > 0) window_size = atoi(argv[i + 1]);
    if ((i = find_arg((char *)"-vocab-file", argc, argv)) > 0) strcpy(vocab_file, argv[i + 1]);
    else strcpy(vocab_file, (char *)"vocab.txt");
    if ((i = find_arg((char *)"-overflow-file", argc, argv)) > 0) strcpy(file_head, argv[i + 1]);
    else strcpy(file_head, (char *)"overflow");
    if ((i = find_arg((char *)"-memory", argc, argv)) > 0) memory_limit = atof(argv[i + 1]);
    
    /* The memory_limit determines a limit on the number of elements in bigram_table and the overflow buffer */
    /* Estimate the maximum value that max_product can take so that this limit is still satisfied */
    rlimit = 0.85 * (real)memory_limit * 1073741824/(sizeof(CREC));
    while (fabs(rlimit - n * (log(n) + 0.1544313298)) > 1e-3) n = rlimit / (log(n) + 0.1544313298);
    max_product = (long long) n;
    overflow_length = (long long) rlimit/6; // 0.85 + 1/6 ~= 1
    
    /* Override estimates by specifying limits explicitly on the command line */
    if ((i = find_arg((char *)"-max-product", argc, argv)) > 0) max_product = atoll(argv[i + 1]);
    if ((i = find_arg((char *)"-overflow-length", argc, argv)) > 0) overflow_length = atoll(argv[i + 1]);
    
    return get_bi_cooccurrence();
}

