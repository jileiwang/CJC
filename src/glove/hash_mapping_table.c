//  The mapping table from Kanji to Hanzi, or reverse.
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
#include "hash_mapping_table.h"

MAPTABREC **k2sc, **sc2k;

/* Create hash table, initialise ptrs to NULL */
MAPTABREC ** init_hash_mapping_table() {
    int i;
    MAPTABREC **ht;
    ht = (MAPTABREC **) malloc( sizeof(MAPTABREC *) * M_SIZE );
    for(i = 0; i < M_SIZE; i++) {
        ht[i] = (MAPTABREC *) NULL;
    }
    return(ht);
}

/* Search hash table for given string, return record if found, else NULL */
MAPTABREC * hash_search_mapping_table(MAPTABREC **ht, char *w) {
    MAPTABREC *htmp, *hprv;
    unsigned int hval = M_HASHFN(w, M_SIZE, M_SEED);

    for( hprv = NULL, htmp=ht[hval]
        ; htmp != NULL && scmp(htmp->original, w) != 0
        ; hprv = htmp, htmp = htmp->next )
    {
    ;
    }

    if( hprv!=NULL ) {
        /* move to front on access */
        hprv->next = htmp->next;
        htmp->next = ht[hval];
        ht[hval] = htmp;
    }

    return(htmp);
}

/* Search hash table for given string, insert if not found */
void hash_insert_mapping_table(MAPTABREC **ht, char *original, char *corresponding) {
    MAPTABREC *htmp, *hprv;
    unsigned int hval = M_HASHFN(original, M_SIZE, M_SEED);

    for( hprv = NULL, htmp=ht[hval]
        ; htmp != NULL && scmp(htmp->original, original) != 0
        ; hprv = htmp, htmp = htmp->next ) 
    {
    ;
    }

    if(htmp == NULL) {
        htmp = (MAPTABREC *) malloc( sizeof(MAPTABREC) );
        strcpy(htmp->original, original);
        strcpy(htmp->corresponding, corresponding);
        htmp->next = NULL;
        if(hprv == NULL) {
            ht[hval] = htmp;
        }
        else {
            hprv->next = htmp;
        }
        /* new records are not moved to front */
    }
    else {
        if(hprv != NULL) {
            /* move to front on access */
            hprv->next = htmp->next;
            htmp->next = ht[hval];
            ht[hval] = htmp;
        }
    }
    return;
}

/**
 * Read one or several Kanji or Hanzi from text file
 */
int get_characters(char *word, FILE *fin) {
    int i = 0, ch;
    while (!feof(fin)) {
        ch = fgetc(fin);
        if ((ch == ' ') || (ch == '\t') || (ch == '\n')) {
            break;
        }
        word[i++] = ch;
    }
    word[i] = 0;
    return i;
}

/**
 * Read mapping table from file, and sort it
 */
MAPTABREC ** load_one_mapping_table(char *filename, int table_size) {
    int i;
    FILE *f;
    char original[4], corresponding[22];
    MAPTABREC **table;
    table = init_hash_mapping_table();
    fprintf(stderr, "start reading %s\n", filename);
    f = fopen(filename, "rb");
    for (i = 0; i < table_size; i++) {
        get_characters(original, f);
        get_characters(corresponding, f);
        hash_insert_mapping_table(table, original, corresponding);
    }
    fclose(f);


    return table;
}

/**
 * Read 2 mapping tables
 */
void load_mapping_tables() {
    sc2k = load_one_mapping_table(SC2K_FILENAME, SC2K_SIZE);
    k2sc = load_one_mapping_table(K2SC_FILENAME, K2SC_SIZE);
}

/**
 * Compare is in word2 there is some characters in corresponding,
 * where corresponding is the characters mapped from a character
 */
int compare_characters(char *corresponding, char *word2) {
    int len1, len2, i, j;
    len1 = strlen(corresponding);
    len2 = strlen(word2);
    for (i = 0; i < len1; i+= 3) {
        for (j = 0; j < len2; j += 3) {
            if (word2[j] == corresponding[i] && word2[j+1] == corresponding[i+1] && word2[j+2] == corresponding[i+2]) {
                return 1;
            }
        }
    }
    return 0;
}

/**
 * return 1 if the 2 words has a common character.
 */
int has_common_character(char *word1, char *word2, int lang_id) {
    char source_ch[4], target_ch[4];
    int len1, len2;
    int i, j, table_size;
    int parallel = 1 - lang_id;
    MAPTABREC **table;
    MAPTABREC *ret;

    if (lang_id == 1) {
        table = sc2k;
        table_size = SC2K_SIZE;
    }
    else {
        table = k2sc;
        table_size = K2SC_SIZE;
    }

    len1 = strlen(word1);
    len2 = strlen(word2);
    if (len1 % 3 != 0 || len2 % 3 != 0) {
        return 0;
    }

    source_ch[3] = 0;
    target_ch[3] = 0;
    for (i = 0; i < len1; i += 3) {
        // get a character from word1
        for (j = 0; j < 3; j++) {
            source_ch[j] = word1[i + j];
        }
        ret = hash_search_mapping_table(table, source_ch);
        if (ret) {
            if (compare_characters(ret->corresponding, word2) > 0) {
                return 1;
            }
        }
    }
    return 0;
}
