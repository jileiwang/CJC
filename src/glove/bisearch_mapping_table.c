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
#include "bisearch_mapping_table.h"


MAPTABREC **k2sc, **sc2k;

/**
 * Compare 2 mapping table record, used in qsort function.
 * notice the type of a and b are (MAPTABREC **)
 */
int compare_original(const void *a, const void *b) {
    return strcmp((*(MAPTABREC **)a)->original, (*(MAPTABREC **)b)->original);
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
    char *original, *corresponding;
    MAPTABREC *table_record;
    MAPTABREC **table;
    table = (MAPTABREC **)malloc( sizeof(MAPTABREC *) * table_size );

    f = fopen(filename, "rb");
    for (i = 0; i < table_size; i++) {      
        table[i] = (MAPTABREC *)malloc( sizeof(MAPTABREC) );
        get_characters(table[i]->original, f);
        get_characters(table[i]->corresponding, f);
    }
    fclose(f);

    qsort(&table[0], table_size, sizeof(MAPTABREC*), compare_original);

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
 * Binary Search a character from the mapping table
 */
int bisearch_table(MAPTABREC **table, char *ch, int left, int right) {
    int mid, cmp;
    if (right < left) {
        return -1;
    }
    mid = left + (right - left) / 2;
    
    cmp = strcmp(table[mid]->original, ch);

    if (cmp == 0) {
        return mid;
    }
    else if (cmp > 0) {
        return bisearch_table(table, ch, left, mid - 1);
    }
    else {
        return bisearch_table(table, ch, mid + 1, right);
    }
}

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
    int i, j, pos, table_size;
    int parallel = 1 - lang_id;
    MAPTABREC **table;

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
        pos = bisearch_table(table, source_ch, 0, table_size - 1);
        if (pos >= 0) {
            if (compare_characters(table[pos]->corresponding, word2) > 0) {
                return 1;
            }
        }
    }
    return 0;
}