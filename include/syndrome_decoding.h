/*  This file is part of FBBL (File-Based BKW for LWE).
 *
 *  FBBL is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  FBBL is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Nome-Programma.  If not, see <http://www.gnu.org/licenses/>
 */

#ifndef SYNDROME_DECODING_H
#define SYNDROME_DECODING_H

#include "bkw_step_parameters.h"
#include "assert_utils.h"
#include "memory_utils.h"
#include "lwe_sorting.h"
#include <stdio.h>

typedef struct
{
    int t1;
    int t2;
} intpair;

typedef struct
{
    int t1;
    int t2;
    int t3;
} inttriplet;

typedef struct
{
    short t1;
    short t2;
    short t3;
    short t4;
} intquad;

//void load_syndrome_table_21(int q);
//void load_syndrome_table_31(int q);
//void load_syndrome_table_41(int q);
int load_syndrome_decoding_table(int q, codingType ct, int generateIfFileDoesNotExist);
void freeSyndromeTables(void);

int generate_syndrome_decoding_table(int q, codingType ct);
int is_syndrome_decoding_table_generated(int q, codingType ct);
int is_syndrome_decoding_table_loaded(int q, codingType ct);
int is_code_word_2_1(int q, int a1, int a2);
int is_code_word_3_1(int q, int a1, int a2, int a3);
int is_code_word_4_1(int q, int a1, int a2, int a3, int a4);
void closest_code_word_2_1(int *c1, int *c2, int q, int a1, int a2);
void closest_code_word_3_1(int *c1, int *c2, int *c3, int q, int a1, int a2, int a3);
void closest_code_word_4_1(int *c1, int *c2, int *c3, int *c4, int q, int a1, int a2, int a3, int a4);

#endif
