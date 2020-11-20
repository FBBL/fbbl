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

#include "bkw_step_parameters.h"
#include "syndrome_decoding.h"
#include "workplace_localization.h"
#include "storage_file_utilities.h"

static struct
{
    int q;
    codingType ct;
    union
    {
        intpair *_2;
        inttriplet *_3;
        intquad *_4;
    } table;
} syndrome_repo = { .q = -1, .ct = blockCode_21, .table._2 = NULL };

/* check if corresponding syndrome table has been already been loaded (read from file into memory buffer) */
int is_syndrome_decoding_table_loaded(int q, codingType ct)
{
    if (syndrome_repo.q != q) return 0;
    if (ct == concatenatedCode_21_21)
    {
        if (syndrome_repo.ct != blockCode_21) return 0; /* concatenatedCode_21_21 uses tables for [2,1] block code */
    }
    else
    {
        if (syndrome_repo.ct != ct) return 0;
    }
    switch (ct)
    {
    case blockCode_21:
        return syndrome_repo.table._2 != NULL;
    case blockCode_31:
        return syndrome_repo.table._3 != NULL;
    case blockCode_41:
        return syndrome_repo.table._4 != NULL;
    case concatenatedCode_21_21:
        return syndrome_repo.table._2 != NULL;
    default: ; /* intentionally left blank */
    }
    ASSERT_ALWAYS("unhandled case in is_syndrome_decoding_table_loaded (all coding types listed?)");
    return 0;
}

static void free_current_syndrome_table(void)
{
    if (syndrome_repo.q == -1) return; /* no syndrome table loaded */
    switch (syndrome_repo.ct)
    {
    case concatenatedCode_21_21: /* intentional fall-through */
    case blockCode_21:
        FREE(syndrome_repo.table._2);
        syndrome_repo.table._2 = NULL;
        break;
    case blockCode_31:
        FREE(syndrome_repo.table._3);
        syndrome_repo.table._3 = NULL;
        break;
    case blockCode_41:
        FREE(syndrome_repo.table._4);
        syndrome_repo.table._4 = NULL;
        break;
    default:
        ASSERT_ALWAYS("unhandled case in is_syndrome_decoding_table_loaded (all coding types listed?)");
    }
    syndrome_repo.q = -1;
}

static void intpairs_from_file_2_buf(FILE *f, intpair *buf, int numTableEntriesToRead)
{
    int rd = 0; /* num table entries actually read */
    while (!feof(f))
    {
        int read = fread(buf, sizeof(intpair), numTableEntriesToRead-rd+1, f);
        rd += read;
        buf += read;
    }
    ASSERT(rd == numTableEntriesToRead, "Read from file not complete!\n");
}

static void inttriplets_from_file_2_buf(FILE *f, inttriplet *buf, int numTableEntriesToRead)
{
    int rd = 0; /* num table entries actually read */
    while (!feof(f))
    {
        int read = fread(buf, sizeof(inttriplet), numTableEntriesToRead-rd+1, f);
        rd += read;
        buf += read;
    }
    ASSERT(rd == numTableEntriesToRead, "Read from file not complete!\n");
}

static void intquads_from_file_2_buf(FILE *f, intquad *buf, int numTableEntriesToRead)
{
    int rd = 0; /* num table entries actually read */
    while (!feof(f))
    {
        int read = fread(buf, sizeof(intquad), numTableEntriesToRead-rd+1, f);
        rd += read;
        buf += read;
    }
    ASSERT(rd == numTableEntriesToRead, "Read from file not complete!\n");
}

int load_syndrome_decoding_table(int q, codingType ct, int generateIfFileDoesNotExist)
{
    char fileName[256];
    FILE *f;
    int tableSize, entrySize;
    int bl; /* block length */
    int ml; /* message length */

    /* check if (correct) syndrome table is already loaded into memory */
    if (is_syndrome_decoding_table_loaded(q, ct))
    {
        return 0; /* this syndrome table is already loaded and ready to use */
    }

    /* free currently loaded syndrome table (in case another one was loaded previously) */
    free_current_syndrome_table();

    /* make sure that there is a suitable syndrome decoding table on file */
    if (!is_syndrome_decoding_table_generated(q, ct))
    {
        if (!generateIfFileDoesNotExist)
        {
            return 1; /* caller specified to not create new file */
        }
        if (generate_syndrome_decoding_table(q, ct))   /* generate a new one if necessary */
        {
            ASSERT_ALWAYS("could not generate syndrome decoding table");
            return 2; /* generation failed */
        }
        ASSERT(is_syndrome_decoding_table_generated(q, ct), "syndrome decoding table not found on file");
    }

    /* open syndrome table from file */
    switch (ct)
    {
    case concatenatedCode_21_21: /* intentional fall-though */
    case blockCode_21:
        bl = 2;
        ml = 1;
        tableSize = q;
        entrySize = sizeof(intpair);
        break;
    case blockCode_31:
        bl = 3;
        ml = 1;
        tableSize = q * q;
        entrySize = sizeof(inttriplet);
        break;
    case blockCode_41:
        bl = 4;
        ml = 1;
        tableSize = q * q * q;
        entrySize = sizeof(intquad);
        break;
    default:
        ASSERT_ALWAYS("unhandled case in load_syndrome_table (all coding types listed?)");
    }
    sprintf(fileName, LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A LOCAL_PATH_DELIMITER "syndrome_decoding_table_%d%d_%d.dat", bl, ml, q);
    f = fopen(fileName, "rb");
    if (!f)
    {
        ASSERT_ALWAYS("could not find or open suitable syndrome decoding table");
        exit(1);
    }

    /* read into memory buffer */
    switch (ct)
    {
    case concatenatedCode_21_21: /* intentional fall-though */
    case blockCode_21:
        syndrome_repo.ct = blockCode_21;
        syndrome_repo.table._2 = MALLOC(tableSize * entrySize); /* allocate memory for syndrome table buffer p */
        ASSERT(syndrome_repo.table._2, "Allocation failed!\n");
        intpairs_from_file_2_buf(f, syndrome_repo.table._2, tableSize); /* read syndrome table from file into buffer */
        break;
    case blockCode_31:
        syndrome_repo.ct = blockCode_31;
        syndrome_repo.table._3 = MALLOC(tableSize * entrySize); /* allocate memory for syndrome table buffer p */
        ASSERT(syndrome_repo.table._3, "Allocation failed!\n");
        inttriplets_from_file_2_buf(f, syndrome_repo.table._3, tableSize); /* read syndrome table from file into buffer */
        break;
    case blockCode_41:
        syndrome_repo.ct = blockCode_41;
        syndrome_repo.table._4 = MALLOC(tableSize * entrySize); /* allocate memory for syndrome table buffer p */
        ASSERT(syndrome_repo.table._4, "Allocation failed!\n");
        intquads_from_file_2_buf(f, syndrome_repo.table._4, tableSize); /* read syndrome table from file into buffer */
        break;
    default:
        ASSERT_ALWAYS("unhandled case in load_syndrome_table (all coding types listed?)");
    }
    syndrome_repo.q = q;
    fclose(f);
    return 0; /* syndrome decoding table successfully loaded into memory */
}

void freeSyndromeTables2(intpair *syndrome_table)
{
    if (syndrome_table)
    {
        FREE(syndrome_table);
        syndrome_table = NULL;
    }
    else
    {
        ASSERT_ALWAYS("Unsupported sorting method, can't free syndrome table!\n");
    }
//  if (syndrome_table_21_2053) { FREE(syndrome_table_21_2053); syndrome_table_21_2053 = NULL; }
//  if (syndrome_table_21_16411) { FREE(syndrome_table_21_16411); syndrome_table_21_16411 = NULL; }
//  if (syndrome_table_31_631) { FREE(syndrome_table_31_631); syndrome_table_31_631 = NULL; }
//  if (syndrome_table_31_2053) { FREE(syndrome_table_31_2053); syndrome_table_31_2053 = NULL; }
//  if (syndrome_table_31_16411) { FREE(syndrome_table_31_16411); syndrome_table_31_16411 = NULL; }
//  if (syndrome_table_41_631) { FREE(syndrome_table_41_631); syndrome_table_41_631 = NULL; }
//  if (syndrome_table_41_2053) { FREE(syndrome_table_41_2053); syndrome_table_41_2053 = NULL; }
//  if (syndrome_table_32_631) { FREE(syndrome_table_32_631); syndrome_table_32_631 = NULL; }
//  if (syndrome_table_32_2053) { FREE(syndrome_table_32_2053); syndrome_table_32_2053 = NULL; }
}

#if 0
q=631
  [2,1] code, min for (1,73), variance 101.255151 (optimal)
      [3,1] code, min for (1,205,303), variance 1277.293512 (optimal)
          [4,1] code, min for (1,126,9,332), variance 4951.529059 (probably not optimal)

              q=1601
                [2,1] code, min for (1,335), variance 256.798251 (optimal)
                    [3,1] code, min for (1,8,118), variance 5333.082583 (probably not optimal)

                        q=2053
                          [2,1] code, min for (1,175), variance 329.241111 (optimal)
                              [3,1] code, min for (1,14,443), variance 6151.251737 (probably not optimal)
                                  [4,1] code, min for (1,1462,1140), variance 29107.733089 (probably not optimal)

                                      q=16411
                                        [2,1] code, min for (1,2584), variance 2631.986716 (optimal)
                                            [3,1] code, min for (1,3872,7445), variance 99166.247216 (probably not optimal)
#endif

                                                /* [2,1]-codes
                                                 * G = [1   30] for q=  101, H = [  -30 1]
                                                 * G = [1   73] for q=  631, H = [  -73 1]
                                                 * G = [1  335] for q= 1601, H = [ -335 1]
                                                 * G = [1  175] for q= 2053, H = [ -175 1]
                                                 * G = [1 2584] for q=16411, H = [-2584 1]
                                                 */

                                                /* variance optimal [2,1]-code has G=[1 73] for q=631 *
                                                 * variance optimal [2,1]-code has G=[1 175] for q=2053 *
                                                 * variance optimal [2,1]-code has G=[1 2584] for q=16411 */

                                                static int g2_2[5] = {30, 73, 335, 175, 2584};

static int g2fromq_2(int q)
{
    switch (q)
    {
    case 101:
        return g2_2[0];
    case 631:
        return g2_2[1];
    case 1601:
        return g2_2[2];
    case 2053:
        return g2_2[3];
    case 16411:
        return g2_2[4];
    default:
        ASSERT_ALWAYS("Unsupported q!\n");
    }
    return -1;
}

/* [3,1]-codes
 * G = [1 66 88] for q=101, H = [35 1 0]
 *                              [13 0 1]
 * G = [1 205 303] for q=631, H = [-205 1 0]
 *                                [-303 0 1]
 * G = [1 8 118] for q=1601, H = [ -8  1 0]
 *                               [-118 0 1]
 * G = [1 14 443] for q=2053, H = [ -14 1 0]
 *                                [-443 0 1]
 * G = [1 3872 7445] for q=16411, H = [-3872 1 0]
 *                                    [-7445 0 1]
 */
static int g2_3[5] = {35, 205, 8, 14, 3872};
static int g3_3[5] = {13, 303, 118, 443, 7445};

static int g2fromq_3(int q)
{
    switch (q)
    {
    case 101:
        return g2_3[0];
    case 631:
        return g2_3[1];
    case 1601:
        return g2_3[2];
    case 2053:
        return g2_3[3];
    case 16411:
        return g2_3[4];
    default:
        ASSERT_ALWAYS("Unsupported q!\n");
    }
    return -1;
}

static int g3fromq_3(int q)
{
    switch (q)
    {
    case 101:
        return g3_3[0];
    case 631:
        return g3_3[1];
    case 1601:
        return g2_3[2];
    case 2053:
        return g3_3[3];
    case 16411:
        return g3_3[4];
    default:
        ASSERT_ALWAYS("Unsupported q!\n");
    }
    return -1;
}


/* [4,1]-code with G=[1 g1 g2 g3] and
 * H = [-g1 1 0 0]
 *     [-g2 0 1 0]
 *     [-g3 0 0 1]
 *
 * G = [1 20 6 44] for q=53, H = [ -20 1 0 0] - var = 121.057423
 *                               [  -6 0 1 0]
 *                               [ -44 0 0 1]
 * G = [1 69 7 91] for q=101, H = [ -69 1 0 0] - var = 315.023699
 *                                [  -7 0 1 0]
 *                                [ -91 0 0 1]
 * G = [1 126 9 332] for q=631, H = [-126 1 0 0]
 *                                  [  -9 0 1 0]
 *                                  [-332 0 0 1]
 * G = [1 123 456 789] for q=2053, H = [-123 1 0 0]
 *                                     [-456 0 1 0]
 *                                     [-789 0 0 1]
 */

static int g2_4[4] = {20, 69, 126, 123};
static int g3_4[4] = {6, 7, 9, 456};
static int g4_4[4] = {44, 91, 332, 789};

static int g2fromq_4(int q)
{
    switch (q)
    {
    case 53:
        return g2_4[0];
    case 101:
        return g2_4[1];
    case 631:
        return g2_4[2];
    case 2053:
        return g2_4[3];
    default:
        ASSERT_ALWAYS("Unsupported q!\n");
    }
    return -1;
}

static int g3fromq_4(int q)
{
    switch (q)
    {
    case 53:
        return g3_4[0];
    case 101:
        return g3_4[1];
    case 631:
        return g3_4[2];
    case 2053:
        return g3_4[3];;
    default:
        ASSERT_ALWAYS("Unsupported q!\n");
    }
    return -1;
}

static int g4fromq_4(int q)
{
    switch (q)
    {
    case 53:
        return g4_4[0];
    case 101:
        return g4_4[1];
    case 631:
        return g4_4[2];
    case 2053:
        return g4_4[3];
    default:
        ASSERT_ALWAYS("Unsupported q!\n");
    }
    return -1;
}

static int syndrome_2_1(int q, int g2, int a1, int a2)
{
    ASSERT(0 <= a1 && a1 < q, "Unexpected value for a1!\n");
    ASSERT(0 <= a2 && a2 < q, "Unexpected value for a2!\n");
    return ((q-g2) * a1 + a2) % q;
}

static int syndrome_3_1(int *s1, int *s2, int q, int g2, int g3, int a1, int a2, int a3)
{
    int syndrome1, syndrome2;
    ASSERT(0 <= a1 && a1 < q, "Unexpected value for a1!\n");
    ASSERT(0 <= a2 && a2 < q, "Unexpected value for a2!\n");
    ASSERT(0 <= a3 && a3 < q, "Unexpected value for a3!\n");
    syndrome1 = ((q-g2) * a1 + a2) % q;
    syndrome2 = ((q-g3) * a1 + a3) % q;
    if (s1)
    {
        *s1 = syndrome1;
    }
    if (s2)
    {
        *s2 = syndrome2;
    }
    return q * syndrome2 + syndrome1;
}

static int syndrome_4_1(int *s1, int *s2, int *s3, int q, int g2, int g3, int g4, int a1, int a2, int a3, int a4)
{
    int syndrome1, syndrome2, syndrome3;
    ASSERT(0 <= a1 && a1 < q, "Unexpected value for a1!\n");
    ASSERT(0 <= a2 && a2 < q, "Unexpected value for a2!\n");
    ASSERT(0 <= a3 && a3 < q, "Unexpected value for a3!\n");
    ASSERT(0 <= a4 && a4 < q, "Unexpected value for a4!\n");
    syndrome1 = ((q-g2) * a1 + a2) % q;
    syndrome2 = ((q-g3) * a1 + a3) % q;
    syndrome3 = ((q-g4) * a1 + a4) % q;
    if (s1)
    {
        *s1 = syndrome1;
    }
    if (s2)
    {
        *s2 = syndrome2;
    }
    if (s3)
    {
        *s3 = syndrome3;
    }
    return syndrome1 + q * syndrome2 + q * q * syndrome3;
}

int is_code_word_2_1(int q, int a1, int a2)
{
    ASSERT(0 <= a1 && a1 < q, "Unexpected value for a1!\n");
    ASSERT(0 <= a2 && a2 < q, "Unexpected value for a2!\n");
    return syndrome_2_1(q, g2fromq_2(q), a1, a2) == 0;
}

int is_code_word_3_1(int q, int a1, int a2, int a3)
{
    ASSERT(0 <= a1 && a1 < q, "Unexpected value for a1!\n");
    ASSERT(0 <= a2 && a2 < q, "Unexpected value for a2!\n");
    ASSERT(0 <= a3 && a3 < q, "Unexpected value for a3!\n");
    return syndrome_3_1(NULL, NULL, q, g2fromq_3(q), g3fromq_3(q), a1, a2, a3) == 0;
}

int is_code_word_4_1(int q, int a1, int a2, int a3, int a4)
{
    ASSERT(0 <= a1 && a1 < q, "Unexpected value for a1!\n");
    ASSERT(0 <= a2 && a2 < q, "Unexpected value for a2!\n");
    ASSERT(0 <= a3 && a3 < q, "Unexpected value for a3!\n");
    ASSERT(0 <= a4 && a4 < q, "Unexpected value for a4!\n");
    return syndrome_4_1(NULL, NULL, NULL, q, g2fromq_4(q), g3fromq_4(q), g4fromq_4(q), a1, a2, a3, a4) == 0;
}

void closest_code_word_2_1(int *c1, int *c2, int q, int a1, int a2)
{
    int syndrome = syndrome_2_1(q, g2fromq_2(q), a1, a2);
    int e1, e2; /* error vector */
    ASSERT(is_syndrome_decoding_table_loaded(q, blockCode_21), "syndrome table not loaded ([21]-code)");
    e1 = syndrome_repo.table._2[syndrome].t1;
    e2 = syndrome_repo.table._2[syndrome].t2;
    *c1 = (a1 - e1 + q) % q;
    *c2 = (a2 - e2 + q) % q;
    ASSERT(is_code_word_2_1(q, *c1, *c2), "Not a code word!\n");
}

void closest_code_word_3_1(int *c1, int *c2, int *c3, int q, int a1, int a2, int a3)
{
    int syndrome = syndrome_3_1(NULL, NULL, q, g2fromq_3(q), g3fromq_3(q), a1, a2, a3);
    int e1, e2, e3; /* error vector */
    ASSERT(is_syndrome_decoding_table_loaded(q, blockCode_31), "syndrome table for 31-code not loaded");
    e1 = syndrome_repo.table._3[syndrome].t1;
    e2 = syndrome_repo.table._3[syndrome].t2;
    e3 = syndrome_repo.table._3[syndrome].t3;
    *c1 = (a1 - e1 + q) % q;
    *c2 = (a2 - e2 + q) % q;
    *c3 = (a3 - e3 + q) % q;
    ASSERT(is_code_word_3_1(q, *c1, *c2, *c3), "Not a code word!\n");
}

void closest_code_word_4_1(int *c1, int *c2, int *c3, int *c4, int q, int a1, int a2, int a3, int a4)
{
    int syndrome = syndrome_4_1(NULL, NULL, NULL, q, g2fromq_4(q), g3fromq_4(q), g4fromq_4(q), a1, a2, a3, a4);
    int e1, e2, e3, e4; /* error vector */
    ASSERT(is_syndrome_decoding_table_loaded(q, blockCode_41), "syndrome table for 41-code not loaded");
    e1 = syndrome_repo.table._4[syndrome].t1;
    e2 = syndrome_repo.table._4[syndrome].t2;
    e3 = syndrome_repo.table._4[syndrome].t3;
    e4 = syndrome_repo.table._4[syndrome].t4;
    *c1 = (a1 - e1 + q) % q;
    *c2 = (a2 - e2 + q) % q;
    *c3 = (a3 - e3 + q) % q;
    *c4 = (a4 - e4 + q) % q;
    //printf("(%5d,%5d,%5d,%5d) [%s] = (%5d,%5d,%5d,%5d) [%s] - (%5d,%5d,%5d,%5d)\n", *c1, *c2, *c3, *c4, is_code_word_4_1(q, *c1, *c2, *c3, *c4) ? "YES" : " NO", a1, a2, a3, a4, is_code_word_4_1(q, a1, a2, a3, a4) ? "YES" : " NO", e1, e2, e3, e4);
    ASSERT(is_code_word_4_1(q, *c1, *c2, *c3, *c4), "Not a code word!\n");
}

int generate_syndrome_decoding_table_2_1_code(int q)
{
    char fileName[256];
    FILE *f;
    int i, e1, e2, g2 = g2fromq_2(q), wr;
    int *num = CALLOC(q, sizeof(int));
    int *minvar = MALLOC(q * sizeof(int));
    intpair *t = MALLOC(q * sizeof(intpair));

    sprintf(fileName, "%s/syndrome_decoding_table_21_%d.dat", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A, q);
    f = fopen(fileName, "wb");
    ASSERT(f, "Could not open output file!\n");
    ASSERT(num, "Could not allocate num!\n");
    ASSERT(minvar, "Could not allocate minvar!\n");
    ASSERT(t, "Could not allocate t!\n");
    for (i=0; i<q; i++)
    {
        minvar[i] = -1;
        t[i].t1 = -1;
        t[i].t2 = -1;
    }

    for (e1=0; e1<q; e1++)
    {
        for (e2=0; e2<q; e2++)   // for all error vectors
        {
            int syndrome = syndrome_2_1(q, g2, e1, e2);
            int u = e1 > q/2 ? q-e1 : e1;
            int v = e2 > q/2 ? q-e2 : e2;
            int var = u*u + v*v; // squared distance
            num[syndrome]++;

            if (minvar[syndrome] == -1 || var < minvar[syndrome])
            {
                t[syndrome].t1 = e1;
                t[syndrome].t2 = e2;
                minvar[syndrome] = var;
            }
        }
    }

#ifdef _DEBUG
    for (i=0; i<q; i++)
    {
        ASSERT(num[i] == q, "Unexpected class size!\n");
        ASSERT(0 <= t[i].t1 && t[i].t1 < q, "Unexpected value for t1!\n");
        ASSERT(0 <= t[i].t2 && t[i].t2 < q, "Unexpected value for t2!\n");
//    printf("syndrome %4d: t1=%4d t2=%4d minvar=%4d num=%4d\n", i, t1[i], t2[i], minvar[i], num[i]);
    }
#endif
    wr = fwrite(t, sizeof(intpair), q, f);
    ASSERT(wr == q, "Write to file not completed!\n");
    fflush(f);
    if (f) fclose(f);

    FREE(t);
    FREE(minvar);
    FREE(num);
    return 0;
}

int generate_syndrome_decoding_table_3_1_code(int q)
{
    char fileName[256];
    FILE *f;
    int i, e1, e2, e3, g2 = g2fromq_3(q), g3 = g3fromq_3(q), wr;
    int numSyndromes = q * q;
    int *num = CALLOC(numSyndromes, sizeof(int));
    int *minvar = MALLOC(numSyndromes * sizeof(int));
    inttriplet *t = MALLOC(numSyndromes * sizeof(inttriplet));

    sprintf(fileName, "%s/syndrome_decoding_table_31_%d.dat", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A, q);
    f = fopen(fileName, "wb");
    ASSERT(f, "Could not open output file!\n");

    ASSERT(num, "Could not allocate num!\n");
    ASSERT(minvar, "Could not allocate minvar!\n");
    ASSERT(t, "Could not allocate t!\n");
    for (i=0; i<numSyndromes; i++)
    {
        minvar[i] = -1;
        t[i].t1 = -1;
        t[i].t2 = -1;
        t[i].t3 = -1;
    }

    for (e1=0; e1<q; e1++)
    {
        for (e2=0; e2<q; e2++)
        {
            for (e3=0; e3<q; e3++)   // for all error vectors
            {
                int syndrome = syndrome_3_1(NULL, NULL, q, g2, g3, e1, e2, e3);
                int u = e1 > q/2 ? q-e1 : e1;
                int v = e2 > q/2 ? q-e2 : e2;
                int w = e3 > q/2 ? q-e3 : e3;
                int var = u*u + v*v + w*w; // squared distance
                num[syndrome]++;

                if (minvar[syndrome] == -1 || var < minvar[syndrome])
                {
                    t[syndrome].t1 = e1;
                    t[syndrome].t2 = e2;
                    t[syndrome].t3 = e3;
                    minvar[syndrome] = var;
                }
            }
        }
    }

#ifdef _DEBUG
    for (i=0; i<numSyndromes; i++)
    {
        ASSERT(num[i] == q, "Unexpected class size!\n");
        ASSERT(0 <= t[i].t1 && t[i].t1 < q, "Unexpected value for t1!\n");
        ASSERT(0 <= t[i].t2 && t[i].t2 < q, "Unexpected value for t2!\n");
        ASSERT(0 <= t[i].t3 && t[i].t3 < q, "Unexpected value for t3!\n");
//    printf("syndrome %4d: t1=%4d t2=%4d t2=%4d minvar=%4d num=%4d\n", i, t1[i], t2[i], t3[i], minvar[i], num[i]);
    }
#endif
    wr = fwrite(t, sizeof(inttriplet), numSyndromes, f);
    ASSERT(wr == numSyndromes, "Write to file not completed!\n");
    fflush(f);
    if (f) fclose(f);

    FREE(t);
    FREE(minvar);
    FREE(num);
    return 0;
}

int generate_syndrome_decoding_table_4_1_code(int q, int maxComponentError)
{
    char fileName[256];
    FILE *f;
    int i, e1, e2, e3, e4, wr, maxLoop = (0 < maxComponentError && maxComponentError < q) ? maxComponentError : q;
    int u, v, w, x, var, syndrome;
    int numSyndromes = q * q * q;
    intquad *t = MALLOC(numSyndromes * sizeof(intquad));
    int *minvar = MALLOC(numSyndromes * sizeof(int));
    int numDeterminedTableEntries = 0;
    int g2 = g2fromq_4(q);
    int g3 = g3fromq_4(q);
    int g4 = g4fromq_4(q);
    long double expected_variance = 0;

    printf("Generating syndrome decoding table for [4,1]-code [1 %d %d %d] for q=%d\n", g2, g3, g4, q);

    sprintf(fileName, "%s/syndrome_decoding_table_41_%d.dat", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A, q);
    f = fopen(fileName, "wb");
    ASSERT(f, "Could not open output file!\n");

    ASSERT(minvar, "Could not allocate minvar!\n");
    ASSERT(t, "Could not allocate t!\n");
    for (i=0; i<numSyndromes; i++)
    {
        minvar[i] = -1;
        t[i].t1 = -1;
        t[i].t2 = -1;
        t[i].t3 = -1;
        t[i].t4 = -1;
    }

    for (e1=0; e1<maxLoop; e1++)
    {
        for (e2=0; e2<maxLoop; e2++)
        {
            for (e3=0; e3<maxLoop; e3++)
            {
                for (e4=0; e4<maxLoop; e4++)   // for all error vectors with components up to maxLoop
                {
                    syndrome = syndrome_4_1(NULL, NULL, NULL, q, g2, g3, g4, e1, e2, e3, e4);
                    u = e1 > q/2 ? q-e1 : e1;
                    v = e2 > q/2 ? q-e2 : e2;
                    w = e3 > q/2 ? q-e3 : e3;
                    x = e4 > q/2 ? q-e4 : e4;
                    var = u*u + v*v + w*w + x*x; // squared distance

                    if (minvar[syndrome] == -1 || var < minvar[syndrome])
                    {
                        t[syndrome].t1 = e1;
                        t[syndrome].t2 = e2;
                        t[syndrome].t3 = e3;
                        t[syndrome].t4 = e4;
                        minvar[syndrome] = var;
                    }
                }
            }
        }
    }

//   check if all entries are determined
    for (i=0; i<numSyndromes; i++)
    {
        if (minvar[i] > 0)   // value assigned
        {
            numDeterminedTableEntries++;
        }
    }
    printf("Num determined entries: %d/%d (%.2f%%)\n", numDeterminedTableEntries, numSyndromes, numDeterminedTableEntries/(double)numSyndromes*100);
    if (numDeterminedTableEntries < numSyndromes)
    {
        printf("Warning: Table incomplete! %d < %d\n", numDeterminedTableEntries, numSyndromes);
    }
    else
    {
        printf("Table complete!\n");
    }

    wr = fwrite(t, sizeof(intquad), numSyndromes, f);
    ASSERT(wr == numSyndromes, "Write to file not completed!\n");
    fflush(f);
    if (f) fclose(f);

    for (i=0; i<numSyndromes; i++)
    {
        if (minvar[i] > 0)
        {
            expected_variance += minvar[i];
        }
    }
    expected_variance /= (long double)numDeterminedTableEntries;
    printf("expected variance = %Lf\n", expected_variance);

    FREE(t);
    FREE(minvar);
    return 0;
}

// Check if the syndrome table has been already generated
int is_syndrome_decoding_table_generated(int q, codingType ct)
{
    char fileName[256];
    int bl; /* block length */
    switch (ct)
    {
    case concatenatedCode_21_21: /* intentional fall-through */
    case blockCode_21:
        bl = 2;
        break;
    case blockCode_31:
        bl = 3;
        break;
    case blockCode_41:
        bl = 4;
        break;
    default:
        ASSERT_ALWAYS("coding type not handled (all types listed?)");
    }
    sprintf(fileName, "%s/syndrome_decoding_table_%d%d_%d.dat", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A, bl, 1, q);
    return fileExists(fileName);
}

int generate_syndrome_decoding_table(int q, codingType ct)
{
    switch (ct)
    {
    case concatenatedCode_21_21: /* intentional fall-through */
    case blockCode_21:
        return generate_syndrome_decoding_table_2_1_code(q);
    case blockCode_31:
        return generate_syndrome_decoding_table_3_1_code(q);
    case blockCode_41:
        return generate_syndrome_decoding_table_4_1_code(q, q);
    default: /* intentionally left empty */
        ;
    }
    return -1;
}
