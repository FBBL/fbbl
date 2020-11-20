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

#include "position_values_2_category_index.h"
#include "assert_utils.h"
#include "memory_utils.h"
#include "bkw_step_parameters.h"
#include "syndrome_decoding.h"
#include <math.h>

/* plain BKW with 2 positions */

static struct
{
    int q;
    u64 **table;
} lookup_plain_bkw_2_positions;

static u64 internal_position_values_2_category_index_plain_bkw(int q, int p1, int p2)
{
    /* case 00 */
    if (p1 == 0 && p2 == 0)
    {
        return 0;
    }

    /* case p1 == 0 */
    if (p1 == 0)
    {
        if (p2 <= (q - 1)/2)
        {
            return 2*p2 - 1;
        }
        else
        {
            return 2*(q - p2);
        }
    }

    /* remaining cases */
    if (p1 <= (q - 1)/2)
    {
        return (2*p1 - 1)*q + 2*p2;
    }
    else
    {
        if (p2 == 0)
        {
            return (2*q - 1 - 2*p1)*q + 1;
        }
        else
        {
            return (2*q - 1 - 2*p1)*q + 2*(q - p2) + 1;
        }
    }

    ASSERT_ALWAYS("Case not handled!\n");
}

static void generate_table_plain_bkw_2_positions(int q)
{
    if (lookup_plain_bkw_2_positions.q == q)
    {
        return; /* table already generated */
    }
    lookup_plain_bkw_2_positions.table = MALLOC(q * sizeof(u64*));
    ASSERT(lookup_plain_bkw_2_positions.table, "memory allocation failed");
    for (int i=0; i<q; i++)
    {
        lookup_plain_bkw_2_positions.table[i] = MALLOC(q * sizeof(u64));
        ASSERT(lookup_plain_bkw_2_positions.table[i], "memory allocation failed");
        for (int j=0; j<q; j++)
        {
            lookup_plain_bkw_2_positions.table[i][j] = internal_position_values_2_category_index_plain_bkw(q, i, j);
        }
    }
    lookup_plain_bkw_2_positions.q = q;
}

void free_table_plain_bkw_2_positions(void)
{
    if (lookup_plain_bkw_2_positions.q == -1)
    {
        return; /* no table */
    }
    for (int i=0; i<lookup_plain_bkw_2_positions.q; i++)
    {
        FREE(lookup_plain_bkw_2_positions.table[i]);
    }
    FREE(lookup_plain_bkw_2_positions.table);
    lookup_plain_bkw_2_positions.table = NULL;
    lookup_plain_bkw_2_positions.q = -1;
}

u64 position_values_2_category_index_plain_bkw(int q, short *a)
{
    int p1 = a[0]; /* not passing sample so not using columnValue() */
    int p2 = a[1]; /* not passing sample so not using columnValue() */
    /* use lookup table if available */
    if (lookup_plain_bkw_2_positions.table && lookup_plain_bkw_2_positions.q == q)   /* remove check for q for increased speed (is check superfluous)? */
    {
        return lookup_plain_bkw_2_positions.table[p1][p2];
    }
    /* if not, try to generate lookup table */
    generate_table_plain_bkw_2_positions(q);
    /* now try to use lookup table again */
    if (lookup_plain_bkw_2_positions.table && lookup_plain_bkw_2_positions.q == q)
    {
        return lookup_plain_bkw_2_positions.table[p1][p2];
    }
    /* if usage and generation of lookup table fails, use old school calculation */
    return internal_position_values_2_category_index_plain_bkw(q, p1, p2);
}

void category_index_2_position_values_plain_bkw(int q, u64 category_index, short *a)
{
    int p1, p2;
    if (category_index == 0)
    {
        p1 = 0;
        p2 = 0;
    }
    else if ((int)category_index < q)
    {
        p1 = 0;
        if (category_index % 2 == 1)
        {
            p2 = (category_index + 1)/2;
        }
        else
        {
            p2 = q - category_index/2;
        }
    }
    else
    {
        int box = (category_index - q)/(2*q) + 1; // Calculates which box we are in.
        int innerIndex = category_index - (2*box - 1)*q; // Calculates the index within the box
        if (innerIndex % 2 == 0)
        {
            p1 = box;
            p2 = innerIndex/2;
        }
        else
        {
            p1 = q - box;
            p2 = (q - (innerIndex + 1)/2 + 1) % q;
        }
    }
    a[0] = p1;
    a[1] = p2;
}

/* LMS */

static struct
{
    int numPositions;
    int c;
    u64 *table;
    /* we need to keep track of which categories that are singletons.
       if c is odd, then the only singleton category is at index 0 (zero).
       if c is even there are 2^(numPositions) singleton categories
    */
    int *singletonCategoryIndices;
} lookup_lms;

/* 1 position */
static u64 *createLMStableLevel1(int c)
{
    u64 *indexTable1 = MALLOC(c * sizeof(u64));
    indexTable1[0] = 0;

    /* The new solution */
    for (int p1 = 1; 2*p1 <= c; p1++)
    {
        indexTable1[p1] = 2*p1 - 1;
    }

    for (int p1 = c/2 + 1; p1 < c; p1++)
    {
        indexTable1[p1] = 2*(c - p1);
    }

    return indexTable1;
}

/* 2 positions */
static u64 *createLMStableLevel2(int c)
{
    u64 *indexTable1 = createLMStableLevel1(c); /* Index table for 1 position */
    u64 *indexTable2 = MALLOC(c*c * sizeof(u64));

    /* The trivial positions with p2 = 0 */
    for (int p1 = 0; p1 < c; p1++)
    {
        indexTable2[p1] = indexTable1[p1];
    }

    /* The positions with small p2 values */
    for (int p2 = 1; p2*2 < c; p2++)
    {
        for (int p1 = 0; p1 < c; p1++)
        {
            int curPos = p1 + p2*c; /* Current position in the index table */
            indexTable2[curPos] = (2*p2 - 1)*c + 2*indexTable1[p1];
        }
    }

    /* Add an extra part in the cases where q is even */
    int qHalf = c/2;
    if (qHalf*2 == c)
    {
        for (int p1 = 0; p1 < c; p1++)
        {
            int curPos = p1 + qHalf*c; /* Current position in the index table */
            indexTable2[curPos] = (2*qHalf - 1)*c + indexTable1[p1];
        }
    }


    /* The positions with large p2 values */
    for (int p2 = c/2 + 1; p2 < c; p2++)
    {
        for (int p1 = 0; p1 < c; p1++)
        {
            int curPos = p1 + p2*c; /* Current position in the index table */
            indexTable2[curPos] = (2*(c - p2) - 1)*c + 1 + 2*indexTable1[(c - p1)%c];
        }
    }

    FREE(indexTable1);
    return indexTable2;
}

/* 3 positions */
static u64 *createLMStableLevel3(int c)
{
    u64 *indexTable2 = createLMStableLevel2(c); /* Index table for 2 positions */
    u64 *indexTable3 = MALLOC(c*c*c * sizeof(u64));

    /* The trivial positions with p3 = 0 */
    for (int p2 = 0; p2 < c; p2++)
    {
        for (int p1 = 0; p1 < c; p1++)
        {
            int curPos = p1 + p2*c;
            indexTable3[curPos] = indexTable2[curPos];
        }
    }

    /* The positions with small p3 values */
    for (int p3 = 1; p3*2 < c; p3++)
    {
        for (int p2 = 0; p2 < c; p2++)
        {
            for (int p1 = 0; p1 < c; p1++)
            {
                int oldPos = p1 + p2*c; /* Position in previous index table */
                int curPos = oldPos + p3*c*c; /* Current position in the index table */
                indexTable3[curPos] = (2*p3 - 1)*c*c + 2*indexTable2[oldPos];
            }
        }
    }

    /* Add an extra part in the cases where q is even */
    int qHalf = c/2;
    if (qHalf*2 == c)
    {
        for (int p2 = 0; p2 < c; p2++)
        {
            for (int p1 = 0; p1 < c; p1++)
            {
                int oldPos = p1 + p2*c; /* Position in previous index table */
                int curPos = oldPos + qHalf*c*c; /* Current position in the index table */
                indexTable3[curPos] = (2*qHalf - 1)*c*c + indexTable2[oldPos];
            }
        }
    }

    /* The positions with large p3 values */
    for (int p3 = c/2 + 1; p3 < c; p3++)
    {
        for (int p2 = 0; p2 < c; p2++)
        {
            for (int p1 = 0; p1 < c; p1++)
            {
                int oldPos = p1 + p2*c; /* Position in previous index table */
                int oldPosInverted = (c - p1)%c + ((c - p2)%c) * c; /* Position in previous index table */
                int curPos = oldPos + p3*c*c; /* Current position in the index table */
                indexTable3[curPos] = (2*(c - p3) - 1)*c*c + 1 + 2*indexTable2[oldPosInverted];
            }
        }
    }

    FREE(indexTable2);
    return indexTable3;
}

/* 4 positions */
static u64 *createLMStableLevel4(int c)
{
    u64 *indexTable3 = createLMStableLevel3(c); /* Index table for 3 positions */
    u64 *indexTable4 = MALLOC(c*c*c*c * sizeof(u64));

    /* The trivial positions with p4 = 0 */
    for (int p3 = 0; p3 < c; p3++)
    {
        for (int p2 = 0; p2 < c; p2++)
        {
            for (int p1 = 0; p1 < c; p1++)
            {
                int curPos = p1 + p2*c + p3*c*c;
                indexTable4[curPos] = indexTable3[curPos];
            }
        }
    }

    /* The positions with small p4 values */
    for (int p4 = 1; p4*2 < c; p4++)
    {
        for (int p3 = 0; p3 < c; p3++)
        {
            for (int p2 = 0; p2 < c; p2++)
            {
                for (int p1 = 0; p1 < c; p1++)
                {
                    int oldPos = p1 + p2*c + p3*c*c; /* Position in previous index table */
                    int curPos = oldPos + p4*c*c*c; /* Current position in the index table */
                    indexTable4[curPos] = (2*p4 - 1)*c*c*c + 2*indexTable3[oldPos];
                }
            }
        }
    }

    /* Add an extra part in the cases where q is even */
    int qHalf = c/2;
    if (qHalf*2 == c)
    {
        for (int p3 = 0; p3 < c; p3++)
        {
            for (int p2 = 0; p2 < c; p2++)
            {
                for (int p1 = 0; p1 < c; p1++)
                {
                    int oldPos = p1 + p2*c + p3*c*c; /* Position in previous index table */
                    int curPos = oldPos + qHalf*c*c*c; /* Current position in the index table */
                    indexTable4[curPos] = (2*qHalf - 1)*c*c*c + indexTable3[oldPos];
                }
            }
        }
    }

    /* The positions with large p4 values */
    for (int p4 = c/2 + 1; p4 < c; p4++)
    {
        for (int p3 = 0; p3 < c; p3++)
        {
            for (int p2 = 0; p2 < c; p2++)
            {
                for (int p1 = 0; p1 < c; p1++)
                {
                    int oldPos = p1 + p2*c + p3*c*c; /* Position in previous index table */
                    int oldPosInverted = (c - p1)%c + ((c - p2)%c) * c + ((c - p3)%c) * c*c; /* Position in previous index table */
                    int curPos = oldPos + p4*c*c*c; /* Current position in the index table */
                    indexTable4[curPos] = (2*(c - p4) - 1)*c*c*c + 1 + 2*indexTable3[oldPosInverted];
                }
            }
        }
    }

    FREE(indexTable3);
    return indexTable4;
}

/* 5 positions */
static u64 *createLMStableLevel5(int c)
{
    u64 *indexTable4 = createLMStableLevel4(c); /* Index table for 4 positions */
    u64 *indexTable5 = MALLOC(c*c*c*c*c * sizeof(u64));

    /* The trivial positions with p5 = 0 */
    for (int p4 = 0; p4 < c; p4++)
    {
        for (int p3 = 0; p3 < c; p3++)
        {
            for (int p2 = 0; p2 < c; p2++)
            {
                for (int p1 = 0; p1 < c; p1++)
                {
                    int curPos = p1 + p2*c + p3*c*c + p4*c*c*c;
                    indexTable5[curPos] = indexTable4[curPos];
                }
            }
        }
    }

    /* The positions with small p5 values */
    for (int p5 = 1; p5*2 < c; p5++)
    {
        for (int p4 = 0; p4 < c; p4++)
        {
            for (int p3 = 0; p3 < c; p3++)
            {
                for (int p2 = 0; p2 < c; p2++)
                {
                    for (int p1 = 0; p1 < c; p1++)
                    {
                        int oldPos = p1 + p2*c + p3*c*c + p4*c*c*c; /* Position in previous index table */
                        int curPos = oldPos + p5*c*c*c*c; /* Current position in the index table */
                        indexTable5[curPos] = (2*p5 - 1)*c*c*c*c + 2*indexTable4[oldPos];
                    }
                }
            }
        }
    }

    /* Add an extra part in the cases where q is even */
    int qHalf = c/2;
    if (qHalf*2 == c)
    {
        for (int p4 = 0; p4 < c; p4++)
        {
            for (int p3 = 0; p3 < c; p3++)
            {
                for (int p2 = 0; p2 < c; p2++)
                {
                    for (int p1 = 0; p1 < c; p1++)
                    {
                        int oldPos = p1 + p2*c + p3*c*c + p4*c*c*c; /* Position in previous index table */
                        int curPos = oldPos + qHalf*c*c*c*c; /* Current position in the index table */
                        indexTable5[curPos] = (2*qHalf - 1)*c*c*c*c + indexTable4[oldPos];
                    }
                }
            }
        }
    }

    /* The positions with large p5 values */
    for (int p5 = c/2 + 1; p5 < c; p5++)
    {
        for (int p4 = 0; p4 < c; p4++)
        {
            for (int p3 = 0; p3 < c; p3++)
            {
                for (int p2 = 0; p2 < c; p2++)
                {
                    for (int p1 = 0; p1 < c; p1++)
                    {
                        int oldPos = p1 + p2*c + p3*c*c + p4*c*c*c; /* Position in previous index table */
                        int oldPosInverted = (c - p1)%c + ((c - p2)%c) * c + ((c - p3)%c) * c*c + ((c - p4)%c) * c*c*c; /* Position in previous index table */
                        int curPos = oldPos + p5*c*c*c*c; /* Current position in the index table */
                        indexTable5[curPos] = (2*(c - p5) - 1)*c*c*c*c + 1 + 2*indexTable4[oldPosInverted];
                    }
                }
            }
        }
    }

    FREE(indexTable4);
    return indexTable5;
}

/* 6 positions */
static u64 *createLMStableLevel6(int c)
{
    u64 *indexTable5 = createLMStableLevel5(c); /* Index table for 5 positions */
    u64 *indexTable6 = MALLOC(c*c*c*c*c*c * sizeof(u64));

    /* The trivial positions with p6 = 0 */
    for (int p5 = 0; p5 < c; p5++)
    {
        for (int p4 = 0; p4 < c; p4++)
        {
            for (int p3 = 0; p3 < c; p3++)
            {
                for (int p2 = 0; p2 < c; p2++)
                {
                    for (int p1 = 0; p1 < c; p1++)
                    {
                        int curPos = p1 + p2*c + p3*c*c + p4*c*c*c + p5*c*c*c*c;
                        indexTable6[curPos] = indexTable5[curPos];
                    }
                }
            }
        }
    }

    /* The positions with small p6 values */
    for (int p6 = 1; p6*2 < c; p6++)
    {
        for (int p5 = 0; p5 < c; p5++)
        {
            for (int p4 = 0; p4 < c; p4++)
            {
                for (int p3 = 0; p3 < c; p3++)
                {
                    for (int p2 = 0; p2 < c; p2++)
                    {
                        for (int p1 = 0; p1 < c; p1++)
                        {
                            int oldPos = p1 + p2*c + p3*c*c + p4*c*c*c + p5*c*c*c*c; /* Position in previous index table */
                            int curPos = oldPos + p6*c*c*c*c*c; /* Current position in the index table */
                            indexTable6[curPos] = (2*p6 - 1)*c*c*c*c*c + 2*indexTable5[oldPos];
                        }
                    }
                }
            }
        }
    }

    /* Add an extra part in the cases where q is even */
    int qHalf = c/2;
    if (qHalf*2 == c)
    {
        for (int p5 = 0; p5 < c; p5++)
        {
            for (int p4 = 0; p4 < c; p4++)
            {
                for (int p3 = 0; p3 < c; p3++)
                {
                    for (int p2 = 0; p2 < c; p2++)
                    {
                        for (int p1 = 0; p1 < c; p1++)
                        {
                            int oldPos = p1 + p2*c + p3*c*c + p4*c*c*c + p5*c*c*c*c; /* Position in previous index table */
                            int curPos = oldPos + qHalf*c*c*c*c*c; /* Current position in the index table */
                            indexTable6[curPos] = (2*qHalf - 1)*c*c*c*c*c + indexTable5[oldPos];
                        }
                    }
                }
            }
        }
    }

    /* The positions with large p6 values */
    for (int p6 = c/2 + 1; p6 < c; p6++)
    {
        for (int p5 = 0; p5 < c; p5++)
        {
            for (int p4 = 0; p4 < c; p4++)
            {
                for (int p3 = 0; p3 < c; p3++)
                {
                    for (int p2 = 0; p2 < c; p2++)
                    {
                        for (int p1 = 0; p1 < c; p1++)
                        {
                            int oldPos = p1 + p2*c + p3*c*c + p4*c*c*c + p5*c*c*c*c; /* Position in previous index table */
                            int oldPosInverted = (c - p1)%c + ((c - p2)%c) * c + ((c - p3)%c) * c*c + ((c - p4)%c) * c*c*c + + ((c - p5)%c) * c*c*c*c; /* Position in previous index table */
                            int curPos = oldPos + p6*c*c*c*c*c; /* Current position in the index table */
                            indexTable6[curPos] = (2*(c - p6) - 1)*c*c*c*c*c + 1 + 2*indexTable5[oldPosInverted];
                        }
                    }
                }
            }
        }
    }

    FREE(indexTable5);
    return indexTable6;
}

/* BKW with LMS steps */
static u64 *createLMStable(int q, int p, int numPositions)
{
    int c = q/p + 1; /* The number of possible values for each position */
    switch (numPositions)
    {
    case 1:
        return createLMStableLevel1(c);
    case 2:
        return createLMStableLevel2(c);
    case 3:
        return createLMStableLevel3(c);
    case 4:
        return createLMStableLevel4(c);
    case 5:
        return createLMStableLevel5(c);
    case 6:
        return createLMStableLevel6(c);
    }
    ASSERT_ALWAYS("unsupported parameter numPositions");
    return NULL; /* Unsupported value for b */
}

static void generate_table_lms(int q, int p, int numPositions)
{
    int c = q/p + 1;
    if (lookup_lms.table && lookup_lms.c == c && lookup_lms.numPositions == numPositions)
    {
        return; /* table already built */
    }
    if (lookup_lms.table)
    {
        FREE(lookup_lms.table);
        lookup_lms.table = NULL;
        lookup_lms.numPositions = 0;
        lookup_lms.c = 0;
        if (lookup_lms.singletonCategoryIndices)
        {
            FREE(lookup_lms.singletonCategoryIndices);
            lookup_lms.singletonCategoryIndices = NULL;
        }
    }
    lookup_lms.table = createLMStable(q, p, numPositions);
    if (lookup_lms.table)
    {
        lookup_lms.numPositions = numPositions;
        lookup_lms.c = c;
        if ((lookup_lms.c & 1) == 0)   // if c is even
        {
            lookup_lms.singletonCategoryIndices = MALLOC((1 << numPositions) * sizeof(int));
            ASSERT(lookup_lms.singletonCategoryIndices, "allocation failed");
            /* usage of locals lwe and bkwStepPar is not pretty,
               so perhaps change the API for position_values_2_category_index_lms?
               */
            /* ugliness start */
            lweInstance lwe;
            bkwStepParameters bkwStepPar;
            lwe.q = q;
            bkwStepPar.sorting = LMS;
            bkwStepPar.sortingPar.LMS.p = p;
            bkwStepPar.numPositions = numPositions;
            /* ugliness end */
            short *pn = MALLOC(numPositions * sizeof(short));
            for (u32 i = 0; i < (u32)(1) << numPositions; i++)
            {
//        short pn[MAX_LMS_POSITIONS];
                ASSERT(numPositions <= MAX_LMS_POSITIONS, "insufficient size for pn");
                for (u32 j = 0; j < (u32)numPositions; j++)
                {
                    pn[j] = ((i >> j) & 1) == 0 ? 0 : q/2;
                }
                lookup_lms.singletonCategoryIndices[i] = position_values_2_category_index_lms(&lwe, &bkwStepPar, pn);
            }
            FREE(pn);
        }
        else     /* c is odd */
        {
            lookup_lms.singletonCategoryIndices = MALLOC((1 << numPositions) * sizeof(int));
            ASSERT(lookup_lms.singletonCategoryIndices, "allocation failed");
            for (u32 i = 0; i < (u32)(1) << numPositions; i++)
            {
                lookup_lms.singletonCategoryIndices[i] = 0;
            }
        }
        return;
    }
    ASSERT_ALWAYS("could not create table");
}

/* the number of singleton categories in [a,b) */
int num_lms_singletons_in_category_interval(int a, int b)
{
    int numSingletonsInInterval = 0;
    for (int i=0; i<1<<lookup_lms.numPositions; i++)
    {
        int singletonIndex = lookup_lms.singletonCategoryIndices[i];
        if (a <= singletonIndex && singletonIndex < b)
        {
            numSingletonsInInterval++;
        }
    }
    return numSingletonsInInterval;
}

int is_lms_singleton(u64 categoryIndex)
{
    if (categoryIndex == 0)
    {
        return 1;
    }
    if (!lookup_lms.singletonCategoryIndices)
    {
        ASSERT_ALWAYS("LMS singleton table not built");
        return 0;
    }
    for (int i=0; i<1<<lookup_lms.numPositions; i++)
    {
        if (categoryIndex == (u64)lookup_lms.singletonCategoryIndices[i])
        {
            return 1;
        }
    }
    return 0;
}

int is_smooth_lms_singleton(u64 categoryIndex, u64 numCategories)
{

    if(numCategories & 1)
    {
        return (categoryIndex == 0);
    }
    return 0;
}


int is_singleton(bkwStepParameters *bkwStepPar, u64 categoryIndex, u64 numCategories)
{
    switch (bkwStepPar->sorting)
    {
    case plainBKW:
        return categoryIndex == 0;
    case LMS:
        return is_lms_singleton(categoryIndex);
    case smoothLMS:
        return is_smooth_lms_singleton(categoryIndex, numCategories);
    case codedBKW:
        switch(bkwStepPar->sortingPar.CodedBKW.ct)
        {
        case blockCode_21:
            return categoryIndex == 0;
        case blockCode_31:
            return categoryIndex == 0;
        case blockCode_41:
            return categoryIndex == 0;
        case concatenatedCode_21_21:
            return categoryIndex == 0;
        default:
            ASSERT_ALWAYS("is_singleton not implemented for this choice of coded bkw");
            break;
        }
        break;
    default:
        break;
    }
    ASSERT_ALWAYS("unsupported sorting");
    return -1;
}

/* ACCESSING INDICES FROM LMS TABLES */

/* The plain cases */
#if 0
/* 1 position */
static u64 positionValuesFromTablePlain1(u64 *table, short p1, int q)
{
    return table[p1];
}

/* 2 positions */
static u64 positionValuesFromTablePlain2(u64 *table, short p1, short p2, int q)
{
    return table[p1 + p2*q];
}

/* 3 positions */
static u64 positionValuesFromTablePlain3(u64 *table, short p1, short p2, short p3, int q)
{
    return table[p1 + p2*q + p3*q*q];
}

/* 4 positions */
static u64 positionValuesFromTablePlain4(u64 *table, short p1, short p2, short p3, short p4, int q)
{
    return table[p1 + p2*q + p3*q*q + p4*q*q*q];
}

/* 5 positions */
static u64 positionValuesFromTablePlain5(u64 *table, short p1, short p2, short p3, short p4, short p5, int q)
{
    return table[p1 + p2*q + p3*q*q + p4*q*q*q + p5*q*q*q*q];
}

/* 6 positions */
static u64 positionValuesFromTablePlain6(u64 *table, short p1, short p2, short p3, short p4, short p5, short p6, int q)
{
    return table[p1 + p2*q + p3*q*q + p4*q*q*q + p5*q*q*q*q + p6*q*q*q*q*q];
}
#endif
static u64 positionValuesFromTablePlain(u64 *table, int numPositions, short *t, int c)
{
    ASSERT(0 < numPositions && numPositions <= MAX_LMS_POSITIONS, "unexpected parameter numPositions");
    u64 index = t[numPositions - 1];
    for (int i=numPositions-2; i>=0; i--)
    {
        index = c * index + t[i];
    }
    return table[index];
}

/* makes an LMS mapping of pi using q and p as moduli */
static int positionLMSMap(short pi, int q, int p)
{
    pi += p/2;
    pi /= p;

    int c = q/p + 1;
    if (pi == c)
    {
        pi -= c;
    }

    return pi;
}
#if 0
/* 1 position */
static u64 positionValuesFromTableLMS1(u64 *table, short p1, int q, int p)
{
    short t1 = (p1, q, p);
    short c = q/p + 1;
    return positionValuesFromTablePlain1(table, t1, c);
}

/* 2 positions */
static u64 positionValuesFromTableLMS2(u64 *table, short p1, short p2, int q, int p)
{
    short t1 = positionLMSMap(p1, q, p);
    short t2 = positionLMSMap(p2, q, p);
    short c = q/p + 1;
    return positionValuesFromTablePlain2(table, t1, t2, c);
}

/* 3 positions */
static u64 positionValuesFromTableLMS3(u64 *table, short p1, short p2, short p3, int q, int p)
{
    short t1 = positionLMSMap(p1, q, p);
    short t2 = positionLMSMap(p2, q, p);
    short t3 = positionLMSMap(p3, q, p);
    short c = q/p + 1;
    return positionValuesFromTablePlain3(table, t1, t2, t3, c);
}

/* 4 positions */
static u64 positionValuesFromTableLMS4(u64 *table, short p1, short p2, short p3, short p4, int q, int p)
{
    short t1 = positionLMSMap(p1, q, p);
    short t2 = positionLMSMap(p2, q, p);
    short t3 = positionLMSMap(p3, q, p);
    short t4 = positionLMSMap(p4, q, p);
    short c = q/p + 1;
    return positionValuesFromTablePlain4(table, t1, t2, t3, t4, c);
}

/* 5 positions */
static u64 positionValuesFromTableLMS5(u64 *table, short p1, short p2, short p3, short p4, short p5, int q, int p)
{
    short t1 = positionLMSMap(p1, q, p);
    short t2 = positionLMSMap(p2, q, p);
    short t3 = positionLMSMap(p3, q, p);
    short t4 = positionLMSMap(p4, q, p);
    short t5 = positionLMSMap(p5, q, p);
    short c = q/p + 1;
    return positionValuesFromTablePlain5(table, t1, t2, t3, t4, t5, c);
}

/* 6 positions */
static u64 positionValuesFromTableLMS6(u64 *table, short p1, short p2, short p3, short p4, short p5, short p6, int q, int p)
{
    short t1 = positionLMSMap(p1, q, p);
    short t2 = positionLMSMap(p2, q, p);
    short t3 = positionLMSMap(p3, q, p);
    short t4 = positionLMSMap(p4, q, p);
    short t5 = positionLMSMap(p5, q, p);
    short t6 = positionLMSMap(p6, q, p);
    short c = q/p + 1;
    return positionValuesFromTablePlain6(table, t1, t2, t3, t4, t5, t6, c);
}
#endif
static u64 positionValuesFromTableLMS(u64 *table, int q, int p, int numPositions, short *pn)
{
    short t[MAX_LMS_POSITIONS];
    ASSERT(0 < numPositions && numPositions <= MAX_LMS_POSITIONS, "unexpected parameter numPositions");
    for (int i=0; i<numPositions; i++)
    {
        t[i] = positionLMSMap(pn[i], q, p);
    }
    int c = q/p + 1;
    return positionValuesFromTablePlain(table, numPositions, t, c);
}

u64 position_values_2_category_index_lms(lweInstance *lwe, bkwStepParameters *dstBkwStepPar, short *pn)
{
    int q = lwe->q;
    int p = dstBkwStepPar->sortingPar.LMS.p;
    int numPositions = dstBkwStepPar->numPositions;

    /* use lookup table if available */
    int c = q/p + 1;
    if (lookup_lms.table && lookup_lms.c == c && lookup_lms.numPositions == numPositions)   /* remove check increased speed (if check can be made superfluous) */
    {
        u64 cat = positionValuesFromTableLMS(lookup_lms.table, q, p, numPositions, pn);
        ASSERT(cat < num_categories(lwe, dstBkwStepPar), "category index calculated incorrectly");
        return cat;
    }
    /* if not, try to generate lookup table */
    generate_table_lms(q, p, numPositions);
    /* now try to use lookup table again */
    if (lookup_lms.table && lookup_lms.c == c && lookup_lms.numPositions == numPositions)
    {
        u64 cat = positionValuesFromTableLMS(lookup_lms.table, q, p, numPositions, pn);
        ASSERT(cat < num_categories(lwe, dstBkwStepPar), "category index calculated incorrectly");
        return cat;
    }
    /* if this did not work, throw failure */
    ASSERT_ALWAYS("index lms error");
    return -1;
}

/* Smooth LMS */

// n is the length, t is a list containing the indices, p is a list containing the reduction factor for each position.
// efficiency may be improved with a look-up table
static u64 positionValuesToCategoryGeneralized(int n, short *t, int *c)
{

    u64 index;

    if (c[n-1] & 1)   // c odd case
    {

        if (n == 1)   // single position case
        {
            if (t[0] == 0)
                return 0;
            else if (2*t[0] < c[0])
                return 2*t[0] -1;
            else
                return 2*(c[0]-t[0]);

        }
        else     // other recursive cases
        {
            if (t[n-1] == 0)
            {
                return positionValuesToCategoryGeneralized(n-1, t, c);
            }
            else if (2*t[n-1] < c[n-1])
            {
                index = (2*t[n-1]-1);
                for (int i = 0; i<n-1; i++)
                    index *= c[i];
                index += 2*positionValuesToCategoryGeneralized(n-1, t, c);
            }
            else
            {
                index = 2*(c[n-1]-t[n-1])-1;
                short nt[n-1];
                for (int i = 0; i<n-1; i++)
                {
                    if (c[i] & 1)   /* Category i is odd */
                    {
                        nt[i] = (c[i] - t[i])%c[i];
                    }
                    else     /* Category i is even */
                    {
                        nt[i] = (c[i] - t[i]-1);
                    }
                    index *= c[i];
                }
                index += 1 + 2*positionValuesToCategoryGeneralized(n-1, nt, c);
            }
        }

    }
    else     // c even  case
    {

        if (n == 1)    // single position case
        {
            if (2*t[0] < c[0])
            {
                return 2*t[0];
            }
            else
            {
                return 2*(c[0]-t[0])-1;
            }
        }
        else     // other recursive cases
        {
            if (t[n-1] < (c[n-1]/2))
            {
                index = 2*t[n-1];
                for (int i = 0; i<n-1; i++)
                    index *= c[i];
                index += 2*positionValuesToCategoryGeneralized(n-1, t, c);
            }
            else
            {
                index = 2*(c[n-1]-t[n-1]-1);
                short nt[n-1];
                for (int i = 0; i<n-1; i++)
                {
                    if (c[i] & 1)   /* Category i is odd */
                    {
                        nt[i] = (c[i] - t[i])%c[i];
                    }
                    else     /* Category i is even */
                    {
                        nt[i] = (c[i] - t[i]-1);
                    }
                    index *= c[i];
                }
                index += 1 + 2*positionValuesToCategoryGeneralized(n-1, nt, c);
            }
        }

    }
    return index;
}

/* makes an LMS mapping of pi using q and p as moduli - version for smooth LMS */
/* here q_ must be ceil(lwe.q/2) if the position is not already reduced, otherwise it must be 2*p from previous step */
short positionSmoothLMSMap(short pn, int q, int q_, int p, int c)
{

    int Delta;

    if (c & 1)   // if c is odd
    {
        Delta = p*((c/2) +1) -q_;
        if (pn < q_)
            return (pn + Delta) / p;
        else
            return (c -(((q - pn) + Delta) / p)) % c;
    }
    else     // if c is even
    {
        Delta = p*(c/2)-q_; // shift
        if (pn < q_)
            return (pn+Delta)/p;
        else
            return c -1 -(((q - pn) + Delta) / p);
    }
}

u64 position_values_2_category_index_smooth_lms(lweInstance *lwe, bkwStepParameters *dstBkwStepPar, short *pn)
{
    int q = lwe->q;
    int Ni = dstBkwStepPar->numPositions;
    short p = dstBkwStepPar->sortingPar.smoothLMS.p;
    short p1 = dstBkwStepPar->sortingPar.smoothLMS.p1;
    short p2;
    u64 index_cat = 0;

    short t[MAX_SMOOTH_LMS_POSITIONS+1];
    int c[MAX_SMOOTH_LMS_POSITIONS+1];
    int q_;

    /* Differentiate the first step for general steps */
    if (dstBkwStepPar->sortingPar.smoothLMS.prev_p1 == -1)   // first step
    {
        q_ = q%2 == 1 ? (q+1)/2 : q/2;
        for(int i = 0; i < Ni; i++)
        {
            c[i] = ((2*q_-1) % p) == 0 ? ((2*q_-1) / p) : ((2*q_-1) / p) + 1;
            t[i] = positionSmoothLMSMap(pn[i], q, q_, p, c[i]);
        }
        c[Ni] = ((2*q_-1) % p1) == 0 ? ((2*q_-1) / p1) : ((2*q_-1) / p1) + 1;
        t[Ni] = positionSmoothLMSMap(pn[Ni], q, q_, p1, c[Ni]);
        index_cat = positionValuesToCategoryGeneralized(Ni+1, t, c);
    }
    else if (dstBkwStepPar->startIndex + Ni == lwe->n)     // last step
    {
        p2 = dstBkwStepPar->sortingPar.smoothLMS.p2;
        q_= dstBkwStepPar->sortingPar.smoothLMS.prev_p1;
        c[0] = ((2*q_-1) % p2) == 0 ? ((2*q_-1) / p2) : ((2*q_-1) / p2) + 1;
        t[0] = positionSmoothLMSMap(pn[0], q, q_, p2, c[0]);
        q_ = q%2 == 1 ? (q+1)/2 : q/2;
        for(int i = 1; i < Ni; i++)
        {
            c[i] = ((2*q_-1) % p) == 0 ? ((2*q_-1) / p) : ((2*q_-1) / p) + 1;
            t[i] = positionSmoothLMSMap(pn[i], q, q_, p, c[i]);
        }
        index_cat = positionValuesToCategoryGeneralized(Ni, t, c);
    }
    else      // middle steps
    {
        p2 = dstBkwStepPar->sortingPar.smoothLMS.p2;
        q_= dstBkwStepPar->sortingPar.smoothLMS.prev_p1;
        c[0] = ((2*q_-1) % p2) == 0 ? ((2*q_-1) / p2) : ((2*q_-1) / p2) + 1;
        t[0] = positionSmoothLMSMap(pn[0], q, q_, p2, c[0]);
        q_ = q%2 == 1 ? (q+1)/2 : q/2;
        for(int i = 1; i < Ni; i++)
        {
            c[i] = ((2*q_-1) % p) == 0 ? ((2*q_-1) / p) : ((2*q_-1) / p) + 1;
            t[i] = positionSmoothLMSMap(pn[i], q, q_, p, c[i]);
        }
        c[Ni] = ((2*q_-1) % p1) == 0 ? ((2*q_-1) / p1) : ((2*q_-1) / p1) + 1;
        t[Ni] = positionSmoothLMSMap(pn[Ni], q, q_, p1, c[Ni]);
        index_cat = positionValuesToCategoryGeneralized(Ni+1, t, c);
    }
    return index_cat;
}

/* TODO Switch notation from Ni to ni */
/* TODO implement lookup table to speedup */
/* TODO simplify to not calculate unnecessary stuff not needed when skipping the last positions*/
/* Smooth LMS mapping where we skip some position(s) at the end */
u64 position_values_2_category_index_smooth_lms_meta(lweInstance *lwe, bkwStepParameters *dstBkwStepPar, short *pn)
{
    int q = lwe->q;
    int Ni = dstBkwStepPar->numPositions;
    short p = dstBkwStepPar->sortingPar.smoothLMS.p;
    short p1 = dstBkwStepPar->sortingPar.smoothLMS.p1;
    short p2;
    short meta_skipped = dstBkwStepPar->sortingPar.smoothLMS.meta_skipped; /* Number of positions to skip when dividing samples into metacategories */
    u64 index_cat = 0;

    short t[MAX_SMOOTH_LMS_POSITIONS+1];
    int c[MAX_SMOOTH_LMS_POSITIONS+1];
    int q_;

    /* Differentiate the first step for general steps */
    if (dstBkwStepPar->sortingPar.smoothLMS.prev_p1 == -1)   // first step
    {
        q_ = q%2 == 1 ? (q+1)/2 : q/2;
        for(int i = 0; i < Ni; i++)
        {
            c[i] = ((2*q_-1) % p) == 0 ? ((2*q_-1) / p) : ((2*q_-1) / p) + 1;
            t[i] = positionSmoothLMSMap(pn[i], q, q_, p, c[i]);
        }
        c[Ni] = ((2*q_-1) % p1) == 0 ? ((2*q_-1) / p1) : ((2*q_-1) / p1) + 1;
        t[Ni] = positionSmoothLMSMap(pn[Ni], q, q_, p1, c[Ni]);
        index_cat = positionValuesToCategoryGeneralized(Ni+1-meta_skipped, t, c);
    }
    else if (dstBkwStepPar->startIndex + Ni == lwe->n)     // last step
    {
        p2 = dstBkwStepPar->sortingPar.smoothLMS.p2;
        q_= dstBkwStepPar->sortingPar.smoothLMS.prev_p1;
        c[0] = ((2*q_-1) % p2) == 0 ? ((2*q_-1) / p2) : ((2*q_-1) / p2) + 1;
        t[0] = positionSmoothLMSMap(pn[0], q, q_, p2, c[0]);
        q_ = q%2 == 1 ? (q+1)/2 : q/2;
        for(int i = 1; i < Ni; i++)
        {
            c[i] = ((2*q_-1) % p) == 0 ? ((2*q_-1) / p) : ((2*q_-1) / p) + 1;
            t[i] = positionSmoothLMSMap(pn[i], q, q_, p, c[i]);
        }
        index_cat = positionValuesToCategoryGeneralized(Ni-meta_skipped, t, c);
    }
    else      // middle steps
    {
        p2 = dstBkwStepPar->sortingPar.smoothLMS.p2;
        q_= dstBkwStepPar->sortingPar.smoothLMS.prev_p1;
        c[0] = ((2*q_-1) % p2) == 0 ? ((2*q_-1) / p2) : ((2*q_-1) / p2) + 1;
        t[0] = positionSmoothLMSMap(pn[0], q, q_, p2, c[0]);
        q_ = q%2 == 1 ? (q+1)/2 : q/2;
        for(int i = 1; i < Ni; i++)
        {
            c[i] = ((2*q_-1) % p) == 0 ? ((2*q_-1) / p) : ((2*q_-1) / p) + 1;
            t[i] = positionSmoothLMSMap(pn[i], q, q_, p, c[i]);
        }
        c[Ni] = ((2*q_-1) % p1) == 0 ? ((2*q_-1) / p1) : ((2*q_-1) / p1) + 1;
        t[Ni] = positionSmoothLMSMap(pn[Ni], q, q_, p1, c[Ni]);
        index_cat = positionValuesToCategoryGeneralized(Ni+1-meta_skipped, t, c);
    }
    return index_cat;
}

/* Coded BKW */

static u64 positionValuesFromTableCodedBKW(int q, codingType ct, short *t)
{
    int c1, c2, c3, c4;
    switch(ct)
    {
    case blockCode_21:
        closest_code_word_2_1(&c1, &c2, q, t[0], t[1]);
        return c1;
    case blockCode_31:
        closest_code_word_3_1(&c1, &c2, &c3, q, t[0], t[1], t[2]);
        return c1;
    case blockCode_41:
        closest_code_word_4_1(&c1, &c2, &c3, &c4, q, t[0], t[1], t[2], t[3]);
        return c1;
    case concatenatedCode_21_21:
        closest_code_word_2_1(&c1, &c2, q, t[0], t[1]);
        closest_code_word_2_1(&c3, &c4, q, t[2], t[3]);
        return c1 + q*c3;
    default:
        ASSERT_ALWAYS("codedBKW parameters set not implemented!");
    }
    return 1;
}

u64 position_values_2_category_index_coded_bkw(lweInstance *lwe, bkwStepParameters *dstBkwStepPar, short *pn)
{
    int q = lwe->q;
    codingType ct = dstBkwStepPar->sortingPar.CodedBKW.ct;

    /* use syndrome lookup table if available */
    if (is_syndrome_decoding_table_loaded(q, ct))   /* remove check increased speed (if check can be made superfluous) */
    {
        u64 cat = positionValuesFromTableCodedBKW(q, ct, pn);
        ASSERT(cat < num_categories(lwe, dstBkwStepPar), "category index calculated incorrectly");
        return cat;
    }

    /* if not, try to load it into memory from existing file */
    load_syndrome_decoding_table(q, ct, 1 /* generate if non-existing */);
    if (is_syndrome_decoding_table_loaded(q, ct))
    {
        u64 cat = positionValuesFromTableCodedBKW(q, ct, pn);
        ASSERT(cat < num_categories(lwe, dstBkwStepPar), "category index calculated incorrectly");
        return cat;
    }

    /* if this did not work, throw failure */
    ASSERT_ALWAYS("index coded bkw error, failed to load/create syndrome decoding table");
    return -1;
}

u64 position_values_2_category_index(lweInstance *lwe, lweSample *sample, bkwStepParameters *bkwStepPar)
{
    switch (bkwStepPar->sorting)
    {
    case plainBKW:
        if (bkwStepPar->numPositions == 2)
        {
            return position_values_2_category_index_plain_bkw(lwe->q, sample->col.a + bkwStepPar->startIndex); /* calculate destination category */
        }
        if (bkwStepPar->numPositions == 3)
        {
            return position_values_2_category_index_plain_bkw(lwe->q, sample->col.a + bkwStepPar->startIndex); /* calculate destination category, same as for 2-position case as third position value is suppressed */
        }
        break;
    case LMS:
        ASSERT(0 < bkwStepPar->numPositions && bkwStepPar->numPositions <= MAX_LMS_POSITIONS, "unsupported parameter set (LMS)");
        return position_values_2_category_index_lms(lwe, bkwStepPar, sample->col.a + bkwStepPar->startIndex);
    case smoothLMS:
        ASSERT(0 < bkwStepPar->numPositions && bkwStepPar->numPositions <= MAX_SMOOTH_LMS_POSITIONS, "unsupported parameter set (smooth LMS)");
        if (bkwStepPar->sortingPar.smoothLMS.meta_skipped == 0)   /* Don't combine categories into metacategories */
        {
            return position_values_2_category_index_smooth_lms(lwe, bkwStepPar, sample->col.a + bkwStepPar->startIndex);
        }
        else   /* Combine categories into metacategories */
        {
            return position_values_2_category_index_smooth_lms_meta(lwe, bkwStepPar, sample->col.a + bkwStepPar->startIndex);
        }
    case codedBKW:
        ASSERT(0 < bkwStepPar->numPositions && bkwStepPar->numPositions <= MAX_CODED_BKW_POSITIONS, "unsupported parameter set (coded BKW)");
        return position_values_2_category_index_coded_bkw(lwe, bkwStepPar, sample->col.a + bkwStepPar->startIndex);
    default:
        ASSERT_ALWAYS("unsupported parameter set (default)");
        return -1;
    }
    ASSERT_ALWAYS("unsupported parameter set (unknown)");
    return -1;
}

u64 position_values_2_category_index_from_partial_sample(lweInstance *lwe, short *a, bkwStepParameters *bkwStepPar)
{
    switch (bkwStepPar->sorting)
    {
    case plainBKW:
        ASSERT(bkwStepPar->numPositions == 2 || bkwStepPar->numPositions == 3, "unsupported parameter set (plain BKW)");
        return position_values_2_category_index_plain_bkw(lwe->q, a); /* calculate destination category */
    case LMS:
        ASSERT(0 < bkwStepPar->numPositions && bkwStepPar->numPositions <= MAX_LMS_POSITIONS, "unsupported parameter set (LMS)");
        return position_values_2_category_index_lms(lwe, bkwStepPar, a);
    case smoothLMS:
        ASSERT(0 < bkwStepPar->numPositions && bkwStepPar->numPositions <= MAX_SMOOTH_LMS_POSITIONS, "unsupported parameter set (smooth LMS)");
        return position_values_2_category_index_smooth_lms(lwe, bkwStepPar, a);
    case codedBKW:
        ASSERT(0 < bkwStepPar->numPositions && bkwStepPar->numPositions <= MAX_CODED_BKW_POSITIONS, "unsupported parameter set (coded BKW)");
        return position_values_2_category_index_coded_bkw(lwe, bkwStepPar, a);
    default:
        ASSERT_ALWAYS("unsupported parameter set (default)");
        return -1;
    }
    ASSERT_ALWAYS("unsupported parameter set (unknown)");
    return -1;
}
