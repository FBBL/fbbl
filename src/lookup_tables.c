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

#include "lookup_tables.h"
#include "memory_utils.h"
#include "assert_utils.h"

int tableQ = 0;
int **sum_table = NULL;
int **diff_table = NULL;

int createSumAndDiffTables(int q)
{
    if ((tableQ == q) && sum_table && diff_table)
    {
        return 0; /* tables already created */
    }
    tableQ = q;
    sum_table = MALLOC(q * sizeof(int*));
    ASSERT(sum_table, "allocation failed");
    diff_table = MALLOC(q * sizeof(int*));
    ASSERT(diff_table, "allocation failed");
    for (int i=0; i<q; i++)
    {
        sum_table[i] = MALLOC(q * sizeof(int));
        ASSERT(sum_table[i], "allocation failed");
        diff_table[i] = MALLOC(q * sizeof(int));
        ASSERT(diff_table[i], "allocation failed");
        for (int j=0; j<q; j++)
        {
            sum_table[i][j] = (i + j) % q;
            ASSERT(sum_table[i][j] < q, "unexpected value!");
            diff_table[i][j] = (i + q - j) % q;
            ASSERT(diff_table[i][j] < q, "unexpected value!");
        }
    }
    return 0;
}

void freeSumAndDiffTables(void)
{
    if ((tableQ == 0) || !sum_table || !diff_table)
    {
        return; /* tables already deleted */
    }
    for (int i=0; i<tableQ; i++)
    {
        FREE(sum_table[i]);
        FREE(diff_table[i]);
    }
    FREE(sum_table);
    FREE(diff_table);
    sum_table = NULL;
    diff_table = NULL;
    tableQ = 0;
}

int sumTable(int a, int b)
{
    ASSERT(sum_table, "sum_table not initialized");
    return sum_table[a][b];
}

int diffTable(int a, int b)
{
    ASSERT(diff_table, "diff_table not initialized");
    return diff_table[a][b];
}
