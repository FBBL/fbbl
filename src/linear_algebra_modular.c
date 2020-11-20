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

#include "linear_algebra_modular.h"
#include "memory_utils.h"
#include <stdio.h>

/* multiplicative inverse of a mod b, from rosetta code */
static int mul_inv(int a, int b)
{
    int t, nt, r, nr, q, tmp;
    if (b < 0) b = -b;
    if (a < 0) a = b - (-a % b);
    t = 0;
    nt = 1;
    r = b;
    nr = a % b;
    while (nr != 0)
    {
        q = r/nr;
        tmp = nt;
        nt = t - q*nt;
        t = tmp;
        tmp = nr;
        nr = r - q*nr;
        r = tmp;
    }
    if (r > 1) return -1;  /* No inverse */
    if (t < 0) t += b;
    return t;
}

/*
static void print_matrix(int n, short **matrix) {
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      printf(" %6d", matrix[i][j]);
    }
    printf("\n");
  }
}
*/
/* input: triangulated T */
static void back_substitute(int n, int q, short **T, short **A_inverse)
{
    int qq = (q-1)*q;
    /*
      printf("T before back sub:\n");
      print_matrix(n, T);
      printf("\n\n");
    */
    for (int i=n-2; i>=0; i--)   /* for each row (from bottom) */
    {
        for (int j=i+1; j<n; j++)
        {
            int lc = T[i][j];
            for (int k=j; k<n; k++)   /* row i = row i - lc * row j (mod q) */
            {
                T[i][k] = (T[i][k] - lc * T[j][k] + qq) % q;
                ASSERT(T[i][k] >= 0, "Unexpected value");
                ASSERT(T[i][k] < q, "Unexpected value");
            }
            for (int k=0; k<n; k++)   /* row i = row i - lc * row j (mod q) */
            {
                A_inverse[i][k] = (A_inverse[i][k] - lc * A_inverse[j][k] + qq) % q;
                ASSERT(A_inverse[i][k] >= 0, "Unexpected value");
                ASSERT(A_inverse[i][k] < q, "Unexpected value");
            }
        }
        /*
            printf("T with %d back subs remaining:\n", i);
            print_matrix(n, T);
            printf("\n\n");
            printf("A_inverse with %d back subs remaining:\n", i);
            print_matrix(n, A_inverse);
            printf("\n\n");
        */
    }
    /*
      printf("T when finished:\n");
      print_matrix(n, T);
      printf("\n\n");
      printf("A_inverse when finished:\n");
      print_matrix(n, A_inverse);
      printf("\n\n");
    */
}

int compute_matrix_inverse_modular(lweInstance *lwe, lweSample *samples, int num_samples, short **A, short **A_inverse, short *b)
{
    int q = lwe->q;
    int n = lwe->n;
    if (num_samples < n)
    {
        return 0;
    }

    /* allocate temp T */
    short **T = MALLOC(n * sizeof(short*));
    for (int i=0; i<n; i++)
    {
        T[i] = MALLOC(n * sizeof(short));
    }

    int dim = 0; /* number of dimensions (rows) successfully added to A */
    lweSample *sample = samples;
    for (int i=0; i<num_samples; i++, sample++)
    {
        /* copy new row to A */
        for (int j=0; j<n; j++)
        {
            A[dim][j] = sample->col.a[j];
            T[dim][j] = sample->col.a[j];
            A_inverse[dim][j] = j == dim ? 1 : 0;
        }
        /* store b-value (used in reverse transform for recovering the secret vector s) */
        b[dim] = sample->sumWithError;
        /*
            printf("T initial at dim=%d:\n", dim);
            print_matrix(n, T);
            printf("\n\n");
            printf("A_inverse initial at dim=%d:\n", dim);
            print_matrix(n, A_inverse);
            printf("\n\n");
        */
        /* triangulate new row */
        int qq = (q-1)*q;
        for (int j=0; j<dim; j++)   /* for each previous row (already triangulated) */
        {
            int lc = T[dim][j];
            for (int k=0; k<n; k++)   /* row dim = row dim - lc * row j (mod q) */
            {
                T[dim][k] = (T[dim][k] - lc * T[j][k] + qq) % q;
                ASSERT(T[dim][k] >= 0, "Unexpected value");
                ASSERT(T[dim][k] < q, "Unexpected value");
                A_inverse[dim][k] = (A_inverse[dim][k] - lc * A_inverse[j][k] + qq) % q;
                ASSERT(A_inverse[dim][k] >= 0, "Unexpected value");
                ASSERT(A_inverse[dim][k] < q, "Unexpected value");
            }
            ASSERT(T[dim][j] == 0, "Unexpected value");
            /*
                  printf("T at dim=%d:\n", dim);
                  print_matrix(n, T);
                  printf("\n\n");
                  printf("A_inverse at dim=%d:\n", dim);
                  print_matrix(n, A_inverse);
                  printf("\n\n");
            */
        }
        /*
            printf("T triangulated at dim=%d:\n", dim);
            print_matrix(n, T);
            printf("\n\n");
            printf("A_inverse at dim=%d:\n", dim);
            print_matrix(n, A_inverse);
            printf("\n\n");
        */
        /* compute pivot inverse */
        int inv = mul_inv(T[dim][dim], q);
        if (inv == -1)   /* inverse does not exist */
        {
            continue; /* try next sample */
        }
        ASSERT(inv >= 0, "Unexpected value");
        ASSERT(inv < q, "Unexpected value");
        for (int j=0; j<n; j++)
        {
            T[dim][j] = (T[dim][j] * inv) % q;
            ASSERT(T[dim][j] >= 0, "Unexpected value");
            ASSERT(T[dim][j] < q, "Unexpected value");
            A_inverse[dim][j] = (A_inverse[dim][j] * inv) % q;
            ASSERT(A_inverse[dim][j] >= 0, "Unexpected value");
            ASSERT(A_inverse[dim][j] < q, "Unexpected value");
        }
        /*
            printf("T completed at dim=%d:\n", dim);
            print_matrix(n, T);
            printf("\n\n");
            printf("A_inverse completed at dim=%d:\n", dim);
            print_matrix(n, A_inverse);
            printf("\n\n");
        */
        dim++;

        if (dim == n)   /* sufficently many linearly independent samples have been found, and T is triangulated */
        {
            back_substitute(n, q, T, A_inverse);
            break;
        }
    }

    /* cleanup */
    for (int i=0; i<n; i++)
    {
        FREE(T[i]);
    }
    FREE(T);

    return dim == n ? sample - samples + 1: 0;
}
