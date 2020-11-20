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

#include "transform_secret.h"
#include "memory_utils.h"

/* transform (known) secret vector s */
/* s -> b - A^{-1}s */
/* \vector{s} -> b - A^{-1}\vector{s} */
void transformSecret(lweInstance *lwe, short *s)
{
    ASSERT(lwe, "no lwe instance");
    int n = lwe->n;
    int q = lwe->q;
    short **A = lwe->A;
    ASSERT(A, "no transformation matrix A");
    short *b = lwe->b;
    ASSERT(b, "no b-vector");

    short *t = MALLOC(n * sizeof(short));
    MEMCPY(t, s, n * sizeof(short)); /* copy s to temp */

    /* compute transformed s-vector */
    for (int i=0; i<n; i++)
    {
        u64 sum = 0;
        for (int j=0; j<n; j++)
        {
            sum += A[i][j] * t[j];
        }
        sum = sum % q;
        s[i] = (q + sum - b[i]) % q;
        ASSERT(0 <= s[i], "unexpected value");
        ASSERT(s[i] < q, "unexpected value");
    }
    FREE(t);
}

/* inverse transform secret */
/* s -> b - A^{-1}s */
/* \vector{s} -> b - A^{-1}\vector{s} */
void inverseTransformSecret(lweInstance *lwe, short *s)
{
    ASSERT(lwe, "no lwe instance");
    int n = lwe->n;
    int q = lwe->q;
    short **Ainv = lwe->A_inverse;
    ASSERT(Ainv, "no transformation matrix Ainv");
    short *b = lwe->b;
    ASSERT(b, "no b-vector");

    short *t = MALLOC(n * sizeof(short));
    for (int i=0; i<n; i++)
    {
        t[i] = (s[i] + b[i]) % q;
        ASSERT(0 <= t[i], "unexpected value");
        ASSERT(t[i] < q, "unexpected value");
    }

    /* compute inverse transform of s-vector */
    for (int i=0; i<n; i++)
    {
        u64 sum = 0;
        for (int j=0; j<n; j++)
        {
            sum += Ainv[i][j] * t[j];
        }
        s[i] = sum % q;
        ASSERT(0 <= s[i], "unexpected value");
        ASSERT(s[i] < q, "unexpected value");
    }
    FREE(t);
}
