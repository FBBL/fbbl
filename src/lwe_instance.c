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

#include <stdio.h>
#include <math.h>
#include "lwe_instance.h"
#include "memory_utils.h"
#include "assert_utils.h"
#include "random_utils.h"
#include "platform_types.h"

#define PI 3.14159265358979

static int roundInt(double d)
{
    return d > 0.0 ? d + 0.5 : d - 0.5;
}

/* adapted malb */
int chi(double sigma, rand_ctx *rnd)
{
#if 1
    const double aa = randomUtilDouble(rnd);
    const double a = aa == 0 ? 0 : sqrt(-2 * log(aa));
    const double b = 2 * PI * randomUtilDouble(rnd);
    const double x = sigma * a * cos(b);
    //const double y = sigma * a * sin(b);
    int ret = roundInt(x);
    return ret;
#else
    const double aa = rand() / (double)RAND_MAX;
    const double a = aa == 0 ? 0 : sqrt(-2 * log(aa));
    const double b = 2 * PI * rand() / (double)RAND_MAX;
    const double x = sigma * a * cos(b);
    //const double y = sigma * a * sin(b);
    int ret = roundInt(x);
    return ret;
#endif
}

/* compute a simple hash value for a column
 * zero column corresponds to hash value zero */
u64 bkwColumnComputeHash(lweSample *sample, int n, int startRow)
{
    u64 h = 0;
    for (int i=startRow; i<n; i++)
    {
        h = (h << 13) + columnValue(sample,i); /* shift up some and add */
        h = (h ^ (h >> 48)) & ((u64)0x00FFFFFFFFFFFFFF); /* add 16 hi bits to lo end */
    }
    return h;
}

void printColumn(bkwColumn *col, int n)
{
    printf("(%d", col->a[0]);
    for (int i=1; i<n; i++)
    {
        printf(",%d", col->a[i]);
    }
    printf(")");
}

void printSample(lweSample *sample, int n)
{
    printf("[");
    printColumn(&sample->col, n);
    printf(", %d]", sample->sumWithError);
}

int columnIsZero(lweSample *sample, int n)
{
    if (columnHash(sample) != 0)
    {
        return 0;
    }
    for (int i=0; i<n; i++)
    {
        if (columnValue(sample, i) != 0)
        {
            return 0;
        }
    }
    return 1;
}

static lweSample *newEmptySample()
{
    return MALLOC(LWE_SAMPLE_SIZE_IN_BYTES);
}

static lweSample *newRandomSample(int n, int q, double sigma, rand_ctx *rnd, short *s)
{
    int sum = 0, err;
    lweSample *sample = newEmptySample();
    ASSERT(sample, "Allocation failed!\n");
    ASSERT(s, "Unexpected parameter s!\n");
    for (int i=0; i<n; i++)
    {
        sample->col.a[i] = randomUtilInt(rnd, q);
        sum = (sum + sample->col.a[i]*s[i]) % q;
    }
    err = (chi(sigma, rnd) + q) % q;
    sample->error = err; // store error only
    sum = (sum + err + q) % q;
    sample->sumWithError = sum; // store sum a_i*s_i + error
    sample->col.hash = bkwColumnComputeHash(sample, n, 0); // compute hash
    return sample;
}

static void newInPlaceRandomSample(lweSample *sample, int n, int q, double sigma, rand_ctx *rnd, short *s)
{
    int sum = 0, err;
    ASSERT(sample, "Allocation failed!\n");
    ASSERT(s, "Unexpected parameter s!\n");
    for (int i=0; i<n; i++)
    {
        sample->col.a[i] = randomUtilInt(rnd, q);
        sum = (sum + sample->col.a[i]*s[i]) % q;
    }
    err = (chi(sigma, rnd) + q) % q;
    sample->error = err; // store error only
    sum = (sum + err + q) % q;
    sample->sumWithError = sum; // store sum a_i*s_i + error
    sample->col.hash = bkwColumnComputeHash(sample, n, 0); // compute hash
}

static void freeSample(lweSample *sample)
{
    FREE(sample);
}

void lweInit(lweInstance *lwe, int n, int q, double alpha)
{
    int i;
    ASSERT(n <= MAX_N, "Invalid n!\n");
    if (n > MAX_N)
    {
        printf("Invalid parameter n (larger than MAX_N)!\n");
        exit(1);
    }
    if (n < MAX_N)
    {
        printf("Warning: Wasting memory with %d = n < MAX_N = %d!\n", n, MAX_N);
    }
    lwe->n = n;
    lwe->q = q;
    lwe->alpha = alpha;
    lwe->sigma = alpha * q;
    lwe->A = NULL;
    lwe->A_inverse = NULL;
    lwe->b = NULL;
    randomUtilInit(&lwe->rnd);
    for (i=0; i<n; i++)
    {
//      lwe->s[i] = randomUtilInt(&lwe->rnd, q); /* original lwe */
        lwe->s[i] = (chi(lwe->sigma, NULL) + q) % q; /* with implicit initial transformation */
    }
    for (i=n; i<MAX_N; i++)
    {
        lwe->s[i] = 0;
    }
    lwe->newEmptySample = newEmptySample;
    lwe->newRandomSample = newRandomSample;
    lwe->newInPlaceRandomSample = newInPlaceRandomSample;
    lwe->freeSample = freeSample;
}

void lweInstanceAllocateLinearTransformationMatrices(lweInstance *lwe)
{
    if (!lwe->A)   /* if not already allocated, allocate A */
    {
        lwe->A = (short**)MALLOC(lwe->n * sizeof(short*));
        for (int i=0; i<lwe->n; i++)
        {
            lwe->A[i] = (short*)MALLOC(lwe->n * sizeof(short));
        }
    }
    if (!lwe->A_inverse)   /* if not already allocated, allocate A_inverse */
    {
        lwe->A_inverse = (short**)MALLOC(lwe->n * sizeof(short*));
        for (int i=0; i<lwe->n; i++)
        {
            lwe->A_inverse[i] = (short*)MALLOC(lwe->n * sizeof(short));
        }
    }
    if (!lwe->b)   /* if not already allocated, allocate b */
    {
        lwe->b = (short*)MALLOC(lwe->n * sizeof(short));
    }
}

static void lweInstanceFreeLinearTransformationMatrices(lweInstance *lwe)
{
    if (lwe->A)
    {
        for (int i=0; i<lwe->n; i++)
        {
            FREE(lwe->A[i]);
        }
        FREE(lwe->A);
        lwe->A = NULL;
    }
    if (lwe->A_inverse)
    {
        for (int i=0; i<lwe->n; i++)
        {
            FREE(lwe->A_inverse[i]);
        }
        FREE(lwe->A_inverse);
        lwe->A_inverse = NULL;
    }
    if (lwe->b)
    {
        FREE(lwe->b);
        lwe->b = NULL;
    }
}

void lweDestroy(lweInstance *lwe)
{
    lweInstanceFreeLinearTransformationMatrices(lwe);
}
