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

#include "solve_fft.h"
#include "lwe_instance.h"
#include "storage_reader.h"
#include "storage_file_utilities.h"
#include "memory_utils.h"
#include "string_utils.h"
#include <math.h> /* for fftw3 */
#include <complex.h> /* for fftw3 */
#include "fftw3.h" /* for FFT */
#include <inttypes.h>

#define MIN(x,y) ((x)<(y)?(x):(y))

#define TAU 6.283185307179586476925286766559005768394338798750211641949

/* double-precision case */
static void processOneCategory(lweInstance *lwe, lweSample *buf, u64 numSamples, short *solution, int numSolvedCoordinates, fftw_complex *in, int fftPositions)
{
    int q = lwe->q;
    int n = lwe->n;
    ASSERT(0 <= numSolvedCoordinates && numSolvedCoordinates <= n, "Bad parameter numSolvedCoordinates!\n");
    int startIndex = n - numSolvedCoordinates - fftPositions;

    /* process all samples */
    for (u64 i=0; i<numSamples; i++)
    {
        lweSample *sample = &buf[i];

        /* find function value (index) to update */
        int p1, p2, p3;
        u64 currentIndex;
        switch (fftPositions)
        {
        case 1:
            p1 = columnValue(sample, startIndex);
            ASSERT(0 <= p1 && p1 < q, "unexpected p1");
            currentIndex = p1;
            ASSERT(0 <= currentIndex && currentIndex < q, "unexpected value");
            break;
        case 2:
            p1 = columnValue(sample, startIndex);
            ASSERT(0 <= p1 && p1 < q, "unexpected p1");
            p2 = columnValue(sample, startIndex + 1);
            ASSERT(0 <= p2 && p2 < q, "unexpected p2");
            currentIndex = q * p1 + p2;
            ASSERT(0 <= currentIndex && currentIndex < q * q, "unexpected value");
            break;
        case 3:
            p1 = columnValue(sample, startIndex);
            ASSERT(0 <= p1 && p1 < q, "unexpected p1");
            p2 = columnValue(sample, startIndex + 1);
            ASSERT(0 <= p2 && p2 < q, "unexpected p2");
            p3 = columnValue(sample, startIndex + 2);
            ASSERT(0 <= p3 && p3 < q, "unexpected p3");
            currentIndex = q * (q * p1 + p2) + p3;
            ASSERT(0 <= currentIndex && currentIndex < q * q * q, "unexpected value");
            break;
        default:
            ASSERT_ALWAYS("More than 3 positions not supported for fft solving");
            exit(1);
        }

        /* compute partial sum for solved part */
        u64 solvedSum = 0;
        for (int j=n-numSolvedCoordinates; j<n; j++)
        {
            solvedSum += columnValue(sample, j) * solution[j];
        }
        solvedSum = solvedSum % q;

        u64 se = sumWithError(sample);
        se = (se - solvedSum) % q;

#if 0
        /* update function value (index) to update */
        int p1 = columnValue(sample, startIndex);
        u64 currentIndex = p1;
        ASSERT(0 <= p1 && p1 < q, "unexpected p1");
        if (fftPositions >= 2)
        {
            int p2 = columnValue(sample, startIndex + 1);
            ASSERT(0 <= p2 && p2 < q, "unexpected p2");
            currentIndex = q * currentIndex + p2;
        }
        if (fftPositions >= 3)
        {
            int p3 = columnValue(sample, startIndex + 2);
            ASSERT(0 <= p3 && p3 < q, "unexpected p3");
            currentIndex = q * currentIndex + p3;
        }
#endif

        /* update fft input function */
        fftw_complex toAdd = cexp(TAU*I*se/q);
        in[currentIndex] += toAdd;
    }
}

/* single-precision case */
static void processOneCategorySingle(lweInstance *lwe, lweSample *buf, u64 numSamples, short *solution, int numSolvedCoordinates, fftwf_complex *in, int fftPositions)
{
    int q = lwe->q;
    int n = lwe->n;
    ASSERT(0 <= numSolvedCoordinates && numSolvedCoordinates <= n, "Bad parameter numSolvedCoordinates!\n");
    int startIndex = n - numSolvedCoordinates - fftPositions;

    /* process all samples */
    for (u64 i=0; i<numSamples; i++)
    {
        lweSample *sample = &buf[i];

        /* find function value (index) to update */
        int p1, p2, p3;
        u64 currentIndex;
        switch (fftPositions)
        {
        case 1:
            p1 = columnValue(sample, startIndex);
            ASSERT(0 <= p1 && p1 < q, "unexpected p1");
            currentIndex = p1;
            ASSERT(0 <= currentIndex && currentIndex < q, "unexpected value");
            break;
        case 2:
            p1 = columnValue(sample, startIndex);
            ASSERT(0 <= p1 && p1 < q, "unexpected p1");
            p2 = columnValue(sample, startIndex + 1);
            ASSERT(0 <= p2 && p2 < q, "unexpected p2");
            currentIndex = q * p1 + p2;
            ASSERT(0 <= currentIndex && currentIndex < q * q, "unexpected value");
            break;
        case 3:
            p1 = columnValue(sample, startIndex);
            ASSERT(0 <= p1 && p1 < q, "unexpected p1");
            p2 = columnValue(sample, startIndex + 1);
            ASSERT(0 <= p2 && p2 < q, "unexpected p2");
            p3 = columnValue(sample, startIndex + 2);
            ASSERT(0 <= p3 && p3 < q, "unexpected p3");
            currentIndex = q * (q * p1 + p2) + p3;
            ASSERT(0 <= currentIndex && currentIndex < q * q * q, "unexpected value");
            break;
        default:
            ASSERT_ALWAYS("More than 3 positions not supported for fft solving");
            exit(1);
        }

        /* compute partial sum for solved part */
        u64 solvedSum = 0;
        for (int j=n-numSolvedCoordinates; j<n; j++)
        {
            solvedSum += columnValue(sample, j) * solution[j];
        }
        solvedSum = solvedSum % q;

        u64 se = sumWithError(sample);
        se = (se - solvedSum) % q;

        fftwf_complex toAdd = cexp(TAU*I*se/q);
        in[currentIndex] += toAdd;
    }
}

/* calculates the fft of f in 1 dimension - double-precision*/
void calculateFft1d(fftw_complex *in, fftw_complex *out, int q)
{
    fftw_plan p = fftw_plan_dft_1d(q, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
}

/* calculates the fft of f in 1 dimension - single-precision*/
void calculateFft1dSingle(fftwf_complex *in, fftwf_complex *out, int q)
{
    fftwf_plan p = fftwf_plan_dft_1d(q, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(p);
    fftwf_destroy_plan(p);
}

/* calculates the fft of f in 2 dimensions - double-precision*/
void calculateFft2d(fftw_complex *in, fftw_complex *out, int q)
{
    fftw_plan p = fftw_plan_dft_2d(q, q, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
}

/* calculates the fft of f in 2 dimensions - single-precision*/
void calculateFft2dSingle(fftwf_complex *in, fftwf_complex *out, int q)
{
    fftwf_plan p = fftwf_plan_dft_2d(q, q, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(p);
    fftwf_destroy_plan(p);
}

/* calculates the fft of f in 3 dimensions - double-precision */
void calculateFft3d(fftw_complex *in, fftw_complex *out, int q)
{
    fftw_plan p = fftw_plan_dft_3d(q, q, q, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
}

/* calculates the fft of f in 3 dimensions - single-precision */
void calculateFft3dSingle(fftwf_complex *in, fftwf_complex *out, int q)
{
    fftwf_plan p = fftwf_plan_dft_3d(q, q, q, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(p);
    fftwf_destroy_plan(p);
}

/* find the solution by looking at the output of the fft - 1 position and double-precision case */
void updateSolution1d(fftw_complex *out, short *solution, int startIndex, int q)
{
    double current_max = -42; // arbitrary low value
    double current_value;
    int p1 = 0;

    for (int i = 0; i < q; i++)
    {
        current_value = creal(out[i]);
        if (current_value > current_max)
        {
            current_max = current_value;
            p1 = i;
        }
    }

//  printf("solution[%d] = %d\n", startIndex, p1);
    solution[startIndex] = p1;
}

/* Find the solution by looking at the output of the fft - 1 position single-precision case */
void updateSolution1dSingle(fftwf_complex *out, short *solution, int startIndex, int q)
{
    float current_max = creal(out[0]);
    int p1 = 0;

    for (int i = 1; i < q; i++)
    {
        float current_value = creal(out[i]);
        if (current_value > current_max)
        {
            current_max = current_value;
            p1 = i;
        }
    }

//  printf("solution[%d] = %d\n", startIndex, p1);
    solution[startIndex] = p1;
}

/* find the solution by looking at the output of the fft - 2 positions and double-precision case */
void updateSolution2d(fftw_complex *out, short *solution, int startIndex, int q)
{
    double current_max = -42; // arbitrary low value
    double current_value;
    int p1 = 0, p2 = 0;

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            current_value = creal(out[j + q*i]);
            if (current_value > current_max)
            {
                current_max = current_value;
                p1 = i;
                p2 = j;
            }
        }
    }

    solution[startIndex] = p1;
    solution[startIndex + 1] = p2;
}

/* Find the solution by looking at the output of the fft - 2 positions single-precision case */
void updateSolution2dSingle(fftwf_complex *out, short *solution, int startIndex, int q)
{
    float current_max = -42; // arbitrary low value
    float current_value;
    int p1 = 0, p2 = 0;

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            current_value = creal(out[j + q*i]);
            if (current_value > current_max)
            {
                current_max = current_value;
                p1 = i;
                p2 = j;
            }
        }
    }

    solution[startIndex] = p1;
    solution[startIndex + 1] = p2;
}

/* Find the solution by looking at the output of the fft - 3 positions and double-precision case */
void updateSolution3d(fftw_complex *out, short *solution, int startIndex, int q)
{
    double current_max = -42; // arbitrary low value
    double current_value;
    int p1 = 0, p2 = 0, p3 = 0;

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            for (int k = 0; k < q; k++)
            {
                current_value = creal(out[k + q*(j + q*i)]);
                if (current_value > current_max)
                {
                    current_max = current_value;
                    p1 = i;
                    p2 = j;
                    p3 = k;
                }
            }
        }
    }

    solution[startIndex] = p1;
    solution[startIndex + 1] = p2;
    solution[startIndex + 2] = p3;
}

/* Find the solution by looking at the output of the fft - 3 positions and single-precision case */
void updateSolution3dSingle(fftwf_complex *out, short *solution, int startIndex, int q)
{
    double current_max = -42; // arbitrary low value
    double current_value;
    int p1 = 0, p2 = 0, p3 = 0;

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            for (int k = 0; k < q; k++)
            {
                current_value = creal(out[k + q*(j + q*i)]);
                if (current_value > current_max)
                {
                    current_max = current_value;
                    p1 = i;
                    p2 = j;
                    p3 = k;
                }
            }
        }
    }

    solution[startIndex] = p1;
    solution[startIndex + 1] = p2;
    solution[startIndex + 2] = p3;
}

int solve_fft_search(const char *srcFolder, short *solution, int numSolvedCoordinates, int fftPositions, int precision)
{
    ASSERT(1 <= fftPositions && fftPositions <= 3, "The number of positions for fft is not supported!\n");

    lweInstance lwe;
    lweParametersFromFile(&lwe, srcFolder);
    int n = lwe.n;
    int q = lwe.q;

    u64 numCategories;
    sampleInfoFromFile(srcFolder, NULL, &numCategories, NULL, NULL, NULL);

    if (numSolvedCoordinates == n)
    {
        lweDestroy(&lwe);
        return 1; /* all positions solved */
    }

    storageReader sr;
    int ret = storageReaderInitialize(&sr, srcFolder);
    if (ret)
    {
        printf("*** solve_fft.c: storage reader returned %d on initialize\n", ret);
        lweDestroy(&lwe);
        return 2;
    }

    u64 numFftSlots;
    switch (fftPositions)
    {
    case 1:
        numFftSlots = q;
        break;
    case 2:
        numFftSlots = q * q;
        break;
    case 3:
        numFftSlots = q * q * q;
        break;
    default:
        ASSERT_ALWAYS("unsupported number of positions");
        exit(1);
    }

    /* create function to be fft input and output */
    fftw_complex *in;
    fftw_complex *out;
    fftwf_complex *inS;
    fftwf_complex *outS;
    fftw_plan p;
    fftwf_plan pS;

    switch (precision)
    {
    case FFT_SOLVER_SINGLE_PRECISION:
        inS = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * numFftSlots);
        outS = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * numFftSlots);
        ASSERT(inS, "allocation failed");
        ASSERT(outS, "allocation failed");
        if (!inS || !outS)
        {
            char s[128];
            printf("*** solve_plain_bkw_sorted: allocation failure (tried to allocate %s bytes)\n", sprintf_u64_delim(s, sizeof(fftwf_complex) * numFftSlots));
            ASSERT_ALWAYS("allocation failure");
            lweDestroy(&lwe);
            return 4;
        }
        MEMSET(inS, 0, sizeof(fftwf_complex) * numFftSlots);
        MEMSET(outS, 0, sizeof(fftwf_complex) * numFftSlots);
        break;
    case FFT_SOLVER_DOUBLE_PRECISION:
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFftSlots);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFftSlots);
        ASSERT(in, "allocation failed");
        ASSERT(out, "allocation failed");
        if (!in || !out)
        {
            char s[128];
            printf("*** solve_plain_bkw_sorted: allocation failure (tried to allocate %s bytes)\n", sprintf_u64_delim(s, sizeof(fftw_complex) * numFftSlots));
            ASSERT_ALWAYS("allocation failure");
            lweDestroy(&lwe);
            return 3;
        }
        MEMSET(in, 0, sizeof(fftw_complex) * numFftSlots);
        MEMSET(out, 0, sizeof(fftw_complex) * numFftSlots);
        break;
    default:
        ASSERT_ALWAYS("Unhandled precision");
        exit(1);
    }

    switch (fftPositions)
    {
    case 1:
        switch (precision)
        {
        case FFT_SOLVER_SINGLE_PRECISION:
            pS = fftwf_plan_dft_1d(q, inS, outS, FFTW_FORWARD, FFTW_ESTIMATE);
            MEMSET(inS, 0, sizeof(fftwf_complex) * numFftSlots);
            MEMSET(outS, 0, sizeof(fftwf_complex) * numFftSlots);
            break;
        case FFT_SOLVER_DOUBLE_PRECISION:
            p = fftw_plan_dft_1d(q, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            MEMSET(in, 0, sizeof(fftw_complex) * numFftSlots);
            MEMSET(out, 0, sizeof(fftw_complex) * numFftSlots);
            break;
        default:
            ASSERT_ALWAYS("Unhandled precision");
            exit(1);
        }
        break;
    case 2:
        break;
    case 3:
        break;
    default:
        ASSERT_ALWAYS("unsupported number of positions");
        exit(1);
    }

    /* collect statistics from samples */
    u64 categoryIndexCounter = 0;
    lweSample *buf1;
    lweSample *buf2;
    u64 numSamplesInBuf1, numSamplesInBuf2;
    int startIndex = n - fftPositions - numSolvedCoordinates;
    ret = storageReaderGetNextAdjacentCategoryPair(&sr, &buf1, &numSamplesInBuf1, &buf2, &numSamplesInBuf2);

    while (ret)
    {

        /* process buf1 */
        if (buf1)
        {
            switch (precision)
            {
            case FFT_SOLVER_SINGLE_PRECISION:
                categoryIndexCounter++;
                processOneCategorySingle(&lwe, buf1, numSamplesInBuf1, solution, numSolvedCoordinates, inS, fftPositions);
                break;
            case FFT_SOLVER_DOUBLE_PRECISION:
                categoryIndexCounter++;
                processOneCategory(&lwe, buf1, numSamplesInBuf1, solution, numSolvedCoordinates, in, fftPositions);
                break;
            default:
                ASSERT_ALWAYS("Unhandled precision");
                exit(1);
            }
        }

        /* process buf2 */
        if (buf2)
        {
            switch (precision)
            {
            case FFT_SOLVER_SINGLE_PRECISION:
                categoryIndexCounter++;
                processOneCategorySingle(&lwe, buf2, numSamplesInBuf2, solution, numSolvedCoordinates, inS, fftPositions);
                break;
            case FFT_SOLVER_DOUBLE_PRECISION:
                categoryIndexCounter++;
                processOneCategory(&lwe, buf2, numSamplesInBuf2, solution, numSolvedCoordinates, in, fftPositions);
                break;
            default:
                ASSERT_ALWAYS("Unhandled precision");
                exit(1);
            }
        }

        ret = storageReaderGetNextAdjacentCategoryPair(&sr, &buf1, &numSamplesInBuf1, &buf2, &numSamplesInBuf2);
    }
    if (categoryIndexCounter != numCategories)
    {
        printf("*** solve_fft_search: %" PRIu64 " categories processed (%" PRIu64 " expected)\n", categoryIndexCounter, numCategories);
    }

    switch (fftPositions)
    {
    case 1:
        switch (precision)
        {
        case FFT_SOLVER_SINGLE_PRECISION:
//      calculateFft1dSingle(inS, outS, q); /* Calculate fft */
//      printf("in[1577] = %10.2f + %10.2fI (before execution)\n", creal(inS[1577]), cimag(inS[1577]));
            fftwf_execute(pS);
//      printf("in[1577] = %10.2f + %10.2fI (after execution)\n", creal(inS[1577]), cimag(inS[1577]));
            updateSolution1dSingle(outS, solution, startIndex, q); /* deduce position values from fft */
//      printf("in[1577] = %10.2f + %10.2fI (after update)\n", creal(inS[1577]), cimag(inS[1577]));
            fftwf_destroy_plan(pS);
            break;
        case FFT_SOLVER_DOUBLE_PRECISION:
//      calculateFft1d(in, out, q); /* Calculate fft */
            fftw_execute(p);
            updateSolution1d(out, solution, startIndex, q); /* deduce position values from fft */
            fftw_destroy_plan(p);
            break;
        default:
            ASSERT_ALWAYS("Unhandled precision");
        }
        break;
    case 2:
        switch (precision)
        {
        case FFT_SOLVER_SINGLE_PRECISION:
            calculateFft2dSingle(inS, outS, q); /* Calculate fft */
            updateSolution2dSingle(outS, solution, startIndex, q); /* deduce position values from fft */
            break;
        case FFT_SOLVER_DOUBLE_PRECISION:
            calculateFft2d(in, out, q); /* Calculate fft */
            updateSolution2d(out, solution, startIndex, q); /* deduce position values from fft */
            break;
        default:
            ASSERT_ALWAYS("Unhandled precision");
        }
        break;
    case 3:
        switch (precision)
        {
        case FFT_SOLVER_SINGLE_PRECISION:
            calculateFft3dSingle(inS, outS, q); /* Calculate fft */
            updateSolution3dSingle(outS, solution, startIndex, q); /* deduce position values from fft */
            break;
        case FFT_SOLVER_DOUBLE_PRECISION:
            calculateFft3d(in, out, q); /* Calculate fft */
            updateSolution3d(out, solution, startIndex, q); /* deduce position values from fft */
            break;
        default:
            ASSERT_ALWAYS("Unhandled precision");
        }
        break;
    default:
        ASSERT_ALWAYS("unsupported number of positions");
    }

    /* cleanup */
    storageReaderFree(&sr);
    switch (precision)
    {
    case FFT_SOLVER_SINGLE_PRECISION:
        fftw_free(inS);
        fftw_free(outS);
        break;
    case FFT_SOLVER_DOUBLE_PRECISION:
        fftw_free(in);
        fftw_free(out);
        break;
    default:
        ASSERT_ALWAYS("Unhandled precision");
    }

    lweDestroy(&lwe);
    return 0;
}

/* HYBRID CONTENT STARTS HERE */

/* Hybrid category processing - double-precision case */
static void processOneCategoryHybrid(lweInstance *lwe, lweSample *buf, u64 numSamples, short *solution, int numSolvedCoordinates, fftw_complex *in, int fftPositions, int *guess, int bruteForcePositions)
{
    int q = lwe->q;
    int n = lwe->n;
    ASSERT(0 <= numSolvedCoordinates && numSolvedCoordinates <= n, "Bad parameter numSolvedCoordinates!\n");
    int startIndex = n - numSolvedCoordinates - fftPositions - bruteForcePositions; /* Subtract positions due to the BF guessing */

    /* process all samples */
    for (u64 i=0; i<numSamples; i++)
    {
        lweSample *sample = &buf[i];

        /* compute partial sum for solved part */
        u64 solvedSum = 0;

        /* Calculates the contribution from the guessing part */
        for (int j = 0; j < bruteForcePositions; j++)
        {
            solvedSum += columnValue(sample, n - numSolvedCoordinates - bruteForcePositions + j) * guess[j];
        }

        for (int j=n-numSolvedCoordinates; j<n; j++)
        {
            solvedSum += columnValue(sample, j) * solution[j];
        }
        solvedSum = solvedSum % q;

        u64 se = sumWithError(sample);
        se = (se - solvedSum) % q;

        /* update function value (index) to update */
        int p1 = columnValue(sample, startIndex);
        u64 currentIndex = p1;
        ASSERT(0 <= p1 && p1 < q, "unexpected p1");
        if (fftPositions >= 2)
        {
            int p2 = columnValue(sample, startIndex + 1);
            ASSERT(0 <= p2 && p2 < q, "unexpected p2");
            currentIndex = q * currentIndex + p2;
        }
        if (fftPositions >= 3)
        {
            int p3 = columnValue(sample, startIndex + 2);
            ASSERT(0 <= p3 && p3 < q, "unexpected p3");
            currentIndex = q * currentIndex + p3;
        }

        /* update fft input function */
        fftw_complex toAdd = cexp(TAU*I*se/q);
        in[currentIndex] = in[currentIndex] + toAdd;
    }
}

/* Hybrid category processing - single-precision case */
static void processOneCategoryHybridSingle(lweInstance *lwe, lweSample *buf, u64 numSamples, short *solution, int numSolvedCoordinates, fftwf_complex *in, int fftPositions, int *guess, int bruteForcePositions)
{
    int q = lwe->q;
    int n = lwe->n;
    ASSERT(0 <= numSolvedCoordinates && numSolvedCoordinates <= n, "Bad parameter numSolvedCoordinates!\n");
    int startIndex = n - numSolvedCoordinates - fftPositions - bruteForcePositions; /* Subtract positions due to the BF guessing */

    /* process all samples */
    for (u64 i=0; i<numSamples; i++)
    {
        lweSample *sample = &buf[i];

        /* compute partial sum for solved part */
        u64 solvedSum = 0;

        /* Calculates the contribution from the guessing part */
        for (int j = 0; j < bruteForcePositions; j++)
        {
            solvedSum += columnValue(sample, n - numSolvedCoordinates - bruteForcePositions + j) * guess[j];
        }

        for (int j=n-numSolvedCoordinates; j<n; j++)
        {
            solvedSum += columnValue(sample, j) * solution[j];
        }
        solvedSum = solvedSum % q;

        u64 se = sumWithError(sample);
        se = (se - solvedSum) % q;

        /* update function value (index) to update */
        int p1 = columnValue(sample, startIndex);
        u64 currentIndex = p1;
        ASSERT(0 <= p1 && p1 < q, "unexpected p1");
        if (fftPositions >= 2)
        {
            int p2 = columnValue(sample, startIndex + 1);
            ASSERT(0 <= p2 && p2 < q, "unexpected p2");
            currentIndex = q * currentIndex + p2;
        }
        if (fftPositions >= 3)
        {
            int p3 = columnValue(sample, startIndex + 2);
            ASSERT(0 <= p3 && p3 < q, "unexpected p3");
            currentIndex = q * currentIndex + p3;
        }

        /* update fft input function */
        fftwf_complex toAdd = cexp(TAU*I*se/q);
        in[currentIndex] = in[currentIndex] + toAdd;
    }
}

/* find the solution by looking at the output of the fft - 1 position and double-precision case */
double updateSolution1dHybrid(fftw_complex *out, short *solution, int startIndex, int q, double current_max)
{
    double current_value;
    double local_max = -42; /* Arbitrary low value */

    for (int i = 0; i < q; i++)
    {
        current_value = creal(out[i]);

        if (current_value > local_max)
        {
            local_max = current_value;
        }

        if (current_value > current_max)
        {
            current_max = current_value;
            solution[startIndex] = i;
        }
    }

    return local_max;
}

/* find the solution by looking at the output of the fft - 1 position and single-precision case */
float updateSolution1dHybridSingle(fftwf_complex *out, short *solution, int startIndex, int q, double current_max)
{
    float current_value;
    float local_max = -42; /* Arbitrary low value */

    for (int i = 0; i < q; i++)
    {
        current_value = creal(out[i]);

        if (current_value > local_max)
        {
            local_max = current_value;
        }

        if (current_value > current_max)
        {
            current_max = current_value;
            solution[startIndex] = i;
        }
    }

    return local_max;
}

/* find the solution by looking at the output of the fft - 2 positions and double-precision case */
double updateSolution2dHybrid(fftw_complex *out, short *solution, int startIndex, int q, double current_max)
{
    double current_value;
    double local_max = -42; /* Arbitrary low value */

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            current_value = creal(out[j + q*i]);

            if (current_value > local_max)
            {
                local_max = current_value;
            }

            if (current_value > current_max)
            {
                current_max = current_value;
                solution[startIndex] = i;
                solution[startIndex + 1] = j;
            }
        }
    }

    return local_max;
}

/* find the solution by looking at the output of the fft - 2 positions and single-precision case */
float updateSolution2dHybridSingle(fftwf_complex *out, short *solution, int startIndex, int q, double current_max)
{
    float current_value;
    float local_max = -42; /* Arbitrary low value */

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            current_value = creal(out[j + q*i]);

            if (current_value > local_max)
            {
                local_max = current_value;
            }

            if (current_value > current_max)
            {
                current_max = current_value;
                solution[startIndex] = i;
                solution[startIndex + 1] = j;
            }
        }
    }

    return local_max;
}

/* find the solution by looking at the output of the fft - 3 positions and double-precision case */
double updateSolution3dHybrid(fftw_complex *out, short *solution, int startIndex, int q, double current_max)
{
    double current_value;
    double local_max = -42; /* Arbitrary low value */

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            for (int k = 0; k < q; k++)
            {
                current_value = creal(out[k + q*(j + q*i)]);

                if (current_value > local_max)
                {
                    local_max = current_value;
                }

                if (current_value > current_max)
                {
                    current_max = current_value;
                    solution[startIndex] = i;
                    solution[startIndex + 1] = j;
                    solution[startIndex + 2] = k;
                }
            }
        }
    }

    return local_max;
}

/* find the solution by looking at the output of the fft - 3 positions and double-precision case */
float updateSolution3dHybridSingle(fftwf_complex *out, short *solution, int startIndex, int q, double current_max)
{
    float current_value;
    float local_max = -42; /* Arbitrary low value */

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            for (int k = 0; k < q; k++)
            {
                current_value = creal(out[k + q*(j + q*i)]);

                if (current_value > local_max)
                {
                    local_max = current_value;
                }

                if (current_value > current_max)
                {
                    current_max = current_value;
                    solution[startIndex] = i;
                    solution[startIndex + 1] = j;
                    solution[startIndex + 2] = k;
                }
            }
        }
    }

    return local_max;
}

/* Hybrid solver that uses sparse brute-force for bruteForcePositions number of positions and FFT for fftPositions number of positions */
int solve_fft_search_hybrid(const char *srcFolder, short *solution, int numSolvedCoordinates, int fftPositions, int bruteForcePositions, int doublePrecision)
{
    ASSERT(1 <= fftPositions && fftPositions <= 3, "The number of positions for fft is not supported!\n");

    lweInstance lwe;
    lweParametersFromFile(&lwe, srcFolder);
    int n = lwe.n;
    int q = lwe.q;
    double sigma = lwe.alpha * q; /* The noise level */

    u64 numCategories;
    sampleInfoFromFile(srcFolder, NULL, &numCategories, NULL, NULL, NULL);

    if (numSolvedCoordinates == n)
    {
        lweDestroy(&lwe);
        return 1; /* all positions solved */
    }

    /* create function to be fft input and output */
    fftw_complex *in, *out;
    fftwf_complex *inS, *outS;
    u64 numFftSlots;
    switch (fftPositions)
    {
    case 1:
        numFftSlots = q;
        break;
    case 2:
        numFftSlots = q * q;
        break;
    case 3:
        numFftSlots = q * q * q;
        break;
    default:
        ASSERT_ALWAYS("unsupported number of positions");
        exit(1);
    }

    if (doublePrecision)   /* double-precision */
    {
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFftSlots);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFftSlots);
        ASSERT(in, "allocation failed");
        ASSERT(out, "allocation failed");
        if (!in || !out)
        {
            char s[128];
            printf("*** solve_plain_bkw_sorted: allocation failure (tried to allocate %s bytes)\n", sprintf_u64_delim(s, sizeof(fftw_complex) * numFftSlots));
            lweDestroy(&lwe);
            return 3;
        }
    }
    else     /* single-precision */
    {
        inS = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * numFftSlots);
        outS = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * numFftSlots);
        ASSERT(inS, "allocation failed");
        ASSERT(outS, "allocation failed");
        if (!inS || !outS)
        {
            char s[128];
            printf("*** solve_plain_bkw_sorted: allocation failure (tried to allocate %s bytes)\n", sprintf_u64_delim(s, sizeof(fftwf_complex) * numFftSlots));
            lweDestroy(&lwe);
            return 4;
        }
    }

    /* Checks that we don't try to solve more positions than there is */
    while (numSolvedCoordinates + fftPositions + bruteForcePositions > n)
    {
        if (bruteForcePositions > 0)
        {
            bruteForcePositions--;
            continue;
        }
        else
        {
            fftPositions--;
        }
    }

    /* Variable declaration related to the storage reader */
    int ret;

    /* Variable declaration related sample statistics */
    u64 categoryIndexCounter;
    lweSample *buf1;
    lweSample *buf2;
    u64 numSamplesInBuf1, numSamplesInBuf2;
    int startIndex = n - fftPositions - numSolvedCoordinates - bruteForcePositions; /* The position where the FFT should start */

    int upperLimitGuess = (int) 3*sigma;
    printf("Upper limit of the guess is: %d\n", upperLimitGuess);
    double current_max = -42; // arbitrary low value for the guesses
    double current_value; // Current value for the guesses

    /* Initial guess is that all positions are equal to the lower limit */
    int guess[bruteForcePositions];
    int j;/* j is a loop variable to be used below */
    for(int i = 0 ; i < bruteForcePositions ; i++)
    {
        guess[i] = -upperLimitGuess;
    }

    while (1)
    {

        /* Initialize the storage reader */
        storageReader sr;
        ret = storageReaderInitialize(&sr, srcFolder);
        if (ret)
        {
            printf("*** solve_plain_bkw_sorted: storage reader returned %d on initialize\n", ret);
            return 2;
        }

        /* collect statistics from samples */
        categoryIndexCounter = 0;
        ret = storageReaderGetNextAdjacentCategoryPair(&sr, &buf1, &numSamplesInBuf1, &buf2, &numSamplesInBuf2);

        /* Resets the values of the FFT input function */
        if (doublePrecision)
        {
            for (j = 0; j < (int)numFftSlots; j++)
            {
                in[j] = 0;
            }
        }
        else
        {
            for (j = 0; j < (int)numFftSlots; j++)
            {
                inS[j] = 0;
            }
        }

        while (ret)
        {

            /* process buf1 */
            if (buf1)
            {
                if (doublePrecision)
                {
                    categoryIndexCounter++;
                    processOneCategoryHybrid(&lwe, buf1, numSamplesInBuf1, solution, numSolvedCoordinates, in, fftPositions, guess, bruteForcePositions); /* Here one position is guessed */
                }
                else
                {
                    categoryIndexCounter++;
                    processOneCategoryHybridSingle(&lwe, buf1, numSamplesInBuf1, solution, numSolvedCoordinates, inS, fftPositions, guess, bruteForcePositions);
                }
            }

            /* process buf1 */
            if (buf2)
            {
                if (doublePrecision)
                {
                    categoryIndexCounter++;
                    processOneCategoryHybrid(&lwe, buf2, numSamplesInBuf2, solution, numSolvedCoordinates, in, fftPositions, guess, bruteForcePositions); /* Here one position is guessed */
                }
                else
                {
                    categoryIndexCounter++;
                    processOneCategoryHybridSingle(&lwe, buf2, numSamplesInBuf2, solution, numSolvedCoordinates, inS, fftPositions, guess, bruteForcePositions);
                }
            }

            ret = storageReaderGetNextAdjacentCategoryPair(&sr, &buf1, &numSamplesInBuf1, &buf2, &numSamplesInBuf2);
        }
        if (categoryIndexCounter != numCategories)
        {
            printf("*** solve_plain_bkw_sorted: %" PRIu64 " categories processed (%" PRIu64 " expected)\n", categoryIndexCounter, numCategories);
        }

        switch (fftPositions)
        {
        case 1:
            if (doublePrecision)   /* double-precision */
            {
                calculateFft1d(in, out, q); /* Calculate fft */
                current_value = updateSolution1dHybrid(out, solution, startIndex, q, current_max); /* deduce position values from fft, given current value of guess */
            }
            else     /* single-precision */
            {
                calculateFft1dSingle(inS, outS, q); /* Calculate fft */
                current_value = updateSolution1dHybridSingle(outS, solution, startIndex, q, current_max); /* deduce position values from fft, given current value of guess */
            }
            break;
        case 2:
            if (doublePrecision)   /* double-precision */
            {
                calculateFft2d(in, out, q); /* Calculate fft */
                current_value = updateSolution2dHybrid(out, solution, startIndex, q, current_max); /* deduce position values from fft, given current value of guess */
            }
            else     /* single-precision */
            {
                calculateFft2dSingle(inS, outS, q); /* Calculate fft */
                current_value = updateSolution2dHybridSingle(outS, solution, startIndex, q, current_max); /* deduce position values from fft, given current value of guess */
            }
            break;
        case 3:
            if (doublePrecision)   /* double-precision */
            {
                calculateFft3d(in, out, q); /* Calculate fft */
                current_value = updateSolution3dHybrid(out, solution, startIndex, q, current_max); /* deduce position values from fft, given current value of guess */

            }
            else     /* single-precision */
            {
                calculateFft3dSingle(inS, outS, q); /* Calculate fft */
                current_value = updateSolution3dHybridSingle(outS, solution, startIndex, q, current_max); /* deduce position values from fft, given current value of guess */
            }
            break;
        default:
            ASSERT_ALWAYS("unsupported number of positions");
            exit(1);
        }

        /* Update the guessing values if we have found a better solution */
        if (current_value > current_max)
        {
            current_max = current_value;
            for (int k = 0; k < bruteForcePositions; k++)
            {
                solution[n - numSolvedCoordinates - bruteForcePositions + k] = (q + guess[k])%q; /* Save guess and make each value positive */
            }
        }

        printf("Current FFT value is: %lf, for guess = (", current_value);
        for (int k = 0; k < bruteForcePositions; k++)
        {
            if (k < bruteForcePositions - 1)
            {
                printf("%d, ", (q + guess[k])%q);
            }
            else
            {
                printf("%d", (q + guess[k])%q);
            }
        }
        printf(")\n");

        storageReaderFree(&sr); /* Free the current storage reader - preferable not to have to do this each lap */

        /* Update the guess */
        for(j = bruteForcePositions - 1 ; j >= 0 ; j--)
        {
            guess[j]++;
            if(guess[j] <= upperLimitGuess)
            {
                break;
            }
            else
            {
                guess[j] = -upperLimitGuess;
            }
        }
        if(j < 0)
        {
            break;
        }

    }

    /* cleanup */
    // storageReaderFree(&sr);
    if (doublePrecision)   /* double-precision */
    {
        fftw_free(in);
        fftw_free(out);
    }
    else     /* single-precision */
    {
        fftw_free(inS);
        fftw_free(outS);
    }

    lweDestroy(&lwe);
    return 0;
}

/* HYBRID CONTENT ENDS HERE */

#if 0
void test_fft_solver(const char *srcFolder)
{
    lweInstance lwe;
    lweParametersFromFile(&lwe, srcFolder);
    int q = lwe.q; // Toy number
    int N = q*q*q; // 3 dimensions

    // Double precision
    /*
    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_3d(q, q, q, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    */

    // Single precision
    fftwf_complex *in, *out;
    fftwf_plan p;

    in = (fftwf_complex*) fftw_malloc(sizeof(fftwf_complex) * N);
    out = (fftwf_complex*) fftw_malloc(sizeof(fftwf_complex) * N);
    p = fftwf_plan_dft_3d(q, q, q, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftwf_execute(p); /* repeat as needed */

    fftwf_destroy_plan(p);
    fftwf_free(out);
    fftwf_free(in);
}

void test_fft_pruned(void)
{
    int i, j;
    int N = 128*128;
    int K = 128;

    fftw_complex in[N], out[N], twids[(K-1)*(N/K-1)];

    /* With malloc
    fftwf_complex *in, *out, *twids;
    in = (fftwf_complex*) fftw_malloc(sizeof(fftwf_complex) * N);
    out = (fftwf_complex*) fftw_malloc(sizeof(fftwf_complex) * K);
    twids = (fftwf_complex*) fftw_malloc(sizeof(fftwf_complex) * (K-1)*(N/K-1));
    */


    fftw_plan plan;

    /* plan N/K FFTs of size K */
    plan = fftw_plan_many_dft(1, &K, N/K,
                              in, NULL, N/K, 1,
                              out, NULL, 1, K,
                              FFTW_FORWARD, FFTW_ESTIMATE);

    /* precompute twiddle factors (since we usually want more than one FFT) */
    for (j = 1; j < N/K; ++j)
        for (i = 1; i < K; ++i)
            twids[(j-1)*(K-1) + (i-1)] = cexp((I * FFTW_FORWARD * TAU/N) * (i*j));

//  ...initialize in[N] data....

    fftw_execute(plan);

    /* set *only* first K outputs, in out[], to values for full size-N DFT: */
    for (j = 1; j < N/K; ++j)
    {
        out[0] += out[j*K];
        for (i = 1; i < K; ++i)
            out[i] += out[i + j*K] * twids[(j-1)*(K-1) + (i-1)];
    }

    fftw_destroy_plan(plan);
}
#endif
