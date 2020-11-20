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

#include "test_functions.h"
#include "memory_utils.h"
#include "log_utils.h"
#include "string_utils.h"
#include "storage_file_utilities.h"
#include "position_values_2_category_index.h"
#include "transition_unsorted_2_sorted.h"
#include "storage_reader.h"
#include "verify_samples.h"
#include "linear_algebra_modular.h"
#include <stdio.h>
#include <inttypes.h>
#include <time.h>

#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))

/*static void print_matrix(time_t start, int n, short **matrix)
{
    for (int i=0; i<n; i++)
    {
        timeStamp(start);
        for (int j=0; j<n; j++)
        {
            printf(" %6d", matrix[i][j]);
        }
        printf("\n");
    }
}
*/

/* A*B = C */
/*static void matrix_mul(int n, int q, short **A, short **B, short **C)
{
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            int t = 0;
            for (int k=0; k<n; k++)
            {
                t += A[i][k] * B[k][j];
            }
            C[i][j] = t % q;
        }
    }
}
*/
/*void test_modular_matrix_inversion(lweInstance *lwe, int num_samples)
{
    time_t start = time(NULL);
    int q = lwe->q;
    int n = lwe->n;
    // generate samples
    lweSample *lweSamples = MALLOC(num_samples * sizeof(lweSample));
    for (int i=0; i<num_samples; i++)
    {
        lwe->newInPlaceRandomSample(&lweSamples[i], n, q, lwe->sigma, &(lwe->rnd), lwe->s);
    }
    // matrix storages
    short **A = MALLOC(n * sizeof(short*));
    short **A_inverse = MALLOC(n * sizeof(short*));
    short **I = MALLOC(n * sizeof(short*));
    for (int i=0; i<n; i++)
    {
        A[i] = MALLOC(n * sizeof(short));
        A_inverse[i] = MALLOC(n * sizeof(short));
        I[i] = MALLOC(n * sizeof(short));
    }
    short *b = MALLOC(n * sizeof(short));
    // compute
    clock_t t1 = clock();
    int num_samples_used = compute_matrix_inverse_modular(lwe, lweSamples, num_samples, A, A_inverse, b);
    clock_t t2 = clock();
    if (num_samples_used == 0)
    {
        timeStamp(start);
        printf("Matrix not computed!\n");
        return;
    }
    // display
    timeStamp(start);
    printf("%d samples used to compute inverse\n", num_samples_used);
    timeStamp(start);
    printf("b = (%d", b[0]);
    for (int i=1; i<n; i++)
    {
        printf(",%d", b[i]);
    }
    printf(")\n");
    if (n <= 40)
    {
        timeStamp(start);
        printf("A =\n");
        print_matrix(start, n, A);
        printf("\n\n");
        timeStamp(start);
        printf("A_inverse =\n");
        print_matrix(start, n, A_inverse);
        printf("\n\n");
        timeStamp(start);
        printf("I =\n");
    }
    matrix_mul(n, q, A, A_inverse, I);
    int error = 0;
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            if (i == j && I[i][j] != 1)
            {
                error = 1;
                printf("***************************** error in I at (%d,%d)\n", i, j);
            }
            if (i != j && I[i][j] != 0)
            {
                error = 1;
                printf("***************************** error in I at (%d,%d)\n", i, j);
            }
        }
    }
    if (!error)
    {
        timeStamp(start);
        printf("Matrix inversion verified to be correct!\n");
    }
    if (n <= 40)
    {
        print_matrix(start, n, I);
        timeStamp(start);
        printf("\n");
    }
    printf("Time to compute modular matrix inverse of size %d x %d: %.3f seconds\n", n, n, (double)(t2 - t1)/CLOCKS_PER_SEC);

    // cleanup
    FREE(b);
    for (int i=0; i<n; i++)
    {
        FREE(I[i]);
        FREE(A_inverse[i]);
        FREE(A[i]);
    }
    FREE(I);
    FREE(A_inverse);
    FREE(A);
    FREE(lweSamples);
}


void readSamplesTest(char *folderName, int n, int start, int numSamples)
{
    lweSample lweSamples[100];
    ASSERT(numSamples <= 100, "Insufficient buffer size!\n");
    printf("Reading samples %d - %d: ", start, start + numSamples - 1);
    int numReadSamples = readSamplesFromSampleFile(lweSamples, folderName, start, numSamples);
    printf(" Read %d samples from file\n", numReadSamples);
    start = start < 0 ? 0 : start;
    for (int i=0; i<numReadSamples; i++)
    {
        printf("%2d: ", i + start);
        printSample(&lweSamples[i], n);
        printf("\n");
    }
}

void bkwFileBasedTest(int n, int q, double alpha)
{
    char folderName[256];
    time_t start = time(NULL);
    lweInstance lwe;

    sprintf(folderName, "Step0");
    if (newStorageFolder(&lwe, folderName, n, q, alpha))
    {
        printf("Could not create LWE instance folder (Directory already exists)\n");
        return;
    }

    printf("n = %d\nq = %d\nalpha = %f\ns = (%d", lwe.n, lwe.q, lwe.alpha, lwe.s[0]);
    for (int i=1; i<lwe.n; i++)
    {
        printf(",%d", lwe.s[i]);
    }
    printf(")\n");

    lweInstance lwe2;
    lweParametersFromFile(&lwe2, folderName);
    printf("\nRetrieved as:\nn = %d\nq = %d\nalpha = %f\ns = (%d", lwe2.n, lwe2.q, lwe2.alpha, lwe2.s[0]);
    for (int i=1; i<lwe2.n; i++)
    {
        printf(",%d", lwe2.s[i]);
    }
    printf(")\n\n");

    // Tries to add samples to the file of samples
    int nbrOfSamplesToAdd = 10;

    int numberOfAddedSamples = addSamplesToSampleFile(folderName, nbrOfSamplesToAdd, 0);
    timeStamp(start);
    printf("Added %d samples to file %s\n", numberOfAddedSamples);
    timeStamp(start);
    printf("%d samples\n\n", numSamplesInSampleFile(folderName));

    numberOfAddedSamples = addSamplesToSampleFile(folderName, nbrOfSamplesToAdd, 0);
    timeStamp(start);
    printf("Added %d samples\n", numberOfAddedSamples);
    timeStamp(start);
    printf("%d samples\n\n", numSamplesInSampleFile(folderName));

    // Tries to read the samples from the file of samples
    readSamplesTest(folderName, lwe.n, 0, numSamplesInSampleFile(folderName));
    printf("\n");
    readSamplesTest(folderName, lwe.n, 7, 18);
    printf("\n");
    readSamplesTest(folderName, lwe.n, 17, 28);
    printf("\n");
    readSamplesTest(folderName, lwe.n, -5, 11);
    printf("\n");
    readSamplesTest(folderName, lwe.n, -3, 59);
    printf("\n");
    readSamplesTest(folderName, lwe.n, 33, 11);
    printf("\n");

    deleteStorageFolder(folderName, 1, 1, 1);
}

void testPositionValue2IndexBWK(int q)
{
    int numTests = 10000000;
    int errors = 0;

    for (int i=0; i<numTests; i++)
    {
        short pn[2];
        int p1, p2;
        p1 = pn[0] = randomUtilInt(NULL, q);
        p2 = pn[1] = randomUtilInt(NULL, q);
        int p1_additive_complement = p1 == 0 ? 0 : q - p1;
        int p2_additive_complement = p2 == 0 ? 0 : q - p2;
        int i = position_values_2_category_index_plain_bkw(q, pn);
        short p_recovered[2];
        category_index_2_position_values_plain_bkw(q, i, p_recovered);
        int p1_recovered = p_recovered[0];
        int p2_recovered = p_recovered[1];
        if (p1_recovered != p1 || p2_recovered != p2)
        {
            printf("Recovery test 1 failed for (p1, p2) = (%4d, %4d) -> index %d -> (%4d, %4d)\n", p1, p2, i, p1_recovered, p2_recovered);
            errors++;
        }
        pn[0] = p1_additive_complement;
        pn[1] = p2_additive_complement;
        int j = position_values_2_category_index_plain_bkw(q, pn); // passing additive complements
        category_index_2_position_values_plain_bkw(q, j, p_recovered);
        p1_recovered = p_recovered[0];
        p2_recovered = p_recovered[1];
        if (p1_recovered != p1_additive_complement || p2_recovered != p2_additive_complement)
        {
            printf("Recovery test 2 failed for (p1, p2) = (%4d, %4d) -> index %d -> (%4d, %4d)\n", p1_additive_complement, p2_additive_complement, j, p1_recovered, p2_recovered);
            errors++;
        }
        if (abs(i - j) > 1)
        {
            printf("Distance test failed for (p1, p2) = (%4d, %4d), (i, j) = (%10d, %10d)\n", p1, p2, i, j);
            errors++;
        }
    }
    printf("Performed %d randomized tests\n", numTests);
    if (errors)
    {
        printf("%d errors encountered!\n", errors);
    }
    else
    {
        printf("Seems to be working fine!\n\n");
    }
}
*/

void testCreateNewInstanceFolder(const char *folderName, int n, int q, double alpha)
{
    lweInstance lwe;
    int ret = newStorageFolder(&lwe, folderName, n, q, alpha);
    switch (ret)
    {
    case 0:
        printf("instance successfully created\n");
        break;
    case 1:
        printf("directory already exists\n");
        break;
    case 2:
        printf("could not create parameter file\n");
        break;
    case 3:
        printf("could not create samples file\n");
        break;
    default:
        printf("unknown fault\n");
    }
}

/*
void testAddSamplesToSampleFile(const char *folderName, u64 numSamplesToAdd)
{
    printf("num samples in sample file before: %" PRIu64 "\n", numSamplesInSampleFile(folderName));
    printf("trying to add %" PRIu64 " samples...", numSamplesToAdd);
    u64 ret = addSamplesToSampleFile(folderName, numSamplesToAdd, 0);
    printf("done\n");
    printf("num samples added                : %" PRIu64 "\n", ret);
    printf("num samples in sample file after : %" PRIu64 "\n", numSamplesInSampleFile(folderName));
}

u64 expectedNumBKWsamplesToGenerate(int numCategories, int *numSamplesPerCategory)
{
    ASSERT(numCategories & 1, "odd number of categories expected");
    u64 t = numSamplesPerCategory[0];
    u64 total = t * (t - 1) / 2;
    for (int i=1; i <numCategories; i+=2)
    {
        u64 t1 = numSamplesPerCategory[i];
        u64 t2 = numSamplesPerCategory[i + 1];
        total += t1 * (t1 - 1) / 2;
        total += t2 * (t2 - 1) / 2;
        total += t1 * t2;
    }
    return total;
}

static void printErrorDistribution(u64 *h, int q)
{
    int previousZero = 0;
    int previousPrinted = 1;
    for (int i=q-1; i>=0; i--)
    {
        if (h[i] == 0)
        {
            if (previousZero)
            {
                if (previousPrinted)
                {
                    printf("...\n");
                }
                previousPrinted = 0;
            }
            else
            {
                printf("%4d: 0\n", i);
                previousZero = 1;
            }
        }
        else
        {
            printf("%4d: %" PRIu64 "\n", i, h[i]);
            previousZero = 0;
            previousPrinted = 1;
        }
    }
}
*/

void printSampleVerificationOfSortedFolder(const char *srcFolderName, time_t start, bkwStepParameters *bkwStepPar)
{
    u64 totalNumSamplesProcessed, numIncorrectSums, numIncorrectHashes, numIncorrectCategoryClassifications;
    int printOnError = 1;
    timeStamp(start);
    printf("verifying (%s-sorted) samples in %s ... ", sortingAsString(bkwStepPar->sorting), srcFolderName);
    verifySortedSamples(srcFolderName, bkwStepPar, &totalNumSamplesProcessed, &numIncorrectSums, &numIncorrectHashes, &numIncorrectCategoryClassifications, printOnError);
    printf("done!\n");
    if (numIncorrectSums == 0 && numIncorrectHashes == 0 && numIncorrectCategoryClassifications == 0)
    {
        timeStamp(start);
        char s[256];
        printf("all %s samples verified ok\n", sprintf_u64_delim(s, totalNumSamplesProcessed));
    }
    else
    {
        char s[256];
        timeStamp(start);
        printf("%s samples processed\n", sprintf_u64_delim(s, totalNumSamplesProcessed));
        printf("\n  num incorrect sums = %s\n", sprintf_u64_delim(s, numIncorrectSums));
        printf("  num incorrect hashes = %s\n", sprintf_u64_delim(s, numIncorrectHashes));
        printf("  num incorrect category classifications = %s\n", sprintf_u64_delim(s, numIncorrectCategoryClassifications));
    }
    printf("\n");
}

/* utility of LPN */
void printBinarySampleVerification(const char *folderName, time_t start)
{
    u64 numIncorrectSums, numIncorrectHashes, numSamplesProcessed;
    int printOnError = 0;

    timeStamp(start);
    printf("verifying (unsorted) samples in %s ... ", folderName);
    verifyUnsortedSamples(folderName, &numSamplesProcessed, &numIncorrectSums, &numIncorrectHashes, printOnError);
    printf("done!\n");

    /* print statistics summary */
    char s1[256], s2[256];
    if (numIncorrectHashes != 0)
    {
        char s[256];
        timeStamp(start);
        printf("  num incorrect hashes = %s\n", sprintf_u64_delim(s, numIncorrectHashes));
    }
    timeStamp(start);
    printf("num incorrect sums = %s over %s processed\n", sprintf_u64_delim(s1, numIncorrectSums), sprintf_u64_delim(s2, numSamplesProcessed));
    timeStamp(start);
    printf("Error rate: %lf \n", (double)numIncorrectSums/(double)numSamplesProcessed);
    printf("\n");
}

void printSampleVerificationOfUnsortedFolder(const char *folderName, time_t start)
{
    u64 numIncorrectSums, numIncorrectHashes, numSamplesProcessed;
    int printOnError = 1;

    timeStamp(start);
    printf("verifying (unsorted) samples in %s ... ", folderName);
    verifyUnsortedSamples(folderName, &numSamplesProcessed, &numIncorrectSums, &numIncorrectHashes, printOnError);
    printf("done!\n");

    /* print statistics summary */
    timeStamp(start);
    char s[256];
    if (numIncorrectSums == 0 && numIncorrectHashes == 0)
    {
        printf("all %s samples verified ok\n", sprintf_u64_delim(s, numSamplesProcessed));
    }
    else
    {
        printf("%s samples processed\n", sprintf_u64_delim(s, numSamplesProcessed));
        printf("\n  num incorrect sums = %s\n", sprintf_u64_delim(s, numIncorrectSums));
        printf("  num incorrect hashes = %s\n", sprintf_u64_delim(s, numIncorrectHashes));
    }
    printf("\n");
}

int tuDarmstadtFileFormatConversionWithErrorChecking(const char *srcFileName, const char *dstFolderName, int useSampleAmplification, u64 totalNumSamples, time_t start)
{
    timeStamp(start);
    printf("src file  : %s\n", srcFileName);
    timeStamp(start);
    printf("dst folder: %s\n", dstFolderName);
    if (folderExists(dstFolderName))
    {
        return 1; /* tu darmstadt file already converted (destination folder already exists) */
    }
    char s[256];
    timeStamp(start);
    printf("amplifying to %s initial samples (random-ish triplet combinations from original samples)\n", sprintf_u64_delim(s, totalNumSamples));
    lweInstance lwe;
    int ret = convertTUDarmstadtProblemInstanceToNativeFormat(&lwe, srcFileName, dstFolderName, useSampleAmplification, totalNumSamples);
    if (ret)
    {
        timeStamp(start);
        printf("conversion from tu darmstadt file did not work (error %d)\n", ret);
        ASSERT_ALWAYS("conversion from tu darmstadt file did not work");
        return 2;
    }
#if 0
    /* print parameters */
    printf("n = %d\n", lwe.n);
    printf("q = %d\n", lwe.q);
    printf("alpha = %lf\n", lwe.alpha);
    printf("sigma = %lf\n", lwe.sigma);

    /* print samples */
    u64 numSamples = numSamplesInSampleFile(dstFolderName);
    char s[256];
    printf("num samples on file = %s\n", sprintf_u64_delim(s, numSamples));

    int numSamplesToPrint = 10;
    lweSample *sampleBuf = (lweSample*)MALLOC(numSamplesToPrint * LWE_SAMPLE_SIZE_IN_BYTES);
    u64 numSamplesRead = readSamplesFromSampleFile(sampleBuf, dstFolderName, 0, numSamplesToPrint);

    printf("%d first samples:\n", numSamplesToPrint);
    for (int i=0; i<numSamplesRead; i++)
    {
        printf("%10d: (%4d, [%d", i, sampleBuf[i].sumWithError, columnValue(&sampleBuf[i], 0));
        for (int j=1; j<lwe.n; j++)
        {
            printf(" %d", columnValue(&sampleBuf[i], j));
        }
        printf("])\n");
    }
    printf("\n");

    printf("%d last samples:\n", numSamplesToPrint);
    numSamplesRead = readSamplesFromSampleFile(sampleBuf, dstFolderName, numSamples - numSamplesToPrint, numSamplesToPrint);
    for (int i=0; i<numSamplesRead; i++)
    {
        char s[256];
        printf("%10s: (%4d, [%d", sprintf_u64_delim(s, numSamples - numSamplesToPrint + i), sampleBuf[i].sumWithError, columnValue(&sampleBuf[i], 0));
        for (int j=1; j<lwe.n; j++)
        {
            printf(" %d", columnValue(&sampleBuf[i], j));
        }
        printf("])\n");
    }

    FREE(sampleBuf);
#endif
    return 0;
}
