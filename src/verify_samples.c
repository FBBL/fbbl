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

#include "verify_samples.h"
#include "lwe_instance.h"
#include "storage_file_utilities.h"
#include "storage_reader.h"
#include "memory_utils.h"
#include "position_values_2_category_index.h"
#include <inttypes.h>

static void verifySum(lweInstance *lwe, lweSample *sample, u64 *numIncorrectSums, int printOnError)
{
    if (error(sample) == -1)
    {
        return; /* error term unknown, cannot verify sum */
    }
    if (lwe->s[0] == -1)
    {
        return; /* secret vector unknown, cannot verify sum */
    }
    int sum = sumWithError(sample);
    int computedSum = 0;
    for (int i=0; i<lwe->n; i++)
    {
        computedSum += columnValue(sample, i) * lwe->s[i];
    }
    computedSum += error(sample);
    computedSum = computedSum % lwe->q;
    if (sum != computedSum)
    {
        *numIncorrectSums = *numIncorrectSums + 1;
        if (printOnError)
        {
            printf("\n  sum         = %" PRIu64 " *****************************************************************\n", sum);
            printf("  computedSum = %" PRIu64 " *****************************************************************\n", computedSum);
        }
    }
}

/* verification of column hash when computed as a hash of the remaining part of the sample column vector */
static void verifyPartialHash(lweInstance *lwe, lweSample *sample, int startIndex, u64 *numIncorrectHashes, int printOnError)
{
    u64 hash = columnHash(sample);
    u64 computedHash = bkwColumnComputeHash(sample, lwe->n, startIndex);
    if (hash != computedHash)
    {
        *numIncorrectHashes = *numIncorrectHashes + 1;
        if (printOnError)
        {
            printf("  hash         = %" PRIu64 " *****************************************************************\n", hash);
            printf("  computedHash = %" PRIu64 " *****************************************************************\n", computedHash);
        }
    }
}

/* verification of column hash when computed as a hash of the entire sample column vector, not just the remaining part */
static void verifyEntireHash(lweInstance *lwe, lweSample *sample, u64 *numIncorrectHashes, int printOnError)
{
    verifyPartialHash(lwe, sample, 0, numIncorrectHashes, printOnError);
}

static void verifyOneSampleUnsorted(lweInstance *lwe, lweSample *sample, u64 *numIncorrectSums, u64 *numIncorrectHashes, int printOnError)
{
    verifySum(lwe, sample, numIncorrectSums, printOnError);
    //verifyPartialHash(lwe, sample, startI*ndex, numIncorrectHashes, printOnError);
    verifyEntireHash(lwe, sample, numIncorrectHashes, printOnError);
}

void verifyOneSampleSorted(lweInstance *lwe, lweSample *sample, bkwStepParameters *bkwStepPar, u64 *numIncorrectSums, u64 *numIncorrectHashes, u64 expectedCategoryIndex, u64 *numIncorrectCategoryClassifications, int printOnError)
{
    verifySum(lwe, sample, numIncorrectSums, printOnError);
    //verifyPartialHash(lwe, sample, bkwStepPar->startIndex, numIncorrectHashes, printOnError);
    verifyEntireHash(lwe, sample, numIncorrectHashes, printOnError);
    u64 index = position_values_2_category_index(lwe, sample, bkwStepPar);
    if (index != expectedCategoryIndex)
    {
        *numIncorrectCategoryClassifications = *numIncorrectCategoryClassifications + 1;
        if (printOnError)
        {
            printf("  category classification error: expectedCategoryIndex = %" PRIu64 ", computed index = %" PRIu64 "\n", expectedCategoryIndex, index);
        }
    }
}

int verifyUnsortedSamples(const char *folderName, u64 *numSamplesProcessed, u64 *numIncorrectSums, u64 *numIncorrectHashes, int printOnError)
{
    /* read lwe parameters */
    lweInstance lwe;
    lweParametersFromFile(&lwe, folderName);

    /* create read buffer */
    int numSamples = 100000; /* read buffer size */
    lweSample *sampleBuf = MALLOC(numSamples * LWE_SAMPLE_SIZE_IN_BYTES);

    /* initialize counters */
    *numIncorrectSums = 0;
    *numIncorrectHashes = 0;
    *numSamplesProcessed = 0;

    /* verify samples */
    FILE *f = fopenSamples(folderName, "rb");
    if (!f)
    {
        return 1;
    }
    while (!feof(f))
    {
        u64 numRead = freadSamples(f, sampleBuf, numSamples);
        for (u64 i=0; i<numRead; i++)
        {
            lweSample *sample = &sampleBuf[i];
            verifyOneSampleUnsorted(&lwe, sample, numIncorrectSums, numIncorrectHashes, printOnError);
        }
        *numSamplesProcessed = *numSamplesProcessed + numRead;
    }
    FREE(sampleBuf);
    fclose(f);
    return 0;
}

int verifySortedSamples(const char *folderName, bkwStepParameters *bkwStepPar, u64 *totalNumSamplesProcessed, u64 *numIncorrectSums, u64 *numIncorrectHashes, u64 *numIncorrectCategoryClassifications, int printOnError)
{
    ASSERT(folderName, "unexpected parameter");
    ASSERT(bkwStepPar, "unexpected parameter");
    lweInstance lwe;
    u64 categoryIndexCounter = 0;
    lweSample *buf1;
    lweSample *buf2;
    u64 numSamplesInBuf1, numSamplesInBuf2;

    /* get parameters */
    lweParametersFromFile(&lwe, folderName);
    u64 numCategories = num_categories(&lwe, bkwStepPar);

    /* initialize parameters */
    *totalNumSamplesProcessed = 0;
    *numIncorrectSums = 0;
    *numIncorrectHashes = 0;
    *numIncorrectCategoryClassifications = 0;

    /* verify samples */
    storageReader sr;
    int ret = storageReaderInitialize(&sr, folderName);
    if (ret)
    {
        if (printOnError)
        {
            printf("storage reader returned %d on initialize\n", ret);
        }
        return 1;
    }
    ret = storageReaderGetNextAdjacentCategoryPair(&sr, &buf1, &numSamplesInBuf1, &buf2, &numSamplesInBuf2);
    while (ret)
    {

        if (buf1)
        {
            for (u64 i=0; i<numSamplesInBuf1; i++)
            {
                lweSample *sample = &buf1[i];
                verifyOneSampleSorted(&lwe, sample, bkwStepPar, numIncorrectSums, numIncorrectHashes, categoryIndexCounter, numIncorrectCategoryClassifications, printOnError);
            }
            *totalNumSamplesProcessed = *totalNumSamplesProcessed + numSamplesInBuf1;
            categoryIndexCounter++;
        }

        if (buf2)
        {
            for (u64 i=0; i<numSamplesInBuf2; i++)
            {
                lweSample *sample = &buf2[i];
                verifyOneSampleSorted(&lwe, sample, bkwStepPar, numIncorrectSums, numIncorrectHashes, categoryIndexCounter, numIncorrectCategoryClassifications, printOnError);
            }
            *totalNumSamplesProcessed = *totalNumSamplesProcessed + numSamplesInBuf2;
            categoryIndexCounter++;
        }

        ret = storageReaderGetNextAdjacentCategoryPair(&sr, &buf1, &numSamplesInBuf1, &buf2, &numSamplesInBuf2);
    }
    if (categoryIndexCounter != numCategories)
    {
        if (printOnError)
        {
            printf("*** getStdBkwStatsForSortedSamples: %" PRIu64 " categories processed (%" PRIu64 " expected)\n", categoryIndexCounter, numCategories);
        }
    }
    storageReaderFree(&sr);
    return 0;
}











#if 0
void printSampleDataForUnsortedSamples(lweInstance *lwe, const char *folderName, time_t start)
{
    u64 num = numSamplesInSampleFile(folderName);
    timeStamp(start);
    printf("num samples (according to file size) = %I64d\n\n", num);
    int *numSamples = MALLOC(lwe->q * lwe->q * sizeof(int));
    getStdBkwStatsForUnsortedSamples(folderName, 0, 2, numSamples, NULL);

    u64 maxInstancesInCategory = 0;
    u64 minInstancesInCategory = num;
    u64 totalNum = 0;
    for (int i=0; i<lwe->q*lwe->q; i++)
    {
        if (numSamples[i] > maxInstancesInCategory)
        {
            maxInstancesInCategory = numSamples[i];
        }
        if (numSamples[i] < minInstancesInCategory)
        {
            minInstancesInCategory = numSamples[i];
        }
        totalNum += numSamples[i];
    }
    int *density = MALLOC((maxInstancesInCategory + 1) * sizeof(int));
    MEMSET(density, 0, (maxInstancesInCategory + 1) * sizeof(int));
    for (int i=0; i<lwe->q*lwe->q; i++)
    {
        density[numSamples[i]]++;
    }
    timeStamp(start);
    printf("Number of categories with given weight:\n");
    for (int i=0; i<maxInstancesInCategory+1; i++)
    {
        timeStamp(start);
        printf("%3d: %d\n", i, density[i]);
    }
    FREE(density);
    timeStamp(start);
    printf("min num items in category = %9I64d\n", minInstancesInCategory);
    timeStamp(start);
    printf("max num items in category = %9I64d\n", maxInstancesInCategory);
    timeStamp(start);
    printf("total num items           = %9I64d\n", totalNum);
    FREE(numSamples);
}
#endif
