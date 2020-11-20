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

#include "transition_bkw_step_final_smooth_lms_meta.h"
#include "storage_file_utilities.h"
#include "memory_utils.h"
#include "log_utils.h"
#include "string_utils.h"
#include "lwe_sorting.h"
#include "storage_reader.h"
#include "storage_writer.h"
#include "position_values_2_category_index.h"
#include "lookup_tables.h"
#include "config_bkw.h"
#include <inttypes.h>
#include <math.h>

static u64 numZeroColumns;
static u64 numZeroColumnsAdd;

static u64 subtractSamples(lweInstance *lwe, lweSample *sample1, lweSample *sample2, bkwStepParameters *srcBkwStepPar, FILE *wf)
{
    int n = lwe->n;

    // printf("We try to subtract the samples\n");

    /* perform Unnatural Selection */
    if(srcBkwStepPar->sorting == smoothLMS && srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts)
    {
        int q_half = (lwe->q)/2;
        short tmp_a;
        for (int i=srcBkwStepPar->startIndex; i<srcBkwStepPar->startIndex + srcBkwStepPar->numPositions; i++)
        {
            tmp_a = diffTable(columnValue(sample1, i), columnValue(sample2, i));
            if( (tmp_a <= q_half && tmp_a >= srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts ) || (tmp_a > q_half && tmp_a <= lwe->q-srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts))
                return 0; // discard the sample
        }
    }

    //  int storageWriterStatus = 0;
    lweSample *newSample = MALLOC(LWE_SAMPLE_SIZE_IN_BYTES);

    if (newSample)
    {
        for (int i=0; i<n; i++)
        {
            newSample->col.a[i] = diffTable(columnValue(sample1, i),columnValue(sample2, i));
        }
        newSample->col.hash = bkwColumnComputeHash(newSample, n, 0 /* startRow */);
        int err1 = error(sample1);
        int err2 = error(sample2);
        newSample->error = (err1 == -1 || err2 == -1) ? -1 : diffTable(err1,err2); /* undefined if either parent error term is undefined */
        newSample->sumWithError = diffTable(sumWithError(sample1),sumWithError(sample2));
    }
    else
    {
        printf("********* ERROR transition_bkw_step_final_smooth_lms_meta.c: allocation for new sample failed\n");
        exit(-1);
    }

    /* discard zero columns (assuming that these are produced by coincidental cancellation due to sample amplification) */
    if (columnIsZero(newSample, n))
    {
        numZeroColumns++;
        numZeroColumnsAdd++;
        free(newSample);
        return 0; /* sample processed but not added */
    }

    int numWritten = fwrite(newSample, LWE_SAMPLE_SIZE_IN_BYTES, 1, wf);
    ASSERT(numWritten == 1, "Error in writing new sample\n");
    if (numWritten != 1)
    {
        printf("numWritten = %d\n", numWritten);
    }

    free(newSample);

    return 1; /* one sample processed (and actually added) */
}

static u64 addSamples(lweInstance *lwe, lweSample *sample1, lweSample *sample2, bkwStepParameters *srcBkwStepPar, FILE *wf)
{
    int n = lwe->n;

    // printf("We try to subtract the samples\n");

    /* perform Unnatural Selection */
    if(srcBkwStepPar->sorting == smoothLMS && srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts)
    {
        int q_half = (lwe->q)/2;
        short tmp_a;
        for (int i=srcBkwStepPar->startIndex; i<srcBkwStepPar->startIndex + srcBkwStepPar->numPositions; i++)
        {
            tmp_a = sumTable(columnValue(sample1, i), columnValue(sample2, i));
            if( (tmp_a <= q_half && tmp_a >= srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts ) || (tmp_a > q_half && tmp_a <= lwe->q-srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts))
                return 0; // discard the sample
        }
    }
    //  int storageWriterStatus = 0;
    lweSample *newSample = MALLOC(LWE_SAMPLE_SIZE_IN_BYTES);

    if (newSample)
    {
        for (int i=0; i<n; i++)
        {
            newSample->col.a[i] = sumTable(columnValue(sample1, i),columnValue(sample2, i));
        }
        newSample->col.hash = bkwColumnComputeHash(newSample, n, 0 /* startRow */);
        int err1 = error(sample1);
        int err2 = error(sample2);
        newSample->error = (err1 == -1 || err2 == -1) ? -1 : sumTable(err1,err2); /* undefined if either parent error term is undefined */
        newSample->sumWithError = sumTable(sumWithError(sample1),sumWithError(sample2));
    }
    else
    {
        printf("********* ERROR transition_bkw_step_final_smooth_lms_meta.c: allocation for new sample failed\n");
        exit(-1);
    }

    /* discard zero columns (assuming that these are produced by coincidental cancellation due to sample amplification) */
    if (columnIsZero(newSample, n))
    {
        numZeroColumns++;
        numZeroColumnsAdd++;
        free(newSample);
        return 0; /* sample processed but not added */
    }

    int numWritten = fwrite(newSample, LWE_SAMPLE_SIZE_IN_BYTES, 1, wf);
    ASSERT(numWritten == 1, "Error in writing new sample\n");
    if (numWritten != 1)
    {
        printf("numWritten = %d\n", numWritten);
    }

    free(newSample);
    return 1; /* one sample processed (and actually added) */
}

static inline int additiveInverse(int c, int val)
{
    if (c % 2 == 0)   /* Even number of categories - all categories have an additive inverse */
    {
        return c - val - 1;
    }
    else     /* Odd number of categories - The zero category is its own additive inverse */
    {
        return val == 0 ? 0 : c - val;
    }
}

static u64 processSingleCategoryLF1(lweInstance *lwe, lweSample **categorySamplePointers, int numSamplesInCategory, bkwStepParameters *srcBkwStepPar, FILE *wf, time_t start)
{
    u64 numAdded = 0;
    lweSample *firstSample;
    lweSample *sample;

    if (numSamplesInCategory < 2)
    {
        return 0;
    }

    /* subtract all samples from the first one (linear) */
    firstSample = categorySamplePointers[0];
    for (int i=1; i<numSamplesInCategory; i++)
    {
        sample = categorySamplePointers[i];
        numAdded += subtractSamples(lwe, firstSample, sample, srcBkwStepPar, wf);
    }
    return numAdded;
}

static u64 processAdjacentCategoriesLF1(lweInstance *lwe, lweSample **categorySamplePointers1, int numSamplesInCategory1, lweSample **categorySamplePointers2, int numSamplesInCategory2, bkwStepParameters *srcBkwStepPar, FILE *wf, time_t start)
{
    u64 numAdded = 0;
    lweSample *firstSample;
    lweSample *sample;

    if (numSamplesInCategory1 > 0)
    {
        /* LF1-process category 1 */
        /* subtract all samples from the first one (linear) */
        firstSample = categorySamplePointers1[0];
        for (int i=1; i<numSamplesInCategory1; i++)
        {
            sample = categorySamplePointers1[i];
            numAdded += subtractSamples(lwe, firstSample, sample, srcBkwStepPar, wf);
        }
        /* add all samples in adjacent category to first (same as above) sample (linear) */
        for (int i=0; i<numSamplesInCategory2; i++)
        {
            sample = categorySamplePointers2[i];
            numAdded += addSamples(lwe, firstSample, sample, srcBkwStepPar, wf);
        }
    }
    else     /* numSamplesInCategory1 == 0 */
    {
        if (numSamplesInCategory2 >= 2)
        {
            /* LF1-process samples in category 2 only */
            firstSample = categorySamplePointers2[0];
            for (int i=1; i<numSamplesInCategory2; i++)
            {
                sample = categorySamplePointers2[i];
                numAdded += subtractSamples(lwe, firstSample, sample, srcBkwStepPar, wf);
            }
        }
    }
    return numAdded;
}

static u64 processSingleCategoryLF2(lweInstance *lwe, lweSample **categorySamplePointers, int numSamplesInCategory, bkwStepParameters *srcBkwStepPar, FILE *wf, time_t start)
{
    u64 numAdded = 0;
    lweSample *sample1;
    lweSample *sample2;

    /* subtract all pairs of samples (quadratic) */
    for (int i=0; i<numSamplesInCategory; i++)
    {
        sample1 = categorySamplePointers[i];
        for (int j=i+1; j<numSamplesInCategory; j++)
        {
            sample2 = categorySamplePointers[j];
            numAdded += subtractSamples(lwe, sample1, sample2, srcBkwStepPar, wf);
        }
    }
    return numAdded;
}

static u64 processAdjacentCategoriesLF2(lweInstance *lwe, lweSample **categorySamplePointers1, int numSamplesInCategory1, lweSample **categorySamplePointers2, int numSamplesInCategory2, bkwStepParameters *srcBkwStepPar, FILE *wf, time_t start)
{
    u64 numAdded = 0;
    lweSample *sample1;
    lweSample *sample2;

    /* Paul's note: sample dependency may be reduced beyond that given by SAMPLE_DEPENDENCY_SMEARING combining samples in a smarter order (all LF1 samples first, then...) */

    /* process all pairs in category 1 (subtract sample pairs) */
    numAdded += processSingleCategoryLF2(lwe, categorySamplePointers1, numSamplesInCategory1, srcBkwStepPar, wf, start);

    /* process all pairs in category 2 (subtract sample pairs) */
    numAdded += processSingleCategoryLF2(lwe, categorySamplePointers2, numSamplesInCategory2, srcBkwStepPar, wf, start);

    /* process all pairs in categories 1 and 2 (add sample pairs) */
    for (int i=0; i<numSamplesInCategory1; i++)
    {
        sample1 = categorySamplePointers1[i];
        for (int j=0; j<numSamplesInCategory2; j++)
        {
            sample2 = categorySamplePointers2[j];
            numAdded += addSamples(lwe, sample1, sample2, srcBkwStepPar, wf);
        }
    }

    return numAdded;
}

static void getLastCategorySizes(int q, int n, bkwStepParameters *srcBkwStepPar, int *cMidPosition, int *cLastPosition)
{
    ASSERT(srcBkwStepPar->sortingPar.smoothLMS.meta_skipped > 0 && srcBkwStepPar->sortingPar.smoothLMS.meta_skipped <= 2, "Unsopported number of skipped positions!");
    int p = srcBkwStepPar->sortingPar.smoothLMS.p;
    int p1 = srcBkwStepPar->sortingPar.smoothLMS.p1;
    int q_ = q%2 == 1 ? (q+1)/2 : q/2;
    *cMidPosition = ((2*q_-1) % p) == 0 ? ((2*q_-1) / p) : ((2*q_-1) / p) + 1;
    if (srcBkwStepPar->startIndex + srcBkwStepPar->numPositions == n)   /* Last step of smooth LMS */
    {
        *cLastPosition = *cMidPosition;
    }
    else
    {
        // q_ = srcBkwStepPar->sortingPar.smoothLMS.prev_p1;
        *cLastPosition = ((2*q_-1) % p1) == 0 ? ((2*q_-1) / p1) : ((2*q_-1) / p1) + 1;
    }
}

static lweSample ***sortMetacategoryTwoPositions(lweSample *category, int numSamplesInCategory, lweInstance *lwe, bkwStepParameters *srcBkwStepPar, int *valueCounter, int n, int cMidPosition, int cLastPosition)
{
    int maxNumRepetitionsOfCategory = 0;
    int q = lwe->q;

    int pLast;
    int pSecondToLast = srcBkwStepPar->sortingPar.smoothLMS.p;

    int metaCategorySize = cMidPosition*cLastPosition;

    int q_ = q%2 == 1 ? (q+1)/2 : q/2;

    int lastPos = srcBkwStepPar->startIndex + srcBkwStepPar->numPositions; /* Index for the last position */
    if (lastPos == n)   /* Last step of smooth LMS */
    {
        lastPos--;
        pLast = srcBkwStepPar->sortingPar.smoothLMS.p;
    }
    else
    {
        pLast = srcBkwStepPar->sortingPar.smoothLMS.p1;
    }

    /* compute distribution of last position value */
    MEMSET(valueCounter, 0, metaCategorySize * sizeof(int));
    /* TODO: Implement this function... */
    for (int i=0; i<numSamplesInCategory; i++)
    {
        int pnLast = columnValue(&category[i], lastPos); /* Last position value */
        int pnSecondToLast = columnValue(&category[i], lastPos - 1); /* Second to last position value */
        int currentLastPositionCategory = positionSmoothLMSMap(pnLast, q, q_, pLast, cLastPosition); /* Last position category */
        int currentSecondToLastPositionCategory = positionSmoothLMSMap(pnSecondToLast, q, q_, pSecondToLast, cMidPosition); /* Second to last position category */
        int currentCategory = currentLastPositionCategory*cMidPosition + currentSecondToLastPositionCategory;
        // printf("The current category for sample %d is %d\n", i, currentCategory);
        valueCounter[currentCategory]++; /* Increase counter for corresponding category */
    }
    /* find maximum */
    for (int i=0; i<metaCategorySize; i++)   /* Loop over all categories */
    {
        if (valueCounter[i] > maxNumRepetitionsOfCategory)
        {
            maxNumRepetitionsOfCategory = valueCounter[i];
        }
    }

    /* allocate storage structure for (properly) sorted samples (no copying, pointers only) */
    lweSample ***metaCategory = (lweSample***)MALLOC(metaCategorySize * sizeof(lweSample**));

    for (int i=0; i<metaCategorySize; i++)
    {
        metaCategory[i] = (lweSample**)MALLOC(maxNumRepetitionsOfCategory * sizeof(lweSample*));
    }



    /* TODO: Finish this function! */
    /* sort */
    MEMSET(valueCounter, 0, metaCategorySize * sizeof(int));
    for (int i=0; i<numSamplesInCategory; i++)
    {
        lweSample *thisSample = &category[i];
        int pnLast = columnValue(thisSample, lastPos);
        int pnSecondToLast = columnValue(thisSample, lastPos - 1);
        int currentLastPositionCategory = positionSmoothLMSMap(pnLast, q, q_, pLast, cLastPosition); /* Last position category */
        int currentSecondToLastPositionCategory = positionSmoothLMSMap(pnSecondToLast, q, q_, pSecondToLast, cMidPosition); /* Second to last position category */
        int currentCategory = currentLastPositionCategory*cMidPosition + currentSecondToLastPositionCategory;
        metaCategory[currentCategory][valueCounter[currentCategory]++] = thisSample;
    }

    return metaCategory;
}

static lweSample ***sortMetacategoryOnePosition(lweSample *category, int numSamplesInCategory, lweInstance *lwe, bkwStepParameters *srcBkwStepPar, int *valueCounter, int n, int metaCategorySize)
{
    int maxNumRepetitionsOfLastPos = 0;
    int q = lwe->q;
    int p;
    int q_ = q%2 == 1 ? (q+1)/2 : q/2;
    int lastPos = srcBkwStepPar->startIndex + srcBkwStepPar->numPositions; /* Index for the last position */
    if (lastPos == n)   /* Last step of smooth LMS */
    {
        lastPos--;
        p = srcBkwStepPar->sortingPar.smoothLMS.p;
    }
    else
    {
        p = srcBkwStepPar->sortingPar.smoothLMS.p1;
    }

//  printf("Meta category size is %d\n", metaCategorySize);
//  printf("The whole category has %d samples\n", numSamplesInCategory);

    /* compute distribution of last position value */
    MEMSET(valueCounter, 0, metaCategorySize * sizeof(int));
    for (int i=0; i<numSamplesInCategory; i++)
    {
        int pn = columnValue(&category[i], lastPos); /* Last position value */
        int currentCategory = positionSmoothLMSMap(pn, q, q_, p, metaCategorySize); /* Corresponding category */
        valueCounter[currentCategory]++; /* Increase counter for corresponding category */
    }
    /* find maximum */
    for (int i=0; i<metaCategorySize; i++)
    {
//    printf("Category at index %d has %d samples\n", i, valueCounter[i]);
        if (valueCounter[i] > maxNumRepetitionsOfLastPos)
        {
            maxNumRepetitionsOfLastPos = valueCounter[i];
        }
    }

    /* allocate storage structure for (properly) sorted samples (no copying, pointers only) */
    lweSample ***metaCategory = (lweSample***)MALLOC(metaCategorySize * sizeof(lweSample**));
    for (int i=0; i<metaCategorySize; i++)
    {
        metaCategory[i] = (lweSample**)MALLOC(maxNumRepetitionsOfLastPos * sizeof(lweSample*));
    }

    /* sort */
    MEMSET(valueCounter, 0, metaCategorySize * sizeof(int));
    for (int i=0; i<numSamplesInCategory; i++)
    {
        lweSample *thisSample = &category[i];
        int pn = columnValue(&category[i], lastPos); /* Last position value */
        int currentCategory = positionSmoothLMSMap(pn, q, q_, p, metaCategorySize); /* Corresponding category */
        metaCategory[currentCategory][valueCounter[currentCategory]++] = thisSample;
    }

    return metaCategory;
}

int transition_bkw_step_final_smooth_lms_meta(const char *srcFolderName, const char *dstFolderName, bkwStepParameters *srcBkwStepPar, u64 *numSamplesStored, time_t start)
{

    if (folderExists(dstFolderName))   /* if destination folder already exists, assume that we have performed this reduction step already */
    {
        return 100; /* reduction step already performed (destination folder already exists) */
    }

    /* get lwe parameters from file */
    lweInstance lwe;
    lweParametersFromFile(&lwe, srcFolderName);

    /* get sample info from file */
    u64 srcNumCategories, srcCategoryCapacity, srcNumTotalSamples;

    if (sampleInfoFromFile(srcFolderName, srcBkwStepPar, &srcNumCategories, &srcCategoryCapacity, &srcNumTotalSamples, NULL))
    {
        return 1; /* error reading from samples info file */
    }

    /* TODO Calculate the number of metacategories when using smooth LMS with metacategories */
    char nc[256], cc[256], ns[256];
    timeStamp(start);
    printf("transition_bkw_step_lms_meta: num src categories is %s, category capacity is %s, total num src samples is %s (%5.2f%% full)\n", sprintf_u64_delim(nc, srcNumCategories), sprintf_u64_delim(cc, srcCategoryCapacity), sprintf_u64_delim(ns, srcNumTotalSamples), 100*srcNumTotalSamples/(double)(srcNumCategories * srcCategoryCapacity));
    if (srcBkwStepPar->sorting != smoothLMS)
    {
        return 2; /* unexpected sample sorting at src folder */
    }

    /* initialize storage reader */
    storageReader sr;
    if (storageReaderInitialize(&sr, srcFolderName))
    {
        ASSERT_ALWAYS("could not initialize storage reader");
        return 3; /* could not initialize storage reader */
    }

    /* Initialize destination folder and file */
    newStorageFolderWithGivenLweInstance(&lwe, dstFolderName);
    newStorageFolder(&lwe, dstFolderName, lwe.n, lwe.q, lwe.alpha);
    FILE *wf = fopenSamples(dstFolderName, "ab");
    if (!wf)
    {
        return -1;
    }

    /* initialize add and diff tables for faster operation */
    /* TODO: move to initialization */
    if (createSumAndDiffTables(lwe.q))
    {
        return 6; /* could not create addition and difference tables */
    }

    /* process samples */
    u64 cat = 0; /* current category index */
    u64 nextPrintLimit = 2;
    lweSample *buf1;
    lweSample *buf2;
    u64 numSamplesInBuf1, numSamplesInBuf2, numSamplesAdded = 0;
    int numReadCategories = storageReaderGetNextAdjacentCategoryPair(&sr, &buf1, &numSamplesInBuf1, &buf2, &numSamplesInBuf2);

    int meta_skipped = srcBkwStepPar->sortingPar.smoothLMS.meta_skipped;
    int cMidPosition;
    int cLastPosition;
    int metaCategorySize;
    getLastCategorySizes(lwe.q, lwe.n, srcBkwStepPar, &cMidPosition, &cLastPosition);

    if (meta_skipped == 1)
    {
        metaCategorySize = cLastPosition;
    }
    else if (meta_skipped == 2)
    {
        metaCategorySize = cMidPosition*cLastPosition;
    }
    else
    {
        ASSERT_ALWAYS("Unsupported number of skipped positions");
        return 7;
    }

    u64 abortSampleLimit = 4*MAX_NUM_SAMPLES/3;

    /* The main while loop */
    while (numReadCategories)
    {

        /* The meta category pair and the corresponding counter for number samples within each category within the meta category */
        int *valueCounter1 = NULL;
        int *valueCounter2 = NULL;
        lweSample ***metaCategory1 = NULL;
        lweSample ***metaCategory2 = NULL;

        /* pre-process metacategory 1 (if available) */
        if (buf1)
        {
            valueCounter1 = MALLOC(metaCategorySize * sizeof(int)); /* We should find out what this number is */
            if (meta_skipped == 1)
            {
                metaCategory1 = sortMetacategoryOnePosition(buf1, numSamplesInBuf1, &lwe, srcBkwStepPar, valueCounter1, lwe.n, metaCategorySize);
            }
            else if (meta_skipped == 2)
            {
                metaCategory1 = sortMetacategoryTwoPositions(buf1, numSamplesInBuf1, &lwe, srcBkwStepPar, valueCounter1, lwe.n, cMidPosition, cLastPosition);
            }
            else
            {
                ASSERT_ALWAYS("Unsupported number of skipped positions!");
            }
        }

        /* pre-process metacategory 2 (if available) */
        if (buf2)
        {
            ASSERT(buf1, "buf2 should only be available if buf1 is");
            valueCounter2 = MALLOC(metaCategorySize * sizeof(int));
            if (meta_skipped == 1)
            {
                metaCategory2 = sortMetacategoryOnePosition(buf2, numSamplesInBuf2, &lwe, srcBkwStepPar, valueCounter2, lwe.n, metaCategorySize);
            }
            else if (meta_skipped == 2)
            {
                metaCategory2 = sortMetacategoryTwoPositions(buf2, numSamplesInBuf2, &lwe, srcBkwStepPar, valueCounter2, lwe.n, cMidPosition, cLastPosition);
            }
            else
            {
                ASSERT_ALWAYS("Unsupported number of skipped positions!");
            }
        }

        /* Process meta categories based on selection method, single/pair and number of positions skipped */
        switch (srcBkwStepPar->selection)
        {
        case LF1:
            switch (numReadCategories)
            {
            case 1: /* single meta category (no corresponding meta category with first two coordinates having (differing) additive inverses) */
                ASSERT(metaCategory1, "unexpected parameter");
                for (int i=0; i<cLastPosition; i++)
                {
                    if (meta_skipped == 1)   /* One position skipped for meta categories */
                    {
                        numSamplesAdded += processSingleCategoryLF1(&lwe, metaCategory1[i], valueCounter1[i], srcBkwStepPar, wf, start);
                    }
                    else     /* Two positions skipped for meta categories */
                    {
                        for (int k = 0; k < cMidPosition; k++)
                        {
                            int index = i*cLastPosition + k; /* Index in the meta category to access when skipping to positions */
                            numSamplesAdded += processSingleCategoryLF1(&lwe, metaCategory1[index], valueCounter1[index], srcBkwStepPar, wf, start);
                        }
                    }
                }
                break;
            case 2: /* two meta categories (first two coordinates are additive inverses) */
                for (int i=0; i<cLastPosition; i++)
                {
                    int j = additiveInverse(cLastPosition, i);
                    if (meta_skipped == 1)
                    {
                        numSamplesAdded +=  processAdjacentCategoriesLF1(&lwe, metaCategory1[i], valueCounter1[i], metaCategory2[j], valueCounter2[j], srcBkwStepPar, wf, start); /* note: does not matter if i == j or not */
                    }
                    else
                    {
                        for (int k = 0; k < cMidPosition; k++)
                        {
                            int l = additiveInverse(cMidPosition, k);
                            int index = i*cLastPosition + k;
                            int additiveInverseIndex = j*cLastPosition + l;
                            numSamplesAdded +=  processAdjacentCategoriesLF1(&lwe, metaCategory1[index], valueCounter1[index], metaCategory2[additiveInverseIndex], valueCounter2[additiveInverseIndex], srcBkwStepPar, wf, start); /* note: does not matter if i == j or not */
                        }
                    }
                }
                break;
            default:
                timeStamp(start);
                printf("*** transition_bkw_step_smooth_lms_meta: Unexpected number of categories\n");
            }
            break;
        case LF2:
            switch (numReadCategories)
            {
            case 1: /* single meta category (no corresponding meta category with first two coordinates having (differing) additive inverses) */
                for (int i=0; i<cLastPosition; i++)
                {
                    if (meta_skipped == 1)   /* One position skipped for meta categories */
                    {
                        numSamplesAdded +=  processSingleCategoryLF2(&lwe, metaCategory1[i], valueCounter1[i], srcBkwStepPar, wf, start);
                    }
                    else     /* Two positions skipped for meta categories */
                    {
                        for (int k = 0; k < cMidPosition; k++)
                        {
                            int index = i*cLastPosition + k; /* Index in the meta category to access when skipping to positions */
                            numSamplesAdded += processSingleCategoryLF2(&lwe, metaCategory1[index], valueCounter1[index], srcBkwStepPar, wf, start);
                        }
                    }
                }
                break;
            case 2:  /* two meta categories (first two coordinates are additive inverses) */
                for (int i=0; i<cLastPosition; i++)
                {
                    int j = additiveInverse(cLastPosition, i);
                    if (meta_skipped == 1)
                    {
                        numSamplesAdded += processAdjacentCategoriesLF2(&lwe, metaCategory1[i], valueCounter1[i], metaCategory2[j], valueCounter2[j], srcBkwStepPar, wf, start); /* note: does not matter if i == j or not */
                    }
                    else
                    {
                        for (int k = 0; k < cMidPosition; k++)
                        {
                            int l = additiveInverse(cMidPosition, k);
                            int index = i*cMidPosition + k;
                            int additiveInverseIndex = j*cMidPosition + l;
                            numSamplesAdded += processAdjacentCategoriesLF2(&lwe, metaCategory1[index], valueCounter1[index], metaCategory2[additiveInverseIndex], valueCounter2[additiveInverseIndex], srcBkwStepPar, wf, start); /* note: does not matter if i == j or not */
                        }
                    }
                }
                break;
            default:
                timeStamp(start);
                printf("*** transition_bkw_step_final_smooth_lms_meta: Unexpected number of categories\n");
            }
            break;
        default:
            ASSERT_ALWAYS("Unsupported selection parameter");
        }

        cat += numReadCategories;

        /* free intermediate meta category information */
        if (metaCategory1)
        {
            for (int i=0; i<metaCategorySize; i++)
            {
                FREE(metaCategory1[i]);
            }
            FREE(metaCategory1);
        }
        if (metaCategory2)
        {
            for (int i=0; i<metaCategorySize; i++)
            {
                FREE(metaCategory2[i]);
            }
            FREE(metaCategory2);
        }
        if (valueCounter1)
        {
            FREE(valueCounter1);
        }
        if (valueCounter2)
        {
            FREE(valueCounter2);
        }

        while (cat > nextPrintLimit)
        {
            char s1[256], s2[256], s3[256];
            nextPrintLimit *= 2;
            timeStamp(start);
            printf("transition_bkw_step_final_smooth_lms_meta: num src categories read so far / all %10s /%10s, num sample added %s \n", sprintf_u64_delim(s1, cat), sprintf_u64_delim(s2, srcNumCategories), sprintf_u64_delim(s3, numSamplesAdded));
        }

        if (numSamplesAdded >= abortSampleLimit)
        {
            break;
        }

        numReadCategories = storageReaderGetNextAdjacentCategoryPair(&sr, &buf1, &numSamplesInBuf1, &buf2, &numSamplesInBuf2);
    }

    char s1[256], s2[256], s3[256];
    timeStamp(start);
    printf("transition_bkw_step_final_smooth_lms_meta: num src categories read so far / all %10s /%10s, num sample added %s \n", sprintf_u64_delim(s1, cat), sprintf_u64_delim(s2, srcNumCategories), sprintf_u64_delim(s3, numSamplesAdded));

    *numSamplesStored = numSamplesAdded;

    /* close storage handlers */
    storageReaderFree(&sr);
    freeSumAndDiffTables();

    fclose(wf);

    return 0;
}
