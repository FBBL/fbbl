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

#include "transition_bkw_step_smooth_lms_meta.h"
#include "storage_file_utilities.h"
#include "memory_utils.h"
#include "log_utils.h"
#include "string_utils.h"
#include "lwe_sorting.h"
#include "storage_reader.h"
#include "storage_writer.h"
#include "position_values_2_category_index.h"
#include "config_bkw.h"
#include <inttypes.h>
#include <math.h>

#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))

static u64 subtractSamples(lweInstance *lwe, lweSample *sample1, lweSample *sample2, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, storageWriter *sw)
{
    int n = lwe->n;
    int q = lwe->q;

    /* perform Unnatural Selection */
    if(srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts)
    {
        double a_norm_squared = 0;
        short tmp_a;
        for (int i=srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_start_index; i<srcBkwStepPar->startIndex + srcBkwStepPar->numPositions; i++)
        {
            tmp_a = (columnValue(sample1, i) - columnValue(sample2, i) + q) % q;
            tmp_a = MIN(tmp_a, q - tmp_a);
            a_norm_squared += tmp_a*tmp_a;
        }
        int numSelectionPositions = srcBkwStepPar->startIndex + srcBkwStepPar->numPositions - srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_start_index; /* Number of positions to apply unnatural selection on */
        double limit = numSelectionPositions * srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts*srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts; /* The threshold below which we accept samples */
        if (a_norm_squared >= limit)
        {
            return 1; /* Sample discarded */
        }
    }

    /* Check that the +1 position behaves correctly! */
    short tmp_a = (columnValue(sample1, srcBkwStepPar->startIndex + srcBkwStepPar->numPositions) - columnValue(sample2, srcBkwStepPar->startIndex + srcBkwStepPar->numPositions) + q) % q;
    if (tmp_a >= srcBkwStepPar->sortingPar.smoothLMS.p1 && tmp_a <= lwe->q - srcBkwStepPar->sortingPar.smoothLMS.p1)
    {
        printf("tmp_a = %d \n p_1 = %d \n", tmp_a, srcBkwStepPar->sortingPar.smoothLMS.p1);
        exit(1);
    }

    // printf("We will not discard the sample\n");

    int startIndex = dstBkwStepPar->startIndex;
    int numPositions = dstBkwStepPar->numPositions;
    short pn[MAX_SMOOTH_LMS_POSITIONS];

    /* compute category index of new sample (without computing entire new sample) */
    int Ni_ = (startIndex + numPositions) == lwe->n ? numPositions : numPositions+1; // differentiate last step
    for (int i=0; i<Ni_; i++)
    {
        pn[i] = (columnValue(sample1, startIndex + i) - columnValue(sample2, startIndex + i) + q) % q;
//    printf("sample1[%d] = %d\n", i, columnValue(sample1, startIndex + i));
//    printf("sample2[%d] = %d\n", i, columnValue(sample2, startIndex + i));
//    printf("pn     [%d] = %d\n", i, pn[i]);
    }

    u64 categoryIndex = position_values_2_category_index_smooth_lms_meta(lwe, dstBkwStepPar, pn);
    // printf("Difference is in category %" PRIu64 " \n", categoryIndex);
    ASSERT(categoryIndex >= 0, "ERROR: invalid category");

    /* retrieve memory area for new sample in destination storage */
    int storageWriterStatus = 0;
    lweSample *newSample = storageWriterAddSample(sw, categoryIndex, &storageWriterStatus);

    // printf("Do we get here?\n");

    /* if no room, exit */
    if (storageWriterStatus >= 2)
    {
        return 1; /* sample processed but not added */
    }
    if (!newSample)
    {
        printf("operation error in subtractSamples!\n");
        return 1;
    }
    ASSERT(newSample, "No sample area returned from storage writer");

    /* compute new sample (subtract), write to reserved sample memory area */
    for (int i=0; i<n; i++)
    {
        newSample->col.a[i] = (columnValue(sample1, i) - columnValue(sample2, i) + q) % q;
    }

    newSample->col.hash = bkwColumnComputeHash(newSample, n, 0 /* startRow */);
    int err1 = error(sample1);
    int err2 = error(sample2);
    newSample->error = (err1 == -1 || err2 == -1) ? -1 : (err1 - err2 + q) % q; /* undefined if either parent error term is undefined */
    newSample->sumWithError = (sumWithError(sample1) - sumWithError(sample2) + q) % q;

    /* discard zero columns (assuming that these are produced by coincidental cancellation due to sample amplification) */
    if (columnIsZero(newSample, n))
    {
        storageWriterUndoAddSample(sw, categoryIndex); /* return memory area to storage writer */
//    numZeroColumns++;
//    numZeroColumnsSub++;
        return 1; /* sample processed but not added */
    }

    return 1; /* one sample processed (and actually added) */
}

static u64 addSamples(lweInstance *lwe, lweSample *sample1, lweSample *sample2, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, storageWriter *sw)
{
    int n = lwe->n;
    int q = lwe->q;

    /* perform Unnatural Selection */
    if(srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts)
    {
        double a_norm_squared = 0;
        short tmp_a;
        for (int i=srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_start_index; i<srcBkwStepPar->startIndex + srcBkwStepPar->numPositions; i++)
        {
            tmp_a = (columnValue(sample1, i) + columnValue(sample2, i)) % q;
            tmp_a = MIN(tmp_a, q - tmp_a);
            a_norm_squared += tmp_a*tmp_a;
        }
        int numSelectionPositions = srcBkwStepPar->startIndex + srcBkwStepPar->numPositions - srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_start_index; /* Number of positions to apply unnatural selection on */
        double limit = numSelectionPositions * srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts*srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts; /* The threshold below which we accept samples */
        if (a_norm_squared >= limit)
        {
            return 1; /* Sample processed but discarded */
        }
    }

    /* Check that the +1 position behaves correctly! */
    short tmp_a = (columnValue(sample1, srcBkwStepPar->startIndex + srcBkwStepPar->numPositions) + columnValue(sample2, srcBkwStepPar->startIndex + srcBkwStepPar->numPositions)) % q;
    if (tmp_a >= srcBkwStepPar->sortingPar.smoothLMS.p1 && tmp_a <= lwe->q - srcBkwStepPar->sortingPar.smoothLMS.p1)
    {
        printf("tmp_a = %d \n p_1 = %d \n", tmp_a, srcBkwStepPar->sortingPar.smoothLMS.p1);
        exit(1);
    }

    int startIndex = dstBkwStepPar->startIndex;
    int numPositions = dstBkwStepPar->numPositions;
    short pn[MAX_LMS_POSITIONS];

    /* compute category index of new sample (without computing entire new sample) */
    int Ni_ = (startIndex + numPositions) == lwe->n ? numPositions : numPositions+1; // differentiate last step
    for (int i=0; i<Ni_; i++)
    {
        pn[i] = (columnValue(sample1, startIndex + i) + columnValue(sample2, startIndex + i)) % q;
    }

    u64 categoryIndex = position_values_2_category_index_smooth_lms_meta(lwe, dstBkwStepPar, pn);
    ASSERT(categoryIndex >= 0, "ERROR: invalid category");

    /* retrieve memory area for new sample in destination storage */
    int storageWriterStatus = 0;
    lweSample *newSample = storageWriterAddSample(sw, categoryIndex, &storageWriterStatus);
    /* if no room, exit */
    if (storageWriterStatus >= 2)
    {
        return 1; /* sample processed but not added */
    }
    if (!newSample)
    {
        printf("operation error in addSamples!\n");
        return 1;
    }
    ASSERT(newSample, "No sample area returned from storage writer");

    /* compute new sample (add), write to reserved sample memory area */
    for (int i=0; i<n; i++)
    {
        newSample->col.a[i] = (columnValue(sample1, i) + columnValue(sample2, i)) % q;
    }
    newSample->col.hash = bkwColumnComputeHash(newSample, n, 0 /* startRow */);
    int err1 = error(sample1);
    int err2 = error(sample2);
    newSample->error = (err1 == -1 || err2 == -1) ? -1 : (err1 + err2) % q; /* undefined if either parent error term is undefined */
    newSample->sumWithError = (sumWithError(sample1) + sumWithError(sample2)) % q;

    /* discard zero columns (assuming that these are produced by coincidental cancellation due to sample amplification) */
    if (columnIsZero(newSample, n))
    {
        storageWriterUndoAddSample(sw, categoryIndex); /* return memory area to storage writer */
//    numZeroColumns++;
//    numZeroColumnsAdd++;
        return 1; /* sample processed but not added */
    }

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

static void flushStorageWriterIfCloseToFull(storageWriter *sw, time_t start)
{
    if (sw->categoryCapacityBuf < sw->categoryCapacityFile)
    {
        /* flush storage writer category if nearly full */
        double cacheLoad = storageWriterCurrentLoadPercentageCache(sw);
        if (cacheLoad >= MIN_STORAGE_WRITER_CACHE_LOAD_BEFORE_FLUSH)
        {
            timeStamp(start);
            printf("flushing storage writer cache at %6.02g%% load\n", cacheLoad);
            int ret = storageWriterFlush(sw);
            if (ret)
            {
                printf("*** Error: storageWriterFlush returned %d\n", ret);
            }
            timeStamp(start);
            printf("flushing finished (destination file storage now at %6.02g%% load)\n", storageWriterCurrentLoadPercentageFile(sw));
        }
    }
}

static u64 processSingleCategoryLF1(lweInstance *lwe, lweSample **categorySamplePointers, int numSamplesInCategory, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, storageWriter *sw, time_t start)
{
    u64 numProcessed = 0;
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
        numProcessed += subtractSamples(lwe, firstSample, sample, srcBkwStepPar, dstBkwStepPar, sw);
    }
    flushStorageWriterIfCloseToFull(sw, start);
    return numProcessed;
}

static u64 processAdjacentCategoriesLF1(lweInstance *lwe, lweSample **categorySamplePointers1, int numSamplesInCategory1, lweSample **categorySamplePointers2, int numSamplesInCategory2, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, storageWriter *sw, time_t start)
{
    u64 numProcessed = 0;
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
            numProcessed += subtractSamples(lwe, firstSample, sample, srcBkwStepPar, dstBkwStepPar, sw);
        }
        /* add all samples in adjacent category to first (same as above) sample (linear) */
        for (int i=0; i<numSamplesInCategory2; i++)
        {
            sample = categorySamplePointers2[i];
            numProcessed += addSamples(lwe, firstSample, sample, srcBkwStepPar, dstBkwStepPar, sw);
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
                numProcessed += subtractSamples(lwe, firstSample, sample, srcBkwStepPar, dstBkwStepPar, sw);
            }
        }
    }
    flushStorageWriterIfCloseToFull(sw, start);
    return numProcessed;
}

static u64 processSingleCategoryLF2(lweInstance *lwe, lweSample **categorySamplePointers, int numSamplesInCategory, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, storageWriter *sw, u64 maxNumSamplesPerCategory, int flush, time_t start)
{
    u64 numProcessed = 0;
    lweSample *sample1;
    lweSample *sample2;



    /* subtract all pairs of samples (quadratic) */
    for (int i=0; i<numSamplesInCategory; i++)
    {
        // printf("Hello there i = %d!\n", i);
        sample1 = categorySamplePointers[i];
        for (int j=i+1; j<numSamplesInCategory; j++)
        {
            // printf("Hello there j = %d!\n", j);
            sample2 = categorySamplePointers[j];
            // printf("Hello there j = %d!\n", j);
            numProcessed += subtractSamples(lwe, sample1, sample2, srcBkwStepPar, dstBkwStepPar, sw);
            if (numProcessed >= maxNumSamplesPerCategory)
            {
                if (flush)
                {
                    flushStorageWriterIfCloseToFull(sw, start);
                }
                return numProcessed;
            }
        }
    }
    if (flush)
    {
        flushStorageWriterIfCloseToFull(sw, start);
    }
    return numProcessed;
}

static u64 processAdjacentCategoriesLF2(lweInstance *lwe, lweSample **categorySamplePointers1, int numSamplesInCategory1, lweSample **categorySamplePointers2, int numSamplesInCategory2, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, storageWriter *sw, u64 maxNumSamplesPerCategory, time_t start)
{
    u64 numProcessed = 0;
    lweSample *sample1;
    lweSample *sample2;

    /* Paul's note: sample dependency may be reduced beyond that given by SAMPLE_DEPENDENCY_SMEARING combining samples in a smarter order (all LF1 samples first, then...) */

    /* process all pairs in category 1 (subtract sample pairs) */
    numProcessed += processSingleCategoryLF2(lwe, categorySamplePointers1, numSamplesInCategory1, srcBkwStepPar, dstBkwStepPar, sw, maxNumSamplesPerCategory, 0, start);
    if (numProcessed >= maxNumSamplesPerCategory)
    {
        flushStorageWriterIfCloseToFull(sw, start);
        return numProcessed;
    }

    /* process all pairs in category 2 (subtract sample pairs) */
    numProcessed += processSingleCategoryLF2(lwe, categorySamplePointers2, numSamplesInCategory2, srcBkwStepPar, dstBkwStepPar, sw, maxNumSamplesPerCategory - numProcessed, 0, start);
    if (numProcessed >= maxNumSamplesPerCategory)
    {
        flushStorageWriterIfCloseToFull(sw, start);
        return numProcessed;
    }

    /* process all pairs in categories 1 and 2 (add sample pairs) */
    for (int i=0; i<numSamplesInCategory1; i++)
    {
        sample1 = categorySamplePointers1[i];
        for (int j=0; j<numSamplesInCategory2; j++)
        {
            sample2 = categorySamplePointers2[j];
            numProcessed += addSamples(lwe, sample1, sample2, srcBkwStepPar, dstBkwStepPar, sw);
            if (numProcessed >= maxNumSamplesPerCategory)
            {
                flushStorageWriterIfCloseToFull(sw, start);
                return numProcessed;
            }
        }
    }

    flushStorageWriterIfCloseToFull(sw, start);

    return numProcessed;
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

int transition_bkw_step_smooth_lms_meta(const char *srcFolderName, const char *dstFolderName, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, u64 *numSamplesStored, time_t start)
{
    /* get lwe parameters from file */
    lweInstance lwe;
    lweParametersFromFile(&lwe, srcFolderName);

    if (folderExists(dstFolderName))   /* if destination folder already exists, assume that we have performed this reduction step already */
    {
        lweDestroy(&lwe);
        return 100; /* reduction step already performed (destination folder already exists) */
    }

    /* get sample info from file */
    u64 srcNumCategories, srcCategoryCapacity, srcNumTotalSamples;

    if (sampleInfoFromFile(srcFolderName, srcBkwStepPar, &srcNumCategories, &srcCategoryCapacity, &srcNumTotalSamples, NULL))
    {
        lweDestroy(&lwe);
        return 1; /* error reading from samples info file */
    }

    /* TODO Calculate the number of metacategories when using smooth LMS with metacategories */
    char nc[256], cc[256], ns[256];
    timeStamp(start);
    printf("transition_bkw_step_lms_meta: num src categories is %s, category capacity is %s, total num src samples is %s (%5.2f%% full)\n", sprintf_u64_delim(nc, srcNumCategories), sprintf_u64_delim(cc, srcCategoryCapacity), sprintf_u64_delim(ns, srcNumTotalSamples), 100*srcNumTotalSamples/(double)(srcNumCategories * srcCategoryCapacity));
    if (srcBkwStepPar->sorting != smoothLMS)
    {
        lweDestroy(&lwe);
        return 2; /* unexpected sample sorting at src folder */
    }

    /* initialize storage reader */
    storageReader sr;
    if (storageReaderInitialize(&sr, srcFolderName))
    {
        lweDestroy(&lwe);
        ASSERT_ALWAYS("could not initialize storage reader");
        return 3; /* could not initialize storage reader */
    }

    /* initialize storage writer */
    u64 dstNumCategories = num_categories(&lwe, dstBkwStepPar);
    u64 minDestinationStorageCapacityInSamples = round((double)(MAX_NUM_SAMPLES*4)/3);
    u64 dstCategoryCapacity = (minDestinationStorageCapacityInSamples + dstNumCategories - 1) / dstNumCategories;
    timeStamp(start);
    printf("transition_bkw_step_lms_meta: num dst categories is %s, category capacity is %s\n", sprintf_u64_delim(nc, dstNumCategories), sprintf_u64_delim(cc, dstCategoryCapacity));

    storageWriter sw;
    if (storageWriterInitialize(&sw, dstFolderName, &lwe, dstBkwStepPar, dstCategoryCapacity))
    {
        lweDestroy(&lwe);
        ASSERT_ALWAYS("could not initialize storage writer");
        return 4; /* could not initialize storage writer */
    }

    /* process samples */
    // u64 maxNumSamplesPerCategory = dstCategoryCapacity * EARLY_ABORT_LOAD_LIMIT_PERCENTAGE / SAMPLE_DEPENDENCY_SMEARING + 1;
    u64 maxNumSamplesPerCategory = 1000000000;
    u64 cat = 0; /* current category index */
    u64 nextPrintLimit = 2;
    lweSample *buf1;
    lweSample *buf2;
    u64 numSamplesInBuf1, numSamplesInBuf2;
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
        lweDestroy(&lwe);
        ASSERT_ALWAYS("Unsupported number of skipped positions");
        return 7;
    }

    /* The main while loop */
    while (numReadCategories && (storageWriterCurrentLoadPercentage(&sw) < EARLY_ABORT_LOAD_LIMIT_PERCENTAGE))
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
        switch (dstBkwStepPar->selection)
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
                        processSingleCategoryLF1(&lwe, metaCategory1[i], valueCounter1[i], srcBkwStepPar, dstBkwStepPar, &sw, start);
                    }
                    else     /* Two positions skipped for meta categories */
                    {
                        for (int k = 0; k < cMidPosition; k++)
                        {
                            int index = i*cLastPosition + k; /* Index in the meta category to access when skipping to positions */
                            processSingleCategoryLF1(&lwe, metaCategory1[index], valueCounter1[index], srcBkwStepPar, dstBkwStepPar, &sw, start);
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
                        processAdjacentCategoriesLF1(&lwe, metaCategory1[i], valueCounter1[i], metaCategory2[j], valueCounter2[j], srcBkwStepPar, dstBkwStepPar, &sw, start); /* note: does not matter if i == j or not */
                    }
                    else
                    {
                        for (int k = 0; k < cMidPosition; k++)
                        {
                            int l = additiveInverse(cMidPosition, k);
                            int index = i*cLastPosition + k;
                            int additiveInverseIndex = j*cLastPosition + l;
                            processAdjacentCategoriesLF1(&lwe, metaCategory1[index], valueCounter1[index], metaCategory2[additiveInverseIndex], valueCounter2[additiveInverseIndex], srcBkwStepPar, dstBkwStepPar, &sw, start); /* note: does not matter if i == j or not */
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
                // printf("Number of categories for last position: %d\n", cLastPosition);
                for (int i=0; i<cLastPosition; i++)
                {
                    // printf("Hello!\n");
                    if (meta_skipped == 1)   /* One position skipped for meta categories */
                    {
                        // printf("i = %d\n", i);
                        processSingleCategoryLF2(&lwe, metaCategory1[i], valueCounter1[i], srcBkwStepPar, dstBkwStepPar, &sw, maxNumSamplesPerCategory / metaCategorySize + 1, 1, start);
                    }
                    else     /* Two positions skipped for meta categories */
                    {
                        for (int k = 0; k < cMidPosition; k++)
                        {
                            int index = i*cLastPosition + k; /* Index in the meta category to access when skipping to positions */
                            processSingleCategoryLF2(&lwe, metaCategory1[index], valueCounter1[index], srcBkwStepPar, dstBkwStepPar, &sw, maxNumSamplesPerCategory / metaCategorySize + 1, 1, start);
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
                        processAdjacentCategoriesLF2(&lwe, metaCategory1[i], valueCounter1[i], metaCategory2[j], valueCounter2[j], srcBkwStepPar, dstBkwStepPar, &sw, maxNumSamplesPerCategory / metaCategorySize + 1, start); /* note: does not matter if i == j or not */
                    }
                    else
                    {
                        for (int k = 0; k < cMidPosition; k++)
                        {
                            int l = additiveInverse(cMidPosition, k);
                            int index = i*cMidPosition + k;
                            int additiveInverseIndex = j*cMidPosition + l;
                            processAdjacentCategoriesLF2(&lwe, metaCategory1[index], valueCounter1[index], metaCategory2[additiveInverseIndex], valueCounter2[additiveInverseIndex], srcBkwStepPar, dstBkwStepPar, &sw, maxNumSamplesPerCategory / metaCategorySize + 1, start); /* note: does not matter if i == j or not */
                        }
                    }
                }
                break;
            default:
                timeStamp(start);
                printf("*** transition_bkw_step_smooth_lms_meta: Unexpected number of categories\n");
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
            printf("transition_bkw_step_smooth_lms_meta: num src categories read so far / all %10s /%10s, dst storage load %5.2f%% (%s samples)\n", sprintf_u64_delim(s1, cat), sprintf_u64_delim(s2, srcNumCategories), storageWriterCurrentLoadPercentage(&sw), sprintf_u64_delim(s3, sw.totalNumSamplesAddedToStorageWriter));
        }
        numReadCategories = storageReaderGetNextAdjacentCategoryPair(&sr, &buf1, &numSamplesInBuf1, &buf2, &numSamplesInBuf2);
    }

    char s1[256], s2[256], s3[256];
    timeStamp(start);
    printf("transition_bkw_step_smooth_lms_meta: num src categories read so far / all %10s /%10s, dst storage load %5.2f%% (%s samples)\n", sprintf_u64_delim(s1, cat), sprintf_u64_delim(s2, srcNumCategories), storageWriterCurrentLoadPercentage(&sw), sprintf_u64_delim(s3, sw.totalNumSamplesAddedToStorageWriter));

    /* close storage handlers */
    storageReaderFree(&sr);
    *numSamplesStored = sw.totalNumSamplesAddedToStorageWriter;
    storageWriterFree(&sw); /* flushes automatically */

    lweDestroy(&lwe);

    return 0;
}
