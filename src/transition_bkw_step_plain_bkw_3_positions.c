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

#include "transition_bkw_step_plain_bkw_3_positions.h"
#include "lookup_tables.h"
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

static u64 subtractSamples(lweInstance *lwe, lweSample *sample1, lweSample *sample2, bkwStepParameters *dstBkwStepPar, storageWriter *sw)
{

    int n = lwe->n;
    int startIndex = dstBkwStepPar->startIndex;
    int Ni_ = dstBkwStepPar->numPositions;

    if (dstBkwStepPar->sorting == smoothLMS && startIndex+Ni_ < lwe->n)
        Ni_++;
    short p01[Ni_];
    /* compute category index of new sample (without computing entire new sample) */
    for (int i=0; i<Ni_; i++)
        p01[i] = diffTable(columnValue(sample1, startIndex + i), columnValue(sample2, startIndex + i));
    u64 categoryIndex = position_values_2_category_index_from_partial_sample(lwe, p01, dstBkwStepPar);

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
        printf("operation error in subtractSamples!\n");
        return 1;
    }
    ASSERT(newSample, "No sample area returned from storage writer");

    /* compute new sample (subtract), write to reserved sample memory area */
    for (int i=0; i<n; i++)
    {
        newSample->col.a[i] = diffTable(columnValue(sample1, i), columnValue(sample2, i));
    }
    newSample->col.hash = bkwColumnComputeHash(newSample, n, 0 /* startRow */);
    int err1 = error(sample1);
    int err2 = error(sample2);
    newSample->error = (err1 == -1 || err2 == -1) ? -1 : diffTable(err1, err2); /* undefined if either parent error term is undefined */
    newSample->sumWithError = diffTable(sumWithError(sample1), sumWithError(sample2));

    /* discard zero columns (assuming that these are produced by coincidental cancellation due to sample amplification) */
    if (columnIsZero(newSample, n))
    {
        storageWriterUndoAddSample(sw, categoryIndex); /* return memory area to storage writer */
        return 1; /* sample processed but not added */
    }

    return 1; /* one sample processed (and actually added) */
}

static u64 addSamples(lweInstance *lwe, lweSample *sample1, lweSample *sample2, bkwStepParameters *dstBkwStepPar, storageWriter *sw)
{

    int n = lwe->n;
    int startIndex = dstBkwStepPar->startIndex;
    int Ni_ = dstBkwStepPar->numPositions;

    if (dstBkwStepPar->sorting == smoothLMS && startIndex+Ni_ < lwe->n)
        Ni_++;
    short p01[Ni_];
    /* compute category index of new sample (without computing entire new sample) */
    for (int i=0; i<Ni_; i++)
    {
        // printf("INDEX %d\n", startIndex + i);
        p01[i] = sumTable(columnValue(sample1, startIndex + i), columnValue(sample2, startIndex + i));
    }
    u64 categoryIndex = position_values_2_category_index_from_partial_sample(lwe, p01, dstBkwStepPar);

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
        printf("operation error in subtractSamples!\n");
        return 1;
    }
    ASSERT(newSample, "No sample area returned from storage writer");

    /* compute new sample (add), write to reserved sample memory area */
    for (int i=0; i<n; i++)
    {
        newSample->col.a[i] = sumTable(columnValue(sample1, i), columnValue(sample2, i));
    }
    newSample->col.hash = bkwColumnComputeHash(newSample, n, 0 /* startRow */);
    int err1 = error(sample1);
    int err2 = error(sample2);
    newSample->error = (err1 == -1 || err2 == -1) ? -1 : sumTable(err1, err2); /* undefined if either parent error term is undefined */
    newSample->sumWithError = sumTable(sumWithError(sample1), sumWithError(sample2));

    /* discard zero columns (assuming that these are produced by coincidental cancellation due to sample amplification) */
    if (columnIsZero(newSample, n))
    {
        storageWriterUndoAddSample(sw, categoryIndex); /* return memory area to storage writer */
        return 1; /* sample processed but not added */
    }

    return 1; /* one sample processed (and actually added) */
}

static inline int additiveInverse(int q, int val)
{
    return val == 0 ? 0 : q - val;
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

static u64 processSingleCategoryLF1(lweInstance *lwe, lweSample **categorySamplePointers, int numSamplesInCategory, bkwStepParameters *dstBkwStepPar, storageWriter *sw, time_t start)
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
        numProcessed += subtractSamples(lwe, firstSample, sample, dstBkwStepPar, sw);
    }
    flushStorageWriterIfCloseToFull(sw, start);
    return numProcessed;
}

static u64 processAdjacentCategoriesLF1(lweInstance *lwe, lweSample **categorySamplePointers1, int numSamplesInCategory1, lweSample **categorySamplePointers2, int numSamplesInCategory2, bkwStepParameters *dstBkwStepPar, storageWriter *sw, time_t start)
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
            numProcessed += subtractSamples(lwe, firstSample, sample, dstBkwStepPar, sw);
        }
        /* add all samples in adjacent category to first (same as above) sample (linear) */
        for (int i=0; i<numSamplesInCategory2; i++)
        {
            sample = categorySamplePointers2[i];
            numProcessed += addSamples(lwe, firstSample, sample, dstBkwStepPar, sw);
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
                numProcessed += subtractSamples(lwe, firstSample, sample, dstBkwStepPar, sw);
            }
        }
    }
    flushStorageWriterIfCloseToFull(sw, start);
    return numProcessed;
}

static u64 processSingleCategoryLF2(lweInstance *lwe, lweSample **categorySamplePointers, int numSamplesInCategory, bkwStepParameters *dstBkwStepPar, storageWriter *sw, u64 maxNumSamplesPerCategory, int flush, time_t start)
{
    u64 numProcessed = 0;
    lweSample *sample1;
    lweSample *sample2;

    /* subtract all pairs of samples (quadratic) */
    for (int i=0; i<numSamplesInCategory; i++)
    {
        sample1 = categorySamplePointers[i];
        for (int j=i+1; j<numSamplesInCategory; j++)
        {
            sample2 = categorySamplePointers[j];
            numProcessed += subtractSamples(lwe, sample1, sample2, dstBkwStepPar, sw);
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
    numProcessed += processSingleCategoryLF2(lwe, categorySamplePointers1, numSamplesInCategory1, dstBkwStepPar, sw, maxNumSamplesPerCategory, 0, start);
    if (numProcessed >= maxNumSamplesPerCategory)
    {
        flushStorageWriterIfCloseToFull(sw, start);
        return numProcessed;
    }

    /* process all pairs in category 2 (subtract sample pairs) */
    numProcessed += processSingleCategoryLF2(lwe, categorySamplePointers2, numSamplesInCategory2, dstBkwStepPar, sw, maxNumSamplesPerCategory - numProcessed, 0, start);
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
            numProcessed += addSamples(lwe, sample1, sample2, dstBkwStepPar, sw);
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

static lweSample ***sortMetacategorySamplesOnThirdCoordinate(lweSample *category, int numSamplesInCategory, lweInstance *lwe, bkwStepParameters *srcBkwStepPar, int *valueCounter)
{
    int q = lwe->q;
    int pos3 = srcBkwStepPar->startIndex + 2;
    int maxNumRepetitionsOfPos3 = 0;

    /* compute distribution of third coordinate value */
    MEMSET(valueCounter, 0, q * sizeof(int));
    for (int i=0; i<numSamplesInCategory; i++)
    {
        valueCounter[columnValue(&category[i], pos3)]++;
    }
    /* find maximum */
    for (int i=0; i<q; i++)
    {
        if (valueCounter[i] > maxNumRepetitionsOfPos3)
        {
            maxNumRepetitionsOfPos3 = valueCounter[i];
        }
    }

    /* allocate storage structure for (properly) sorted samples (no copying, pointers only) */
    lweSample ***metaCategory = (lweSample***)MALLOC(q * sizeof(lweSample**));
    for (int i=0; i<q; i++)
    {
        metaCategory[i] = (lweSample**)MALLOC(maxNumRepetitionsOfPos3 * sizeof(lweSample*));
    }

    /* sort */
    MEMSET(valueCounter, 0, q * sizeof(int));
    for (int i=0; i<numSamplesInCategory; i++)
    {
        lweSample *thisSample = &category[i];
        int val = columnValue(thisSample, pos3); /* get value of third coordinate from sample */
        metaCategory[val][valueCounter[val]++] = thisSample;
    }

    return metaCategory;
}

int transition_bkw_step_plain_bkw_3_positions(const char *srcFolderName, const char *dstFolderName, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, u64 *numSamplesStored, time_t start)
{

    /* get lwe parameters from file */
    lweInstance lwe;
    lweParametersFromFile(&lwe, srcFolderName);

    /* get sample info from file */
    u64 srcNumCategories, srcCategoryCapacity, srcNumTotalSamples;
    if (sampleInfoFromFile(srcFolderName, srcBkwStepPar, &srcNumCategories, &srcCategoryCapacity, &srcNumTotalSamples, NULL))
    {
        return 1; /* error reading from samples info file */
    }

    char nc[256], cc[256], ns[256];
    timeStamp(start);
    printf("transition_bkw_step_plain_bkw_3_positions: num src categories is %s, category capacity is %s, total num src samples is %s (%5.2f%% full)\n", sprintf_u64_delim(nc, srcNumCategories), sprintf_u64_delim(cc, srcCategoryCapacity), sprintf_u64_delim(ns, srcNumTotalSamples), 100*srcNumTotalSamples/(double)(srcNumCategories * srcCategoryCapacity));
    if (srcBkwStepPar->sorting != plainBKW || srcBkwStepPar->numPositions != 3)
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

    /* initialize storage writer */
    u64 dstNumCategories = num_categories(&lwe, dstBkwStepPar);
    u64 minDestinationStorageCapacityInSamples = round((double)(MAX_NUM_SAMPLES*4)/3);
    u64 dstCategoryCapacity = (minDestinationStorageCapacityInSamples + dstNumCategories - 1) / dstNumCategories;
    timeStamp(start);
    printf("transition_bkw_step_plain_bkw_3_positions: num dst categories is %s, category capacity is %s\n", sprintf_u64_delim(nc, dstNumCategories), sprintf_u64_delim(cc, dstCategoryCapacity));


    /* initialize storage writer */
    storageWriter sw;
    if (storageWriterInitialize(&sw, dstFolderName, &lwe, dstBkwStepPar, dstCategoryCapacity /* categoryCapacityFile same as in src file */))
    {
        ASSERT_ALWAYS("could not initialize storage writer");
        return 4; /* could not initialize storage writer */
    }

    /* initialize add and diff tables for faster operation */
    if (createSumAndDiffTables(lwe.q))
    {
        return 5; /* could not create addition and difference tables */
    }

    /* process samples */
    u64 maxNumSamplesPerCategory = dstCategoryCapacity * EARLY_ABORT_LOAD_LIMIT_PERCENTAGE / SAMPLE_DEPENDENCY_SMEARING + 1;
    u64 cat = 0; /* current category index */
    u64 nextPrintLimit = 2;
    lweSample *buf1;
    lweSample *buf2;
    u64 numSamplesInBuf1, numSamplesInBuf2;
    int numReadCategories = storageReaderGetNextAdjacentCategoryPair(&sr, &buf1, &numSamplesInBuf1, &buf2, &numSamplesInBuf2);
    while (numReadCategories && (storageWriterCurrentLoadPercentage(&sw) < EARLY_ABORT_LOAD_LIMIT_PERCENTAGE))
    {

        int *valueCounter1 = NULL;
        int *valueCounter2 = NULL;
        lweSample ***metaCategory1 = NULL;
        lweSample ***metaCategory2 = NULL;

        /* pre-process metacategory 1 (if available) */
        if (buf1)
        {
            valueCounter1 = MALLOC(lwe.q * sizeof(int));
            metaCategory1 = sortMetacategorySamplesOnThirdCoordinate(buf1, numSamplesInBuf1, &lwe, srcBkwStepPar, valueCounter1);
        }

        /* pre-process metacategory 2 (if available) */
        if (buf2)
        {
            ASSERT(buf1, "buf2 should only be available if buf1 is");
            valueCounter2 = MALLOC(lwe.q * sizeof(int));
            metaCategory2 = sortMetacategorySamplesOnThirdCoordinate(buf2, numSamplesInBuf2, &lwe, srcBkwStepPar, valueCounter2);
        }

        switch (dstBkwStepPar->selection)
        {
        case LF1:
            switch (numReadCategories)
            {
            case 1: /* single meta category (no corresponding meta category with first two coordinates having (differing) additive inverses) */
                ASSERT(metaCategory1, "unexpected parameter");
                for (int i=0; i<lwe.q; i++)
                {
                    processSingleCategoryLF1(&lwe, metaCategory1[i], valueCounter1[i], dstBkwStepPar, &sw, start);
                }
                break;
            case 2: /* two meta categories (first two coordinates are additive inverses) */
                for (int i=0; i<lwe.q; i++)
                {
                    int j = additiveInverse(lwe.q, i);
                    processAdjacentCategoriesLF1(&lwe, metaCategory1[i], valueCounter1[i], metaCategory2[j], valueCounter2[j], dstBkwStepPar, &sw, start); /* note: does not matter if i == j or not */
                }
                break;
            default:
                timeStamp(start);
                printf("*** transition_bkw_step_plain_bkw_3_positions: Unexpected number of categories\n");
            }
            break;
        case LF2:
            switch (numReadCategories)
            {
            case 1: /* single meta category (no corresponding meta category with first two coordinates having (differing) additive inverses) */
                for (int i=0; i<lwe.q; i++)
                {
                    processSingleCategoryLF2(&lwe, metaCategory1[i], valueCounter1[i], dstBkwStepPar, &sw, maxNumSamplesPerCategory / lwe.q + 1, 1, start);
                }
                break;
            case 2:  /* two meta categories (first two coordinates are additive inverses) */
                for (int i=0; i<lwe.q; i++)
                {
                    int j = additiveInverse(lwe.q, i);
                    processAdjacentCategoriesLF2(&lwe, metaCategory1[i], valueCounter1[i], metaCategory2[j], valueCounter2[j], srcBkwStepPar, dstBkwStepPar, &sw, 2*maxNumSamplesPerCategory / lwe.q + 1, start);
                }
                break;
            default:
                timeStamp(start);
                printf("*** transition_bkw_step_plain_bkw_3_positions: Unexpected number of categories\n");
            }
            break;
        default:
            ASSERT_ALWAYS("Unsupported selection parameter");
        }

        cat += numReadCategories;

        /* free intermediate meta category information */
        if (metaCategory1)
        {
            for (int i=0; i<lwe.q; i++)
            {
                FREE(metaCategory1[i]);
            }
            FREE(metaCategory1);
        }
        if (metaCategory2)
        {
            for (int i=0; i<lwe.q; i++)
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
            printf("transition_bkw_step_plain_bkw_3_positions: num src categories read so far / all %10s /%10s, dst storage load %5.2f%% (%s samples)\n", sprintf_u64_delim(s1, cat), sprintf_u64_delim(s2, srcNumCategories), storageWriterCurrentLoadPercentage(&sw), sprintf_u64_delim(s3, sw.totalNumSamplesAddedToStorageWriter));
        }
        numReadCategories = storageReaderGetNextAdjacentCategoryPair(&sr, &buf1, &numSamplesInBuf1, &buf2, &numSamplesInBuf2);
    }
    char s1[256], s2[256], s3[256];
    timeStamp(start);
    printf("transition_bkw_step_plain_bkw_3_positions: num src categories read so far / all %10s /%10s, dst storage load %5.2f%% (%s samples)\n", sprintf_u64_delim(s1, cat), sprintf_u64_delim(s2, srcNumCategories), storageWriterCurrentLoadPercentage(&sw), sprintf_u64_delim(s3, sw.totalNumSamplesAddedToStorageWriter));

    /* close storage handlers */
    storageReaderFree(&sr);
    *numSamplesStored = sw.totalNumSamplesAddedToStorageWriter;
    storageWriterFree(&sw); /* flushes automatically */

    freeSumAndDiffTables();
    return 0;
}
