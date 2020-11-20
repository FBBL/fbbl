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

#include "storage_file_utilities.h"
#include "memory_utils.h"
#include "lookup_tables.h"
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
            tmp_a = diffTable(columnValue(sample1, i), columnValue(sample2, i));
            tmp_a = MIN(tmp_a, q - tmp_a);
            a_norm_squared += tmp_a*tmp_a;
        }
        int numSelectionPositions = srcBkwStepPar->startIndex + srcBkwStepPar->numPositions - srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_start_index; /* Number of positions to apply unnatural selection on */
        double limit = numSelectionPositions * srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts*srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts; /* The threshold below which we accept samples */
        if (a_norm_squared >= limit)
            return 1; /* Sample processed but discarded */
    }

    int startIndex = dstBkwStepPar->startIndex;
    int numPositions = dstBkwStepPar->numPositions;
    short pn[MAX_SMOOTH_LMS_POSITIONS];

    /* compute category index of new sample (without computing entire new sample) */
    int Ni_ = (startIndex + numPositions) == lwe->n ? numPositions : numPositions+1; // differentiate last step
    for (int i=0; i<Ni_; i++)
    {
        pn[i] = diffTable(columnValue(sample1, startIndex + i), columnValue(sample2, startIndex + i));
    }

    u64 categoryIndex = position_values_2_category_index_from_partial_sample(lwe, pn, dstBkwStepPar);
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
            tmp_a = sumTable(columnValue(sample1, i), columnValue(sample2, i));
            tmp_a = MIN(tmp_a, q - tmp_a);
            a_norm_squared += tmp_a*tmp_a;
        }
        int numSelectionPositions = srcBkwStepPar->startIndex + srcBkwStepPar->numPositions - srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_start_index; /* Number of positions to apply unnatural selection on */
        double limit = numSelectionPositions * srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts*srcBkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts; /* The threshold below which we accept samples */
        if (a_norm_squared >= limit)
            return 1; /* Sample processed but discarded */
    }

    int startIndex = dstBkwStepPar->startIndex;
    int numPositions = dstBkwStepPar->numPositions;
    short pn[MAX_LMS_POSITIONS];

    /* compute category index of new sample (without computing entire new sample) */
    int Ni_ = (startIndex + numPositions) == lwe->n ? numPositions : numPositions+1; // differentiate last step
    for (int i=0; i<Ni_; i++)
    {
        pn[i] = sumTable(columnValue(sample1, startIndex + i), columnValue(sample2, startIndex + i));
    }

    u64 categoryIndex = position_values_2_category_index_from_partial_sample(lwe, pn, dstBkwStepPar);
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
//    numZeroColumns++;
//    numZeroColumnsAdd++;
        return 1; /* sample processed but not added */
    }
    return 1; /* one sample processed (and actually added) */
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

static u64 processSingleCategoryLF1(lweInstance *lwe, lweSample *category, int numSamplesInCategory, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, storageWriter *sw, time_t start)
{
    ASSERT(dstBkwStepPar->selection == LF1, "unexpected selection type");
    if (numSamplesInCategory < 2)
    {
        return 0;
    }
    u64 numProcessed = 0;
    lweSample *firstSample = &category[0];
    for (int j=1; j<numSamplesInCategory; j++)
    {
        lweSample *thisSample = &category[j];
        numProcessed += subtractSamples(lwe, firstSample, thisSample, srcBkwStepPar, dstBkwStepPar, sw);
    }
    flushStorageWriterIfCloseToFull(sw, start);
    return numProcessed;
}

static u64 processSingleCategoryLF2(lweInstance *lwe, lweSample *category, int numSamplesInCategory, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, storageWriter *sw, u64 maxNumSamplesPerCategory, time_t start)
{
    ASSERT(dstBkwStepPar->selection == LF2, "unexpected selection type");
    u64 numProcessed = 0;
    for (int i=0; i<numSamplesInCategory; i++)
    {
        lweSample *sample1 = &category[i];
        for (int j=i+1; j<numSamplesInCategory; j++)
        {
            lweSample *sample2 = &category[j];
            numProcessed += subtractSamples(lwe, sample1, sample2, srcBkwStepPar, dstBkwStepPar, sw);
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

static u64 processAdjacentCategoriesLF1(lweInstance *lwe, lweSample *category1, int numSamplesInCategory1, lweSample *category2, int numSamplesInCategory2, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, storageWriter *sw, time_t start)
{
    ASSERT(dstBkwStepPar->selection == LF1, "unexpected selection type");
    u64 numProcessed = 0;
    lweSample *firstSample;
    lweSample *sample;

    if (numSamplesInCategory1 > 0)
    {
        /* LF1-process category 1 */
        /* subtract all samples from the first one (linear) */
        firstSample = &category1[0];
        for (int i=1; i<numSamplesInCategory1; i++)
        {
            sample = &category1[i];
            numProcessed += subtractSamples(lwe, firstSample, sample, srcBkwStepPar, dstBkwStepPar, sw);
        }
        /* add all samples in adjacent category to first (same as above) sample (linear) */
        for (int i=0; i<numSamplesInCategory2; i++)
        {
            sample = &category2[i];
            numProcessed += addSamples(lwe, firstSample, sample, srcBkwStepPar, dstBkwStepPar, sw);
        }
    }
    else     /* numSamplesInCategory1 == 0 */
    {
        if (numSamplesInCategory2 >= 2)
        {
            /* LF1-process samples in category 2 only */
            firstSample = &category2[0];
            for (int i=1; i<numSamplesInCategory2; i++)
            {
                sample = &category2[i];
                numProcessed += subtractSamples(lwe, firstSample, sample, srcBkwStepPar, dstBkwStepPar, sw);
            }
        }
    }

    flushStorageWriterIfCloseToFull(sw, start);
    return numProcessed;
}

static u64 processAdjacentCategoriesLF2(lweInstance *lwe, lweSample *category1, int numSamplesInCategory1, lweSample *category2, int numSamplesInCategory2, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, storageWriter *sw, u64 maxNumSamplesPerCategory, time_t start)
{
    ASSERT(dstBkwStepPar->selection == LF2, "unexpected selection type");
    u64 numProcessed = 0;

    /* Paul's note: sample dependency may be reduced beyond that given by SAMPLE_DEPENDENCY_SMEARING combining samples in a smarter order (all LF1 samples first, then...) */

    /* process all pairs in category 1 (subtract sample pairs) */
    numProcessed += processSingleCategoryLF2(lwe, category1, numSamplesInCategory1, srcBkwStepPar, dstBkwStepPar, sw, maxNumSamplesPerCategory, start);
    if (numProcessed >= maxNumSamplesPerCategory)
    {
        flushStorageWriterIfCloseToFull(sw, start);
        return numProcessed;
    }

    /* process all pairs in category 2 (subtract sample pairs) */
    numProcessed += processSingleCategoryLF2(lwe, category2, numSamplesInCategory2, srcBkwStepPar, dstBkwStepPar, sw, maxNumSamplesPerCategory - numProcessed, start);
    if (numProcessed >= maxNumSamplesPerCategory)
    {
        flushStorageWriterIfCloseToFull(sw, start);
        return numProcessed;
    }

    /* process all pairs in categories 1 and 2 (add sample pairs) */
    for (int i=0; i<numSamplesInCategory1; i++)
    {
        for (int j=0; j<numSamplesInCategory2; j++)
        {
            numProcessed += addSamples(lwe, &category1[i], &category2[j], srcBkwStepPar, dstBkwStepPar, sw);
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

int transition_bkw_step_smooth_lms(const char *srcFolderName, const char *dstFolderName, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, u64 *numSamplesStored, time_t start)
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

    char nc[256], cc[256], ns[256], mc[256];
    timeStamp(start);
    printf("transition_bkw_step_smooth_lms: num src categories is %s, category capacity is %s, total num src samples is %s (%5.2f%% full)\n", sprintf_u64_delim(nc, srcNumCategories), sprintf_u64_delim(cc, srcCategoryCapacity), sprintf_u64_delim(ns, srcNumTotalSamples), 100*srcNumTotalSamples/(double)(srcNumCategories * srcCategoryCapacity));
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

    /* initialize storage writer */
    u64 dstNumCategories = num_categories(&lwe, dstBkwStepPar);
    u64 minDestinationStorageCapacityInSamples = round(((double)MAX_NUM_SAMPLES*5)/4);
    u64 dstCategoryCapacity = (minDestinationStorageCapacityInSamples + dstNumCategories - 1) / dstNumCategories;
    timeStamp(start);
    printf("transition_bkw_step_smooth_lms: num dst categories is %s, category capacity is %s, minDestinationSamplesCapacity %s\n", sprintf_u64_delim(nc, dstNumCategories), sprintf_u64_delim(cc, dstCategoryCapacity), sprintf_u64_delim(mc, minDestinationStorageCapacityInSamples));

    storageWriter sw;
    if (storageWriterInitialize(&sw, dstFolderName, &lwe, dstBkwStepPar, dstCategoryCapacity))
    {
        ASSERT_ALWAYS("could not initialize storage writer");
        return 4; /* could not initialize storage writer */
    }

    /* initialize add and diff tables for faster operation */
    /* TODO: move to initialization */
    if (createSumAndDiffTables(lwe.q))
    {
        return 6; /* could not create addition and difference tables */
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

        switch (dstBkwStepPar->selection)
        {
        case LF1:
            switch (numReadCategories)
            {
            case 1: /* single meta category (no corresponding meta category with first two coordinates having (differing) additive inverses) */
                processSingleCategoryLF1(&lwe, buf1, numSamplesInBuf1, srcBkwStepPar, dstBkwStepPar, &sw, start);
                break;
            case 2: /* two meta categories (first two coordinates are additive inverses) */
                processAdjacentCategoriesLF1(&lwe, buf1, numSamplesInBuf1, buf2, numSamplesInBuf2, srcBkwStepPar, dstBkwStepPar, &sw, start);
                break;
            default:
                timeStamp(start);
                printf("*** transition_bkw_step_smooth_lms: Unexpected number of categories\n");
            }
            break;
        case LF2:
            switch (numReadCategories)
            {
            case 1: /* single meta category (no corresponding meta category with first two coordinates having (differing) additive inverses) */
                processSingleCategoryLF2(&lwe, buf1, numSamplesInBuf1, srcBkwStepPar, dstBkwStepPar, &sw, maxNumSamplesPerCategory, start);
                break;
            case 2:  /* two meta categories (first two coordinates are additive inverses) */
                processAdjacentCategoriesLF2(&lwe, buf1, numSamplesInBuf1, buf2, numSamplesInBuf2, srcBkwStepPar, dstBkwStepPar, &sw, 2*maxNumSamplesPerCategory, start);
                break;
            default:
                timeStamp(start);
                printf("*** transition_bkw_step_smooth_lms: Unexpected number of categories\n");
            }
            break;
        default:
            ASSERT_ALWAYS("Unsupported selection parameter");
        }

        cat += numReadCategories;
        while (cat > nextPrintLimit)
        {
            char s1[256], s2[256], s3[256];
            nextPrintLimit *= 2;
            timeStamp(start);
            printf("transition_bkw_step_lms: num src categories read so far / all %10s /%10s, dst storage load %5.2f%% (%s samples)\n", sprintf_u64_delim(s1, cat), sprintf_u64_delim(s2, srcNumCategories), storageWriterCurrentLoadPercentage(&sw), sprintf_u64_delim(s3, sw.totalNumSamplesAddedToStorageWriter));
        }
        numReadCategories = storageReaderGetNextAdjacentCategoryPair(&sr, &buf1, &numSamplesInBuf1, &buf2, &numSamplesInBuf2);
    }
    char s1[256], s2[256], s3[256];
    timeStamp(start);
    printf("transition_bkw_step_lms: num src categories read so far / all %10s /%10s, dst storage load %5.2f%% (%s samples)\n", sprintf_u64_delim(s1, cat), sprintf_u64_delim(s2, srcNumCategories), storageWriterCurrentLoadPercentage(&sw), sprintf_u64_delim(s3, sw.totalNumSamplesAddedToStorageWriter));

    /* close storage handlers */
    storageReaderFree(&sr);
//  u64 totNumInserted = sw.totalNumSamplesAddedToStorageWriter;
//  u64 dstNumCategories = sw.numCategories;
//  u64 dstCategoryCapacityFile = sw.categoryCapacityFile;
    *numSamplesStored = sw.totalNumSamplesAddedToStorageWriter;
    storageWriterFree(&sw); /* flushes automatically */
    freeSumAndDiffTables();

    return 0;
}

