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

#include "transition_bkw_step_plain_bkw_2_positions.h"
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

#if 0 /* for testing, counting (over-)production of zero columns */
static u64 numZeroColumns;
static u64 numZeroColumnsAdd;
static u64 numZeroColumnsSub;
#endif

static u64 subtractSamples(lweInstance *lwe, lweSample *sample1, lweSample *sample2, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, storageWriter *sw)
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
    if (storageWriterStatus >= 2 || !newSample)   /* if no room, exit */
    {
        return 1; /* sample processed but not added */
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
        storageWriterUndoAddSample(sw, categoryIndex); /* return reserved memory area to storage writer */
//    numZeroColumns++;
//    numZeroColumnsSub++;
        return 1; /* sample processed but not added */
    }

    return 1; /* one sample processed (and actually added) */
}

static u64 addSamples(lweInstance *lwe, lweSample *sample1, lweSample *sample2, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, storageWriter *sw)
{
    int n = lwe->n;
    int startIndex = dstBkwStepPar->startIndex;
    int Ni_ = dstBkwStepPar->numPositions;

    if (dstBkwStepPar->sorting == smoothLMS && startIndex+Ni_ < lwe->n)
        Ni_++;
    short p01[Ni_];
    /* compute category index of new sample (without computing entire new sample) */
    for (int i=0; i<Ni_; i++)
        p01[i] = sumTable(columnValue(sample1, startIndex + i), columnValue(sample2, startIndex + i));
    u64 categoryIndex = position_values_2_category_index_from_partial_sample(lwe, p01, dstBkwStepPar);

    /* retrieve memory area for new sample in destination storage */
    int storageWriterStatus = 0;
    lweSample *newSample = storageWriterAddSample(sw, categoryIndex, &storageWriterStatus);
    if (storageWriterStatus >= 2 || !newSample)   /* if no room, exit */
    {
        return 1; /* sample processed but not added */
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
        storageWriterUndoAddSample(sw, categoryIndex); /* return reserved memory area to storage writer */
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
        for (int j=i+1; j<numSamplesInCategory; j++)
        {
            numProcessed += subtractSamples(lwe, &category[i], &category[j], srcBkwStepPar, dstBkwStepPar, sw);
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

int transition_bkw_step_plain_bkw_2_positions(const char *srcFolderName, const char *dstFolderName, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, u64 *numSamplesStored, time_t start)
{

#if 0 /* for testing only */
    numZeroColumns = 0;
    numZeroColumnsAdd = 0;
    numZeroColumnsSub = 0;
#endif

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
    printf("transition_bkw_step: num src categories is %s, category capacity is %s, total num src samples is %s (%5.2f%% full)\n", sprintf_u64_delim(nc, srcNumCategories), sprintf_u64_delim(cc, srcCategoryCapacity), sprintf_u64_delim(ns, srcNumTotalSamples), 100*srcNumTotalSamples/(double)(srcNumCategories * srcCategoryCapacity));
    if (srcBkwStepPar->sorting != plainBKW && srcBkwStepPar->sorting != LMS)
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
    printf("transition_bkw_step_plain_bkw_2_positions: num dst categories is %s, category capacity is %s\n", sprintf_u64_delim(nc, dstNumCategories), sprintf_u64_delim(cc, dstCategoryCapacity));

    /* initialize storage writer */
    storageWriter sw;
    if (storageWriterInitialize(&sw, dstFolderName, &lwe, dstBkwStepPar, dstCategoryCapacity))
    {
        ASSERT_ALWAYS("could not initialize storage writer");
        return 4; /* could not initialize storage writer */
    }

    /* initialize add and diff tables for faster operation */
    if (createSumAndDiffTables(lwe.q))
    {
        return 6; /* could not create addition and difference tables */
    }

#if 0
    /* for testing: statistics containers */
    storageWriterReturnCodes = MALLOC(4 * sizeof(u64));
    numAddedPerCategory = MALLOC(sw.numCategories * sizeof(u64));
    MEMSET(storageWriterReturnCodes, 0, 4 * sizeof(u64));
    MEMSET(numAddedPerCategory, 0, sw.numCategories * sizeof(u64));
#endif

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
                printf("*** transition_bkw_step_plain_bkw_3_positions: Unexpected number of categories\n");
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
                printf("*** transition_bkw_step_plain_bkw_3_positions: Unexpected number of categories\n");
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
            printf("transition_bkw_step: num src categories read so far / all %10s /%10s, dst storage load %5.2f%% (%s samples)\n", sprintf_u64_delim(s1, cat), sprintf_u64_delim(s2, srcNumCategories), storageWriterCurrentLoadPercentage(&sw), sprintf_u64_delim(s3, sw.totalNumSamplesAddedToStorageWriter));
        }
        numReadCategories = storageReaderGetNextAdjacentCategoryPair(&sr, &buf1, &numSamplesInBuf1, &buf2, &numSamplesInBuf2);
    }
    char s1[256], s2[256], s3[256];
    timeStamp(start);
    printf("transition_bkw_step: num src categories read so far / all %10s /%10s, dst storage load %5.2f%% (%s samples)\n", sprintf_u64_delim(s1, cat), sprintf_u64_delim(s2, srcNumCategories), storageWriterCurrentLoadPercentage(&sw), sprintf_u64_delim(s3, sw.totalNumSamplesAddedToStorageWriter));

    /* close storage handlers */
    storageReaderFree(&sr);
//  u64 totNumInserted = sw.totalNumSamplesAddedToStorageWriter;
//  u64 dstNumCategories = sw.numCategories;
//  u64 dstCategoryCapacityFile = sw.categoryCapacityFile;
    *numSamplesStored = sw.totalNumSamplesAddedToStorageWriter;
    storageWriterFree(&sw); /* flushes automatically */

#if 0 /* for testing purposes only, counting number of zero columns */
    u64 sum = 0;
    for (int i=0; i<dstNumCategories; i++)
    {
        sum += numAddedPerCategory[i] > dstCategoryCapacityFile ? dstCategoryCapacityFile : numAddedPerCategory[i];
    }
    char s1[256], s2[256], s3[256];
    printf("sum of all entries actually inserted should be %s\n", sprintf_u64_delim(s1, sum));
    printf("sum of all entries actually inserted is        %s (not counting zero columns, so could differ by %d)\n\n", sprintf_u64_delim(s1, totNumInserted), sw.categoryCapacityFile);

    printf("testing: return codes from storage writer\n");
    sum = 0;
    for (int i=0; i<4; i++)
    {
        printf("     %d : %12s\n", i, sprintf_u64_delim(s1, storageWriterReturnCodes[i]));
        sum += storageWriterReturnCodes[i];
    }
    printf("  zcol : %12s (%s in sub + %s in add)\n", sprintf_u64_delim(s1, numZeroColumns), sprintf_u64_delim(s2, numZeroColumnsSub), sprintf_u64_delim(s3, numZeroColumnsAdd));
    printf("-------------------------\n");
    printf("   sum : %12s\n", sprintf_u64_delim(s1, sum + numZeroColumns));

    FREE(storageWriterReturnCodes);
    FREE(numAddedPerCategory);
    storageWriterReturnCodes = NULL;
    numAddedPerCategory = NULL;
#endif

    freeSumAndDiffTables();

    return 0;
}
