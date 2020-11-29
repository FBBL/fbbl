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

#include "transition_bkw_step_final.h"
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

#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))

static u64 numZeroColumns;
static u64 numZeroColumnsAdd;

static u64 subtractSamples(lweInstance *lwe, lweSample *sample1, lweSample *sample2, bkwStepParameters *srcBkwStepPar, FILE *wf)
{

    int n = lwe->n;
    int q = lwe->q;

    /* perform Unnatural Selection */
    if (srcBkwStepPar->sorting == smoothLMS)
    {
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
            {
                return 1; /* Sample discarded */
            }
        }
    }

    lweSample *newSample = MALLOC(LWE_SAMPLE_SIZE_IN_BYTES);

    /* compute new sample (subtract), write to reserved sample memory area */
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
        printf("********* ERROR transition_bkw_step_final.c: allocation for new sample failed\n");
        exit(-1);
    }


    /* discard zero columns (assuming that these are produced by coincidental cancellation due to sample amplification) */
    if (columnIsZero(newSample, n))
    {
        numZeroColumns++;
        numZeroColumnsAdd++;
        free(newSample);
        return 1; /* sample processed but not added */
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

static int addSamples(lweInstance *lwe, lweSample *sample1, lweSample *sample2, bkwStepParameters *srcBkwStepPar, FILE *wf)
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
        {
            return 1; /* Sample processed but discarded */
        }
    }

    lweSample *newSample = MALLOC(LWE_SAMPLE_SIZE_IN_BYTES);

    /* compute new sample (add), write to reserved sample memory area */
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
        printf("********* ERROR transition_bkw_step_final.c: allocation for new sample failed\n");
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

static u64 processSingleCategoryLF1(lweInstance *lwe, lweSample *category, int numSamplesInCategory, bkwStepParameters *srcBkwStepPar, FILE *wf, time_t start)
{
    if (numSamplesInCategory < 2)
    {
        return 0;
    }
    u64 numAdded = 0;
    lweSample *firstSample = &category[0];
    for (int j=1; j<numSamplesInCategory; j++)
    {
        lweSample *thisSample = &category[j];
        numAdded += subtractSamples(lwe, firstSample, thisSample, srcBkwStepPar, wf);
    }
    return numAdded;
}

static u64 processSingleCategoryLF2(lweInstance *lwe, lweSample *category, int numSamplesInCategory, bkwStepParameters *srcBkwStepPar, FILE *wf, u64 maxNewSamples, time_t start)
{
    u64 numAdded = 0;
    for (int i=0; i<numSamplesInCategory; i++)
    {
        lweSample *sample1 = &category[i];
        for (int j=i+1; j<numSamplesInCategory; j++)
        {
            lweSample *sample2 = &category[j];
            numAdded += subtractSamples(lwe, sample1, sample2, srcBkwStepPar, wf);
            if (numAdded >= maxNewSamples)
            {
                return numAdded;
            }
        }
    }
    return numAdded;
}

static u64 processAdjacentCategoriesLF1(lweInstance *lwe, lweSample *category1, int numSamplesInCategory1, lweSample *category2, int numSamplesInCategory2, bkwStepParameters *srcBkwStepPar, FILE *wf, time_t start)
{
    u64 numAdded = 0;
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
            numAdded += subtractSamples(lwe, firstSample, sample, srcBkwStepPar, wf);
        }
        /* add all samples in adjacent category to first (same as above) sample (linear) */
        for (int i=0; i<numSamplesInCategory2; i++)
        {
            sample = &category2[i];
            numAdded += addSamples(lwe, firstSample, sample, srcBkwStepPar, wf);
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
                numAdded += subtractSamples(lwe, firstSample, sample, srcBkwStepPar, wf);
            }
        }
    }
    return numAdded;
}

static u64 processAdjacentCategoriesLF2(lweInstance *lwe, lweSample *category1, int numSamplesInCategory1, lweSample *category2, int numSamplesInCategory2, bkwStepParameters *srcBkwStepPar, FILE *wf, u64 maxNewSamples, time_t start)
{
    u64 numAdded = 0;

    /* Paul's note: sample dependency may be reduced beyond that given by SAMPLE_DEPENDENCY_SMEARING combining samples in a smarter order (all LF1 samples first, then...) */

    /* process all pairs in category 1 (subtract sample pairs) */
    numAdded += processSingleCategoryLF2(lwe, category1, numSamplesInCategory1, srcBkwStepPar, wf, maxNewSamples, start);
    if (numAdded >= maxNewSamples)
    {
        return numAdded;
    }

    /* process all pairs in category 2 (subtract sample pairs) */
    numAdded += processSingleCategoryLF2(lwe, category2, numSamplesInCategory2, srcBkwStepPar, wf, maxNewSamples - numAdded, start);
    if (numAdded >= maxNewSamples)
    {
        return numAdded;
    }

    /* process all pairs in categories 1 and 2 (add sample pairs) */
    for (int i=0; i<numSamplesInCategory1; i++)
    {
        for (int j=0; j<numSamplesInCategory2; j++)
        {
            numAdded += addSamples(lwe, &category1[i], &category2[j], srcBkwStepPar, wf);
            if (numAdded >= maxNewSamples)
            {
                return numAdded;
            }
        }
    }
    return numAdded;
}

/* perform a bkw step when reducing last b coefficients */
int transition_bkw_step_final(const char *srcFolderName, const char *dstFolderName, bkwStepParameters *srcBkwStepPar, u64 *numSamplesStored, time_t start)
{

    if (folderExists(dstFolderName))   /* if destination folder already exists, assume that we have performed this reduction step already */
    {
        return 100; /* reduction step already performed (destination folder already exists) */
    }

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
        lweDestroy(&lwe);
        return 1; /* error reading from samples info file */
    }

    u64 maxNewSamplesPerCategory = srcCategoryCapacity; // we set the same maximum as the source

    char nc[256], cc[256], ns[256];
    timeStamp(start);
    printf("transition_bkw_step_final: num src categories is %s, category capacity is %s, total num src samples is %s (%5.2f%% full)\n", sprintf_u64_delim(nc, srcNumCategories), sprintf_u64_delim(cc, srcCategoryCapacity), sprintf_u64_delim(ns, srcNumTotalSamples), 100*srcNumTotalSamples/(double)(srcNumCategories * srcCategoryCapacity));

    /* initialize storage reader */
    storageReader sr;
    if (storageReaderInitialize(&sr, srcFolderName))
    {
        lweDestroy(&lwe);
        ASSERT_ALWAYS("could not initialize storage reader");
        return 3; /* could not initialize storage reader */
    }

    /* Initialize destination folder and file */
    newStorageFolderWithGivenLweInstance(&lwe, dstFolderName);
    FILE *wf = fopenSamples(dstFolderName, "ab");
    if (!wf)
    {
        lweDestroy(&lwe);
        return -1;
    }

    /* initialize add and diff tables for faster operation */
    if (createSumAndDiffTables(lwe.q))
    {
        lweDestroy(&lwe);
        return 6; /* could not create addition and difference tables */
    }

    /* process samples */
    u64 cat = 0; /* current category index */
    u64 nextPrintLimit = 2;
    lweSample *buf1;
    lweSample *buf2;
    u64 numSamplesInBuf1, numSamplesInBuf2, numSamplesAdded = 0;
    int numReadCategories = storageReaderGetNextAdjacentCategoryPair(&sr, &buf1, &numSamplesInBuf1, &buf2, &numSamplesInBuf2);
    while (numReadCategories)
    {

        switch (srcBkwStepPar->selection)
        {
        case LF1:
            switch (numReadCategories)
            {
            case 1: /* single meta category (no corresponding meta category with first two coordinates having (differing) additive inverses) */
                numSamplesAdded += processSingleCategoryLF1(&lwe, buf1, numSamplesInBuf1, srcBkwStepPar, wf, start);
                break;
            case 2: /* two meta categories (first two coordinates are additive inverses) */
                numSamplesAdded += processAdjacentCategoriesLF1(&lwe, buf1, numSamplesInBuf1, buf2, numSamplesInBuf2, srcBkwStepPar, wf, start);
                break;
            default:
                timeStamp(start);
                printf("*** transition_bkw_step_final: Unexpected number of categories\n");
            }
            break;
        case LF2:
            switch (numReadCategories)
            {
            case 1: /* single meta category (no corresponding meta category with first two coordinates having (differing) additive inverses) */
                numSamplesAdded += processSingleCategoryLF2(&lwe, buf1, numSamplesInBuf1, srcBkwStepPar, wf, maxNewSamplesPerCategory, start);
                break;
            case 2:  /* two meta categories (first two coordinates are additive inverses) */
                numSamplesAdded += processAdjacentCategoriesLF2(&lwe, buf1, numSamplesInBuf1, buf2, numSamplesInBuf2, srcBkwStepPar, wf, 2*maxNewSamplesPerCategory, start);
                break;
            default:
                timeStamp(start);
                printf("*** transition_bkw_step_final: Unexpected number of categories\n");
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
            printf("transition_bkw_step_final: num src categories read so far / all %10s /%10s, num sample added %s \n", sprintf_u64_delim(s1, cat), sprintf_u64_delim(s2, srcNumCategories), sprintf_u64_delim(s3, numSamplesAdded));
        }
        numReadCategories = storageReaderGetNextAdjacentCategoryPair(&sr, &buf1, &numSamplesInBuf1, &buf2, &numSamplesInBuf2);
    }
    char s1[256], s2[256], s3[256];
    timeStamp(start);
    printf("transition_bkw_step_final: num src categories read so far / all %10s /%10s, num sample added %s \n", sprintf_u64_delim(s1, cat), sprintf_u64_delim(s2, srcNumCategories), sprintf_u64_delim(s3, numSamplesAdded));

    *numSamplesStored = numSamplesAdded;

    /* close storage handlers */
    storageReaderFree(&sr);
    fclose(wf);
    freeSumAndDiffTables();
    lweDestroy(&lwe);

    return 0;
}
