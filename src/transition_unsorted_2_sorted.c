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

#include "transition_unsorted_2_sorted.h"
#include "position_values_2_category_index.h"
#include "storage_file_utilities.h"
#include "memory_utils.h"
#include "log_utils.h"
#include "string_utils.h"
#include "storage_writer.h"
#include "bkw_step_parameters.h"
#include "verify_samples.h"
#include <inttypes.h>

int transition_unsorted_2_sorted(const char *srcFolderName, const char *dstFolderName, u64 minDestinationStorageCapacityInSamples, bkwStepParameters *bkwStepPar, time_t start)
{
    if (folderExists(dstFolderName))
    {
        return 1; /* output folder already exists */
    }

    lweInstance lwe;
    lweParametersFromFile(&lwe, srcFolderName); /* read lwe parameters from source folder */

    /* get number of samples in source file */
    u64 totNumUnsortedSamples = numSamplesInSampleFile(srcFolderName);
    if (!totNumUnsortedSamples)
    {
        return 2; /* no samples in source file */
    }
    char str[256];
    timeStamp(start);
    printf("src folder: %s (contains %s samples)\n", srcFolderName, sprintf_u64_delim(str, totNumUnsortedSamples));

    /* initiate storage writer (creates destination folder) */
    u64 numCategories = num_categories(&lwe, bkwStepPar);
    u64 categoryCapacityFile = (minDestinationStorageCapacityInSamples + numCategories - 1) / numCategories;
    storageWriter sw;
    int ret = storageWriterInitialize(&sw, dstFolderName, &lwe, bkwStepPar, categoryCapacityFile);
    if (ret)
    {
        return 100 + ret; /* could not initialize storage writer */
    }
    timeStamp(start);
    printf("dst folder: %s (has room for %s samples)\n", dstFolderName, sprintf_u64_delim(str, numCategories * categoryCapacityFile));

    /* open source sample file */
    FILE *f_src = fopenSamples(srcFolderName, "rb");
    if (!f_src)
    {
        if (storageWriterFree(&sw))
        {
            return 3; /* failed to free storage writer */
        }
        return 4; /* could not open samples file */
    }
//  printf("The size of the destination file before adding the samples is %" PRIu64 "\n", fileSize(sw.f));

    /* allocate sample read buffer */
    lweSample *sampleReadBuf = MALLOC(READ_BUFFER_CAPACITY_IN_SAMPLES * LWE_SAMPLE_SIZE_IN_BYTES);
    if (!sampleReadBuf)
    {
        fclose(f_src);
        if (storageWriterFree(&sw))
        {
            return 5; /* failed to free storage writer */
        }
        return 6; /* could not allocate sample read buffer */
    }
//  printf("Read buffer allocated (%d bytes / %d samples)\n", READ_BUFFER_CAPACITY_IN_SAMPLES * LWE_SAMPLE_SIZE_IN_BYTES, READ_BUFFER_CAPACITY_IN_SAMPLES);

    /* process all samples in source file */
    u64 nextPrintLimit = 100000000;
    while (!feof(f_src))
    {
        /* read chunk of samples from source sample file into read buffer */
        u64 numRead = freadSamples(f_src, sampleReadBuf, READ_BUFFER_CAPACITY_IN_SAMPLES);

        /* add samples to storage writer */
        for (u64 i=0; i<numRead; i++)
        {
            lweSample *s = &sampleReadBuf[i];
            u64 categoryIndex = position_values_2_category_index(&lwe, s, bkwStepPar);
            ASSERT(categoryIndex<numCategories, "ERROR *** invalid category index");

            u64 numIncorrectCategoryClassifications = 0;
            verifyOneSampleSorted(&lwe, s, bkwStepPar, NULL, NULL, categoryIndex, &numIncorrectCategoryClassifications, 0);
            if (numIncorrectCategoryClassifications)
            {
                printf("****** classification error detected\n");
            }
            int storageWriterStatus = 0;
            lweSample *d = storageWriterAddSample(&sw, categoryIndex, &storageWriterStatus); /* reserve memory area in storage writer */
            if (d)
            {
                MEMCPY(d, s, sizeof(lweSample)); /* copy sample to reserved memory area */
            }
            if (sw.totalNumSamplesProcessedByStorageWriter >= nextPrintLimit)
            {
                nextPrintLimit += 100000000;
                timeStamp(start);
                printf("%s samples processed by storage writer so far\n", sprintf_u64_delim(str, sw.totalNumSamplesProcessedByStorageWriter));
            }
            switch (storageWriterStatus)
            {
            case 0: /* insertion succeeded */
                verifyOneSampleSorted(&lwe, d, bkwStepPar, NULL, NULL, categoryIndex, &numIncorrectCategoryClassifications, 0);
                if (numIncorrectCategoryClassifications)
                {
                    printf("****** classification error detected (2)");
                }
                break;
            case 1: /* insertion succeeded but this was the last available slot in the cache */
                /* so it makes sense to flush the cache to file here, but only if the cache is actually smaller than the storage on file. */
                if (sw.categoryCapacityBuf < sw.categoryCapacityFile)   /* if buffer storage is smaller than file storage */
                {
                    double cacheLoad = storageWriterCurrentLoadPercentageCache(&sw);
                    if (cacheLoad >= MIN_STORAGE_WRITER_CACHE_LOAD_PERCENTAGE_BEFORE_FLUSH)   /* if flush threshold has been reached */
                    {
                        timeStamp(start);
                        printf("flushing storage writer at %6.02g%% load\n", cacheLoad);
                        storageWriterFlush(&sw);
                        timeStamp(start);
                        printf("flushing finished (file storage now at %6.02g%% load)\n", storageWriterCurrentLoadPercentageFile(&sw));
                    }
                }
                break;
            case 2: /* this sample was discarded, but could have been added if cache was flushed */
                break;
            case 3: /* this sample could not be added to storage writer, category on file + cache is full */
                break;
            }
        }
    }

    char s1[256], s2[256];
    timeStamp(start);
    printf("%s samples added to storage writer (%s samples discarded because some categories were overfull)\n", sprintf_u64_delim(s1, sw.totalNumSamplesAddedToStorageWriter), sprintf_u64_delim(s2, sw.totalNumSamplesProcessedByStorageWriter - sw.totalNumSamplesAddedToStorageWriter));

    /* cleanup */
    FREE(sampleReadBuf);
    fclose(f_src);
    ret = storageWriterFree(&sw); /* flushes automatically */
    if (ret)
    {
        return 200 + ret; /* could not free storage writer */
    }


//  printf("The size of the destination file after adding the samples is %" PRIu64 "\n", fileSize(sw.f));

    return 0;
}
