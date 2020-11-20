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

#include "storage_writer.h"
#include "memory_utils.h"
#include "storage_file_utilities.h"
#if !defined(STORAGE_WRITER_CACHE_SIZE_IN_BYTES)
#include "physicalmemorysize.h"
#endif
#include <inttypes.h>

int storageWriterInitialize(storageWriter *sw, const char *dstFolderName, lweInstance *lwe, bkwStepParameters *bkwStepPar, u64 categoryCapacityFile)
{
    strncpy(sw->dstFolderName, dstFolderName, 512);
    sw->f = NULL; // handle to samples file
    sw->bkwStepPar = bkwStepPar;
    sw->numCategories = num_categories(lwe, sw->bkwStepPar); /* number of destination categories */
    sw->categoryCapacityBuf = 2 * categoryCapacityFile;
    sw->categoryCapacityFile = categoryCapacityFile;
    sw->totalNumSamplesProcessedByStorageWriter = 0;
    sw->totalNumSamplesCurrentlyInStorageWriter = 0;
    sw->totalNumSamplesAddedToStorageWriter = 0;
    sw->totalNumSamplesWrittenToFile = 0;

    /* allocate container for sample counter (per category) for buffer */
    sw->numStoredBuf = CALLOC(sw->numCategories, sizeof(u64)); /* CALLOC sets counters to zero */
    ASSERT(sw->numStoredBuf, "Allocation failed");
    if (!sw->numStoredBuf)   /* failed to allocate sample counter vector for storage buffer */
    {
        return 1;
    }

    /* allocate container for sample counter (per category) for sample file */
    sw->numStoredFile = CALLOC(sw->numCategories, sizeof(u64)); /* CALLOC sets counters to zero */
    ASSERT(sw->numStoredFile, "Allocation failed");
    if (!sw->numStoredFile)
    {
        FREE(sw->numStoredBuf);
        return 2; /* failed to allocate sample counter vector for file storage */
    }

    /* allocate space for file writing buffer */
    u64 categorySizeInBytes = categoryCapacityFile * LWE_SAMPLE_SIZE_IN_BYTES;
    sw->numCategoriesInFileWritingBuffer = APPROXIMATE_SIZE_IN_BYTES_OF_FILE_WRITER_BUFFER / categorySizeInBytes;
    if (sw->numCategoriesInFileWritingBuffer < 1)
    {
        FREE(sw->numStoredBuf);
        FREE(sw->numStoredFile);
        printf("ERROR: not enough space in FileWritingBuffer\n");
        return 3; /* probably a configuration error, file writing buffer holds very few categories */
    }
    sw->fileWritingBuffer = MALLOC(sw->numCategoriesInFileWritingBuffer * categoryCapacityFile * LWE_SAMPLE_SIZE_IN_BYTES);
    if (!sw->fileWritingBuffer)
    {
        FREE(sw->numStoredBuf);
        FREE(sw->numStoredFile);
        return 4; /* failed to allocate file writing buffer */
    }

    /* allocate (temp) buf */
#if defined(STORAGE_WRITER_CACHE_SIZE_IN_BYTES)
    u64 numBytes = (u64)STORAGE_WRITER_CACHE_SIZE_IN_BYTES;
#else
    u64 numBytes = (u64)getPhysicalMemorySize() / 2; /* half of physical memory*/
#endif
    sw->categoryCapacityBuf = numBytes / sw->numCategories / LWE_SAMPLE_SIZE_IN_BYTES;
    if (sw->categoryCapacityBuf > sw->categoryCapacityFile)
    {
        sw->categoryCapacityBuf = sw->categoryCapacityFile;
    }
    sw->buf = MALLOC(sw->numCategories * sw->categoryCapacityBuf * LWE_SAMPLE_SIZE_IN_BYTES);
    ASSERT(sw->buf, "Allocation failed");
//  double memUsed = sw->numCategories * sw->categoryCapacityBuf * LWE_SAMPLE_SIZE_IN_BYTES / 1024 / 1024 / (double)1024;
//  printf("%.2f GB used for storage writer cache\n", memUsed);

    /* create destination folder with params file and empty sample file */
    int ret = newStorageFolderWithGivenLweInstance(lwe, sw->dstFolderName);
    if (ret)
    {
        FREE(sw->numStoredBuf);
        FREE(sw->numStoredFile);
        FREE(sw->fileWritingBuffer);
        FREE(sw->buf);
        return 100 + ret; /* could not create destination folder */
    }

    /* create sample file of correct size */
    sw->f = fopenSamples(sw->dstFolderName, "wb+");
    if (!sw->f)
    {
        FREE(sw->numStoredBuf);
        FREE(sw->numStoredFile);
        FREE(sw->fileWritingBuffer);
        FREE(sw->buf);
        /* destination folder intentionally not deleted */
        /* lwe params intentionally not deleted */
        return 5; /* could not create destination sample file */
    }
    fileExtend(sw->f, categoryCapacityFile * sw->numCategories * LWE_SAMPLE_SIZE_IN_BYTES);

    /* create sample info file */
    ret = sampleInfoToFile(sw->dstFolderName, sw->bkwStepPar, sw->numCategories, sw->categoryCapacityFile, sw->totalNumSamplesWrittenToFile, sw->numStoredFile);
    if (ret)
    {
        FREE(sw->numStoredBuf);
        FREE(sw->numStoredFile);
        FREE(sw->fileWritingBuffer);
        FREE(sw->buf);
        /* destination folder intentionally not deleted */
        /* lwe params intentionally not deleted */
        /* samples file intentionally not deleted */
        return 6; /* could not create sample info file */
    }

    return 0;
}

inline int storageWriterFlush(storageWriter *sw)
{
//  printf("storageWriterFlush called at %6.02g%% load\n", storageWriterCurrentLoadPercentageCache(sw));
    if (sw->totalNumSamplesCurrentlyInStorageWriter == 0)
    {
        return 0; /* nothing to do, skip flushing */
    }

//  printf("tot %d %d\n", sw->totalNumSamplesCurrentlyInStorageWriter, sw->numCategoriesInFileWritingBuffer);
    ASSERT(sw->numCategoriesInFileWritingBuffer > 0, "numCategoriesInFileWritingBuffer must be > 0");

    u64 szInBytes = fileSize(sw->f);
    u64 szInSamples = szInBytes / LWE_SAMPLE_SIZE_IN_BYTES;
    u64 szInCategories = szInSamples / sw->categoryCapacityFile;

    /* flush to storage writer cache to file */

    /* set file position to beginning of file */
    fseeko64(sw->f, 0L, SEEK_SET);

    int oddBytes = szInBytes - (szInCategories * sw->categoryCapacityFile * LWE_SAMPLE_SIZE_IN_BYTES);
    if (oddBytes)
    {
        printf("*** destinaton file seems to hold %" PRIu64 " categories (and %d additional bytes)\n", szInCategories, oddBytes);
    }

    u64 currentDestinationCategory = 0;
    size_t numReadDestinationCategories;
    u64 totalNumDestinationCategoriesRead = 0;
    u64 numIterationsInReadLoop = 0;
    lweSample *s = sw->buf;
    _off64_t nextLastReadPosition = 0;
    _off64_t lastReadPosition = 0;
    do
    {
        /* read categories from destination file */
        nextLastReadPosition = lastReadPosition;
        lastReadPosition = ftello64(sw->f);
        if ((lastReadPosition - nextLastReadPosition) % (LWE_SAMPLE_SIZE_IN_BYTES * sw->categoryCapacityFile) != 0)
        {
            printf("*** non-integral number of categories (div=%d, mod=%d)\n", (lastReadPosition - nextLastReadPosition) / (LWE_SAMPLE_SIZE_IN_BYTES * sw->categoryCapacityFile), (lastReadPosition - nextLastReadPosition) % (LWE_SAMPLE_SIZE_IN_BYTES * sw->categoryCapacityFile));
        }
        if (lastReadPosition + LWE_SAMPLE_SIZE_IN_BYTES * sw->categoryCapacityFile * sw->numCategoriesInFileWritingBuffer > szInBytes)
        {
//      printf("*** about to read too far ***\n");
        }

        numReadDestinationCategories = fread(sw->fileWritingBuffer, LWE_SAMPLE_SIZE_IN_BYTES * sw->categoryCapacityFile, sw->numCategoriesInFileWritingBuffer, sw->f);
        if (ferror(sw->f))
        {
            perror("error on read");
            clearerr(sw->f);
        }

        totalNumDestinationCategoriesRead += numReadDestinationCategories;
        numIterationsInReadLoop++;
//    printf("%" PRIu64 " read (total %" PRIu64 " read, tried to get %" PRIu64 ", iteration %" PRIu64 ")\n", numReadDestinationCategories, totalNumDestinationCategoriesRead, sw->numCategoriesInFileWritingBuffer, numIterationsInReadLoop);

        u64 pos = ftello64(sw->f);
        if (pos > szInBytes)
        {
            printf("*** read too far (%" PRIu64 " bytes past end)\n", pos - szInBytes);
        }

        /* copy samples from cache buffer to flush buffer */
        lweSample *d = sw->fileWritingBuffer;
        for (u64 i=0; i<numReadDestinationCategories; i++, currentDestinationCategory++)
        {

            /* copy samples from  */
            u64 numSamplesToCopy = sw->numStoredBuf[currentDestinationCategory];
            u64 numInCurrentCategoryFile = sw->numStoredFile[currentDestinationCategory];
            if (numSamplesToCopy + numInCurrentCategoryFile > sw->categoryCapacityFile)   /* all samples do not fit on file */
            {
                numSamplesToCopy = sw->categoryCapacityFile - numInCurrentCategoryFile; /* copy only the ones that fit */
                printf("reduced from %" PRIu64 " to %" PRIu64 "\n", sw->numStoredBuf[currentDestinationCategory], numSamplesToCopy);
            }
            if (numSamplesToCopy > 0)
            {
                MEMCPY(d + sw->numStoredFile[currentDestinationCategory], s, numSamplesToCopy * LWE_SAMPLE_SIZE_IN_BYTES); /* copy samples from current category*/
            }
            sw->totalNumSamplesCurrentlyInStorageWriter -= sw->numStoredBuf[currentDestinationCategory];
            sw->numStoredBuf[currentDestinationCategory] = 0;
            sw->numStoredFile[currentDestinationCategory] += numSamplesToCopy;
            sw->totalNumSamplesWrittenToFile += numSamplesToCopy;

            s = s + sw->categoryCapacityBuf;
            d = d + sw->categoryCapacityFile;
        }

        /* back up to previous read position */
        fseeko64(sw->f, lastReadPosition, SEEK_SET);

        /* write adjusted buffer back to file (same position it was read from, but now with additional samples added) */
        size_t written = fwrite(sw->fileWritingBuffer, sw->categoryCapacityFile * LWE_SAMPLE_SIZE_IN_BYTES, numReadDestinationCategories, sw->f);
        if (ferror(sw->f))
        {
//      perror("error on write");
            clearerr(sw->f);
        }
        fseeko64(sw->f, 0L, SEEK_CUR); /* necessary */
        if (written != numReadDestinationCategories)
        {
//      printf("%d categories written to file (expected %d)\n", written, numReadDestinationCategories);
        }
        pos = ftello64(sw->f);
        if (pos > szInBytes)
        {
            printf("*** wrote too far (%" PRIu64 " bytes past end)\n", pos - szInBytes);
        }
//    printf("%d < %d - %d - %d \n", currentDestinationCategory, sw->numCategories, sw->totalNumSamplesCurrentlyInStorageWriter, sw->numCategoriesInFileWritingBuffer);
    }
    while (currentDestinationCategory < sw->numCategories);

    /* (over)write sample info file */
    int ret = sampleInfoToFile(sw->dstFolderName, sw->bkwStepPar, sw->numCategories, sw->categoryCapacityFile, sw->totalNumSamplesWrittenToFile, sw->numStoredFile);
    if (ret)
    {
        return 1; /* could not overwrite sample info file */
    }

//  printf(" (after flush, file storage at %6.02g%% load)\n", storageWriterCurrentLoadPercentageFile(sw));
    return 0;
}

int storageWriterFree(storageWriter *sw)
{
    int ret = storageWriterFlush(sw);
    if (ret)
    {
        return ret;
    }
    FREE(sw->buf);
    FREE(sw->numStoredBuf);
    FREE(sw->numStoredFile);
    fclose(sw->f);
    return 0;
}

inline int storageWriterHasRoom(storageWriter *sw, u64 categoryIndex)
{
    /* check if there is room for this sample on file and in cache */
    if (sw->numStoredBuf[categoryIndex] + sw->numStoredFile[categoryIndex] < sw->categoryCapacityFile)   /* if there is room on disk (after flush) */
    {

        if (sw->numStoredBuf[categoryIndex] < sw->categoryCapacityBuf - 1)   /* if there is room in storage writer cache */
        {
            return 0; /* sample can be added */
        }

        if (sw->numStoredBuf[categoryIndex] == sw->categoryCapacityBuf - 1)   /* if there is room in storage writer cache */
        {
            return 1; /* sample can be added, and this sample fills last available cache slot */
        }

        /* cache is full, so sample cannot be added */
        /* however, there will definitely be room for this sample if the cache is flushed to disk first */
        return 2; /* discarded, but can be added if cache is flushed (room on file but not in cache) */
    }

    /* category is full, no more samples can be added to this category */
    /* here, full means that the category on file is full, or at least will be when the storage writer is flushed */
    return 3;
}

/* reserves a memory area for a sample in the storage writer cache */
/* the caller is expected to copy the sample into the returned memory area */
/* if used correctly, a slot cache slot will always be available */
/* the variable storageWriterCategoryIsFull will be set when the added sample fills the last spot in the given category */
/* the caller is then required to call the flush function to write the currently cached contents to file */
inline lweSample *storageWriterAddSample(storageWriter *sw, u64 categoryIndex, int *storageWriterCategoryIsFull)
{
    sw->totalNumSamplesProcessedByStorageWriter++;
    *storageWriterCategoryIsFull = storageWriterHasRoom(sw, categoryIndex);
    switch (*storageWriterCategoryIsFull)
    {
    case 0: /* there is room */
    case 1: /* there is room, but this is the last slot in the cache so caller may want to flush */
        /* increase sample counters and return memory area for next sample */
        sw->totalNumSamplesCurrentlyInStorageWriter++;
        sw->totalNumSamplesAddedToStorageWriter++;
        lweSample *s = sw->buf + (categoryIndex * sw->categoryCapacityBuf) + sw->numStoredBuf[categoryIndex];
        sw->numStoredBuf[categoryIndex] = sw->numStoredBuf[categoryIndex] + 1;
        return s; /* return reserved memory area for sample */
    }
    return NULL; /* no room for sample */
}

inline void storageWriterUndoAddSample(storageWriter *sw, u64 categoryIndex)
{
    sw->totalNumSamplesCurrentlyInStorageWriter--;
    sw->totalNumSamplesAddedToStorageWriter--;
    sw->numStoredBuf[categoryIndex] = sw->numStoredBuf[categoryIndex] - 1;
}

double storageWriterCurrentLoadPercentageCache(storageWriter *sw)
{
    return 100 * sw->totalNumSamplesCurrentlyInStorageWriter / (double)(sw->categoryCapacityBuf * sw->numCategories);
}

double storageWriterCurrentLoadPercentageFile(storageWriter *sw)
{
    return 100 * sw->totalNumSamplesWrittenToFile / (double)(sw->categoryCapacityFile * sw->numCategories);
}

double storageWriterCurrentLoadPercentage(storageWriter *sw)
{
    return 100 * (sw->totalNumSamplesAddedToStorageWriter) / (double)(sw->categoryCapacityFile * sw->numCategories);
}

// void storageWriterPrint(storageWriter *sw)
// {
//     printf("sw.totalNumSamplesAddedToStorageWriter     = %12" PRIu64 "\n", sw->totalNumSamplesAddedToStorageWriter);
//     printf("sw.totalNumSamplesCurrentlyInStorageWriter = %12" PRIu64 "\n", sw->totalNumSamplesCurrentlyInStorageWriter);
//     printf("sw.totalNumSamplesProcessedByStorageWriter = %12" PRIu64 "\n", sw->totalNumSamplesProcessedByStorageWriter);
//     printf("sw.totalNumSamplesWrittenToFile            = %12" PRIu64 "\n", sw->totalNumSamplesWrittenToFile);
// }
