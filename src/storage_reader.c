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

#include "storage_reader.h"
#include "memory_utils.h"
#include "storage_file_utilities.h"
#include "position_values_2_category_index.h"
#include <inttypes.h>

int storageReaderInitialize(storageReader *sr, const char *srcFolderName)
{
    ASSERT(sr, "unexpected parameter");
    ASSERT(srcFolderName, "unexpected parameter");
    strncpy(sr->srcFolderName, srcFolderName, 512);
    u64 numTotalSamples;
    int ret = sampleInfoFromFile(sr->srcFolderName, &sr->srcBkwStepPar, &sr->numCategories, &sr->categoryCapacity, &numTotalSamples, NULL);
    if (ret)
    {
        return 1; /* could not read sample info file */
    }
    u64 categorySizeInBytes = sr->categoryCapacity * LWE_SAMPLE_SIZE_IN_BYTES;
    sr->bufferCapacityNumCategories = APPROXIMATE_SIZE_IN_BYTES_OF_FILE_READER_BUFFER / categorySizeInBytes;
    if (sr->bufferCapacityNumCategories < 3)
    {
        sr->bufferCapacityNumCategories = 3;
    }
    if (sr->srcBkwStepPar.sorting == LMS)
    {
        lweInstance lwe;
        lweParametersFromFile(&lwe, srcFolderName);
        bkwStepParameters bkwStepPar;
        if (sampleInfoFromFile(srcFolderName, &bkwStepPar, NULL, NULL, NULL, NULL))
        {
            lweDestroy(&lwe);
            return 2; /* could not get sample info */
        }
        /* when reading lms-sorted samples we care about c = q/p + 1. */
        /* whether c is even or odd determines the number and location of singleton categories. */
        sr->c = lwe.q / bkwStepPar.sortingPar.LMS.p + 1;
    }

    /* allocate container for sample counter (per category) */
    sr->numSamplesPerCategory = CALLOC(sr->numCategories, sizeof(u64)); /* calloc used to zero-initialize */
    if (!sr->numSamplesPerCategory)
    {
        return 3; /* could not allocate num */
    }
    ret = sampleInfoFromFile(sr->srcFolderName, NULL, NULL, NULL, NULL, sr->numSamplesPerCategory);
    if (ret)
    {
        return 4; /* could not read sample counts per category */
    }

    /* allocate read buffer */
    u64 bufferSizeInBytes = sr->bufferCapacityNumCategories * categorySizeInBytes;
    sr->buf = MALLOC(bufferSizeInBytes);
    if (!sr->buf)
    {
        FREE(sr->numSamplesPerCategory);
        return 5; /* could not allocate buf */
    }
    sr->indexOfFirstCategoryInBuffer = 0;
    sr->numCategoriesInBuffer = 0;
    sr->currentCategoryIndex = 0;
    sr->totalNumCategoriesReadFromFile = 0;

    /* allocate mini buffer */
    sr->minibuf = MALLOC(categorySizeInBytes);
    if (!sr->minibuf)
    {
        FREE(sr->numSamplesPerCategory);
        FREE(sr->buf);
        return 6; /* could not allocate minibuf */
    }

    /* open source file */
    sr->f = fopenSamples(sr->srcFolderName, "rb");
    if (!sr->f)
    {
        FREE(sr->numSamplesPerCategory);
        FREE(sr->buf);
        FREE(sr->minibuf);
        return 7; /* could not open source file */
    }

    return 0;
}

void storageReaderFree(storageReader *sr)
{
    FREE(sr->numSamplesPerCategory);
    FREE(sr->buf);
    FREE(sr->minibuf);
    fclose(sr->f);
}

static size_t fillBuf(storageReader *sr, int numCategoriesToRead)
{
    if (feof(sr->f))
    {
        return 0; /* end of file reached */
    }
    size_t numReadDestinationCategories = fread(sr->buf, LWE_SAMPLE_SIZE_IN_BYTES * sr->categoryCapacity, numCategoriesToRead, sr->f);
    if ((numReadDestinationCategories != numCategoriesToRead) || ferror(sr->f))
    {
        clearerr(sr->f);
    }
    if (sr->numCategoriesInBuffer > 0)
    {
        sr->indexOfFirstCategoryInBuffer += sr->numCategoriesInBuffer;
    }
    sr->numCategoriesInBuffer = numReadDestinationCategories;
    return numReadDestinationCategories;
}

int storageReaderGetNextAdjacentCategoryPair(storageReader *sr, lweSample **buf1, u64 *numSamplesInBuf1, lweSample **buf2, u64 *numSamplesInBuf2)
{
    lweSample *cat1;
    lweSample *cat2;
    int firstTimeReadingFromFile = sr->numCategoriesInBuffer == 0;
    int exhaustedCategoriesInReadBuffer = sr->currentCategoryIndex >= (sr->indexOfFirstCategoryInBuffer + sr->numCategoriesInBuffer);

    /* read from file if there are no categories available in the buffer */
    if (firstTimeReadingFromFile || exhaustedCategoriesInReadBuffer)
    {
        u64 numCategoriesReadFromFile = fillBuf(sr, sr->bufferCapacityNumCategories);
        if (numCategoriesReadFromFile == 0)
        {
            *buf1 = *buf2 = NULL;
            *numSamplesInBuf1 = *numSamplesInBuf2 = 0;
            return 0; /* could not read from file, no category returned */
        }
    }

    /* there is now at least one category available in the buffer (could be just one) */

    if (is_singleton(&sr->srcBkwStepPar, sr->currentCategoryIndex, sr->numCategories))   /* singleton category */
    {
        u64 offsetInCategories = sr->currentCategoryIndex - sr->indexOfFirstCategoryInBuffer;
        *buf1 = sr->buf + (offsetInCategories * sr->categoryCapacity); /* category 1 */
        *buf2 = NULL;                                               /* no category 2 */
        *numSamplesInBuf1 = sr->numSamplesPerCategory[sr->currentCategoryIndex];
        *numSamplesInBuf2 = 0;
        sr->currentCategoryIndex = sr->currentCategoryIndex + 1;
        return 1; /* singleton category returned */
    }

    u64 numCategoriesAvailableInBuf = sr->indexOfFirstCategoryInBuffer + sr->numCategoriesInBuffer - sr->currentCategoryIndex;

    /* return adjacent categories */
    u64 offsetInCategories = sr->currentCategoryIndex - sr->indexOfFirstCategoryInBuffer;
    if (numCategoriesAvailableInBuf >= 2)   /* both adjacent categories available in buf */
    {
        /* two categories available in buf */
        cat1 = sr->buf + (offsetInCategories * sr->categoryCapacity);
        cat2 = cat1 + sr->categoryCapacity;
    }
    else     /* only one category available in buf */
    {
        ASSERT(numCategoriesAvailableInBuf == 1, "unexpected number of available categories");
        /* copy last remaining category to minibuf */
        u64 categorySizeInBytes = sr->categoryCapacity * LWE_SAMPLE_SIZE_IN_BYTES;
        cat1 = sr->buf + (offsetInCategories * sr->categoryCapacity);
        MEMCPY(sr->minibuf, cat1, categorySizeInBytes);
        /* reload buffer */
        u64 numCategoriesReadFromFile = fillBuf(sr, sr->bufferCapacityNumCategories);
        if (numCategoriesReadFromFile == 0)
        {
            printf("*** storageReaderGetNextAdjacentCategoryPair_lms: error reading categories");
        }
        ASSERT(numCategoriesReadFromFile > 0, "could not categories from file read as expected");
        /* return buf1 from minibuf, buf2 from buf */
        u64 offsetInCategories = sr->currentCategoryIndex + 1 - sr->indexOfFirstCategoryInBuffer;
        cat1 = sr->minibuf;
        cat2 = sr->buf + (offsetInCategories * sr->categoryCapacity);
    }
    *buf1 = cat1;
    *buf2 = cat2;
    *numSamplesInBuf1 = sr->numSamplesPerCategory[sr->currentCategoryIndex];
    *numSamplesInBuf2 = sr->numSamplesPerCategory[sr->currentCategoryIndex + 1];
    sr->currentCategoryIndex = sr->currentCategoryIndex + 2;
    return 2; /* adjacent categories returned */
}

