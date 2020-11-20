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

#ifndef STORAGE_READER_H
#define STORAGE_READER_H
#include <stdio.h>
#include "bkw_step_parameters.h"

/* a buffer is used when reading the content of the storage to file */
#define APPROXIMATE_SIZE_IN_BYTES_OF_FILE_READER_BUFFER (250 * 1024 * 1024)

typedef struct
{
    char srcFolderName[512];
    FILE *f; /* file handle to sample file */
    lweSample *buf; /* big sample buffer, stores samples read from file */
    lweSample *minibuf; /* mini buffer for one category, used when a category pair is split between two (big) buffer reads */
    bkwStepParameters srcBkwStepPar; /* bkw step parameters */
    u64 numCategories; /* total number of categories */
    u64 categoryCapacity; /* how many samples a category can contain (fixed size for each folder) */
    u64 *numSamplesPerCategory; /* counter array to keep track of the number of stored samples in each category */
    u64 indexOfFirstCategoryInBuffer; /* index of the first category that currently resides in the sample buffer */
    u64 numCategoriesInBuffer; /* number of categories that have currently been read into the sample buffer */
    u64 bufferCapacityNumCategories; /* maximum number of categories that the sample buffer can hold */
    u64 currentCategoryIndex; /* state of the storage reader, indicates which category (index) that is next to be output (sequentially, starting at zero) */
    /* used for LMS only */
    int c;
//  u64 *numSingletonCategories;
//  u64 *singletonCategoryIndices;
    /* stats for testing purposes only */
    u64 totalNumCategoriesReadFromFile;
} storageReader;

int storageReaderInitialize(storageReader *sr, const char *srcFolderName);
void storageReaderFree(storageReader *sr);

int storageReaderGetNextAdjacentCategoryPair(storageReader *sr, lweSample **buf1, u64 *numSamplesInBuf1, lweSample **buf2, u64 *numSamplesInBuf2);

#endif
