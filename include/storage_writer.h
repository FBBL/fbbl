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

#ifndef STORAGE_WRITER_H
#define STORAGE_WRITER_H
#include "bkw_step_parameters.h"
#include "config_compiler.h"
#include<stdio.h>

/*
  this macro is used to set the size of the storage writer cache.
  this is the main storage container for the storage writer.
  flushing to file is slow, and you generally want to make your storage writer cache as large as possible.
  if the storage writer can hold more samples, you do not need to flush as often.

  if the macro is not defined, the code will default to using half of the available physical memory
  (if available, otherwise half of that if available, otherwise half of that if available,...).
 */
//#define STORAGE_WRITER_CACHE_SIZE_IN_BYTES (2 * 1024 * 1024 * (u64)1024)
#define STORAGE_WRITER_CACHE_SIZE_IN_BYTES (100 * 1024 * 1024 * (u64)1024)

/* a file writer buffer is temporarily used when flushing the content of the storage writer to file */
#define APPROXIMATE_SIZE_IN_BYTES_OF_FILE_WRITER_BUFFER (512 * 1024 * 1024)

typedef struct
{
    char dstFolderName[512];
    FILE *f;
    lweSample *buf;
    bkwStepParameters *bkwStepPar;
    u64 numCategories;
    u64 categoryCapacityBuf;
    u64 *numStoredBuf;
    u64 categoryCapacityFile;
    u64 *numStoredFile;
    u64 numCategoriesInFileWritingBuffer;
    lweSample *fileWritingBuffer;
    /* stats for testing purposes only */
    u64 totalNumSamplesProcessedByStorageWriter; /* num items added to storage writer, including those that were discarded for lack of room */
    u64 totalNumSamplesCurrentlyInStorageWriter; /* num items currently in storage writer cache (in memory) */
    u64 totalNumSamplesWrittenToFile; /* num items currently written to file */
    u64 totalNumSamplesAddedToStorageWriter; /* num items in storage writer, counting both cache and on file */
} storageWriter;

int storageWriterInitialize(storageWriter *dsh, const char *dstFolderName, lweInstance *lwe, bkwStepParameters *bkwStepPar, u64 categoryCapacityFile);
int storageWriterFree(storageWriter *dsh);

int storageWriterHasRoom(storageWriter *dsh, u64 categoryIndex);
lweSample *storageWriterAddSample(storageWriter *dsh, u64 categoryIndex, int *storageWriterCategoryIsFull);
void storageWriterUndoAddSample(storageWriter *sw, u64 categoryIndex);
int storageWriterFlush(storageWriter *sw);

/* helper functions */
double storageWriterCurrentLoadPercentage(storageWriter *sw);
double storageWriterCurrentLoadPercentageCache(storageWriter *sw);
double storageWriterCurrentLoadPercentageFile(storageWriter *sw);
void storageWriterPrint(storageWriter *sw);

#endif
