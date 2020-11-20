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

#ifndef STORAGE_FILE_UTILITIES_H
#define STORAGE_FILE_UTILITIES_H
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "lwe_instance.h"
#include "bkw_step_parameters.h"

/* utility functions */
void parameterFileName(char *paramFileName, const char *folderName); /* parameter file name from folder name */
void samplesFileName(char *samplesFileName, const char *folderName); /* samples file name from folder name */
void samplesInfoFileName(char *samplesInfoFileName, const char *folderName); /* samples info file name from folder name */

/* lwe instance folders */
int newStorageFolderWithGivenLweInstance(lweInstance *lwe, const char *folderName); /* creates folder with parameter file (with given parameters) and an empty samples and samples info file */
int newStorageFolder(lweInstance *lwe, const char *folderName, int n, int q, double alpha); /* generates new lwe parameters, creates folder with parameter file and an empty samples and samples info file */
int deleteStorageFolder(const char *folderName, int deleteParFile, int deleteSampleInfoFile, int deleteSamples); /* deletes folder, including parameter and samples file */

/* parameter file */
int parametersToFile(lweInstance *lwe, const char *folderName);
int lweParametersFromFile(lweInstance *lwe, const char *folderName);

/* sample info file */
int sampleInfoToFile(const char *folderName, bkwStepParameters *bkwStepPar, u64 numCategories, u64 categoryCapacity, u64 numTotalSamples, u64 *numSamplesPerCategory);
int sampleInfoFromFile(const char *folderName, bkwStepParameters *bkwStepPar, u64 *numCategories, u64 *categoryCapacity, u64 *numTotalSamples, u64 *numSamplesPerCategory);

/* samples file */
FILE *fopenSamples(const char *folderName, const char *mode);
u64 freadSamples(FILE *f, lweSample *sampleBuf, u64 numSamples); /* read sample range from current position into buffer (does not close file) */
u64 freadCategories(FILE *f, lweSample *sampleBuf, u64 numCategories, u64 categoryCapacityInSamples); /* read category range from current position into buffer (does not close file) */
u64 fsetAndReadSamples(FILE *f, lweSample *sampleBuf, u64 startingSample, u64 nbrOfSamples); /* read specific sample range into buffer (does not close file) */
u64 numSamplesInSampleFile(const char *folderName); /* number of samples in sample file */
u64 addSamplesToSampleFile(const char *folderName, u64 nbrOfSamples, time_t start); /* add samples to sample file */
u64 readSamplesFromSampleFile(lweSample *sampleBuf, const char *folderName, u64 startingSample, u64 nbrOfSamples); /* read sample range into buffer */

/* conversion (import) of TU Darmstadt (lwe) problem instances to local format, with optional sample amplification */
int convertTUDarmstadtProblemInstanceToNativeFormat(lweInstance *lwe, const char *srcFileName, const char *dstfolderName, int useSampleAmplification, u64 totalNumSamples);

/* erik */
void destinationFileName(char *destFileName, const char *folderName);/* destination numbers file name from folder name */
u64 fileSize(FILE *fp);
void fileExtend(FILE *fp, u64 size);

/* possibly superfluous utilities */
int folderExists(const char *folderPath);

int fileExists (const char *Filename);

#endif
