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

#ifndef TEST_FUNCTIONS
#define TEST_FUNCTIONS
#include <time.h>
#include "lwe_instance.h"
#include "bkw_step_parameters.h"

#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))

void readSamplesTest(char *folderName, int n, int start, int numSamples);
void bkwFileBasedTest(int n, int q, double alpha);
void printSampleDataForUnsortedSamples(lweInstance *lwe, const char *folderName, time_t start);
void testPositionValue2IndexBWK(int q);

void testCreateNewInstanceFolder(const char *folderName, int n, int q, double alpha);
void testAddSamplesToSampleFile(const char *folderName, u64 numSamplesToAdd);

u64 expectedNumBKWsamplesToGenerate(int numCategories, int *numSamplesPerCategory);

void printSampleVerificationOfUnsortedFolder(const char *folderName, time_t start);
void printSampleVerificationOfSortedFolder(const char *srcFolderName, time_t start, bkwStepParameters *bkwStepPar);
void printBinarySampleVerification(const char *folderName, time_t start);

void test_modular_matrix_inversion(lweInstance *lwe, int num_samples);

int tuDarmstadtFileFormatConversionWithErrorChecking(const char *srcFileName, const char *dstFolderName, int useSampleAmplification, u64 totalNumSamples, time_t start);

#endif
