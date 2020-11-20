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

#ifndef VERIFY_SAMPLES_H
#define VERIFY_SAMPLES_H
#include "bkw_step_parameters.h"

void verifyOneSampleSorted(lweInstance *lwe, lweSample *sample, bkwStepParameters *bkwStepPar, u64 *numIncorrectSums, u64 *numIncorrectHashes, u64 expectedCategoryIndex, u64 *numIncorrectCategoryClassifications, int printOnError);

int verifyUnsortedSamples(const char *folderName, u64 *numSamplesProcessed, u64 *numIncorrectSums, u64 *numIncorrectHashes, int printOnError);
int verifySortedSamples(const char *folderName, bkwStepParameters *bkwStepPar, u64 *totalNumSamplesProcessed, u64 *numIncorrectSums, u64 *numIncorrectHashes, u64 *numIncorrectCategoryClassifications, int printOnError);

#endif
