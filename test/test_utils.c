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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <dirent.h>
#include <inttypes.h>
#include <sys/stat.h>

#include "lwe_sorting.h"
#include "memory_utils.h"
#include "assert_utils.h"
#include "lwe_instance.h"
#include "log_utils.h"
#include "string_utils.h"
#include "storage_file_utilities.h"
#include "test_functions.h"
#include "transition_bkw_step.h"
#include "workplace_localization.h"
#include "verify_samples.h"
#include "bkw_step_parameters.h"
#include "random_utils.h"
#include "transform_secret.h"
#include "position_values_2_category_index.h"

#define NUM_REDUCTION_STEPS 5
#define BRUTE_FORCE_POSITIONS 0

int main()
{

    time_t start = time(NULL);

    // TEST 1

    if(unordered == sortingFromString("unordered"))
    {
        timeStamp(start);
        printf("Test on sortingFromString(unordered): success\n");
    }
    else
    {
        timeStamp(start);
        printf("Error in sortingFromString(unordered)\n");
        return -1;
    }

    if(plainBKW == sortingFromString("Plain BKW"))
    {
        timeStamp(start);
        printf("Test on sortingFromString(plainBKW): success\n");
    }
    else
    {
        timeStamp(start);
        printf("Error in sortingFromString(Plain BKW)\n");
        return -1;
    }
    if(LMS == sortingFromString("Lazy Modulus Switching (LMS)"))
    {
        timeStamp(start);
        printf("Test on sortingFromString(Lazy Modulus Switching (LMS)): success\n");
    }
    else
    {
        timeStamp(start);
        printf("Error in sortingFromString(LMS)\n");
        return -1;
    }
    if(smoothLMS == sortingFromString("Smooth LMS"))
    {
        timeStamp(start);
        printf("Test on sortingFromString(Smooth LMS): success\n");
    }
    else
    {
        timeStamp(start);
        printf("Error in sortingFromString(smoothLMS)\n");
        return -1;
    }
    if(codedBKW == sortingFromString("Coded BKW"))
    {
        timeStamp(start);
        printf("Test on sortingFromString(Coded BKW): success\n");
    }
    else
    {
        timeStamp(start);
        printf("Error in sortingFromString(Coded BKW)\n");
        return -1;
    }

    // TEST 2

    timeStamp(start);
    printf("Testing LWE transformations\n");

    lweInstance lwe;
    int n, q, ret;
    n = 10;
    q = 101;
    float alpha = 0.01;
    lweInit(&lwe, n, q, alpha);

    char outputfolder[128];
    char challengeFileName[256];
    sprintf(outputfolder, "%s/test_utils", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A);
    mkdir(outputfolder, 0777);
    sprintf(challengeFileName, "%s/LWE_10_101_01.txt", LOCAL_TEST_VECTORS_PATH_PREFIX);
    char convertedFolderName[256];
    sprintf(convertedFolderName, "%s/converted", outputfolder);

    u64 totalNumInitialSamples = 100;

    ret = tuDarmstadtFileFormatConversionWithErrorChecking(challengeFileName, convertedFolderName, 1, totalNumInitialSamples, start);
    if (ret == 1)
    {
        timeStamp(start);
        printf("File already exists: %s\n", challengeFileName);
    }
    else if (ret != 0 && ret != 1)
    {
        timeStamp(start);
        printf("Error converting filie\n");
    }

    int paramsReadError = lweParametersFromFile(&lwe, convertedFolderName);
    if (paramsReadError)
    {
        timeStamp(start);
        printf("Exiting, could not read lwe parameters from %s\n", convertedFolderName);
        exit(1);
    }

    /* known solution */
    short original_s[] = {64, 36, 0, 23, 66, 62, 40, 15, 3, 43};
    short s[] = {64, 36, 0, 23, 66, 62, 40, 15, 3, 43};

    /* apply initial transformation to s */
    transformSecret(&lwe, s);
    for (int i=0; i<n; i++)
        lwe.s[i] = s[i];

    timeStamp(start);
    printf("Transformed secret s\n(");
    for (int i = 0; i < n; ++i)
    {
        printf("%hi ", lwe.s[i]);
    }
    printf(")\n");

    /* apply inverse transformation */
    inverseTransformSecret(&lwe, s);

    timeStamp(start);
    printf("Original s\n(");
    for (int i = 0; i < n; ++i)
    {
        printf("%hi ", original_s[i]);
    }
    printf(")\n");

    timeStamp(start);
    printf("Computed s\n(");
    for (int i = 0; i < n; ++i)
    {
        printf("%hi ", s[i]);
    }
    printf(")\n");

    for (int i = 0; i < n; ++i)
    {
        if (original_s[i] != s[i])
        {
            timeStamp(start);
            printf("Error in transformation\n");
            return 1;
        }
    }

    // TEST 3 - just to increase coverage...
    printf("Start test random functions\n");
    u8 fakeRandBuf[] = {1,2};
    randomUtilAppendRandomness(&lwe.rnd, fakeRandBuf, 2);
    long double randld = randomUtilLongDouble(&lwe.rnd);
    timeStamp(start);
    printf("Long double generated %Lf\n", randld);

    // TEST 4 - just to increase coverage...

    lweSample *emptySample, *randomSample;
    srand(time(NULL));
    randomUtilRandomize();

    emptySample = lwe.newEmptySample();
    randomSample = lwe.newRandomSample(n, q, alpha*q, &lwe.rnd, s);


    // TEST 5 - just to increase coverage...
    srand(time(NULL));
    u64 category_index;
    short a[2], b[2];

    for (int i = 0; i < 1000; ++i)
    {
        a[0] = rand()%q;
        a[1] = rand()%q;
        b[0] = a[0];
        b[1] = a[1];
        category_index = position_values_2_category_index_plain_bkw(q, a);
        category_index_2_position_values_plain_bkw(q, category_index, b);
        if (a[0] != b[0] || a[1] != b[1])
        {
            timeStamp(start);
            printf("Error in position_values_2_category_index_plain_bkw or category_index_2_position_values_plain_bkw\n");
            return 1;
        }

    }
    free_table_plain_bkw_2_positions();
    lwe.freeSample(emptySample);
    lwe.freeSample(randomSample);

    lweDestroy(&lwe);
    timeStamp(start);
    printf("Test passed\n");

    return 0;

}