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

#include "memory_utils.h"
#include "assert_utils.h"
#include "lwe_instance.h"
#include "log_utils.h"
#include "string_utils.h"
#include "transition_unsorted_2_sorted.h"
#include "transition_times2_modq.h"
#include "transition_bkw_step_final.h"
#include "transition_bkw_step_final_smooth_lms_meta.h"
#include "transition_mod2.h"
#include "solve_fwht.h"
#include "storage_file_utilities.h"
#include "test_functions.h"
#include "transition_bkw_step.h"
#include "workplace_localization.h"
#include "verify_samples.h"
#include "config_bkw.h"
// #include "Python.h" // needs to be included before any standard headers, also make sure to add Python37.dll (or whichever version you need to use) to project
#include "bkw_step_parameters.h"
#include "transform_secret.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
// #include <dir.h>
#include <unistd.h>
#include <stdio.h>
#include <dirent.h>
#include <inttypes.h>
#include <sys/stat.h>

#define ZERO_COORDINATES 0

#define NUM_LMZ_STEPS 5

#define NUM_REDUCTION_STEPS 5

#define BRUTE_FORCE_POSITIONS 0

//#define DELETE_FOLDERS_FOR_PREVIOUS_BKW_STEPS

int main()
{
//  u64 totalNumInitialSamples = 1000000000; /* 1 billion */
//  u64 totalNumInitialSamples = 100000000; /* 100 million */
//  u64 totalNumInitialSamples = 1000000; /* 1 million */
//  u64 totalNumInitialSamples = 30000;

    u64 totalNumInitialSamples = 100000;

    time_t start = time(NULL);
    srand(time(NULL));
    randomUtilRandomize();

    lweInstance lwe;
    int n, ret, q;
    n = 10;
    q = 101;
    float alpha = 0.01;
    lweInit(&lwe, n, q, alpha);

    char outputfolder[128];
    char challengeFileName[256];
    char sortedFolderName[256];
    char srcFolderName[256];
    char dstFolderName[256];
    char convertedFolderName[256];

    sprintf(outputfolder, "%s/test_smooth_lms_sample_ampli_fwht_10_101_01", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A);
    mkdir(outputfolder, 0777);

    sprintf(challengeFileName, "%s/LWE_10_101_01.txt", LOCAL_TEST_VECTORS_PATH_PREFIX);
    sprintf(convertedFolderName, "%s/converted", outputfolder);
    sprintf(sortedFolderName, "%s/step_0", outputfolder);

    /* convert (lwe) problem instance to local format */
    timeStamp(start);
    printf("converting instance to local format\n");
    ret = tuDarmstadtFileFormatConversionWithErrorChecking(challengeFileName, convertedFolderName, 1, totalNumInitialSamples, start);
    switch (ret)
    {
    case 0: /* file conversion ok */
        printSampleVerificationOfUnsortedFolder(convertedFolderName, start); /* verify unsorted samples */
        break;
    case 1: /* file conversion unnecessary (destination folder already exists) */
        timeStamp(start);
        printf("skipping, destination folder %s already exists\n\n", convertedFolderName);
        break;
    default:
        timeStamp(start);
        printf("error %d when converting TU Darmstadt instance\n", ret);
        exit(1);
    }
    u64 minDestinationStorageCapacityInSamples = totalNumInitialSamples * 5 / 4; /* add about 25% storage room for sorted samples */

    /* set bkw step parameters */
    bkwStepParameters bkwStepPar[NUM_REDUCTION_STEPS];

//     n = 40, alpha = 0.005*sqrt(3) - (only 1600 samples available), 8 smoothplainBKW + 4 smoothLMS
    // int start_index[NUM_REDUCTION_STEPS] =                     {0,     2,    5,    7,  10,  13,    15,  18,    21,  25,   29,    34};
    // int len_step[NUM_REDUCTION_STEPS] =                        {2,     3,    2,    3,   3,   2,     3,   2,     4,   5,    5,     6};
    // int p_step[NUM_REDUCTION_STEPS] =                          {1,     1,    1,    1,   1,   1,     1,   1,    20,  28,   37,    53};
    // int p1_step[NUM_REDUCTION_STEPS] =                         {10,  180,    2,   27, 540,   6,    97,   1,   170, 601,  401,     1};
    // int prev_p1_step[NUM_REDUCTION_STEPS] =                    {-1,   10,  180,    2,  27, 540,     6,  97,    -1, 170,  601,   401};
    // int unnatural_selection_ts[NUM_REDUCTION_STEPS] =          {0,     0,    0,    0,   0,   0,     0,   0,     0,   0,    0,     0};
    // int unnatural_selection_start_index[NUM_REDUCTION_STEPS] = {0,     0,    0,    0,   0,   0,     0,   0,     0,   0,    0,     0};
    // int meta_skipped[NUM_REDUCTION_STEPS] =                    {0,     0,    0,    0,   0,   0,     0,   0,     0,   0,    0,     0};
    // short knownSecret[] = {451, 410, 1055, 253, 1038, 456, 1300, 406, 1586, 134, 1260, 5, 454, 1247, 1069, 800, 444, 199, 1093, 135, 1175, 676, 1374, 69, 596, 871, 770, 261, 658, 1575, 657, 1370, 601, 1492, 913, 516, 138, 623, 35, 1550};


    // n = 10, alpha = 0.01, 5 smoothLMS
    int start_index[NUM_REDUCTION_STEPS] =  {0,    2,     4,    6,    8};
    int len_step[NUM_REDUCTION_STEPS] =     {2,    2,     2,    2,    2};
    int p_step[NUM_REDUCTION_STEPS] =       {2,    3,     4,    6,    8};
    int p1_step[NUM_REDUCTION_STEPS] =      {21,  21,    21,   21,   21};
    int prev_p1_step[NUM_REDUCTION_STEPS] = {-1,  21,    21,   21,   21};
    int unnatural_selection_ts[NUM_REDUCTION_STEPS] = {0, 0, 0, 0, 0};
    short knownSecret[] = {64, 36, 101, 23, 66, 62, 40, 15, 3, 43};

    /* Set parameters */
    for (int i=0; i<NUM_REDUCTION_STEPS; i++)
    {
        bkwStepPar[i].sorting = smoothLMS;
        bkwStepPar[i].startIndex = start_index[i];
        bkwStepPar[i].numPositions = len_step[i];
        bkwStepPar[i].selection = LF2;
        bkwStepPar[i].sortingPar.smoothLMS.p = p_step[i]; // test
        bkwStepPar[i].sortingPar.smoothLMS.p1 = p1_step[i]; // test
        bkwStepPar[i].sortingPar.smoothLMS.p2 = bkwStepPar[i].sortingPar.smoothLMS.p; // must be set to p
        bkwStepPar[i].sortingPar.smoothLMS.prev_p1 = prev_p1_step[i];
        bkwStepPar[i].sortingPar.smoothLMS.meta_skipped = 0;//meta_skipped[i];
        bkwStepPar[i].sortingPar.smoothLMS.unnatural_selection_ts = unnatural_selection_ts[i];
        bkwStepPar[i].sortingPar.smoothLMS.unnatural_selection_start_index = 0;//unnatural_selection_start_index[i];
        // char ns[256];
        // sprintf_u64_delim(ns, num_categories(&lwe, &bkwStepPar[i]));
        // printf("levl %d: %d %d num Categories  %s \n", i+1, bkwStepPar[i].startIndex, bkwStepPar[i].numPositions, ns);
    }

//  exit(0);

    /* Multiply times 2 and sort samples */
    timeStamp(start);
    printf("Multiply times 2 (mod q) and sort initial samples\n");
    ret = transition_times2_modq(convertedFolderName, sortedFolderName, minDestinationStorageCapacityInSamples, &bkwStepPar[0], start);
    switch (ret)
    {
    case 0: /* transition computed ok */
        printSampleVerificationOfSortedFolder(sortedFolderName, start, &bkwStepPar[0]); /* verify sorted samples */
        break;
    case 1: /* sorting unnecessary (destination folder already exists) */
        timeStamp(start);
        printf("skipping, destination folder %s already exists\n\n", sortedFolderName);
        break;
    default:
        timeStamp(start);
        printf("error %d when sorting initial samples\n", ret);
        exit(1);
    }

#ifdef DELETE_FOLDERS_FOR_PREVIOUS_BKW_STEPS
    /* delete old folder */
    ret = deleteStorageFolder(originalFolderName, 0, 0, 1); /* deletes folder, including parameter and samples file */
    timeStamp(start);
    printf("Deleted folder %s\n", originalFolderName);
#endif


    /* perform all but last smooth LMS BKW reduction steps */
    int numReductionSteps = NUM_REDUCTION_STEPS;

    for (int i=0; i<NUM_REDUCTION_STEPS-1; i++)
    {
        /* process plain BKW step */
        timeStamp(start);
        printf("Reduction step %02d -> %02d, %s reduction at positions %d to %d (destination sorting using positions %d to %d)\n", i, i+1, sortingAsString(bkwStepPar[i+1].sorting), bkwStepPar[i].startIndex, bkwStepPar[i].startIndex + bkwStepPar[i].numPositions - 1, bkwStepPar[i+1].startIndex, bkwStepPar[i+1].startIndex + bkwStepPar[i+1].numPositions);
        int ret;
        sprintf(srcFolderName, "%s/step_%d", outputfolder, i);
        sprintf(dstFolderName, "%s/step_%d", outputfolder, i+1);

        u64 numSamplesStored;
        ret = transition_bkw_step(srcFolderName, dstFolderName, &bkwStepPar[i], &bkwStepPar[i+1], &numSamplesStored, start);
        switch (ret)
        {
        case 0: /* reduction computed ok */
            printSampleVerificationOfSortedFolder(dstFolderName, start, &bkwStepPar[i+1]); /* verify samples in destination folder */
            break;
        case 100: /* reduction step unnecessary (destination folder already exists) */
            timeStamp(start);
            printf("skipping, destination folder %s already exists\n\n", dstFolderName);
            break;
        default:
            timeStamp(start);
            printf("error %d in reduction step %d\n", ret, i);
            timeStamp(start);
            printf("  src folder: %s\n", srcFolderName);
            timeStamp(start);
            printf("  dst folder: %s\n", dstFolderName);
            exit(1);
        }
#ifdef DELETE_FOLDERS_FOR_PREVIOUS_BKW_STEPS
        /* delete old folder */

        if (i == 0)
        {
            ret = deleteStorageFolder(sortedFolderName, 0, 0, 1); /* deletes folder, including parameter and samples file */
            timeStamp(start);
            printf("Deleted folder %s\n", sortedFolderName);
        }
        else
        {
            char deleteFolderName[256];
            sprintf(deleteFolderName, "%s/step_%d", outputfolder, i);
            ret = deleteStorageFolder(deleteFolderName, 0, 0, 1); /* deletes folder, including parameter and samples file */
            timeStamp(start);
            printf("Deleted folder %s\n", deleteFolderName);
        }
        if (ret)
        {
            timeStamp(start);
            printf("*** Error code %d returned from deleteStorageFolder\n", ret);
        }
#endif
    }

    /* perform last reduction step */
    int i = numReductionSteps-1;
    timeStamp(start);
    printf("Last reduction step %02d -> %02d, %s reduction at positions %d to %d\n", i, i+1, sortingAsString(bkwStepPar[i].sorting), bkwStepPar[i].startIndex, bkwStepPar[i].startIndex + bkwStepPar[i].numPositions - 1);
    sprintf(srcFolderName, "%s/step_%d", outputfolder, i);
    sprintf(dstFolderName, "%s/step_final", outputfolder);
    timeStamp(start);
    printf("  src folder: %s\n", srcFolderName);
    timeStamp(start);
    printf("  dst folder: %s\n", dstFolderName);

    u64 numSamplesStored;
    if(bkwStepPar[i].sortingPar.smoothLMS.meta_skipped)
        ret = transition_bkw_step_final_smooth_lms_meta(srcFolderName, dstFolderName, &bkwStepPar[i], &numSamplesStored, start);
    else
        ret = transition_bkw_step_final(srcFolderName, dstFolderName, &bkwStepPar[i], &numSamplesStored, start);
    switch (ret)
    {
    case 0: /* reduction computed ok */
        printSampleVerificationOfUnsortedFolder(dstFolderName, start); /* verify samples in destination folder */
        break;
    case 100: /* reduction step unnecessary (destination folder already exists) */
        timeStamp(start);
        printf("skipping, destination folder %s already exists\n\n", dstFolderName);
        break;
    default:
        timeStamp(start);
        printf("error %d in reduction step %d\n", ret, i);
        timeStamp(start);
        printf("src folder: %s\n", srcFolderName);
        timeStamp(start);
        printf("dst folder: %s\n", dstFolderName);
        exit(1);
    }

    /* reduce all the system modulo 2 - to compute error rate - used only for testing */
    timeStamp(start);
    printf("Start binary sample verification for computing error rate\n");
    sprintf(srcFolderName, "%s/step_final", outputfolder);
    sprintf(dstFolderName, "%s/step_binary", outputfolder);

    ret = transition_mod2(srcFolderName, dstFolderName, start);
    switch (ret)
    {
    case 0: /* transition computed ok */
        printBinarySampleVerification(dstFolderName, start);
        break;
    case 100: /* mod2 unnecessary (destination folder already exists) */
        timeStamp(start);
        printf("skipping, destination folder %s already exists\n\n", dstFolderName);
        break;
    default:
        timeStamp(start);
        printf("error %d when reducing modulo 2 samples. Were there enough initial samples?\n", ret);
        exit(1);
    }


#ifdef DELETE_FOLDERS_FOR_PREVIOUS_BKW_STEPS
    /* delete old folder */
    ret = deleteStorageFolder(dstFolderName, 0, 0, 1); /* deletes folder, including parameter and samples file */
    timeStamp(start);
    printf("Deleted folder %s\n", dstFolderName);
#endif

    lweParametersFromFile(&lwe, convertedFolderName);
    transformSecret(&lwe, knownSecret);

    /* Solving phase - using Fast Walsh Hadamard Tranform */

    int zeroPositions = ZERO_COORDINATES;
    int fwhtPositions = lwe.n - zeroPositions;

    u8 *binary_solution = CALLOC(fwhtPositions, sizeof(u8));

    timeStamp(start);
    printf("Solving phase - Fast Walsh Hadamard Transform from %d to %d\n", zeroPositions, zeroPositions+fwhtPositions-1);

    ret = solve_fwht_search(srcFolderName, binary_solution, zeroPositions, fwhtPositions, start);
    if(ret)
    {
        printf("error %d in solve_fwht_search_hybrid\n", ret);
        exit(-1);
    }

    printf("\nFound Solution\n");
    for(int i = 0; i<fwhtPositions; i++)
        printf("%d ",binary_solution[i]);
    printf("\n");

    printf("\nReal Solution\n");
    int bin;
    for(int i = zeroPositions; i<zeroPositions+fwhtPositions; i++)
    {
        bin = knownSecret[i] > lwe.q/2 ? (knownSecret[i]%2 +1)%2 : knownSecret[i]%2;
        printf("%d ", bin);
    }
    printf("\n");

    for(int i = zeroPositions; i<zeroPositions+fwhtPositions; i++)
    {
        bin = knownSecret[i] > lwe.q/2 ? (knownSecret[i]%2 +1)%2 : knownSecret[i]%2;
        if (binary_solution[i] != bin)
        {
            printf("WRONG retrieved solution!\n");
            return 1;
        }
    }

    printf("Test passed\n");

    return 0;
}
