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

#include "press_any_key.h"
#include "memory_utils.h"
#include "assert_utils.h"
#include "num_cores.h"
#include "lwe_instance.h"
#include "statistics_utils.h"
#include "log_utils.h"
#include "string_utils.h"
#include "transition_unsorted_2_sorted.h"
#include "transition_bkw_step_final.h"
#include "storage_file_utilities.h"
#include "test_functions.h"
#include "transition_bkw_step.h"
#include "workplace_localization.h"
#include "iterator_samples.h"
#include "verify_samples.h"
#include "bkw_step_parameters.h"
#include "random_utils.h"
#include "transition_times2_modq.h"
#include "transition_mod2.h"
#include "solve_fwht.h"

#define NUM_REDUCTION_STEPS 5

int main(int argc, char* argv[])
{
//  u64 totalNumInitialSamples = 1000000000; /* 1 billion */
//  u64 totalNumInitialSamples = 100000000; /* 100 million */
//  u64 totalNumInitialSamples = 1000000; /* 1 million */
    u64 totalNumInitialSamples = 1000;

    time_t start = time(NULL);
    srand(time(NULL));
    randomUtilRandomize();

    lweInstance lwe, lpn;
    int ret;
    int n = 3*NUM_REDUCTION_STEPS;
    int q = 101;
    double alpha = 0.005;
    // double alpha = 0.01;
    // double alpha = 0.02;
    //double alpha = 0.0265;

    lweInit(&lwe, n, q, alpha);

    char originalFolderName[256];
    char sortedFolderName[256];
    sprintf(originalFolderName, "%s/tu_darmstadt_%d_%d_%f_original", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A, n, q, alpha);
    sprintf(sortedFolderName, "%s/tu_darmstadt_%d_%d_%f_step_00", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A, n, q, alpha);

    testCreateNewInstanceFolder(originalFolderName, n, q, alpha);
    newStorageFolder(&lwe, originalFolderName, n, q, alpha);
    addSamplesToSampleFile(originalFolderName, totalNumInitialSamples, start);

    u64 minDestinationStorageCapacityInSamples = round((double)(totalNumInitialSamples*4)/3); /* add about 25% storage room for sorted samples */

    /* set bkw step parameters */
    bkwStepParameters bkwStepPar[NUM_REDUCTION_STEPS];

    /* Set steps: smooth LMS */
    for (int i=0; i<NUM_REDUCTION_STEPS; i++)
    {
        bkwStepPar[i].sorting = smoothLMS;
        bkwStepPar[i].startIndex = i == 0 ? 0 : bkwStepPar[i-1].startIndex + bkwStepPar[i-1].numPositions;
        bkwStepPar[i].numPositions = 3;
        bkwStepPar[i].selection = LF2;
        bkwStepPar[i].sortingPar.smoothLMS.p = 19; // test
        bkwStepPar[i].sortingPar.smoothLMS.p1 = 35; // test
        bkwStepPar[i].sortingPar.smoothLMS.p2 = 19; // Should be equal to p by construction
        bkwStepPar[i].sortingPar.smoothLMS.prev_p1 = i == 0 ? -1 : bkwStepPar[i-1].sortingPar.smoothLMS.p1;
        bkwStepPar[i].sortingPar.smoothLMS.meta_skipped = 2;
        ASSERT(bkwStepPar[i].sortingPar.smoothLMS.p2 != 0, "smooth-LMS p2 parameter not valid");
    }

    /* Perform multiplication by 2 to each sample -  sort (unsorted) samples */
    timeStamp(start);
    printf("Performing multiplication times 2 mod q and sorting samples\n");
    const char *src = originalFolderName;
    const char *dst = sortedFolderName;
    ret = transition_times2_modq(src, dst, minDestinationStorageCapacityInSamples, &bkwStepPar[0], start);
    switch (ret)
    {
    case 0: /* transition computed ok */
        printSampleVerificationOfSortedFolder(dst, start, &bkwStepPar[0]); /* verify sorted samples */
        break;
    case 1: /* sorting unnecessary (destination folder already exists) */
        timeStamp(start);
        printf("skipping, destination folder %s already exists\n\n", dst);
        break;
    default:
        timeStamp(start);
        printf("error %d when sorting initial samples\n", ret);
        exit(1);
    }

    /* perform all but last smooth LMS BKW reduction steps */
    char srcFolderName[256];
    char dstFolderName[256];
    int numReductionSteps = NUM_REDUCTION_STEPS;

    for (int i=0; i<numReductionSteps-1; i++)
    {
        /* process smooth LMS BKW step */
        timeStamp(start);
        printf("Reduction step %02d -> %02d, %s reduction at positions %d to %d (destination sorting using positions %d to %d)\n", i, i+1, sortingAsString(bkwStepPar[i+1].sorting), bkwStepPar[i].startIndex, bkwStepPar[i].startIndex + bkwStepPar[i].numPositions - 1, bkwStepPar[i+1].startIndex, bkwStepPar[i+1].startIndex + bkwStepPar[i+1].numPositions);
        int ret;
        sprintf(srcFolderName, "%s/tu_darmstadt_%d_%d_%f_step_%02d", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A, n, q, alpha,i);
        sprintf(dstFolderName, "%s/tu_darmstadt_%d_%d_%f_step_%02d", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A, n, q, alpha,i+1);

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
    }

    /* perform last reduction step */
    int i = numReductionSteps-1;
    timeStamp(start);
    printf("Last reduction step %02d -> %02d, %s reduction at positions %d to %d\n", i, i+1, sortingAsString(bkwStepPar[i].sorting), bkwStepPar[i].startIndex, bkwStepPar[i].startIndex + bkwStepPar[i].numPositions - 1);
    sprintf(srcFolderName, "%s/tu_darmstadt_%d_%d_%f_step_%02d", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A, n, q, alpha,i);
    sprintf(dstFolderName, "%s/tu_darmstadt_%d_%d_%f_step_final", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A, n, q, alpha);
    timeStamp(start);
    printf("  src folder: %s\n", srcFolderName);
    timeStamp(start);
    printf("  dst folder: %s\n", dstFolderName);

    u64 numSamplesStored;
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
    sprintf(srcFolderName, "%s/tu_darmstadt_%d_%d_%f_step_final", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A, n, q, alpha);
    sprintf(dstFolderName, "%s/tu_darmstadt_%d_%d_%f_binary", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A,  n, q, alpha);
    ret = transition_mod2(srcFolderName, dstFolderName, start);
    switch (ret)
    {
    case 0: /* transition computed ok */
        timeStamp(start);
        printf("Start binary sample verification for computing error rate\n");
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

    lweParametersFromFile(&lpn, dstFolderName);

    /* Solving phase - using Fast Walsh Hadamard Tranform */

    int bruteForcePositions = 4;
    int fwht_positions = lwe.n - bruteForcePositions;

    timeStamp(start);
    printf("Solving phase - brute force from positions %d to %d - Fast Walsh Hadamard Transform from positions %d to %d\n", 0, bruteForcePositions-1, bruteForcePositions, bruteForcePositions+fwht_positions-1);

    u8 *binary_solution = CALLOC(fwht_positions, sizeof(u8));

//  ret = solve_fwht_search(srcFolderName, binary_solution, bruteForcePositions, fwht_positions, start);
    ret = solve_fwht_search_hybrid(srcFolderName, binary_solution, bruteForcePositions, fwht_positions, start);
    if(ret)
    {
        printf("error %d in solve_fwht_search_hybrid\n", ret);
        exit(-1);
    }

    printf("\nFound Solution\n");
    for(int i = 0; i<lpn.n; i++)
        printf("%d ",binary_solution[i]);
    printf("\n");

    printf("\nReal Solution\n");
    for(int i = 0; i<lpn.n; i++)
        printf("%d ",lpn.s[i]);
    printf("\n");

    return 0;
}















