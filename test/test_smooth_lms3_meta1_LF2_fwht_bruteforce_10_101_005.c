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

#include "memory_utils.h"
#include "assert_utils.h"
#include "lwe_instance.h"
#include "log_utils.h"
#include "string_utils.h"
#include "transition_unsorted_2_sorted.h"
#include "transition_bkw_step_final_smooth_lms_meta.h"
#include "storage_file_utilities.h"
#include "test_functions.h"
#include "transition_bkw_step.h"
#include "workplace_localization.h"
#include "verify_samples.h"
#include "bkw_step_parameters.h"
#include "random_utils.h"
#include "transition_times2_modq.h"
#include "transition_mod2.h"
#include "solve_fwht.h"

#define NUM_REDUCTION_STEPS 3
#define BRUTE_FORCE_POSITIONS 1
#define NUMTHREADS 2

int main()
{
//  u64 totalNumInitialSamples = 1000000000; /* 1 billion */
//  u64 totalNumInitialSamples = 100000000; /* 100 million */
//  u64 totalNumInitialSamples = 1000000; /* 1 million */
    u64 totalNumInitialSamples = 50000;

    time_t start = time(NULL);
    srand(time(NULL));
    randomUtilRandomize();

    lweInstance lwe;
    int ret;
    int n = 10;
    int q = 101;
    double alpha = 0.005;
    // double alpha = 0.01;
    // double alpha = 0.02;
    //double alpha = 0.0265;

    lweInit(&lwe, n, q, alpha);

    char outputfolder[128];
    char originalFolderName[256];
    char sortedFolderName[256];
    char srcFolderName[256];
    char dstFolderName[256];

    sprintf(outputfolder, "%s/test_smooth_lms3_meta1_LF2_fwht_bruteforce_10_101_005", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A);
    mkdir(outputfolder, 0777);

    sprintf(originalFolderName, "%s/original", outputfolder);

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
        bkwStepPar[i].sortingPar.smoothLMS.meta_skipped = 1;
        bkwStepPar[i].sortingPar.smoothLMS.unnatural_selection_ts = 0;
        ASSERT(bkwStepPar[i].sortingPar.smoothLMS.p2 != 0, "smooth-LMS p2 parameter not valid");
    }

    int bruteForcePositions = BRUTE_FORCE_POSITIONS;
    int fwht_positions = lwe.n - bruteForcePositions;
    int zeropositions = 0;

    u8 binary_solution[fwht_positions];
    short bf_solution[bruteForcePositions];

    /* sort (unsorted) samples */
    sprintf(originalFolderName, "%s/original", outputfolder);
    sprintf(sortedFolderName, "%s/step_00", outputfolder);

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
    int numReductionSteps = NUM_REDUCTION_STEPS;

    timeStamp(start);
    printf("Start reduction phase\n");

    for (int i=0; i<numReductionSteps-1; i++)
    {
        /* process smooth LMS BKW step */
        timeStamp(start);
        printf("Reduction step %02d -> %02d, %s reduction at positions %d to %d (destination sorting using positions %d to %d)\n", i, i+1, sortingAsString(bkwStepPar[i+1].sorting), bkwStepPar[i].startIndex, bkwStepPar[i].startIndex + bkwStepPar[i].numPositions - 1, bkwStepPar[i+1].startIndex, bkwStepPar[i+1].startIndex + bkwStepPar[i+1].numPositions);
        int ret;
        sprintf(srcFolderName, "%s/step_%02d", outputfolder, i);
        sprintf(dstFolderName, "%s/step_%02d", outputfolder, i+1);

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
    printf("Last reduction step %02d -> %02d, %s reduction at positions %d to %d\n", i, i+1, sortingAsString(bkwStepPar[i].sorting), bkwStepPar[i].startIndex, bkwStepPar[i].startIndex + bkwStepPar[i].numPositions);
    sprintf(srcFolderName, "%s/step_%02d", outputfolder, i);
    sprintf(dstFolderName, "%s/step_final", outputfolder);
    timeStamp(start);
    printf("  src folder: %s\n", srcFolderName);
    timeStamp(start);
    printf("  dst folder: %s\n", dstFolderName);

    u64 numSamplesStored;
    ret = transition_bkw_step_final_smooth_lms_meta(srcFolderName, dstFolderName, &bkwStepPar[i], &numSamplesStored, start);
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

    sprintf(srcFolderName, "%s/step_final", outputfolder);
    lweDestroy(&lwe);
    lweParametersFromFile(&lwe, originalFolderName);

    /* compute binary secret */
    u8 real_binary_secret[lwe.n];
    for (int i = 0; i < lwe.n; ++i)
    {
        if (lwe.s[i] < q/2)
            real_binary_secret[i] = lwe.s[i] % 2;
        else
            real_binary_secret[i] = (lwe.s[i]+1) % 2;
    }

    /* Solving phase - using Fast Walsh Hadamard Tranform */
    timeStamp(start);
    printf("Solving phase - brute force from positions %d to %d - Fast Walsh Hadamard Transform from positions %d to %d\n", 0, bruteForcePositions-1, bruteForcePositions, bruteForcePositions+fwht_positions-1);

    ret = solve_fwht_search_bruteforce(srcFolderName, binary_solution, bf_solution, zeropositions, bruteForcePositions, fwht_positions, start, NUMTHREADS);
    if(ret)
    {
        printf("error %d in solve_fwht_search_hybrid\n", ret);
        exit(-1);
    }

    printf("\n");
    timeStamp(start);
    printf("Binary Solution Found (");
    for(int i = 0; i<lwe.n-bruteForcePositions; i++)
        printf("%hhu ",binary_solution[i]);
    printf("- ");
    for(int i = 0; i<bruteForcePositions; i++)
        printf("%hi ",bf_solution[i]);
    printf(")\n");

    timeStamp(start);
    printf("Real Binary Solution  (");
    for(int i = 0; i<lwe.n-bruteForcePositions; i++)
        printf("%hhu ",real_binary_secret[i]);
    printf("- ");
    for(int i = 0; i<bruteForcePositions; i++)
        printf("%hi ",lwe.s[i+zeropositions+fwht_positions]);
    printf(")\n");

    for(int i = 0; i<lwe.n-bruteForcePositions; i++)
    {
        if (binary_solution[i] != real_binary_secret[i])
        {
            printf("WRONG retrieved solution!\n");
            return 1;
        }
    }
    for(int i = 0; i<bruteForcePositions; i++)
    {
        if (bf_solution[i] != lwe.s[i+zeropositions+fwht_positions])
        {
            printf("WRONG retrieved solution!\n");
            return 1;
        }
    }

    lweDestroy(&lwe);

    printf("Test passed\n");
    return 0;
}




