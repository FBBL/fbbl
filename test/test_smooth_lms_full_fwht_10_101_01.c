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
#include "transition_reduce_secret.h"
#include "transition_unsorted_2_sorted.h"
#include "transition_bkw_step_final.h"
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

#define NUM_REDUCTION_STEPS 5
#define BRUTE_FORCE_POSITIONS 0

int main()
{
//  u64 totalNumInitialSamples = 1000000000; /* 1 billion */
//  u64 totalNumInitialSamples = 100000000; /* 100 million */
//  u64 totalNumInitialSamples = 1000000; /* 1 million */
    u64 totalNumInitialSamples = 10000;

    time_t start = time(NULL);
    srand(time(NULL));
    randomUtilRandomize();

    lweInstance lwe, lpn;
    int ret;
    int n = 10;
    int q = 101;
    double alpha = 0.01;

    lweInit(&lwe, n, q, alpha);

    char outputfolder[128];
    char originalFolderName[256];
    char sortedFolderName[256];
    char srcFolderName[256];
    char dstFolderName[256];

    sprintf(outputfolder, "%s/test_smooth_lms_full_fwht_10_101_01", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A);
    mkdir(outputfolder, 0777);


    sprintf(originalFolderName, "%s/original_0", outputfolder);

    testCreateNewInstanceFolder(originalFolderName, n, q, alpha);
    newStorageFolder(&lwe, originalFolderName, n, q, alpha);
    addSamplesToSampleFile(originalFolderName, totalNumInitialSamples, start);

    u64 minDestinationStorageCapacityInSamples = round((double)(totalNumInitialSamples*4)/3); /* add about 25% storage room for sorted samples */

    /* set bkw step parameters */
    bkwStepParameters bkwStepPar[NUM_REDUCTION_STEPS];

#ifdef CODEDBKW
    /* Set steps: codedBKW(2,1), 2 positions */
    for (int i=0; i<NUM_REDUCTION_STEPS; i++)
    {
        bkwStepPar[i].sorting = codedBKW;
        bkwStepPar[i].startIndex = i == 0 ? 0 : bkwStepPar[i-1].startIndex + bkwStepPar[i-1].numPositions;
        bkwStepPar[i].numPositions = 2;
        bkwStepPar[i].selection = LF2;
        bkwStepPar[i].sortingPar.CodedBKW.ct = blockCode_21;
    }
#else
    /* Set steps: smooth LMS */
    for (int i=0; i<NUM_REDUCTION_STEPS; i++)
    {
        bkwStepPar[i].sorting = smoothLMS;
        bkwStepPar[i].startIndex = i == 0 ? 0 : bkwStepPar[i-1].startIndex + bkwStepPar[i-1].numPositions;
        bkwStepPar[i].numPositions = 2;
        bkwStepPar[i].selection = LF2;
        bkwStepPar[i].sortingPar.smoothLMS.p = 9; // test
        bkwStepPar[i].sortingPar.smoothLMS.p1 = 18; // test
        bkwStepPar[i].sortingPar.smoothLMS.p2 = bkwStepPar[i].sortingPar.smoothLMS.p;
        bkwStepPar[i].sortingPar.smoothLMS.prev_p1 = i == 0 ? -1 : bkwStepPar[i-1].sortingPar.smoothLMS.p1;
        bkwStepPar[i].sortingPar.smoothLMS.meta_skipped = 0;
        bkwStepPar[i].sortingPar.smoothLMS.unnatural_selection_ts = 0;
        ASSERT(bkwStepPar[i].sortingPar.smoothLMS.p2 != 0, "smooth-LMS p2 parameter not valid");
    }
#endif

    int bruteForcePositions = BRUTE_FORCE_POSITIONS;
    int fwht_positions = lwe.n - bruteForcePositions;

    int MAX_digits = ceil(log2(4*alpha*q));

    u8 binary_solution[MAX_digits][fwht_positions]; //CALLOC(fwht_positions*MAX_digits, sizeof(u8));
    timeStamp(start);
    printf("Start reduction phase - MAX Number of Iterations %d\n", MAX_digits);

    for(int j = 0; j < MAX_digits; j++)
    {

        timeStamp(start);
        printf("Start Iteration %d\n", j);

        if (j > 0)
        {
            sprintf(srcFolderName, "%s/original_%d", outputfolder, j-1);
            sprintf(dstFolderName, "%s/original_%d", outputfolder, j);

            /* Get new samples with reduced secret */
            timeStamp(start);
            printf("Get new set of initial samples with reduced secret\n");
            ret =  transition_reduce_secret(srcFolderName, dstFolderName, binary_solution[j-1], start);
            switch (ret)
            {
            case 0: /* transition computed ok */
                printSampleVerificationOfUnsortedFolder(dstFolderName, start);
                break;
            case 1: /* times2_modq unnecessary (destination folder already exists) */
                timeStamp(start);
                printf("skipping, destination folder %s already exists\n\n", dstFolderName);
                break;
            default:
                timeStamp(start);
                printf("error %d when generating new samples with reduced secret\n", ret);
                exit(1);
            }
        }

        sprintf(originalFolderName, "%s/original_%d", outputfolder, j);
        sprintf(sortedFolderName, "%s/step_0_%d", outputfolder, j);

        /* sort (unsorted) samples */
        timeStamp(start);
        printf("Perform Multiplication times 2 mod q and sort samples\n");
        ret = transition_times2_modq(originalFolderName, sortedFolderName, minDestinationStorageCapacityInSamples, &bkwStepPar[0], start);
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

        /* perform all but last smooth LMS BKW reduction steps */
        int numReductionSteps = NUM_REDUCTION_STEPS;

        for (int i=0; i<numReductionSteps-1; i++)
        {
            /* process smooth LMS BKW step */
            timeStamp(start);
            printf("Reduction step %02d -> %02d, %s reduction at positions %d to %d (destination sorting using positions %d to %d)\n", i, i+1, sortingAsString(bkwStepPar[i+1].sorting), bkwStepPar[i].startIndex, bkwStepPar[i].startIndex + bkwStepPar[i].numPositions - 1, bkwStepPar[i+1].startIndex, bkwStepPar[i+1].startIndex + bkwStepPar[i+1].numPositions);
            sprintf(srcFolderName, "%s/step_%d_%d", outputfolder, i, j);
            sprintf(dstFolderName, "%s/step_%d_%d", outputfolder, i+1, j);

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
        sprintf(srcFolderName, "%s/step_%d_%d", outputfolder, i, j);
        sprintf(dstFolderName, "%s/step_final_%d", outputfolder, j);
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
        sprintf(srcFolderName, "%s/step_final_%d", outputfolder, j);
        sprintf(dstFolderName, "%s/binary_%d", outputfolder, j);

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

        timeStamp(start);
        printf("Solving phase - Fast Walsh Hadamard Transform from position %d to %d\n", 0, fwht_positions-1);

        ret = solve_fwht_search(srcFolderName, binary_solution[j], 0, fwht_positions, start);
        //  ret = solve_fwht_search_hybrid(srcFolderName, binary_solution[j], bruteForcePositions, fwht_positions, start);
        if(ret)
        {
            printf("error %d in solve_fwht_search_hybrid\n", ret);
            exit(-1);
        }

        printf("\n");
        timeStamp(start);
        printf("Binary Solution Found (");
        for(int i = 0; i<lpn.n; i++)
            printf("%d ",binary_solution[j][i]);
        printf(")\n");

        timeStamp(start);
        printf("Real Binary Solution  (");
        for(int i = 0; i<lpn.n; i++)
            printf("%d ",lpn.s[i]);
        printf(")\n");

        // TERMINATE if recovered secret is not correct - for TESTING
        for(int i = 0; i<lpn.n; i++)
            if (binary_solution[j][i] != lpn.s[i])
            {
                printf("Recovered binary secret is not correct - Terminating the program.\n");
                exit(0);
            }
    }

    printf("\n");
    timeStamp(start);
    printf("Sequence of binary solutions found\n");
    for(int j=0; j<MAX_digits; j++)
    {
        printf("(");
        for(int i = 0; i<lwe.n; i++)
        {
            printf("%d ", binary_solution[j][i]);
        }
        printf(")\n");
    }

    short full_secret[lwe.n];

    retrieve_full_secret(full_secret, binary_solution, MAX_digits, lwe.n, lwe.q);

    printf("\n");
    timeStamp(start);
    printf("Full Secret Found    (");
    for(int i = 0; i<lwe.n; i++)
        printf("%d ", full_secret[i]);
    printf(")\n");

    sprintf(originalFolderName, "%s/original_0", outputfolder);
    lweParametersFromFile(&lwe, originalFolderName);

    timeStamp(start);
    printf("Original Full Secret (");
    for(int i = 0; i<lpn.n; i++)
        printf("%d ",lwe.s[i]);
    printf(")\n\n");

    for(int i = 0; i<lwe.n; i++)
    {
        if (full_secret[i] != lwe.s[i])
        {
            printf("WRONG retrieved solution!\n");
            return 1;
        }
    }

    printf("Test passed\n");

    return 0;
}

