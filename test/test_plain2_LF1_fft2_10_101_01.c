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
#include "storage_file_utilities.h"
#include "test_functions.h"
#include "transition_bkw_step.h"
#include "workplace_localization.h"
#include "verify_samples.h"
#include "solve_fft.h"
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

/* The solver uses fft for guessing the remaining (after reduction) positions.
   If the amount of memory in your system is limited, you can probably solve for 2 positions.
   If you have more memory (about 32 Gb maybe), then you can probably solve for 3 positions.
   Solving for 4 positions is typically out of scope for a straightforward fft implementation,
   at least for the parameter sets that we are concerned with. */
#define MAX_NUM_FFT_POSITIONS 2
//#define MAX_NUM_FFT_POSITIONS 3

/* for testing */
#define NUM_REDUCTION_STEPS 4

static void printSolutionVector(short *s, int n, time_t start)
{
    timeStamp(start);
    printf("s = (%4d", s[0]);
    for (int i=1; i<n; i++)
    {
        if (i % 10 == 0)
        {
            printf(",\n");
            timeStamp(start);
            printf("     %4d", s[i]);
        }
        else
        {
            printf(", %4d", s[i]);
        }
    }
    printf(")\n");
}

/**************************************************************************
 * Main
 **************************************************************************/
int main()
{
//  u64 totalNumInitialSamples = 1000000000; /* 1 billion */
//  u64 totalNumInitialSamples = 100000000; /* 100 million */
//  u64 totalNumInitialSamples = 1000000; /* 1 million */
//  u64 totalNumInitialSamples = 10000;

    u64 totalNumInitialSamples = 100000;

    time_t start = time(NULL);
    srand(time(NULL));
    randomUtilRandomize();

    lweInstance lwe;
    int n, ret;

    n = 10;
    int q = 101;
    double alpha = 0.01;

    lweInit(&lwe, n, q, alpha);

    char outputfolder[128];
    char unsortedFolderName[256];
    char sortedFolderName[256];
    char srcFolderName[256];
    char dstFolderName[256];

    sprintf(outputfolder, "%s/test_plain2_LF1_fft_10_101_01", LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A);
    mkdir(outputfolder, 0777);

    sprintf(unsortedFolderName, "%s/original", outputfolder);
    sprintf(sortedFolderName, "%s/step_0", outputfolder);
    testCreateNewInstanceFolder(unsortedFolderName, n, q, alpha);
    newStorageFolder(&lwe, unsortedFolderName, n, q, alpha);
    addSamplesToSampleFile(unsortedFolderName, totalNumInitialSamples, start);

    /* set bkw step parameters */
    /* the number of reduction steps is set above in the macro NUM_REDUCTION_STEPS.
       reduction step i uses bkwStepPar[i] to define how to produce new samples (plain bkw, lms,...).
       in addition to this, bkwStepPar[i+1] is used to determine how to sort the resulting samples.
       different reduction types may be mixed arbitrarily.
       for example, while one reduction step may use plain bkw, the next may use lms. */
    bkwStepParameters bkwStepPar[NUM_REDUCTION_STEPS + 1];

    for (int i=0; i<NUM_REDUCTION_STEPS+1; i++)
    {
        bkwStepPar[i].sorting = plainBKW;
        bkwStepPar[i].startIndex = i == 0 ? 0 : bkwStepPar[i-1].startIndex + bkwStepPar[i-1].numPositions;
        bkwStepPar[i].numPositions = 2;
        bkwStepPar[i].selection = LF1;
    }

    /* sort (unsorted) samples */
    timeStamp(start);
    printf("sorting initial samples\n");
    u64 minDestinationStorageCapacityInSamples = totalNumInitialSamples * 4 / 3; /* add about 25% storage room for sorted samples */

    ret = transition_unsorted_2_sorted(unsortedFolderName, sortedFolderName, minDestinationStorageCapacityInSamples, &bkwStepPar[0], start);
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


    /* execute several consecutive bkw steps */
    int numReductionSteps = NUM_REDUCTION_STEPS;
    for (int i=0; i<numReductionSteps; i++)
    {

        /* process BKW step */
        timeStamp(start);
        printf("Reduction step %02d -> %02d, %s reduction at positions %d to %d (destination sorting using positions %d to %d)\n", i, i+1, sortingAsString(bkwStepPar[i+1].sorting), bkwStepPar[i].startIndex, bkwStepPar[i].startIndex + bkwStepPar[i].numPositions - 1, bkwStepPar[i+1].startIndex, bkwStepPar[i+1].startIndex + bkwStepPar[i+1].numPositions - 1);

        int ret;

        /* define src and dst folders */
        sprintf(srcFolderName, "%s/step_%d", outputfolder, i);
        sprintf(dstFolderName, "%s/step_%d", outputfolder, i+1);

        /* perform reduction step */
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

    /* solving phase (using fft to guess) */
    int numSolvedCoordinates = 0;
    sprintf(srcFolderName, "%s/step_%d", outputfolder, numReductionSteps);
    int paramsReadError = lweParametersFromFile(&lwe, srcFolderName);
    if (paramsReadError)
    {
        timeStamp(start);
        printf("Exiting, could not read lwe parameters from %s\n", srcFolderName);
        exit(1);
    }
    n = lwe.n;

    short solution[n];
    for (int i = 0; i < n; ++i)
        solution[i] = 0;

    const int numPositionsToGuessAtEnd = MAX_NUM_FFT_POSITIONS;
    numSolvedCoordinates = n - (bkwStepPar[NUM_REDUCTION_STEPS-1].startIndex + bkwStepPar[NUM_REDUCTION_STEPS-1].numPositions + numPositionsToGuessAtEnd);

    /* provide partial (transformed) solution */
    for (int i=0; i<numSolvedCoordinates; i++)
    {
        solution[n - 1 - i] = lwe.s[n - 1 - i];
    }
    short s[n];
    for (int i = 0; i < n; i++)
    {
        s[i] = lwe.s[i];
    }

    timeStamp(start);
    printf("Correct s\n");
    printSolutionVector(s, n, start);

    /* print initial (partial) solution */
    timeStamp(start);
    printf("initial s\n");
    printSolutionVector(solution, n, start);

    /* solve for secret (generic implementation, handles arbitrary reductions steps, solves as much as possible) */
    for (int i=0; i<numReductionSteps+1; i++)
    {

        sprintf(srcFolderName, "%s/step_%d", outputfolder, numReductionSteps - i);

        if (!folderExists(srcFolderName))
        {
            timeStamp(start);
            printf("Halting solver, cannot find folder %s\n", srcFolderName);
            break;
        }
        timeStamp(start);
        printf("using samples at step %02d (folder %s)\n", numReductionSteps - i, srcFolderName);
        int numUnsolvedPositionsForThisReductionStep = i == 0 ? numPositionsToGuessAtEnd : bkwStepPar[numReductionSteps - i].numPositions;
        int numSolvablePositions = i == 0 ? numPositionsToGuessAtEnd : MIN(numUnsolvedPositionsForThisReductionStep, MAX_NUM_FFT_POSITIONS);
        if (numSolvablePositions < numUnsolvedPositionsForThisReductionStep)
        {
            timeStamp(start);
            printf("not able to solve below positions %d (at folder %s)\n", n - numSolvedCoordinates, srcFolderName);
            break;
        }

        int startIndexForSolving = n - numSolvedCoordinates - numSolvablePositions;
        timeStamp(start);
        if (numSolvablePositions == 1)
        {
            printf("solving for position %d...", startIndexForSolving);
        }
        else
        {
            printf("solving for positions %d to %d...", startIndexForSolving, startIndexForSolving + numSolvablePositions - 1);
        }
//    int ret = solve_fft_search(srcFolderName, solution, numSolvedCoordinates, numSolvablePositions, FFT_SOLVER_SINGLE_PRECISION);
        int ret = solve_fft_search(srcFolderName, solution, numSolvedCoordinates, numSolvablePositions, FFT_SOLVER_DOUBLE_PRECISION);
        printf("done\n");
        if (ret)
        {
            printf("*** main: solve_fft_search returned %d\n", ret);
        }
        printSolutionVector(solution, n, start);
        numSolvedCoordinates += numSolvablePositions;
        numUnsolvedPositionsForThisReductionStep -= numSolvablePositions;
    }

    int correctly_solved = 1;
    for (int i=0; i<n; i++)
    {
        if (solution[i] != lwe.s[i])
        {
            correctly_solved = 0;
        }
    }
    timeStamp(start);
    printf("recovered s\n");
    printSolutionVector(solution, n, start);

    timeStamp(start);
    printf("correct s\n");
    timeStamp(start);
    printf("s = (%4d", lwe.s[0]);
    for (int i=1; i<n; i++)
    {
        if (i % 10 == 0)
        {
            printf(",\n");
            timeStamp(start);
            printf("    %4d", lwe.s[i]);
        }
        else
        {
            printf(", %4d", lwe.s[i]);
        }
    }
    printf(")\n");

    if (correctly_solved)
    {
        timeStamp(start);
        printf("Test passed!\n");
    }
    else
    {
        timeStamp(start);
        printf("Test failed!\n");
    }

    return 0;
}

