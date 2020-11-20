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

#include "press_any_key.h"
#include "memory_utils.h"
#include "assert_utils.h"
#include "num_cores.h"
#include "lwe_instance.h"
#include "statistics_utils.h"
#include "log_utils.h"
#include "string_utils.h"
#include "transition_unsorted_2_sorted.h"
#include "storage_file_utilities.h"
#include "test_functions.h"
#include "transition_bkw_step.h"
#include "workplace_localization.h"
#include "iterator_samples.h"
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

/* The solver uses fft for guessing the remaining (after reduction) positions.
   If the amount of memory in your system is limited, you can probably solve for 2 positions.
   If you have more memory (about 32 Gb maybe), then you can probably solve for 3 positions.
   Solving for 4 positions is typically out of scope for a straightforward fft implementation,
   at least for the parameter sets that we are concerned with. */
#define MAX_NUM_FFT_POSITIONS 2
//#define MAX_NUM_FFT_POSITIONS 3

/* for testing */
#define NUM_REDUCTION_STEPS 3

/* to maximize the number of samples we can use */
// #define DELETE_FOLDERS_FOR_PREVIOUS_BKW_STEPS

/* create samples ourselves */
//#define UNLIMITED_SAMPLES


static void printSolutionVector(int *s, int n, time_t start)
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
int main(int argc, char* argv[])
{
//  u64 totalNumInitialSamples = 1000000000; /* 1 billion */
//  u64 totalNumInitialSamples = 100000000; /* 100 million */
//  u64 totalNumInitialSamples = 1000000; /* 1 million */
//  u64 totalNumInitialSamples = 10000;

    u64 totalNumInitialSamples = 5000000;

    time_t start = time(NULL);
    srand(time(NULL));
    randomUtilRandomize();

    lweInstance lwe;
    int n, ret;

#ifdef UNLIMITED_SAMPLES
    n = 40;
    int q = 1601;
//  int q = 101;
    double alpha = 0.02;

    lweInit(&lwe, n, q, alpha);
    const char *unsortedFolderName = LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A "/tu_darmstadt_40_1601_005_unsorted";
    const char *sortedFolderName = LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_B "/tu_darmstadt_40_1601_005_step_00";
    testCreateNewInstanceFolder(unsortedFolderName, n, q, alpha);
    newStorageFolder(&lwe, unsortedFolderName, n, q, alpha);
    addSamplesToSampleFile(unsortedFolderName, totalNumInitialSamples, start);
#else
    const char *srcFileName = "../tu_darmstadt_instances/LWE_40_005.txt";
    const char *unsortedTransformedFolderName = LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A "/tu_darmstadt_40_1601_005_unsorted_transformed";
    const char *sortedTransformedFolderName = LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_B "/tu_darmstadt_40_1601_005_step_00";

    /* convert TU Darmstadt (lwe) problem instance to local format */
    timeStamp(start);
    printf("converting tu darmstadt instance to local format\n");
    ret = tuDarmstadtFileFormatConversionWithErrorChecking(srcFileName, unsortedTransformedFolderName, 1, totalNumInitialSamples, start);

    switch (ret)
    {
    case 0: /* file conversion ok */
        printSampleVerificationOfUnsortedFolder(unsortedTransformedFolderName, start); /* verify unsorted samples */
        break;
    case 1: /* file conversion unnecessary (destination folder already exists) */
        timeStamp(start);
        printf("skipping, destination folder %s already exists\n\n", unsortedTransformedFolderName);
        break;
    default:
        timeStamp(start);
        printf("error %d when converting TU Darmstadt instance\n", ret);
        exit(1);
    }
#endif

    /* set bkw step parameters */
    /* the number of reduction steps is set above in the macro NUM_REDUCTION_STEPS.
       reduction step i uses bkwStepPar[i] to define how to produce new samples (plain bkw, lms,...).
       in addition to this, bkwStepPar[i+1] is used to determine how to sort the resulting samples.
       different reduction types may be mixed arbitrarily.
       for example, while one reduction step may use plain bkw, the next may use lms. */
    bkwStepParameters bkwStepPar[NUM_REDUCTION_STEPS + 1];

#if 1 /* plain BKW, 2 positions */
    for (int i=0; i<NUM_REDUCTION_STEPS+1; i++)
    {
        bkwStepPar[i].sorting = plainBKW;
        bkwStepPar[i].startIndex = i == 0 ? 0 : bkwStepPar[i-1].startIndex + bkwStepPar[i-1].numPositions;
        bkwStepPar[i].numPositions = 2;
        bkwStepPar[i].selection = LF2;
    }
#endif

#if 0 /* plain BKW, 3 positions */
    for (int i=0; i<NUM_REDUCTION_STEPS+1; i++)
    {
        bkwStepPar[i].sorting = plainBKW;
        bkwStepPar[i].startIndex = i == 0 ? 0 : bkwStepPar[i-1].startIndex + bkwStepPar[i-1].numPositions;
        bkwStepPar[i].numPositions = 3;
        bkwStepPar[i].selection = LF2;
    }
#endif
#if 0 /* LMS */
    for (int i=0; i<NUM_REDUCTION_STEPS+1; i++)
    {
        bkwStepPar[i].sorting = LMS;
        bkwStepPar[i].startIndex = i == 0 ? 0 : bkwStepPar[i-1].startIndex + bkwStepPar[i-1].numPositions;
        bkwStepPar[i].numPositions = 3;
        bkwStepPar[i].selection = LF2;
        bkwStepPar[i].sortingPar.LMS.p = 5;
    }
#endif
#if 0 /* codedBKW(2,1), 2 positions */
    for (int i=0; i<NUM_REDUCTION_STEPS+1; i++)
    {
        bkwStepPar[i].sorting = codedBKW;
        bkwStepPar[i].startIndex = i == 0 ? 0 : bkwStepPar[i-1].startIndex + bkwStepPar[i-1].numPositions;
        bkwStepPar[i].numPositions = 2;
        bkwStepPar[i].selection = LF2;
        bkwStepPar[i].sortingPar.CodedBKW.ct = blockCode_21;
    }
#elif 0 /* codedBKW(3,1), 3 positions */
    for (int i=0; i<NUM_REDUCTION_STEPS+1; i++)
    {
        bkwStepPar[i].sorting = codedBKW;
        bkwStepPar[i].startIndex = i == 0 ? 0 : bkwStepPar[i-1].startIndex + bkwStepPar[i-1].numPositions;
        bkwStepPar[i].numPositions = 3;
        bkwStepPar[i].selection = LF2;
        bkwStepPar[i].sortingPar.CodedBKW.ct = blockCode_31;
    }
#elif 0 /* codedBKW(4,1), 4 positions */
    for (int i=0; i<NUM_REDUCTION_STEPS+1; i++)
    {
        bkwStepPar[i].sorting = codedBKW;
        bkwStepPar[i].startIndex = i == 0 ? 0 : bkwStepPar[i-1].startIndex + bkwStepPar[i-1].numPositions;
        bkwStepPar[i].numPositions = 4;
        bkwStepPar[i].selection = LF2;
        bkwStepPar[i].sortingPar.CodedBKW.ct = blockCode_41;
    }
#elif 0 /* concatenated [2,1][2,1] code, 4 positions, concatenated */
    for (int i=0; i<NUM_REDUCTION_STEPS+1; i++)
    {
        bkwStepPar[i].sorting = codedBKW;
        bkwStepPar[i].startIndex = i == 0 ? 0 : bkwStepPar[i-1].startIndex + bkwStepPar[i-1].numPositions;
        bkwStepPar[i].numPositions = 4;
        bkwStepPar[i].selection = LF2;
        bkwStepPar[i].sortingPar.CodedBKW.ct = concatenatedCode_21_21;
    }
#endif

    /* sort (unsorted) samples */
    timeStamp(start);
    printf("sorting initial samples\n");
    u64 minDestinationStorageCapacityInSamples = totalNumInitialSamples * 4 / 3; /* add about 25% storage room for sorted samples */

#ifdef UNLIMITED_SAMPLES
    const char *src = unsortedFolderName;
    const char *dst = sortedFolderName;
#else
    const char *src = unsortedTransformedFolderName;
    const char *dst = sortedTransformedFolderName;
#endif
    ret = transition_unsorted_2_sorted(src, dst, minDestinationStorageCapacityInSamples, &bkwStepPar[0], start);
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


    /* execute several consecutive bkw steps */
    char srcFolderName[256];
    char dstFolderName[256];
    int numReductionSteps = NUM_REDUCTION_STEPS;
    for (int i=0; i<numReductionSteps; i++)
    {

        /* process BKW step */
        timeStamp(start);
        printf("Reduction step %02d -> %02d, %s reduction at positions %d to %d (destination sorting using positions %d to %d)\n", i, i+1, sortingAsString(bkwStepPar[i+1].sorting), bkwStepPar[i].startIndex, bkwStepPar[i].startIndex + bkwStepPar[i].numPositions - 1, bkwStepPar[i+1].startIndex, bkwStepPar[i+1].startIndex + bkwStepPar[i+1].numPositions - 1);

        int ret;
#ifdef DELETE_FOLDERS_FOR_PREVIOUS_BKW_STEPS
        /* delete old folder */

        if (i == 0)
        {
#ifdef UNLIMITED_SAMPLES
            ret = deleteStorageFolder(unsortedFolderName); /* deletes folder, including parameter and samples file */
            timeStamp(start);
            printf("Deleted folder %s\n", unsortedFolderName);
#else
            ret = deleteStorageFolder(unsortedTransformedFolderName); /* deletes folder, including parameter and samples file */
            timeStamp(start);
            printf("Deleted folder %s\n", unsortedTransformedFolderName);
#endif

        }
        else
        {
            char deleteFolderName[256];
            sprintf(deleteFolderName, i & 1 ? (LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_B "/tu_darmstadt_40_1601_005_step_%02d") : (LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A "/tu_darmstadt_40_1601_005_step_%02d"), i-1);
            ret = deleteStorageFolder(deleteFolderName); /* deletes folder, including parameter and samples file */
            timeStamp(start);
            printf("Deleted folder %s\n", deleteFolderName);
        }
        if (ret)
        {
            timeStamp(start);
            printf("*** Error code %d returned from deleteStorageFolder\n", ret);
        }
#endif

        /* define src and dst folders */
        sprintf(srcFolderName, i & 1 ? (LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A "/tu_darmstadt_40_1601_005_step_%02d") : (LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_B "/tu_darmstadt_40_1601_005_step_%02d"), i);
        sprintf(dstFolderName, i & 1 ? (LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_B "/tu_darmstadt_40_1601_005_step_%02d") : (LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A "/tu_darmstadt_40_1601_005_step_%02d"), i+1);

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
    sprintf(srcFolderName, numReductionSteps & 1 ? (LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A "/tu_darmstadt_40_1601_005_step_%02d") : (LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_B "/tu_darmstadt_40_1601_005_step_%02d"), numReductionSteps);
    int paramsReadError = lweParametersFromFile(&lwe, srcFolderName);
    if (paramsReadError)
    {
        timeStamp(start);
        printf("Exiting, could not read lwe parameters from %s\n", srcFolderName);
        exit(1);
    }
    n = lwe.n;
    int *solution = MALLOC(n * sizeof(int));
    MEMSET(solution, 0, n * sizeof(int));

#ifndef UNLIMITED_SAMPLES
    /* known solution for first test instance */
    short s[] = {451, 410, 1055, 253, 1038, 456, 1300, 406, 1586, 134, 1260, 5, 454, 1247, 1069, 800, 444, 199, 1093, 135, 1175, 676, 1374, 69, 596, 871, 770, 261, 658, 1575, 657, 1370, 601, 1492, 913, 516, 138, 623, 35, 1550};
    for (int i=0; i<n; i++)
    {
        lwe.s[i] = s[i];
    }

    /* apply initial transformation to known s */
    transformSecret(&lwe, s);
    timeStamp(start);
    printf("transformed s\n");
    printSolutionVector(s, n, start);

#endif // UNLIMITED_SAMPLES

    const int numPositionsToGuessAtEnd = MAX_NUM_FFT_POSITIONS;
    numSolvedCoordinates = n - (bkwStepPar[NUM_REDUCTION_STEPS-1].startIndex + bkwStepPar[NUM_REDUCTION_STEPS-1].numPositions + numPositionsToGuessAtEnd);

    /* provide partial (transformed) solution */
#ifdef UNLIMITED_SAMPLES
    for (int i=0; i<numSolvedCoordinates; i++)
    {
        solution[n - 1 - i] = lwe.s[n - 1 - i];
    }
    int s[n];
    for (int i = 0; i < n; i++)
    {
        s[i] = lwe.s[i];
    }
#else
    for (int i=0; i<numSolvedCoordinates; i++)
    {
        solution[n - 1 - i] = s[n - 1 - i];
    }
#endif // UNLIMITED_SAMPLES

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
        sprintf(srcFolderName, (numReductionSteps - i) & 1 ? LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_A "/tu_darmstadt_40_1601_005_step_%02d" : LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX_B "/tu_darmstadt_40_1601_005_step_%02d", numReductionSteps - i);
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



#ifndef UNLIMITED_SAMPLES
    timeStamp(start);
    printf("recovered transformed s\n");
    printSolutionVector(solution, n, start);
    inverseTransformSecret(&lwe, solution);
#endif

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
    timeStamp(start);
    printf("Recovered solution is %s!\n", correctly_solved ? "correct" : "INCORRECT");

    FREE(solution);

    printf("\n");
    timeStamp(start);
    printf("Done!\n");
//  pressAnyKey();
    return 0;
}








#if 0 /* create new instance */
const char *folderName = LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX "darmstadt_40_1601_005_unsorted";
testCreateNewInstanceFolder(folderName, n, q, alpha);
#endif
#if 0 /* add random samples to unsorted sample file */
u64 numSamplesToAdd = 10000000;
//  u64 numSamplesToAdd = 20000000;
//  u64 numSamplesToAdd = 100000000;
const char *folderName = LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX "darmstadt_40_1601_005_unsorted";
testAddSamplesToSampleFile(folderName, numSamplesToAdd, start);
#endif
#if 0 /* check error distribution */
char srcFolderName[256];
sprintf(srcFolderName, LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX "darmstadt_40_1601_005_sorted_plain_bkw_2_step_%02d_transformed", 19);
u64 *h = (u64*)MALLOC(q * sizeof(u64));
MEMSET(h, 0, q * sizeof(u64));
iterate_over_samples(srcFolderName, accumulateSampleErrors, h);
printErrorDistribution(h, q);
FREE(h);
#endif

#if 0
lweInstance lwe;
#if 0
int n = 40;
int q = 1601;
double alpha = 0.005;
lweInit(&lwe, n, q, alpha);
newLweInstanceFolder(&lwe, folderName, n, q, alpha);
#else
if (lweParametersFromFile(&lwe, folderName))
{
    printf("Could not read parameters from file!\n");
    return 1;
}
#endif
//  addSamples(folderName, 60000000);
//  printSampleData(&lwe, folderName, start);
testPositionValue2IndexBWK(lwe.q);
lweDestroy(&lwe);
#endif
#if 0
//  exec_python_code("print('x')\nprint('y')\nprint(1+1)\nprint(13+19)");
n = 15;
char code[1024*10];
const char *code_format =
    "import numpy as np\n"
    "import time\n"
    "from sympy import Matrix\n"
    "q=%d\n"
    "n=%d\n"
    "print('q =',q)\n"
    "print('n =',n)\n"
    "M = Matrix(%d,%d,lambda i,j: np.random.randint(0,%d-1))\n"
    "print('M =')\n"
    "print(np.matrix(M))\n"
    "Mdet = M.det()\n"
    "print('det(M) =', Mdet)\n"
    "t1=time.time()\n"
    "Minv = M.inv_mod(%d)\n"
    "t2=time.time()\n"
    "print('Minv =')\n"
    "print(np.matrix(Minv))\n"
    "MM = (M * Minv).applyfunc(lambda x: x %% %d)\n"
    "print('MM =')\n"
    "print(np.matrix(MM))\n"
    "print()\n"
    "if (Matrix.eye(%d) == MM):\n"
    "  print('YES - Inverse found')\n"
    "else:\n"
    "  print('NO - No inverse found')\n"
    "print(t2-t1,'seconds to invert')\n"
    ;
sprintf(code, code_format, q, n, n, n, q, q, q, n);
exec_python_code(code);
#endif
#if 0
n = 1000; /* about 3.8 seconds for n = 1000 on Ericsson laptop */
lweInstance lwe;
lweInit(&lwe, n, q, alpha);
test_modular_matrix_inversion(&lwe, lwe.n + 10);
#endif
#if 0
char srcFolderName[256];
sprintf(srcFolderName, LOCAL_SIMULATION_DIRECTORY_PATH_PREFIX "darmstadt_20_1601_005_sorted_plain_bkw_2_step_09");
test_fft_pruned();
#endif
#if 0
void accumulateSampleErrors(void *p, lweSample *sample, u64 categoryIndex, u64 sampleIndexInCategory, u64 numSamplesInCategory)
{
    u64 *h = (u64*)p;
    short err = error(sample);
    ASSERT(err < 1601, "*** error too large");
    h[err]++;
}

/*
void exec_python_code(const char *code) {
  Py_Initialize();
  PyRun_SimpleString(code);
  Py_Finalize();
}
*/
#endif

