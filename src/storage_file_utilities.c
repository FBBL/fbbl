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

#include "storage_file_utilities.h"

#include "config_compiler.h"
#if defined(GCC)
#include <unistd.h>
#include <stdio.h>
#include <dirent.h>
#else
#include <dir.h>
#endif

#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include <inttypes.h>
#include "log_utils.h"
#include "memory_utils.h"
#include "string_utils.h"
#include "bkw_step_parameters.h"
#include "linear_algebra_modular.h"

#define _FILE_OFFSET_BITS 64

#ifndef _FILE_OFFSET_BITS
#error _FILE_OFFSET_BITS not defined
/* compile project with -D_FILE_OFFSET_BITS=64 */
#if (_FILE_OFFSET_BITS != 64)
#error _FILE_OFFSET_BITS not set to 64
#endif
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

/* name of parameter file */
static const char *par_file_name = "params.txt";

/* name of samples file */
static const char *sam_file_name = "samples.dat";

/* name of samples info file */
static const char *sam_info_file_name = "samples_info.txt";

void parameterFileName(char *paramFileName, const char *folderName)
{
    sprintf(paramFileName, "%s/%s", folderName, par_file_name);
}

void samplesFileName(char *samplesFileName, const char *folderName)
{
    sprintf(samplesFileName, "%s/%s", folderName, sam_file_name);
}

void samplesInfoFileName(char *samplesInfoFileName, const char *folderName)
{
    sprintf(samplesInfoFileName, "%s/%s", folderName, sam_info_file_name);
}

/* writes (lwe) problem parameters to file */
int parametersToFile(lweInstance *lwe, const char *folderName)
{
    char fileName[512];
    parameterFileName(fileName, folderName);

    FILE *f = fopen(fileName, "w");
    if (!f)
    {
        return 1;
    }
    fprintf(f, "n = %d\n", lwe->n);
    fprintf(f, "q = %d\n", lwe->q);
    fprintf(f, "alpha = %lf\n", lwe->alpha);
    fprintf(f, "sigma = %lf\n", lwe->sigma);
    fprintf(f, "rnd_ctx = (%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%d)\n", lwe->rnd.A1, lwe->rnd.A2, lwe->rnd.B1, lwe->rnd.B2, lwe->rnd.C1, lwe->rnd.C2, lwe->rnd.initialized);
    fprintf(f, "s = (%hi", lwe->s[0]);
    for (int i=1; i<lwe->n; i++)
    {
        fprintf(f, ",%hi", lwe->s[i]);
    }
    fprintf(f, ")\n");
    if (lwe->A && lwe->A_inverse && lwe->b)
    {
        fprintf(f, "A =\n");
        for (int i=0; i<lwe->n; i++)
        {
            fprintf(f, "%hi", lwe->A[i][0]);
            for (int j=1; j<lwe->n; j++)
            {
                fprintf(f, " %hi", lwe->A[i][j]);
            }
            fprintf(f, "\n");
        }
        fprintf(f, "A_inverse =\n");
        for (int i=0; i<lwe->n; i++)
        {
            fprintf(f, "%hi", lwe->A_inverse[i][0]);
            for (int j=1; j<lwe->n; j++)
            {
                fprintf(f, " %hi", lwe->A_inverse[i][j]);
            }
            fprintf(f, "\n");
        }
        fprintf(f, "b = (%hi", lwe->b[0]);
        for (int i=1; i<lwe->n; i++)
        {
            fprintf(f, ",%hi", lwe->b[i]);
        }
        fprintf(f, ")\n");
    }
    fclose(f);
    return 0;
}

/* reads (lwe) problem parameters from file */
int lweParametersFromFile(lweInstance *lwe, const char *folderName)
{
    int n, q;
    double alpha, sigma;
    char fileName[512];
    parameterFileName(fileName, folderName);
    FILE *f = fopen(fileName, "r");
    if (!f)
    {
        return 1;
    }
    if(!fscanf(f, "n = %d\n", &n))
        return 1;
    if(!fscanf(f, "q = %d\n", &q))
        return 1;
    if(!fscanf(f, "alpha = %lf\n", &alpha))
        return 1;
    if(!fscanf(f, "sigma = %lf\n", &sigma))
        return 1;
    lweInit(lwe, n, q, alpha);
    if(!fscanf(f, "rnd_ctx = (%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%d)\n", &lwe->rnd.A1, &lwe->rnd.A2, &lwe->rnd.B1, &lwe->rnd.B2, &lwe->rnd.C1, &lwe->rnd.C2, &lwe->rnd.initialized))
        return 1;
    if(!fscanf(f, "s = (%hi", &lwe->s[0]))
        return 1;
    for (int i=1; i<lwe->n; i++)
    {
        if(!fscanf(f, ",%hi", &lwe->s[i]))
            return 1;
    }
    if(fscanf(f, ")\n"))
        return 1;
    /* check if initial transformation matrix has been written to file */
//  int ret = fscanf(f, "A =\n");
    char line[100];
    if(fgets(line, 100, f))
    {
        if (line[0] != 'A')
        {
            //  if (ret || feof(f)) {
            fclose(f);
            return 0; /* no A and A_inverse not present in parameter file */
        }
    }
    else
    {
        return 0; /* no A and A_inverse not present in parameter file */
    }

    lweInstanceAllocateLinearTransformationMatrices(lwe);
    /* read A */
    for (int i=0; i<lwe->n; i++)
    {
        if(!fscanf(f, "%hi", &lwe->A[i][0]))
            return 1;
        for (int j=1; j<lwe->n; j++)
        {
            if(!fscanf(f, " %hi", &lwe->A[i][j]))
                return 1;
        }
        if(fscanf(f, "\n"))
            return 1;
    }
    /* read A_inverse */
    if(fscanf(f, "A_inverse =\n"))
        return 1;
    for (int i=0; i<lwe->n; i++)
    {
        if(!fscanf(f, "%hi", &lwe->A_inverse[i][0]))
            return 1;
        for (int j=1; j<lwe->n; j++)
        {
            if(!fscanf(f, " %hi", &lwe->A_inverse[i][j]))
                return 1;
        }
        if(fscanf(f, "\n"))
            return 1;
    }
    /* read b */
    if(!fscanf(f, "b = (%hi", &lwe->b[0]))
        return 1;
    for (int i=1; i<lwe->n; i++)
    {
        if(!fscanf(f, ",%hi", &lwe->b[i]))
            return 1;
    }
    if(fscanf(f, "\n"))
        return 1;
    fclose(f);
    return 0;
}

/* write sample information to file */
int sampleInfoToFile(const char *folderName, bkwStepParameters *bkwStepPar, u64 numCategories, u64 categoryCapacity, u64 numTotalSamples, u64 *numSamplesPerCategory)
{
    char fileName[512];
    samplesInfoFileName(fileName, folderName);
    FILE *f = fopen(fileName, "w");
    if (!f)
    {
        return 1; /* could not open sample info file */
    }
    char str[512];
    fprintf(f, "sorting = %s\n", bkwStepParametersAsString(str, bkwStepPar)); /* write sorting and additional parameters */
    fprintf(f, "num categories = %" PRIu64 "\n", numCategories);
    fprintf(f, "category capacity (num samples) = %" PRIu64 "\n", categoryCapacity);
    fprintf(f, "total num samples stored = %" PRIu64 "\n", numTotalSamples);
    fprintf(f, "num samples per category = (%" PRIu64, numSamplesPerCategory[0]);
    for (int i=1; i<numCategories; i++)
    {
        fprintf(f, ",%" PRIu64, numSamplesPerCategory[i]);
    }
    fprintf(f, ")\n");
    fclose(f);
    return 0;
}

/* read sample information from file */
/* note: last parameter is optional by passing NULL */
int sampleInfoFromFile(const char *folderName, bkwStepParameters *bkwStepPar, u64 *numCategories, u64 *categoryCapacity, u64 *numTotalSamples, u64 *numSamplesPerCategory)
{
    ASSERT(folderName, "unexpected parameter");
    /* get q from lwe parameters from file */
    lweInstance lwe;
    lweParametersFromFile(&lwe, folderName);
    /* open sample info file */
    char fileName[512];
    samplesInfoFileName(fileName, folderName);
    FILE *f = fopen(fileName, "r");
    if (!f)
    {
        lweDestroy(&lwe);
        return 1; /* could not open sample info file */
    }

    /* read sample info from file */
    u64 dummy, numCat;

    /* sorting and bkw step parameters */
    char sortingString[512];
    if(!fscanf(f, "sorting = %[^\n]\n", sortingString))  /* read text to end of line to include spaces */
    {
        lweDestroy(&lwe);
        return 2; /* could not determine sorting method */
    }

    if (bkwStepPar && !bkwStepParametersFromString(sortingString, bkwStepPar))
    {
        lweDestroy(&lwe);
        ASSERT_ALWAYS("sorting not determined");
        return 2; /* could not determine sorting method */
    }

    /* number of categories */
    if(!fscanf(f, "num categories = %" PRIu64 "\n", &numCat))
    {
        lweDestroy(&lwe);
        return 3;
    }
    if (numCategories)
    {
        *numCategories = numCat;
    }
    if (bkwStepPar && numCategories)
    {
        // u64 nc = num_categories(&lwe, bkwStepPar);
//    if (*numCategories != nc) {
//      printf("Unexpected number of categories, numCategories=%" PRIu64 ", num_categories=%" PRIu64 "\n", *numCategories, nc);
//    }
        ASSERT(*numCategories == nc, "unexpected number of categories");
    }

    /* category capacity */
    if(!fscanf(f, "category capacity (num samples) = %" PRIu64 "\n", categoryCapacity ? categoryCapacity : &dummy))
    {
        lweDestroy(&lwe);
        return 3;
    }

    /* total number of samples stored */
    if(!fscanf(f, "total num samples stored = %" PRIu64 "\n", numTotalSamples ? numTotalSamples : &dummy))
    {
        lweDestroy(&lwe);
        return 3;
    }

    /* number of samples per category */
    if (numSamplesPerCategory)
    {
        if(!fscanf(f, "num samples per category = (%" PRIu64 "", &numSamplesPerCategory[0]))
        {
            lweDestroy(&lwe);
            return 3;
        }
        for (int i=1; i<numCat; i++)
        {
            if(!fscanf(f, ",%" PRIu64 "", &numSamplesPerCategory[i]))
            {
                lweDestroy(&lwe);
                return 3;
            }
        }
    }

    fclose(f);
    lweDestroy(&lwe);
    return 0;
}

int newStorageFolder(lweInstance *lwe, const char *folderName, int n, int q, double alpha)
{
    lweInit(lwe, n, q, alpha); /* create lwe instance with given parameters */
    return newStorageFolderWithGivenLweInstance(lwe, folderName);
}

int newStorageFolderWithGivenLweInstance(lweInstance *lwe, const char *folderName)
{
#if defined(GCC)
    if (mkdir(folderName, 0777) == -1)   /* create lwe instance folder */
    {
#else
    if (mkdir(folderName) == -1)   /* create lwe instance folder */
    {
#endif
        return 1; /* directory already exists */
    }
    if (parametersToFile(lwe, folderName))   /* create parameter file */
    {
        return 2; /* could not create parameter file */
    }
    /* create empty samples file */
    char sFileName[512];
    samplesFileName(sFileName, folderName);
    FILE *f = fopen(sFileName, "wb"); /* create file */
    if (!f)
    {
        return 3; /* could not create samples file */
    }
    fclose(f);
    return 0;
}

int deleteStorageFolder(const char *folderName, int deleteParFile, int deleteSampleInfoFile, int deleteSamples)
{
    int ret = 0;
    char fileName[512];
    if(deleteParFile)
    {
        parameterFileName(fileName, folderName);
        if (remove(fileName))   /* delete parameter file */
        {
            ret |= 1;
        }
    }
    if(deleteSampleInfoFile)
    {
        samplesInfoFileName(fileName, folderName);
        if (remove(fileName))   /* delete samples info file */
        {
            /* folder with unsorted samples has no samples info file, so do not report failure to delete as error */
            //    ret |= 2;
        }
    }
    if(deleteSamples)
    {
        samplesFileName(fileName, folderName);
        if (remove(fileName))   /* delete samples file */
        {
            ret |= 4;
        }
    }
    if(deleteParFile && deleteSampleInfoFile && deleteSamples)
    {
        if (rmdir(folderName))   /* delete folder */
        {
            ret |= 8;
        }
    }
    return ret;
}

/* opens samples file */
FILE *fopenSamples(const char *folderName, const char *mode)
{
    char sFileName[512];
    samplesFileName(sFileName, folderName);
    return fopen(sFileName, mode);
}

/* read sample range from current position into buffer (does not close file) */
u64 freadSamples(FILE *f, lweSample *sampleBuf, u64 numSamples)
{
    u64 ret = fread(sampleBuf, LWE_SAMPLE_SIZE_IN_BYTES, numSamples, f);
    return ret >= 0 ? ret : 0;
}

/* read category range from current position into buffer (does not close file) */
u64 freadCategories(FILE *f, lweSample *sampleBuf, u64 numCategories, u64 categoryCapacityInSamples)
{
    return freadSamples(f, sampleBuf, numCategories * categoryCapacityInSamples);
}

/* read sample range into buffer (does not close file) */
u64 fsetAndReadSamples(FILE *f, lweSample *sampleBuf, u64 startingSample, u64 nbrOfSamples)
{
    if (startingSample < 0)
    {
        return 0;
    }
    fseeko(f, startingSample * LWE_SAMPLE_SIZE_IN_BYTES, SEEK_SET);
    return freadSamples(f, sampleBuf, nbrOfSamples);
}

u64 sampleFileSizeInBytes(const char *folderName)
{
#if 1
    FILE *f = fopenSamples(folderName, "rb");
    if (!f)
    {
        return 0;
    }
    fseeko(f, 0L, SEEK_END); /* not necessarily supported for binary files or large files, but seems to work with MinGW64 and Windows */
    u64 sz = ftello(f);
    fclose(f);
    return sz;
#else
    struct stat stbuf;
    char sFileName[512];

    samplesFileName(sFileName, folderName);
    int fd = open(sFileName, O_RDONLY);
    if (fd == -1)
    {
        return 0;
    }
    if ((fstat(fd, &stbuf) != 0) || (!S_ISREG(stbuf.st_mode)))
    {
        return 0;
    }
    return stbuf.st_size;
#endif
}

u64 numSamplesInSampleFile(const char *folderName)
{
    u64 sz = sampleFileSizeInBytes(folderName);
    ASSERT(sz % LWE_SAMPLE_SIZE_IN_BYTES == 0, "Non-integral number of samples detected!\n");
    return sz / (u64)LWE_SAMPLE_SIZE_IN_BYTES;
}

/* add samples to samples file */
u64 addSamplesToSampleFile(const char *folderName, u64 nbrOfSamples, time_t start)
{
    u64 n = nbrOfSamples;
    int blockSize = 1000000;
    lweInstance lwe;
    if (lweParametersFromFile(&lwe, folderName))
    {
        return 0;
    }

    FILE *f = fopenSamples(folderName, "ab");
    if (!f)
    {
        return 0;
    }

    lweSample *sampleBuf = MALLOC(blockSize * LWE_SAMPLE_SIZE_IN_BYTES); // A buffer for LWE samples
    ASSERT(sampleBuf, "Allocation failed!\n");

    //print every 1/10
    u64 level = (nbrOfSamples*9)/10;
    while (n)
    {
        int numSamplesThisRound = MIN(n, blockSize);
        for (int i=0; i<numSamplesThisRound; i++)
        {
            lwe.newInPlaceRandomSample(&sampleBuf[i], lwe.n, lwe.q, lwe.sigma, &lwe.rnd, lwe.s);
        }
        int numWritten = fwrite(sampleBuf, LWE_SAMPLE_SIZE_IN_BYTES, numSamplesThisRound, f);
        n -= numWritten;
        if(n < level)
        {
            timeStamp(start);
            printf("Written %5.2f%% of samples\n", 100*(nbrOfSamples-n)/(double)nbrOfSamples);
            level -= nbrOfSamples/10;
        }

    }
    timeStamp(start);
    printf("Written %5.2f%% of samples\n", 100*(double)(nbrOfSamples-n)/(double)nbrOfSamples);
    fclose(f);
    FREE(sampleBuf);
    return nbrOfSamples;
}

/* read sample range into buffer (closes file) */
u64 readSamplesFromSampleFile(lweSample *sampleBuf, const char *folderName, u64 startingSample, u64 nbrOfSamples)
{
    if (startingSample < 0)
    {
        return 0;
    }
    FILE *f = fopenSamples(folderName, "rb");
    if (!f)
    {
        return 0;
    }
    int ret = fsetAndReadSamples(f, sampleBuf, startingSample, nbrOfSamples);
    fclose(f);
    return ret;
}

/* transform sample in-place */
/* (a,b) -> (a',b') = ( -A^{-1}a, b - <A^{-1}a,b> ) */
/* (a,b) = (\vector{a},\scalar{b}) -> ( -A^{-1}\vector{a}, \scalar{b} - <A^{-1}\vector{a},\vector{b}> ), where \vector{b} is the b-vector from the initial transform lwe->b */
static void transformSampleInPlace(lweInstance *lwe, lweSample *sample)
{
    int n = lwe->n;
    int q = lwe->q;
    short **Ainv = lwe->A_inverse;
    lweSample s;
    MEMCPY(s.col.a, sample->col.a, n * sizeof(short)); /* copy a to temp */
    s.sumWithError = sample->sumWithError;
    /* compute transformed a-vector and b-value */
    u64 sumb = 0;
    for (int i=0; i<n; i++)
    {
        u64 suma = 0;
        for (int j=0; j<n; j++)
        {
//      suma += Ainv[i][j] * columnValue(&s, j);
            suma += Ainv[j][i] * columnValue(&s, j); /* transpose */
        }
        suma = suma % q;
//    sample->col.a[i] = suma == 0 ? 0 : q - suma;
        sample->col.a[i] = suma;
        ASSERT(0 <= sample->col.a[i], "unexpected value");
        ASSERT(sample->col.a[i] < q, "unexpected value");
        sumb += sample->col.a[i] * lwe->b[i];
    }
    sumb = sumb % q;
    sample->sumWithError = (q + s.sumWithError - sumb) % q;
    ASSERT(0 <= sample->sumWithError, "unexpected value");
    ASSERT(sample->sumWithError < q, "unexpected value");
    sample->error = -1; /* unknown error */
    sample->col.hash = bkwColumnComputeHash(sample, n, 0); /* compute column hash */
}

/*
static void combineTwoSamplesAdd(lweInstance *lwe, lweSample *dst, lweSample *sample1, lweSample *sample2)
{
    int q = lwe->q;
    for (int i=0; i<lwe->n; i++)
    {
        dst->col.a[i] = (columnValue(sample1, i) + columnValue(sample2, i)) % q;
    }
    dst->col.hash = bkwColumnComputeHash(dst, lwe->n, 0);
    if (sample1->error == -1 || sample2->error == -1)   // if either error term is undefined
    {
        dst->error = -1; // resulting sum of error terms is also undefined
    }
    else
    {
        dst->error = (sample1->error + sample2->error) % q;
    }
    dst->sumWithError = (sample1->sumWithError + sample2->sumWithError) % q;
}
*/

static void combineTwoSamplesSub(lweInstance *lwe, lweSample *dst, lweSample *sample1, lweSample *sample2)
{
    int q = lwe->q;
    for (int i=0; i<lwe->n; i++)
    {
        dst->col.a[i] = (q + columnValue(sample1, i) - columnValue(sample2, i)) % q;
    }
    dst->col.hash = bkwColumnComputeHash(dst, lwe->n, 0);
    if (sample1->error == -1 || sample2->error == -1)   /* if either error term is undefined */
    {
        dst->error = -1; /* resulting sum of error terms is also undefined */
    }
    else
    {
        dst->error = (q + sample1->error - sample2->error) % q;
    }
    dst->sumWithError = (q + sample1->sumWithError - sample2->sumWithError) % q;
}

static void combineThreeSamplesAddAdd(lweInstance *lwe, lweSample *dst, lweSample *sample1, lweSample *sample2, lweSample *sample3)
{
    int q = lwe->q;
    for (int i=0; i<lwe->n; i++)
    {
        dst->col.a[i] = (columnValue(sample1, i) + columnValue(sample2, i) + columnValue(sample3, i)) % q;
    }
    dst->col.hash = bkwColumnComputeHash(dst, lwe->n, 0);
    if (sample1->error == -1 || sample2->error == -1 || sample3->error == -1)   /* if either error term is undefined */
    {
        dst->error = -1; /* resulting sum of error terms is also undefined */
    }
    else
    {
        dst->error = (sample1->error + sample2->error + sample3->error) % q;
    }
    dst->sumWithError = (sample1->sumWithError + sample2->sumWithError + sample3->sumWithError) % q;
}

static void combineThreeSamplesAddSub(lweInstance *lwe, lweSample *dst, lweSample *sample1, lweSample *sample2, lweSample *sample3)
{
    int q = lwe->q;
    for (int i=0; i<lwe->n; i++)
    {
        dst->col.a[i] = (q + columnValue(sample1, i) + columnValue(sample2, i) - columnValue(sample3, i)) % q;
    }
    dst->col.hash = bkwColumnComputeHash(dst, lwe->n, 0);
    if (sample1->error == -1 || sample2->error == -1 || sample3->error == -1)   /* if either error term is undefined */
    {
        dst->error = -1; /* resulting sum of error terms is also undefined */
    }
    else
    {
        dst->error = (q + sample1->error + sample2->error - sample3->error) % q;
    }
    dst->sumWithError = (q + sample1->sumWithError + sample2->sumWithError - sample3->sumWithError) % q;
}

static void combineThreeSamplesSubAdd(lweInstance *lwe, lweSample *dst, lweSample *sample1, lweSample *sample2, lweSample *sample3)
{
    int q = lwe->q;
    for (int i=0; i<lwe->n; i++)
    {
        dst->col.a[i] = (q + columnValue(sample1, i) - columnValue(sample2, i) + columnValue(sample3, i)) % q;
    }
    dst->col.hash = bkwColumnComputeHash(dst, lwe->n, 0);
    if (sample1->error == -1 || sample2->error == -1 || sample3->error == -1)   /* if either error term is undefined */
    {
        dst->error = -1; /* resulting sum of error terms is also undefined */
    }
    else
    {
        dst->error = (q + sample1->error - sample2->error + sample3->error) % q;
    }
    dst->sumWithError = (q + sample1->sumWithError - sample2->sumWithError + sample3->sumWithError) % q;
}

static void combineThreeSamplesSubSub(lweInstance *lwe, lweSample *dst, lweSample *sample1, lweSample *sample2, lweSample *sample3)
{
    int q = lwe->q;
    for (int i=0; i<lwe->n; i++)
    {
        dst->col.a[i] = (q + q + columnValue(sample1, i) - columnValue(sample2, i) - columnValue(sample3, i)) % q;
    }
    dst->col.hash = bkwColumnComputeHash(dst, lwe->n, 0);
    if (sample1->error == -1 || sample2->error == -1 || sample3->error == -1)   /* if either error term is undefined */
    {
        dst->error = -1; /* resulting sum of error terms is also undefined */
    }
    else
    {
        dst->error = (q + q + sample1->error - sample2->error - sample3->error) % q;
    }
    dst->sumWithError = (q + q + sample1->sumWithError - sample2->sumWithError - sample3->sumWithError) % q;
}

static void flushUnsortedSampleBuf(FILE *f, lweSample *sampleBuf, int numSamples)
{
    while (numSamples)
    {
        int numWritten = fwrite(sampleBuf, LWE_SAMPLE_SIZE_IN_BYTES, numSamples, f);
        numSamples -= numWritten;
        sampleBuf += numWritten;
    }
}

/* conversion (import) of TU Darmstadt (lwe) problem instances to local format, with optional sample amplification */
int convertTUDarmstadtProblemInstanceToNativeFormat(lweInstance *lwe, const char *srcFileName, const char *dstFolderName, int useSampleAmplification, u64 totalNumSamples)
{
    int n, q, numSamples;
    double alpha;
    FILE *f = fopen(srcFileName, "r");
    if (!f)
    {
        return 1; /* could not open source file */
    }
    if(!fscanf(f, "%d\n", &n))
        return 1;
    if(!fscanf(f, "%d\n", &numSamples))
        return 1;
    if(!fscanf(f, "%d\n", &q))
        return 1;
    if(!fscanf(f, "%lf\n", &alpha))
        return 1;
    lweInit(lwe, n, q, alpha); /* initializes lwe instance */
    for (int i=0; i<lwe->n; i++)
    {
        lwe->s[i] = -1; /* s is unknown */
    }

    lweSample *sampleBuf = MALLOC(numSamples * sizeof(lweSample));
    ASSERT(sampleBuf, "allocation error");

    /* read sample b-values */
    if(!fscanf(f, "[%hi ", &(sampleBuf[0].sumWithError)))
        return 1;
    sampleBuf[0].error = 0; /* error unknown */
    for (int i=1; i<numSamples; i++)
    {
        if(!fscanf(f, "%hi ", &(sampleBuf[i].sumWithError)))
            return 1;
        sampleBuf[i].error = 0; /* error unknown */
    }
    if(fscanf(f, "]\n"))
        return 1;

    /* read sample vectors a */
    if(fscanf(f, "["))
        return 1;
    for (int i=0; i<numSamples; i++)
    {
        if(!fscanf(f, "[%hi", &(sampleBuf[i].col.a[0])))
            return 1;
        for (int j=1; j<lwe->n; j++)
        {
            if(!fscanf(f, " %hi", &(sampleBuf[i].col.a[j])))
                return 1;
        }
        if(fscanf(f, "]\n"))
            return 1;
        sampleBuf[i].col.hash = bkwColumnComputeHash(&sampleBuf[i], lwe->n, 0);
    }
    if(fscanf(f, "]\n"))
        return 1;

    fclose(f);

#if 0 /* for testing purposes only, enable this to test overwriting with random samples */
    for (int i=0; i<lwe->n; i++)
    {
        lwe->s[i] = randomUtilInt(&lwe->rnd, lwe->q);
    }
    for (int i=0; i<numSamples; i++)
    {
        lwe->newInPlaceRandomSample(&sampleBuf[i], lwe->n, lwe->q, lwe->sigma, &lwe->rnd, lwe->s);
    }
#endif

    /* compute initial linear transformation (uses O(n\log\log{q}) samples) */
    lweInstanceAllocateLinearTransformationMatrices(lwe);
    int numSamplesUsedToComputeInitialTransformation = compute_matrix_inverse_modular(lwe, sampleBuf, numSamples, lwe->A, lwe->A_inverse, lwe->b);

    /* from here on, work with the remaining samples */
    int numRemainingSamples = numSamples - numSamplesUsedToComputeInitialTransformation;
    lweSample *remainingSampleBuf = sampleBuf + numSamplesUsedToComputeInitialTransformation;

    /* transform remaining samples */
    for (int i=0; i<numRemainingSamples; i++)
    {
        transformSampleInPlace(lwe, &remainingSampleBuf[i]);
    }

    /* create instance folder with parameter file */
    int ret = newStorageFolderWithGivenLweInstance(lwe, dstFolderName); /* creates folder with parameter file (with given parameters) and an empty samples and samples info file */
    if (ret != 0)
    {
        printf("Error in parametersFromTUDarmstadtFile, could not create instance folder (error %d returned from newStorageFolderWithGivenLweInstance)\n", ret);
        FREE(sampleBuf);
        return 2;
    }

    /* open samples file */
    FILE *f2 = fopenSamples(dstFolderName, "ab");
    if (!f2)
    {
        printf("failed to open samples file");
        FREE(sampleBuf);
        return 3;
    }

    if (!useSampleAmplification)
    {
        fclose(f2); /* close (unsorted) sample file */
        FREE(sampleBuf);
        return 0;
    }

    /* amplify samples */
    u64 numGeneratedSamples = 0;

    /* exhaust all sample triplets */
    /* quadruples or higher not implemented */
    int amplifiedSampleBufSize = 1000;
    lweSample *amplifiedSampleBuf = MALLOC(amplifiedSampleBufSize * LWE_SAMPLE_SIZE_IN_BYTES);
    int numSamplesInAmplifiedSampleBuf = 0;

    /* loop through all triplets */
    /* for dependency reasons, it is not advisable to implement triplet generation naively with a triple loop */
    /* a few vectors will then appear much more often in the resulting samples, which causes many zero vectors to be generated */
    /* in order to reduce this dependency, we spread the indices chosen for the triplets by using a maximum length lfsr */
    /* quads could be implemented in the same way (using sage to generate a primitive polynomial), but we have settled for triplets for now */
    int ii = 0, jj = 0, kk = 1; /* any starting point except all zeros */
    for (;;)
    {

        /* choose indices for the next triplet */
        /* clock lfsr until indices are good for a triplet */
        /* this could, of course, be implemented more efficiently, but it is fast enough for our purposes */
        do
        {
            int temp = (1594*ii + 1600*jj + 1600*kk) % q;
            ii = jj;
            jj = kk;
            kk = temp;
        }
        while (ii >= jj || jj >= kk || ii >= numRemainingSamples || jj >= numRemainingSamples || kk >= numRemainingSamples);

        /* now we have indices ii < jj < kk for a triplet */
        /* randomizing the signs in a_{ii} +/- a_{jj} +/- a_{kk} is also good for reducing dependencies (verified by simulating and counting number of generated zero vectors) */
        switch (randomUtilInt(&lwe->rnd, 4))
        {
        case 0:
            combineThreeSamplesAddAdd(lwe, &amplifiedSampleBuf[numSamplesInAmplifiedSampleBuf++], &remainingSampleBuf[ii], &remainingSampleBuf[jj], &remainingSampleBuf[kk]);
            break;
        case 1:
            combineThreeSamplesAddSub(lwe, &amplifiedSampleBuf[numSamplesInAmplifiedSampleBuf++], &remainingSampleBuf[ii], &remainingSampleBuf[jj], &remainingSampleBuf[kk]);
            break;
        case 2:
            combineThreeSamplesSubAdd(lwe, &amplifiedSampleBuf[numSamplesInAmplifiedSampleBuf++], &remainingSampleBuf[ii], &remainingSampleBuf[jj], &remainingSampleBuf[kk]);
            break;
        case 3:
            combineThreeSamplesSubSub(lwe, &amplifiedSampleBuf[numSamplesInAmplifiedSampleBuf++], &remainingSampleBuf[ii], &remainingSampleBuf[jj], &remainingSampleBuf[kk]);
            break;
        default:
            ASSERT_ALWAYS("*** ERROR, unhandled case");
        }
        numGeneratedSamples += 1;

        /* flush samples when buffer is full */
        if (numSamplesInAmplifiedSampleBuf >= amplifiedSampleBufSize)
        {
            flushUnsortedSampleBuf(f2, amplifiedSampleBuf, numSamplesInAmplifiedSampleBuf); /* flush samples */
            numSamplesInAmplifiedSampleBuf = 0;
        }

        /* exit if enough samples have been generated */
        if (numGeneratedSamples >= totalNumSamples)
        {
            break;
        }
    }

    /* flush buffered samples */
    flushUnsortedSampleBuf(f2, amplifiedSampleBuf, numSamplesInAmplifiedSampleBuf);

//  printf("sample amplification from quadruples not implemented, so quitting amplification here\n");

    FREE(sampleBuf);
    FREE(amplifiedSampleBuf);
    fclose(f2);

    return 0;
}

/* erik */
// Calculates the size of a file
u64 fileSize(FILE *fp)
{
    u64 cur = ftello64(fp);
    fseeko64(fp, 0L, SEEK_END);
    u64 sz = ftello64(fp);
    fseeko64(fp, cur, SEEK_SET);
    return sz;
}

// Extends a file to a specified file size
void fileExtend(FILE *fp, u64 size)
{
    fseeko64(fp, size - 1, SEEK_SET);
    fputc('\0', fp);
    fseeko64(fp, 0L, SEEK_SET);
    // rewind(fp);
}

int folderExists(const char *folderPath)
{
    struct stat info;
    if( stat( folderPath, &info ) != 0 )
    {
        return 0; /* cannot access folder */
    }
    if( info.st_mode & S_IFDIR )
    {
        return 1; /* path is a directory */
    }
    return 0; /* not a directory */
}

// check if a file exists
int fileExists (const char *Filename)
{
    if(access( Filename, F_OK ) != -1)
        return 1;
    return 0;
}
