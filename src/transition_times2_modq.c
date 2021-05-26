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

#include "transition_times2_modq.h"
#include "position_values_2_category_index.h"
#include "storage_file_utilities.h"
#include "memory_utils.h"
#include "log_utils.h"
#include "string_utils.h"
#include "storage_writer.h"
#include "bkw_step_parameters.h"
#include "verify_samples.h"
#include "test_functions.h"
#include <inttypes.h>
#include <pthread.h>

// global variables
typedef struct {
	int index;
    lweInstance *lwe;
    bkwStepParameters *bkwStepPar;
    lweSample *sampleReadBuf;
    storageWriter *sw;
    u64 min;
    u64 max;
    time_t start;
} Params;

/* define mutexes to protect common resources from concurrent access */
static pthread_mutex_t screen_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t flush_sw_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t resumeCond;
static pthread_mutex_t *storage_mutex;
static int n_storage_mutex;
static u64 max_num_categories_per_mutex;
static u64 nextPrintLimit;
static int m_SuspendFlag = 0;

void suspendWriting()
{ // tell the thread to suspend
    pthread_mutex_lock(&flush_sw_mutex);
    m_SuspendFlag = 1;
    pthread_mutex_unlock(&flush_sw_mutex);
}
void resumeWriting()
{ // tell the thread to resume
    pthread_mutex_lock(&flush_sw_mutex);
    m_SuspendFlag = 0;
    pthread_cond_broadcast(&resumeCond);
    pthread_mutex_unlock(&flush_sw_mutex);
}
void checkSuspend()
{ // if suspended, suspend until resumed
    pthread_mutex_lock(&flush_sw_mutex);
    while (m_SuspendFlag != 0) pthread_cond_wait(&resumeCond, &flush_sw_mutex);
    pthread_mutex_unlock(&flush_sw_mutex);
}

short multiply_time2_modq(short a, int q) {
    a = a << 1;
    if (a >= q)
        a -= q;
    return a;
}

int sample_times2_modq(lweInstance *lwe, lweSample *sample)
{
    int n = lwe->n;
    int q = lwe->q;
    for (int i=0; i<n; i++)
        // sample->col.a[i] = (2*sample->col.a[i]) % q;
        sample->col.a[i] = multiply_time2_modq(sample->col.a[i], q);
    sample->sumWithError = multiply_time2_modq(sample->sumWithError, q);
    sample->col.hash = bkwColumnComputeHash(sample, n, 0);
    sample->error = sample->error < 0 ? -1 : multiply_time2_modq(error(sample), q);
    return 0;
}


void *single_thread_work(void *params){

    Params *p = (Params*)params;

    char str[256];
    int mutex_index;

    pthread_mutex_lock(&screen_mutex);
    printf("screen min %d max %d\n", p->min, p->max);
    pthread_mutex_unlock(&screen_mutex);

    /* add samples to storage writer */
    for (u64 i=p->min; i<p->max; i++)
    {
        lweSample *s = &p->sampleReadBuf[i];
        sample_times2_modq(p->lwe, s);
        u64 categoryIndex = position_values_2_category_index(p->lwe, s, p->bkwStepPar);
        ASSERT(categoryIndex<numCategories, "ERROR *** invalid category index");

        u64 numIncorrectCategoryClassifications = 0;
        verifyOneSampleSorted(p->lwe, s, p->bkwStepPar, NULL, NULL, categoryIndex, &numIncorrectCategoryClassifications, 0);
        if (numIncorrectCategoryClassifications)
        {
        	pthread_mutex_lock(&screen_mutex);
            printf("****** classification error detected\n");
            pthread_mutex_unlock(&screen_mutex);
            exit(-1);
        }
        int storageWriterStatus = 0;
        checkSuspend();
        mutex_index = categoryIndex / max_num_categories_per_mutex;
        pthread_mutex_lock(&storage_mutex[mutex_index]);
        lweSample *d = storageWriterAddSample(p->sw, categoryIndex, &storageWriterStatus); /* reserve memory area in storage writer */
        pthread_mutex_unlock(&storage_mutex[mutex_index]);
        if (d)
        {
            MEMCPY(d, s, sizeof(lweSample)); /* copy sample to reserved memory area */
        }
        if (p->sw->totalNumSamplesProcessedByStorageWriter >= nextPrintLimit)
        {
        	pthread_mutex_lock(&screen_mutex);
            nextPrintLimit += 100000000;
            timeStamp(p->start);
            printf("%s samples processed by storage writer so far\n", sprintf_u64_delim(str, p->sw->totalNumSamplesProcessedByStorageWriter));
        	pthread_mutex_unlock(&screen_mutex);
        }
        switch (storageWriterStatus)
        {
        case 0: /* insertion succeeded */
            verifyOneSampleSorted(p->lwe, d, p->bkwStepPar, NULL, NULL, categoryIndex, &numIncorrectCategoryClassifications, 0);
            if (numIncorrectCategoryClassifications)
            {
            	pthread_mutex_lock(&screen_mutex);
                printf("****** classification error detected (2)");
            	pthread_mutex_unlock(&screen_mutex);
            	exit(-1);
            }
            break;
        case 1: /* insertion succeeded but this was the last available slot in the cache */
            /* so it makes sense to flush the cache to file here, but only if the cache is actually smaller than the storage on file. */
            suspendWriting();
            if (p->sw->categoryCapacityBuf < p->sw->categoryCapacityFile)   /* if buffer storage is smaller than file storage */
            {
                double cacheLoad = storageWriterCurrentLoadPercentageCache(p->sw);
                if (cacheLoad >= MIN_STORAGE_WRITER_CACHE_LOAD_PERCENTAGE_BEFORE_FLUSH)   /* if flush threshold has been reached */
                {
                    pthread_mutex_lock(&screen_mutex);
                    timeStamp(p->start);
                    printf("flushing storage writer at %6.02g%% load\n", cacheLoad);
                    storageWriterFlush(p->sw);
                    timeStamp(p->start);
                    printf("flushing finished (file storage now at %6.02g%% load)\n", storageWriterCurrentLoadPercentageFile(p->sw));
                    pthread_mutex_unlock(&screen_mutex);
                }
            }
            resumeWriting();
            break;
        case 2: /* this sample was discarded, but could have been added if cache was flushed */
            break;
        case 3: /* this sample could not be added to storage writer, category on file + cache is full */
            break;
        }
    }

}



int transition_times2_modq(const char *srcFolderName, const char *dstFolderName, u64 minDestinationStorageCapacityInSamples, bkwStepParameters *bkwStepPar, time_t start)
{
    if (folderExists(dstFolderName))
    {
        return 1; /* output folder already exists */
    }

    lweInstance lwe;
    lweParametersFromFile(&lwe, srcFolderName); /* read lwe parameters from source folder */

    ASSERT(N_THREADS >= 1, "Unexpected number of threads!");

    /* get number of samples in source file */
    u64 totNumUnsortedSamples = numSamplesInSampleFile(srcFolderName);
    if (!totNumUnsortedSamples)
    {
        lweDestroy(&lwe);
        return 2; /* no samples in source file */
    }
    char str[256];
    timeStamp(start);
    printf("src folder: %s (contains %s samples)\n", srcFolderName, sprintf_u64_delim(str, totNumUnsortedSamples));

    /* initiate storage writer (creates destination folder) */
    u64 numCategories = num_categories(&lwe, bkwStepPar);
    u64 categoryCapacityFile = (minDestinationStorageCapacityInSamples + numCategories - 1) / numCategories;
    storageWriter sw;
    int ret = storageWriterInitialize(&sw, dstFolderName, &lwe, bkwStepPar, categoryCapacityFile);
    if (ret)
    {
        lweDestroy(&lwe);
        return 100 + ret; /* could not initialize storage writer */
    }
    timeStamp(start);
    printf("dst folder: %s (has room for %s samples)\n", dstFolderName, sprintf_u64_delim(str, numCategories * categoryCapacityFile));

    pthread_t thread[N_THREADS];
    Params param[N_THREADS]; /* one set of in-/output paramaters per thread, so no need to lock these */

    /* allocate memory for an optimal number of mutex */
    n_storage_mutex = numCategories < MAX_NUM_STORAGE_MUTEXES ? numCategories : MAX_NUM_STORAGE_MUTEXES;
    max_num_categories_per_mutex = (numCategories + n_storage_mutex - 1) / n_storage_mutex;
    storage_mutex = malloc(n_storage_mutex * sizeof(pthread_mutex_t));
    for (int i=0; i<n_storage_mutex; i++) { pthread_mutex_init(&storage_mutex[i], NULL); }

    /* initialize pthread condition */
    pthread_cond_init(&resumeCond, NULL);


    /* open source sample file */
    FILE *f_src = fopenSamples(srcFolderName, "rb");
    if (!f_src)
    {
        lweDestroy(&lwe);
        if (storageWriterFree(&sw))
        {
            return 3; /* failed to free storage writer */
        }
        return 4; /* could not open samples file */
    }
//  printf("The size of the destination file before adding the samples is %" PRIu64 "\n", fileSize(sw.f));

    /* allocate sample read buffer */
    lweSample *sampleReadBuf = MALLOC(READ_BUFFER_CAPACITY_IN_SAMPLES * LWE_SAMPLE_SIZE_IN_BYTES);
    if (!sampleReadBuf)
    {
        lweDestroy(&lwe);
        fclose(f_src);
        if (storageWriterFree(&sw))
        {
            return 5; /* failed to free storage writer */
        }
        return 6; /* could not allocate sample read buffer */
    }
//  printf("Read buffer allocated (%d bytes / %d samples)\n", READ_BUFFER_CAPACITY_IN_SAMPLES * LWE_SAMPLE_SIZE_IN_BYTES, READ_BUFFER_CAPACITY_IN_SAMPLES);

    /* load input parameters */
    for (int i=0; i<N_THREADS; i++) {
    	param[i].index = i;
        param[i].lwe = &lwe; /* set input parameter to thread number */
        param[i].bkwStepPar = bkwStepPar;
    	param[i].sw = &sw;
        param[i].start = start;
    }

    nextPrintLimit = 100000000;

    /* process all samples in source file */
    while (!feof(f_src))
    {

        /* read chunk of samples from source sample file into read buffer */
        u64 numRead = freadSamples(f_src, sampleReadBuf, READ_BUFFER_CAPACITY_IN_SAMPLES);

	    /* load input parameters */
	    for (int i=0; i<N_THREADS; i++) {
	    	param[i].sampleReadBuf = sampleReadBuf;
	        param[i].min = i*(numRead/N_THREADS);
	        param[i].max = (i+1)*(numRead/N_THREADS);
	    }
	    param[N_THREADS-1].max = numRead;

	    /* start threads */
	    for (int i = 0; i < N_THREADS; ++i)
	    {
	        if (!pthread_create(&thread[i], NULL, single_thread_work, (void*)&param[i])) {
	            // pthread_mutex_lock(&screen_mutex);
	            // printf("Thread %d created!\n", i+1);
	            // pthread_mutex_unlock(&screen_mutex);
	        } else {
	            // pthread_mutex_lock(&screen_mutex);
	            // printf("Error creating thread %d!\n", i+1);
	            // pthread_mutex_unlock(&screen_mutex);
	        }
	    }

	    /* wait until all threads have completed */
	    for (int i = 0; i < N_THREADS; i++) {
	        if (!pthread_join(thread[i], NULL)) {
	            // pthread_mutex_lock(&screen_mutex);
	            // printf("Thread %d joined!\n", i+1);
	            // pthread_mutex_unlock(&screen_mutex);
	        } else {
	            // pthread_mutex_lock(&screen_mutex);
	            // printf("Error joining thread %d!\n", i+1);
	            // pthread_mutex_unlock(&screen_mutex);
	        }
	    }

    }

    char s1[256], s2[256];
    timeStamp(start);
    printf("%s samples added to storage writer (%s samples discarded because some categories were overfull)\n", sprintf_u64_delim(s1, sw.totalNumSamplesAddedToStorageWriter), sprintf_u64_delim(s2, sw.totalNumSamplesProcessedByStorageWriter - sw.totalNumSamplesAddedToStorageWriter));

    /* cleanup */
    FREE(sampleReadBuf);
    fclose(f_src);
    lweDestroy(&lwe);
    ret = storageWriterFree(&sw); /* flushes automatically */
    if (ret)
    {
        return 200 + ret; /* could not free storage writer */
    }
    return 0;

    for (int i=0; i<n_storage_mutex; i++) { pthread_mutex_destroy(&storage_mutex[i]); }

    pthread_cond_destroy(&resumeCond);

    free(&storage_mutex);

}

