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

#include "transition_reduce_secret.h"
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

/* a-> 2a
 * b -> b-a*lsb_secret
 */
int transform_sample(lweInstance *lwe, lweSample *sample, u8 * lsb_secret)
{
    int n = lwe->n;
    int q = lwe->q;
    // compute a*lsb_secret
    int subt = 0;
    int sgn = 0;

    for (int i=0; i<n; i++)
    {
        if (!lsb_secret[i])
            sgn = 0;
        else
            sgn = lsb_secret[i] <= q/2 ? 1 : q-1;
        subt = (subt+sample->col.a[i]*sgn) % q;
    }
    for (int i=0; i<n; i++)
        sample->col.a[i] = (2*sample->col.a[i]) % q;
    sample->sumWithError = (sample->sumWithError -subt) >= 0 ? (sample->sumWithError -subt) : (sample->sumWithError -subt) +q;
    sample->col.hash = bkwColumnComputeHash(sample, n, 0);
    return 0;
}

/* given lsb of the secret, get a reduced secret
 */
int reduce_secret(lweInstance *lwe, u8 * lsb_secret)
{
    for(int i = 0; i<lwe->n; i++)
    {
        if(lwe->s[i] <= lwe->q/2)
            lwe->s[i] = (lwe->s[i]-lsb_secret[i])>>1;
        else
            lwe->s[i] = ((short)lwe->q -(((short)lwe->q-lwe->s[i]+lsb_secret[i])>>1))%lwe->q;
    }
    return 0;
}


/*
 * Given the LSBs of the secret, create new samples with a reduced version of the secret
 */
int transition_reduce_secret(const char *srcFolderName, const char *dstFolderName, u8 * lsb_secret, time_t start)
{

    if (folderExists(dstFolderName))   /* if destination folder already exists, assume that we have performed this reduction step already */
    {
        return 1; /* reduction step already performed (destination folder already exists) */
    }
    timeStamp(start);
    printf("src folder: %s\n", srcFolderName);
    timeStamp(start);
    printf("dst folder: %s\n", dstFolderName);

    lweInstance lwe;
    lweParametersFromFile(&lwe, srcFolderName); /* read lwe parameters from source folder */

    /* get number of samples in source file */
    u64 totNumdSamples = numSamplesInSampleFile(srcFolderName);
    if (!totNumdSamples)
    {
        return 2; /* no samples in source file */
    }
    char str[256];
    timeStamp(start);
    printf("src folder: %s (contains %s samples)\n", srcFolderName, sprintf_u64_delim(str, totNumdSamples));

    /* Reduce secret */
    reduce_secret(&lwe, lsb_secret);

    /* open source sample file */
    FILE *f_src = fopenSamples(srcFolderName, "rb");
    if (!f_src)
    {
        return 4; /* could not open samples file */
    }

    /* allocate sample read buffer */
    lweSample *sampleReadBuf = MALLOC(READ_BUFFER_CAPACITY_IN_SAMPLES * LWE_SAMPLE_SIZE_IN_BYTES);
    if (!sampleReadBuf)
    {
        fclose(f_src);
        return 6; /* could not allocate sample read buffer */
    }

    /* process all samples in source file */
    lweSample *s;

    while (!feof(f_src))
    {
        /* read chunk of samples from source sample file into read buffer */
        u64 numRead = freadSamples(f_src, sampleReadBuf, READ_BUFFER_CAPACITY_IN_SAMPLES);

        /* transform sample by subtracting a*s in both sides of the equation */
        for (u64 i=0; i<numRead; i++)
        {
            s = &sampleReadBuf[i];
            transform_sample(&lwe, s, lsb_secret);
        }

        /* add samples to storage writer */
        u64 n = numRead;
        int blockSize = 1000000;

        newStorageFolderWithGivenLweInstance(&lwe, dstFolderName);
        newStorageFolder(&lwe, dstFolderName, lwe.n, lwe.q, lwe.alpha);
        FILE *f = fopenSamples(dstFolderName, "ab");
        if (!f)
        {
            return -1;
        }

        while (n)
        {
            int numSamplesThisRound = MIN(n, blockSize);
            int numWritten = fwrite(sampleReadBuf, LWE_SAMPLE_SIZE_IN_BYTES, numSamplesThisRound, f);
            n -= numWritten;
            // printf("n = %" PRIu64 "\n", n);
        }
        fclose(f);
    }

    /* cleanup */
    FREE(sampleReadBuf);
    fclose(f_src);

    return 0;
}



