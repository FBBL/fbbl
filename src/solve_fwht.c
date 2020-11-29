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

#include "solve_fwht.h"
#include "lwe_instance.h"
#include "storage_reader.h"
#include "storage_file_utilities.h"
#include "memory_utils.h"
#include "string_utils.h"
#include <math.h>
#include <inttypes.h>

/*
 * integer to binary sequence
 */
void int_to_bin(u64 input, u8 *binary, int binlen)
{
    // zero all entries
    for (int i = 0; i<binlen; i++)
        binary[i] = 0;

    // counter for binary array
    int i = 0;
    while (input > 0)
    {
        // storing remainder in binary array
        binary[i] = input % 2;
        input = input / 2;
        i++;
    }
}

/*
 * Convert sample to binary sequence, then to integer (u64)
 */
u64 sample_to_int(short *input, int len, int q)
{

    u64 output = 0;
    short a;

    for (int i = 0; i<len; i++)
    {
        a = input[i] <= q/2 ? input[i] : abs(input[i]-q);
        if (a % 2 != 0)
            output += ((u64)1)<<i;
    }

    return output;
}

/*
 * Fast in-place Walsh-Hadamard Transform. Created by Florian Tramer.
 * (adapted from http://www.musicdsp.org/showone.php?id=18)
 *
 * @data    the data to transform the transform over
 * @size    the length of the data
 */
#ifdef USE_SOFT_INFORMATION
void FWHT (double* data, int size)
#else
void FWHT (long* data, int size)
#endif
{
    int n = log2(size);
#ifdef USE_SOFT_INFORMATION
    double tmp;
#else
    long tmp;
#endif
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < (1 << n); j += 1 << (i+1))
        {
            for (int k = 0; k < (1<<i); ++k)
            {
                int a = j + k;
                int b = j + k + (1<<i);

                tmp = data[a];
                data[a] += data[b];
                data[b] = tmp - data[b];
            }
        }
    }
}

#ifdef USE_SOFT_INFORMATION

#define PI 3.14159265358979323846
#define BOUND 101
double *bias_table;

double normal_pdf(double x, double mu, double sigma)
{
    return (1/(sigma*sqrt(2*PI)))*exp(-((x-mu)*(x-mu))/(2*sigma*sigma));
}

/* NOTE: we assume q is odd */
void initialize_bias_table(int q, float sigma)
{

    bias_table = CALLOC(q, sizeof(double));
    ASSERT(bias_table, "Error allocating memory for bias_table");
    double even, odd, sum;

    for(int z = 0; z<q; z++)
    {
        odd = 0;
        int y = z-(q-1)/2;
        even = normal_pdf(y, 0, sigma);
        for(int j = 0; j<BOUND; j++)
        {
            if(j%2)
                odd += normal_pdf(y+j*q, 0, sigma) + normal_pdf(y-j*q, 0, sigma);
            else
                even += normal_pdf(y+j*q, 0, sigma) + normal_pdf(y-j*q, 0, sigma);
        }
        sum = even + odd;
        even = even/sum;
        odd = odd/sum;
        bias_table[z] = 2*even-1;
    }
}

void free_bias_table()
{
    if(bias_table)
        FREE(bias_table);
}
#endif

/* Retrieve binary secret using Fast Walsh Hadamard Transform */
#ifdef USE_SOFT_INFORMATION
int solve_fwht_search(const char *srcFolder, u8 *binary_solution, int zeroPositions, int fwht_positions, double sigma, time_t start)
{
#else
int solve_fwht_search(const char *srcFolder, u8 *binary_solution, int zeroPositions, int fwht_positions, time_t start)
{
#endif

    lweInstance lwe;
    lweParametersFromFile(&lwe, srcFolder);
    int n = lwe.n;
    int q = lwe.q;

    u64 numCategories;
    sampleInfoFromFile(srcFolder, NULL, &numCategories, NULL, NULL, NULL);

    ASSERT(fwht_positions <= MAX_FWHT, "The number of positions for fwht is not supported in this implementation!\n");

    if(zeroPositions == n)
    {
        printf("Coordinates all solved\n");
        return 1;
    }

    /* create initial list */
    u64 N = (u64)1<<fwht_positions; // N = 2^fwht_positions
#ifdef USE_SOFT_INFORMATION
    double* list = CALLOC(N,sizeof(double));
#else
    long* list = CALLOC(N,sizeof(long));
#endif
    if (!list)
    {
        printf("*** solve_fwht_search: failed to allocate memory for initial list\n");
        exit(-1);
    }

    /* open source sample file */
    FILE *f_src = fopenSamples(srcFolder, "rb");
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

#ifdef USE_SOFT_INFORMATION
    /* initialize bias_table */
    initialize_bias_table(q, sigma);
#endif

    /* process all samples in source file */
    lweSample *sample;
    short z, lsb_z;
    u64 intsample;

    while (!feof(f_src))
    {
        /* read chunk of samples from source sample file into read buffer */
        u64 numRead = freadSamples(f_src, sampleReadBuf, READ_BUFFER_CAPACITY_IN_SAMPLES);

        for (u64 i=0; i<numRead; i++)
        {
            sample = &sampleReadBuf[i];
            intsample = sample_to_int(sample->col.a+zeroPositions, fwht_positions, q);
            z = sample->sumWithError > (q-1)/2 ? (sample->sumWithError -q) : (sample->sumWithError);
            lsb_z = z%2 == 0 ? 0 : 1;
#ifdef USE_SOFT_INFORMATION
            if (lsb_z == 0)
                list[intsample] += bias_table[z+(q-1)/2];
            else
                list[intsample] -= bias_table[z+(q-1)/2];
#else
            if (lsb_z == 0)
                list[intsample] += 1;
            else
                list[intsample] -= 1;
#endif
        }
    }
    fclose(f_src);
    FREE(sampleReadBuf);

#ifdef USE_SOFT_INFORMATION
    /* free bias_table */
    free_bias_table();
#endif


    /* Apply Fast Walsh Hadamard Tranform */
    timeStamp(start);
    printf("Start FWHT\n");
    FWHT(list, N);

    // find maximum
    u64 max_pos = -1;
    double max = 0;
    double tot = 0;
    for (int i = 0; i<N; i++)
    {
#ifdef USE_SOFT_INFORMATION
        tot += fabs(list[i]);
        if (max < fabs(list[i]))
        {
            max = fabs(list[i]);
#else
        tot += labs(list[i]);
        if (max < labs(list[i]))
        {
            max = labs(list[i]);
#endif
            max_pos = i;
        }
    }
    timeStamp(start);
    printf("Index found %lld - max %f - probability %f \n", max_pos, max, (max)/tot);

    // Convert solution into binary
    int_to_bin(max_pos, binary_solution, fwht_positions);



    FREE(list);
    return 0;
}

/* return 0 if guesses are terminated, otherwise 1 */
int nextBruteForceGuess(int ratio, int *BFguess, int lenght)
{

    for(int i = lenght-1; i>=0; i--)
    {
        if(BFguess[i] < ratio)
        {
            BFguess[i]++;
            for(int j = lenght-1; j>i; j--)
                BFguess[j] = -ratio;
            return 1;
        }
    }
    return 0;
}

/* Hybrid solver that uses brute-force for bruteForcePositions number of positions and Fast Walsh Hadamard Transform for fftPositions number of positions
 */
#ifdef USE_SOFT_INFORMATION
int solve_fwht_search_bruteforce(const char *srcFolder, u8 *binary_solution, short *bf_solution, int zeroPositions, int bruteForcePositions, int fwhtPositions, double sigma, time_t start)
{
#else
int solve_fwht_search_bruteforce(const char *srcFolder, u8 *binary_solution, short *bf_solution, int zeroPositions, int bruteForcePositions, int fwhtPositions, time_t start)
{
#endif
    lweInstance lwe;
    lweParametersFromFile(&lwe, srcFolder);
    int q = lwe.q;

    u64 numCategories;
    sampleInfoFromFile(srcFolder, NULL, &numCategories, NULL, NULL, NULL);

    ASSERT(1 <= fwhtPositions && fwhtPositions <= MAX_FWHT, "The number of positions for fwht is not supported in this implementation!\n");
    ASSERT(1 <= bruteForcePositions && bruteForcePositions <= MAX_BRUTE_FORCE, "The number of positions for bruteforce guessing is not supported in this implementation!\n");
    ASSERT(fwhtPositions + bruteForcePositions + zeroPositions == n, "The number of positions for bruteforce and fwht is => n!\n");

    /* create initial list */
    u64 N = (u64)1<<fwhtPositions; // N = 2^fwht_positions
#ifdef USE_SOFT_INFORMATION
    double* list = CALLOC(N,sizeof(double));
#else
    long* list = CALLOC(N,sizeof(long));
#endif
    if (!list)
    {
        printf("*** solve_fwht_search: failed to allocate memory for initial list\n");
        exit(-1);
    }

#ifdef USE_SOFT_INFORMATION
    /* initialize bias_table */
    initialize_bias_table(q, sigma);
#endif

    FILE *f_src;

    /* allocate sample read buffer */
    lweSample *sampleReadBuf;

    lweSample *sample;
    long max_pos = -1;
    double max = 0, global_max = 0;
    u64 intsample;
    short z, lsb_z;

    u8 bin_guess[fwhtPositions];

    timeStamp(start);
    printf("Start FWHT: brute force %d positions, guess %d positions\n", bruteForcePositions, fwhtPositions);

    int ratio = round(lwe.alpha*lwe.q*3); // 2*3*standard_deviation is the interval length where to search
    int BFguess[bruteForcePositions];
    int tmp_si;

    // initialize brute force guess
    for(int i = 0; i<bruteForcePositions; i++)
        BFguess[i] = -ratio;

    do
    {
#ifdef USE_SOFT_INFORMATION
        memset(list, 0, N*sizeof(double));
#else
        memset(list, 0, N*sizeof(long));
#endif
        sampleReadBuf = MALLOC(READ_BUFFER_CAPACITY_IN_SAMPLES * LWE_SAMPLE_SIZE_IN_BYTES);
        if (!sampleReadBuf)
        {
            return 6; /* could not allocate sample read buffer */
        }

        /* process all samples in source file */
        /* open source sample file */
        f_src = fopenSamples(srcFolder, "rb");
        if (!f_src)
        {
            return 4; /* could not open samples file */
        }

        while (!feof(f_src))
        {
            /* read chunk of samples from source sample file into read buffer */
            u64 numRead = freadSamples(f_src, sampleReadBuf, READ_BUFFER_CAPACITY_IN_SAMPLES);
            for (u64 i=0; i<numRead; i++)
            {
                sample = &sampleReadBuf[i];
                intsample = sample_to_int(sample->col.a+zeroPositions, fwhtPositions, q);
                z = sample->sumWithError;
                // update z with the bruteforce-guessed positions
                for(int j = 0; j<bruteForcePositions; j++)
                {
                    if(BFguess[j] < 0)
                        tmp_si = BFguess[j]+q;
                    else
                        tmp_si = BFguess[j];
                    z = (z - (columnValue(sample, zeroPositions+fwhtPositions+j)*tmp_si)%q );
                    if (z < 0)
                        z += q;
                }
                z = z > (q-1)/2 ? z -q : z;
                intsample = sample_to_int(sample->col.a+zeroPositions, fwhtPositions, q);
                lsb_z = z%2 == 0 ? 0 : 1;
#ifdef USE_SOFT_INFORMATION
                if (lsb_z == 0)
                    list[intsample] += bias_table[z+(q-1)/2];
                else
                    list[intsample] -= bias_table[z+(q-1)/2];
#else
                if (lsb_z == 0)
                    list[intsample] += 1;
                else
                    list[intsample] -= 1;
#endif
            }
        }
        fclose(f_src);
        FREE(sampleReadBuf);

        /* Apply Fast Walsh Hadamard Tranform */
        FWHT(list, N);

        // find maximum
        max_pos = -1;
        max = 0;
        for (int i = 0; i<N; i++)
        {
#ifdef USE_SOFT_INFORMATION
            if (max < fabs(list[i]))
            {
                max = fabs(list[i]);
#else
            if (max < labs(list[i]))
            {
                max = labs(list[i]);
#endif
                max_pos = i;
            }
        }

        // Convert solution into binary
        int_to_bin(max_pos, bin_guess, fwhtPositions);

#ifdef PRINT_INTERMEDIATE_SOLUTIONS_BRUTEFORCE
        timeStamp(start);
        printf("Index found %ld - max %f \n(",max_pos, max);
        for(int j = 0; j<bruteForcePositions; j++)
            printf("%d ", BFguess[j]);
        printf(") - (");
        for(int j = 0; j<fwhtPositions; j++)
            printf("%hu ", bin_guess[j]);
        printf(")\n");
#endif

        if(max > global_max)
        {
            global_max = max;
            for(int j = 0; j<fwhtPositions; j++)
                binary_solution[j] = bin_guess[j];
            for(int j = 0; j<bruteForcePositions; j++)
                bf_solution[j] = BFguess[j] >= 0 ? BFguess[j] : BFguess[j] +q;
        }

    }
    while(nextBruteForceGuess(ratio, BFguess, bruteForcePositions));

#ifdef USE_SOFT_INFORMATION
    /* free bias_table */
    free_bias_table();
#endif
    FREE(list);
    return 0;
}

#ifndef USE_SOFT_INFORMATION
/* Hybrid solver that uses sparse brute-force for bruteForcePositions number of positions and Fast Walsh Hadamard Transform for fftPositions number of positions */
/* NOTE: if the brute-force positions are not reduced, the success probability goes down quite a lot, this is more intended to be used for guessing a larger number
 * of positions than the FWHT can handle. It does not support soft-information.
 */
int solve_fwht_search_hybrid(const char *srcFolder, u8 *binary_solution, int zeroPositions, int bruteForcePositions, int fwhtPositions, time_t start)
{

    lweInstance lwe;
    lweParametersFromFile(&lwe, srcFolder);
    int n = lwe.n;
    int q = lwe.q;

    u64 numCategories;
    sampleInfoFromFile(srcFolder, NULL, &numCategories, NULL, NULL, NULL);

    ASSERT(1 <= fwhtPositions && fwhtPositions <= MAX_FWHT, "The number of positions for fwht is not supported in this implementation!\n");
    ASSERT(1 <= bruteForcePositions && bruteForcePositions <= MAX_BRUTE_FORCE, "The number of positions for bruteforce guessing is not supported in this implementation!\n");
    ASSERT(fwhtPositions + bruteForcePositions + zeroPositions == n, "The number of positions for bruteforce and fwht is => n!\n");

    /* create initial list */
    u64 N = (u64)1<<fwhtPositions; // N = 2^fwht_positions
    long *list = CALLOC(N,sizeof(long));  // calloc every time so that we start with all zeros
    if (!list)
    {
        printf("*** solve_fwht_search: failed to allocate memory for initial list\n");
        exit(-1);
    }

    FILE *f_src;

    /* allocate sample read buffer */
    lweSample *sampleReadBuf;

    lweSample *sample;
    short z;
    long max_pos = -1, global_max = 0, max = 0;
    u64 intsample;

    timeStamp(start);
    printf("Start FWHT: brute force %d positions, guess %d positions\n", bruteForcePositions, fwhtPositions);

    u64 MAX = ((int)1)<<bruteForcePositions;
    u8 bin_guess[n-zeroPositions];

    /* bruteforce the last bruteForcePositions positions */
    for(u64 guess = 0; guess<MAX; guess++)
    {
        int_to_bin(guess, bin_guess+fwhtPositions, bruteForcePositions);

        memset(list, 0, N*sizeof(long));
        sampleReadBuf = MALLOC(READ_BUFFER_CAPACITY_IN_SAMPLES * LWE_SAMPLE_SIZE_IN_BYTES);
        if (!sampleReadBuf)
        {
            return 6; /* could not allocate sample read buffer */
        }

        /* process all samples in source file */
        /* open source sample file */
        f_src = fopenSamples(srcFolder, "rb");
        if (!f_src)
        {
            return 4; /* could not open samples file */
        }
        while (!feof(f_src))
        {
            /* read chunk of samples from source sample file into read buffer */
            u64 numRead = freadSamples(f_src, sampleReadBuf, READ_BUFFER_CAPACITY_IN_SAMPLES);
            for (u64 i=0; i<numRead; i++)
            {
                sample = &sampleReadBuf[i];
                intsample = sample_to_int(sample->col.a+zeroPositions, fwhtPositions, q);
                z = sample->sumWithError <= q/2 ? (sample->sumWithError%2) : (sample->sumWithError -q)%2;
                // update z with the bruteforce-guessed positions
                for(int j = fwhtPositions; j<fwhtPositions+bruteForcePositions; j++)
                    z = (z + abs(columnValueSigned(sample, zeroPositions+j, q)%2)*(bin_guess[j]))%2 ;
                if (z == 0)
                    list[intsample] += 1;
                else
                    list[intsample] -= 1;
            }
        }
        fclose(f_src);
        FREE(sampleReadBuf);

        /* Apply Fast Walsh Hadamard Tranform */
        FWHT(list, N);

        // find maximum
        max_pos = -1;
        max = 0;
        for (int i = 0; i<N; i++)
        {
            if (max < labs(list[i]))
            {
                max = labs(list[i]);
                max_pos = i;
            }
        }

        // Convert solution into binary
        int_to_bin(max_pos, bin_guess, fwhtPositions);

#if 1 // print intermediate solutions
        timeStamp(start);
        printf("Index found %ld - max %ld \n(",max_pos, max);
        for(int j = 0; j<n-zeroPositions; j++)
            printf("%hu ", bin_guess[j]);
        printf(")\n");
#endif

        if(max > global_max)
        {
            global_max = max;
            for(int j = 0; j<n; j++)
                binary_solution[j] = bin_guess[j];
        }
    }

    FREE(list);
    return 0;
}

static short power_of_2[14] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192}; // it should be enough

int retrieve_full_secret(short *full_secret, u8 binary_secret[][MAX_N], int n_iterations, int n, int q)
{

    int one_zero, max_abs;

    for(int i=0; i<n; i++)
    {

        full_secret[i] = 0;
        // save last resulting bits, it tells us the sign of the secret entry
        one_zero = binary_secret[n_iterations-1][i];

        if(one_zero)   // negative secret entry case
        {

            max_abs = n_iterations-1;
            while(binary_secret[max_abs][i] == one_zero && max_abs >= 0)
                max_abs--;

            // retrieve correct binary representation
            for (int j = 0; j<n_iterations; j++)
            {
                if (binary_secret[j][i] != 0)
                {
                    for (int k = j+1; k < n_iterations; ++k)
                        binary_secret[k][i] = !binary_secret[k][i];
                    break;
                }
            }

            for(int j = 0; j<max_abs+2; j++)
                full_secret[i] += power_of_2[j]*binary_secret[j][i];
            full_secret[i] = (q-full_secret[i]) % q;

        }
        else     // positive secret entry case
        {

            for(int j = 0; j<n_iterations; j++)
                full_secret[i] += power_of_2[j]*binary_secret[j][i];
        }
    }

    return 0;
}
#endif



