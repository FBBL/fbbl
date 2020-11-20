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

#ifndef SRC_FILE_BASED_LWE_SOLVE_FWHT_H_
#define SRC_FILE_BASED_LWE_SOLVE_FWHT_H_

#include "log_utils.h"
#include "lwe_instance.h"

#define MAX_FWHT 35

#define MAX_BRUTE_FORCE 60

//#define USE_SOFT_INFORMATION
#define PRINT_INTERMEDIATE_SOLUTIONS_BRUTEFORCE

#define APPROXIMATE_SIZE_IN_BYTES_OF_READ_BUFFER (250*1024*1024)
#define MIN_STORAGE_WRITER_CACHE_LOAD_PERCENTAGE_BEFORE_FLUSH 25
#define READ_BUFFER_CAPACITY_IN_SAMPLES (APPROXIMATE_SIZE_IN_BYTES_OF_READ_BUFFER / LWE_SAMPLE_SIZE_IN_BYTES)

int retrieve_full_secret(short *full_secret, u8 binary_secret[][MAX_N], int n_iterations, int n, int q);

#ifdef USE_SOFT_INFORMATION

int solve_fwht_search(const char *srcFolder, u8 *binary_solution, int zeroPositions, int fwht_positions, double sigma, time_t start);
int solve_fwht_search_bruteforce(const char *srcFolder, u8 *binary_solution, short *bf_solution, int zeroPositions, int bruteForcePositions, int fwhtPositions, double sigma, time_t start);

#else

int solve_fwht_search(const char *srcFolder, u8 *binary_solution, int zeroPositions, int fwht_positions, time_t start);
int solve_fwht_search_bruteforce(const char *srcFolder, u8 *binary_solution, short *bf_solution, int zeroPositions, int bruteForcePositions, int fwhtPositions, time_t start);
int solve_fwht_search_hybrid(const char *srcFolder, u8 *binary_solution, int zeroPositions, int bruteForcePositions, int fwht_positions, time_t start);

#endif

#endif /* SRC_FILE_BASED_LWE_SOLVE_FWHT_H_ */
