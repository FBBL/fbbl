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

#ifndef SOLVE_FFT_H
#define SOLVE_FFT_H

#define FFT_SOLVER_SINGLE_PRECISION 0
#define FFT_SOLVER_DOUBLE_PRECISION 1

//void test_fft_solver(const char *srcFolder);
int solve_fft_search(const char *srcFolder, short *solution, int numSolvedCoordinates, int fftPositions, int doublePrecision);
int solve_fft_search_hybrid(const char *srcFolder, short *solution, int numSolvedCoordinates, int fftPositions, int bruteForcePositions, int doublePrecision);

#endif
