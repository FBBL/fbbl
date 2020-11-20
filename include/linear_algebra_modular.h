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

#ifndef LINEAR_ALGEBRA_MODULAR_H
#define LINEAR_ALGEBRA_MODULAR_H
#include "lwe_instance.h"

/* Given m = n + k (where k >= 0) bkw samples,
   select the first n linearly independent samples to form matrix A.
   List of sample indices returned for A. (samples represent rows in A)
   Also computes and returns inverse of A. (samples represent rows here too)
   Returns number of samples used on success, returns 0 on failure.
   Failure probability is (for prime q) at most 1/q^{k+1}
   */
int compute_matrix_inverse_modular(lweInstance *lwe, lweSample *samples, int num_samples, short **A, short **A_inverse, short *);

#endif
