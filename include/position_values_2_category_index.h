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

#ifndef POSITION_VALUES_2_CATEGORY_INDEX
#define POSITION_VALUES_2_CATEGORY_INDEX
#include "bkw_step_parameters.h"

/* position values to index value */
u64 position_values_2_category_index(lweInstance *lwe, lweSample *sample, bkwStepParameters *bkwStepPar);
u64 position_values_2_category_index_from_partial_sample(lweInstance *lwe, short *a, bkwStepParameters *bkwStepPar);


/* for slightly more efficient loop unrolling and such, these functions can be used directly */

/* plain BKW 2 positions */
u64 position_values_2_category_index_plain_bkw(int q, short *a);
void category_index_2_position_values_plain_bkw(int q, u64 category_index, short *a);

/* LMS */
u64 position_values_2_category_index_lms(lweInstance *lwe, bkwStepParameters *dstBkwStepPar, short *pn);
int num_lms_singletons_in_category_interval(int a, int b);
int is_lms_singleton(u64 categoryIndex);

/* smooth LMS */
u64 position_values_2_category_index_smooth_lms(lweInstance *lwe, bkwStepParameters *dstBkwStepPar, short *pn);

/* smooth LMS with meta categories*/
u64 position_values_2_category_index_smooth_lms_meta(lweInstance *lwe, bkwStepParameters *dstBkwStepPar, short *pn);

/* For convenience in smooth LMS with meta categories */
short positionSmoothLMSMap(short pn, int q, int q_, int p, int c);

/* coded BKW 2 */
u64 position_values_2_category_index_coded_bkw(lweInstance *lwe, bkwStepParameters *dstBkwStepPar, short *pn);

/* singleton categories */
int is_singleton(bkwStepParameters *bkwStepPar, u64 categoryIndex, u64 numCategories);

/* free table */
void free_table_plain_bkw_2_positions(void);

#endif
