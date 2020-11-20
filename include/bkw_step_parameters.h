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

#ifndef BKW_STEP_PARAMETERS_H
#define BKW_STEP_PARAMETERS_H
#include "lwe_instance.h"
#include "lwe_sorting.h"
#include "lwe_sample_selection.h"

#define MAX_LMS_POSITIONS 6
#define MAX_CODED_BKW_POSITIONS 4
#define MAX_SMOOTH_LMS_POSITIONS 10


typedef enum
{
    blockCode_21,
    blockCode_31,
    blockCode_41,
    concatenatedCode_21_21,
    numCodingTypes
} codingType;

typedef struct
{
    sortingMethod sorting;
    int startIndex;
    int numPositions; /* = Ni, in case of codedBKW */
    selectionMethod selection;
    union
    {
        struct
        {
            short p; /* reduction factor, same over all positions */
        } LMS;
        struct
        {
            short p; /* general reduction factor */
            short p1; /* additional reduction factor for position Ni+1 */
            short p2; /* additional reduction factor for first position after the first step */
            short prev_p1; /* reduction factor for position Ni+1 in the previous step - used on all but the first steps */
            short meta_skipped;
            short unnatural_selection_ts; /* unnatural selection threshold */
            short unnatural_selection_start_index; /* Starting index for the unnatural selection */
        } smoothLMS;
        struct
        {
            short b; /* dimension of the linear code (Ni,b) */
            codingType ct;
        } CodedBKW;
    } sortingPar;
} bkwStepParameters;

const char *codingTypeAsString(codingType ct);
codingType codingTypeFromString(const char *str);

char *bkwStepParametersAsString(char *str, bkwStepParameters *bkwStepPar);
int bkwStepParametersFromString(const char *str, bkwStepParameters *bkwStepPar);

u64 num_categories(lweInstance *lwe, bkwStepParameters *bkwStepPar);

#endif
