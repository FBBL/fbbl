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

#ifndef TRANSITION_BKW_STEP_SMOOTH_LMS_H
#define TRANSITION_BKW_STEP_SMOOTH_LMS_H

#include "bkw_step_parameters.h"
#include <time.h>

int transition_bkw_step_smooth_lms(const char *srcFolderName, const char *dstFolderName, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, u64 *numSamplesStored, time_t start);

#endif /* TRANSITION_BKW_STEP_SMOOTH_LMS_H */
