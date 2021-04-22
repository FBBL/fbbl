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

#ifndef LWE_TRANSITION_MOD2_H
#define LWE_TRANSITION_MOD2_H

#include "bkw_step_parameters.h"
#include <time.h>

#define APPROXIMATE_SIZE_IN_BYTES_OF_READ_BUFFER (250*1024*1024)
#define MIN_STORAGE_WRITER_CACHE_LOAD_PERCENTAGE_BEFORE_FLUSH 25
#define READ_BUFFER_CAPACITY_IN_SAMPLES (APPROXIMATE_SIZE_IN_BYTES_OF_READ_BUFFER / LWE_SAMPLE_SIZE_IN_BYTES)

int transition_mod2(const char *srcFolderName, const char *dstFolderName, time_t start);

#endif /* LWE_TRANSITION_MOD2_H */
