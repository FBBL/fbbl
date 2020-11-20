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

#ifndef TRANSITION_UNSORTED_2_SORTED_H
#define TRANSITION_UNSORTED_2_SORTED_H
#include "bkw_step_parameters.h"
#include <time.h>

/* 250 MB default size for read buffer */
#define APPROXIMATE_SIZE_IN_BYTES_OF_READ_BUFFER (250*1024*1024)

/* Short explanation:
   The storage writer cache is flushed to file when it is "full enough",
   as defined by the macro below.

   Long explanation:
   The unsorted samples are read sequentially from the source file
   and are then written to a destination file using a storage writer,
   which sorts the samples into their corresponding category. The precise
   destination category for a given sample depends on the sample itself
   and sorting parameters (how sorting is performed).
   If this sorting behaves well, the process can be modeled as a coupon
   collector's problem.
   The storage writer has a cache buffer (RAM) into which the inserted
   samples are written. The buffer is typically much smaller than
   the storage available on file. That is, categories on file typically
   hold more samples than the corresponding category storage in the cache
   buffer. Every time the storage writer is flushed, the cached samples
   are written into (added) to file. This entails one sequential read/write
   across the entire full-size sample file. If you have many samples, this
   will take some time.
   At the expense of discarding a few samples (that over-fill the cache
   buffer categories, but would potentially fit within the category storage
   on file), we may improve running time by reducing the amount of flushing
   that we perform.
   If sorting behaves reasonably well, the first cache category may be filled
   at a storage writer cache load of about 5% (heavily depends on problem
   parameters, of course).
   */
#define MIN_STORAGE_WRITER_CACHE_LOAD_PERCENTAGE_BEFORE_FLUSH 25

int transition_unsorted_2_sorted(const char *srcFolderName, const char *dstFolderName, u64 minDestinationStorageCapacityInSamples, bkwStepParameters *bkwStepPar, time_t start);


/* automatically defined */
#define READ_BUFFER_CAPACITY_IN_SAMPLES (APPROXIMATE_SIZE_IN_BYTES_OF_READ_BUFFER / LWE_SAMPLE_SIZE_IN_BYTES)

#endif
