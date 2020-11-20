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

#include "transition_bkw_step.h"
#include "transition_bkw_step_plain_bkw_2_positions.h"
#include "transition_bkw_step_plain_bkw_3_positions.h"
#include "transition_bkw_step_lms.h"
#include "transition_bkw_step_smooth_lms.h"
#include "transition_bkw_step_smooth_lms_meta.h"
#include "transition_bkw_step_coded_bkw.h"
#include "log_utils.h"
#include "string_utils.h"
#include "storage_file_utilities.h"

int transition_bkw_step(const char *srcFolderName, const char *dstFolderName, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, u64 *numSamplesStored, time_t start)
{
    int ret = -1;
    timeStamp(start);
    printf("src folder: %s\n", srcFolderName);
    timeStamp(start);
    printf("dst folder: %s\n", dstFolderName);
    if (folderExists(dstFolderName))   /* if destination folder already exists, assume that we have performed this reduction step already */
    {
        return 100; /* reduction step already performed (destination folder already exists) */
    }

    switch (srcBkwStepPar->sorting)
    {

    case plainBKW:

        switch (srcBkwStepPar->numPositions)
        {
        case 2:
            ret = transition_bkw_step_plain_bkw_2_positions(srcFolderName, dstFolderName, srcBkwStepPar, dstBkwStepPar, numSamplesStored, start);
            break;
        case 3:
            ret = transition_bkw_step_plain_bkw_3_positions(srcFolderName, dstFolderName, srcBkwStepPar, dstBkwStepPar, numSamplesStored, start);
            break;
        default:
            ASSERT_ALWAYS("unsupported number of bkw positions detected in transition_bkw_step");
        }
        break;

    case LMS:
        ret = transition_bkw_step_lms(srcFolderName, dstFolderName, srcBkwStepPar, dstBkwStepPar, numSamplesStored, start);
        break;

    case smoothLMS:
        if (srcBkwStepPar->sortingPar.smoothLMS.meta_skipped == 0)
            ret = transition_bkw_step_smooth_lms(srcFolderName, dstFolderName, srcBkwStepPar, dstBkwStepPar, numSamplesStored, start);
        else
            ret = transition_bkw_step_smooth_lms_meta(srcFolderName, dstFolderName, srcBkwStepPar, dstBkwStepPar, numSamplesStored, start);
        break;

    case codedBKW:
        ret = transition_bkw_step_coded_bkw(srcFolderName, dstFolderName, srcBkwStepPar, dstBkwStepPar, numSamplesStored, start);
        break;

    default:
        ASSERT_ALWAYS("unsupported sorting detected in transition_bkw_step");
    }

    if (ret)   /* catch errors thrown in bkw reduction step */
    {
        timeStamp(start);
        printf("*** transition_bkw_step returned %d\n", ret);
    }
    else
    {
        char s[256];
        timeStamp(start);
        printf("%s reduction completed\n", sortingAsString(dstBkwStepPar->sorting));
        timeStamp(start);
        printf("%s samples stored\n", sprintf_u64_delim(s, *numSamplesStored));
    }

    return ret;
}
