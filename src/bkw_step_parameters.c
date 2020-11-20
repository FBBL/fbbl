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

#include "bkw_step_parameters.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MIN(a,b) (((a)<(b))?(a):(b))

static const char *coding_type_label[numCodingTypes] =
{
    "[2,1] block code",
    "[3,1] block code",
    "[4,1] block code",
    "concatenated [2,1][2,1] block code"
};

const char *codingTypeAsString(codingType ct)
{
    return ct < numCodingTypes ? coding_type_label[ct] : NULL;
}

codingType codingTypeFromString(const char *str)
{
    for (int i=0; i<numCodingTypes; i++)
    {
        const char *codingTypeName = codingTypeAsString(i);
        int len = strlen(codingTypeName);
        if (!strncmp(codingTypeName, str, len))
        {
            return i;
        }
    }
    ASSERT_ALWAYS("coding type could not be determined");
    return 0;
}

char *bkwStepParametersAsString(char *str, bkwStepParameters *bkwStepPar)
{
    ASSERT(str, "no parameter str");
    ASSERT(bkwStepPar, "no parameter bkwStepPar");

    const char *sortingName = sortingAsString(bkwStepPar->sorting);

    switch (bkwStepPar->sorting)
    {

    case plainBKW:
        if (bkwStepPar->numPositions != 2 &&
                bkwStepPar->numPositions != 3)
        {
            ASSERT_ALWAYS("unspported number of positions for plain BKW");
            return NULL;
        }
        sprintf(str, "%s [%d positions, start index=%d]", sortingName, bkwStepPar->numPositions, bkwStepPar->startIndex);
        return str;

    case LMS:
        ASSERT(2 <= bkwStepPar->numPositions && bkwStepPar->numPositions <= MAX_LMS_POSITIONS, "unspported number of positions for LMS");
        sprintf(str, "%s [%d positions, start index=%d, p=%d]", sortingName, bkwStepPar->numPositions, bkwStepPar->startIndex, bkwStepPar->sortingPar.LMS.p);
        return str;

    case smoothLMS:
        ASSERT(2 <= bkwStepPar->numPositions && bkwStepPar->numPositions <= MAX_SMOOTH_LMS_POSITIONS, "unspported number of positions for smooth LMS");
        sprintf(str, "%s [%d positions, start index=%d, p=%d, p1=%d, p2=%d, prev_p1=%d, meta_skipped=%d, unnatural_selection=%d, unnatural_selection_start_index=%d]", sortingName, bkwStepPar->numPositions, bkwStepPar->startIndex, bkwStepPar->sortingPar.smoothLMS.p, bkwStepPar->sortingPar.smoothLMS.p1, bkwStepPar->sortingPar.smoothLMS.p2, bkwStepPar->sortingPar.smoothLMS.prev_p1, bkwStepPar->sortingPar.smoothLMS.meta_skipped, bkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts, bkwStepPar->sortingPar.smoothLMS.unnatural_selection_start_index);
        return str;

    case codedBKW:
        ASSERT(2 <= bkwStepPar->numPositions && bkwStepPar->numPositions <= MAX_CODED_BKW_POSITIONS, "unspported number of positions for coded BKW");
        const char *codingTypeString = codingTypeAsString(bkwStepPar->sortingPar.CodedBKW.ct);
        sprintf(str, "%s [%d positions, start index=%d, code=%s]", sortingName, bkwStepPar->numPositions, bkwStepPar->startIndex, codingTypeString);
        return str;

    default:
        ASSERT_ALWAYS("unhandled sorting in parameter set");
        return NULL;
    }

    ASSERT_ALWAYS("unhandled parameter set");
    return NULL;
}

int bkwStepParametersFromString(const char *str, bkwStepParameters *bkwStepPar)
{
    ASSERT(str, "no parameter str");
    ASSERT(bkwStepPar, "no parameter bkwStepPar");

    const char *p = str;
    char codingTypeName[64];

    int sortingDetermined = 0;
    for (int i=0; i<numSortingMethods; i++)
    {
        const char *sortingName = sortingAsString(i);
        int len = strlen(sortingName);
        if (!strncmp(sortingName, str, len))
        {
            bkwStepPar->sorting = i;
            sortingDetermined = 1;
            p += len + 1;
            break;
        }
    }
    if (!sortingDetermined)
    {
        ASSERT_ALWAYS("sorting could not be determined");
        return 0; /* sorting could not be determined */
    }

    /* extract additional parameters */
    switch (bkwStepPar->sorting)
    {

    case plainBKW:
        sscanf(p, "[%d positions, start index=%d]", &bkwStepPar->numPositions, &bkwStepPar->startIndex);
        if (bkwStepPar->numPositions != 2 &&
                bkwStepPar->numPositions != 3)
        {
            ASSERT_ALWAYS("unspported number of positions for plain BKW");
            return 0; /* error in additional parameters */
        }
        return 1; /* Plain BKW-sorting fully determined */

    case LMS:
        sscanf(p, "[%d positions, start index=%d, p=%hu]", &bkwStepPar->numPositions, &bkwStepPar->startIndex, &bkwStepPar->sortingPar.LMS.p);
        if (2 > bkwStepPar->numPositions || bkwStepPar->numPositions > MAX_LMS_POSITIONS)
        {
            ASSERT_ALWAYS("unsupported number of positions for LMS");
            return 0; /* error in additional parameters */
        }
        return 1; /* LMS-sorting fully determined */

    case smoothLMS:
        sscanf(p, "[%d positions, start index=%d, p=%hi, p1=%hi, p2=%hi, prev_p1=%hi, meta_skipped=%hi, unnatural_selection=%hi, unnatural_selection_start_index=%hi]", &(bkwStepPar->numPositions), &(bkwStepPar->startIndex), &(bkwStepPar->sortingPar.smoothLMS.p), &(bkwStepPar->sortingPar.smoothLMS.p1), &(bkwStepPar->sortingPar.smoothLMS.p2), &(bkwStepPar->sortingPar.smoothLMS.prev_p1), &(bkwStepPar->sortingPar.smoothLMS.meta_skipped), &(bkwStepPar->sortingPar.smoothLMS.unnatural_selection_ts), &(bkwStepPar->sortingPar.smoothLMS.unnatural_selection_start_index));
        if (2 > bkwStepPar->numPositions || bkwStepPar->numPositions > MAX_SMOOTH_LMS_POSITIONS)
        {
            ASSERT_ALWAYS("unsupported number of positions for smooth LMS");
            return 0; /* error in additional parameters */
        }
        return 1; /* smooth LMS-sorting fully determined */

    case codedBKW:
        sscanf(p, "[%d positions, start index=%d, code=%[^\t\n]]", &bkwStepPar->numPositions, &bkwStepPar->startIndex, codingTypeName);
        codingTypeName[strlen(codingTypeName)-1] = 0; /* sscanf reads to end of line (to include spaces), so discard last character (a right bracket) */
        bkwStepPar->sortingPar.CodedBKW.ct = codingTypeFromString(codingTypeName);
        if (2 > bkwStepPar->numPositions || bkwStepPar->numPositions > MAX_CODED_BKW_POSITIONS)
        {
            ASSERT_ALWAYS("unspported number of positions for coded BKW");
            return 0; /* error in additional parameters */
        }
        return 1; /* Coded BKW-sorting fully determined */

    default:
        ASSERT_ALWAYS("error in format for additional sorting parameters");
        return 0; /* error in additional parameters */
    }

    return 0; /* error in additional parameters */
}

u64 num_categories(lweInstance *lwe, bkwStepParameters *bkwStepPar)
{
    ASSERT(lwe, "no parameter lwe");
    ASSERT(bkwStepPar, "no parameter bkwStepPar");
    u64 numCategories = 0;
    u64 c;
    switch (bkwStepPar->sorting)
    {
    case plainBKW:
        if (bkwStepPar->numPositions == 2)
        {
            return lwe->q * lwe->q;
        }
        if (bkwStepPar->numPositions == 3)
        {
            return lwe->q * lwe->q;    /* for 3-position plain BKW we use meta categories, which discards the last position value */
        }
        ASSERT_ALWAYS("unspported number of positions for plain BKW");
        break;
    case LMS:
        ASSERT(2 <= bkwStepPar->numPositions && bkwStepPar->numPositions <= MAX_LMS_POSITIONS, "unsupported number of positions for LMS");
        if (bkwStepPar->numPositions < 2 || bkwStepPar->numPositions > MAX_LMS_POSITIONS)
        {
            ASSERT_ALWAYS("unhandled parameter set for LMS");
            return 0;
        }
        c = lwe->q / bkwStepPar->sortingPar.LMS.p + 1;
        numCategories = c;
        for (int i=1; i<bkwStepPar->numPositions; i++)
        {
            numCategories *= c;
        }
        return numCategories;
        break;
    case smoothLMS:
        ASSERT(2 <= bkwStepPar->numPositions && bkwStepPar->numPositions <= MAX_SMOOTH_LMS_POSITIONS, "unsupported number of positions for smooth LMS");
        /* TODO - add suitable assert for meta_skipped */
        int meta_skipped = bkwStepPar->sortingPar.smoothLMS.meta_skipped; /* Number of positions to skip when dividing into meta categories */
        int p = bkwStepPar->sortingPar.smoothLMS.p;
        int p1 = bkwStepPar->sortingPar.smoothLMS.p1;
        int p2 = bkwStepPar->sortingPar.smoothLMS.p2;
        ASSERT(lwe->q & 1, "Modulo power of 2 not handled!");
        int q_ = lwe->q%2 == 1 ? (lwe->q+1)/2 : lwe->q/2;
        c = ((2*q_-1) % p) == 0 ? ((2*q_-1) / p) : ((2*q_-1) / p) + 1;
        u64 c1 = ((2*q_-1) % p1) == 0 ? ((2*q_-1) / p1) : ((2*q_-1) / p1) + 1;
        if (bkwStepPar->sortingPar.smoothLMS.prev_p1 == -1)   // first step
        {
            numCategories = 1;
            int lastPosition = MIN(bkwStepPar->numPositions, bkwStepPar->numPositions - meta_skipped + 1);
            for (int i=0; i<lastPosition; i++)
            {
                numCategories *= c;
            }
            if (meta_skipped == 0)
            {
                numCategories *= c1;
            }
        }
        else if (bkwStepPar->startIndex + bkwStepPar->numPositions == lwe->n)     // last step
        {
            q_ = bkwStepPar->sortingPar.smoothLMS.prev_p1;
            u64 c2 = ((2*q_-1) % p2) == 0 ? ((2*q_-1) / p2) : ((2*q_-1) / p2) + 1;
            numCategories = c2;
            for (int i=1; i<bkwStepPar->numPositions - meta_skipped; i++)
            {
                numCategories *= c;
            }
        }
        else     // other steps
        {
            q_ = bkwStepPar->sortingPar.smoothLMS.prev_p1;
            u64 c2 = ((2*q_-1) % p2) == 0 ? ((2*q_-1) / p2) : ((2*q_-1) / p2) + 1;
//        printf("%lld ", c2);
            numCategories = c2;
            int lastPosition = MIN(bkwStepPar->numPositions, bkwStepPar->numPositions - meta_skipped + 1);
            for (int i=1; i<lastPosition; i++)
            {
                numCategories *= c;
//          printf("%lld ", c);
            }
            if (meta_skipped == 0)
            {
                numCategories *= c1;
//          printf("%lld\n", c1);
            }
        }
        return numCategories;
        break;
    case codedBKW:
        switch(bkwStepPar->sortingPar.CodedBKW.ct)
        {
        case blockCode_21:
            ASSERT(2 == bkwStepPar->numPositions && bkwStepPar->numPositions <= MAX_CODED_BKW_POSITIONS, "unsupported number of positions for CodedBKW");
            return lwe->q;
            break;
        case blockCode_31:
            ASSERT(3 == bkwStepPar->numPositions && bkwStepPar->numPositions <= MAX_CODED_BKW_POSITIONS, "unsupported number of positions for CodedBKW");
            return lwe->q;
            break;
        case blockCode_41:
            ASSERT(4 == bkwStepPar->numPositions && bkwStepPar->numPositions <= MAX_CODED_BKW_POSITIONS, "unsupported number of positions for CodedBKW");
            return lwe->q;
            break;
        case concatenatedCode_21_21:
            ASSERT(4 == bkwStepPar->numPositions && bkwStepPar->numPositions <= MAX_CODED_BKW_POSITIONS, "unsupported number of positions for CodedBKW");
            return lwe->q * lwe->q;
            break;
        default:
            ASSERT_ALWAYS("codedBKW parameters to be implemented!");
        }
        return 0;
    default:
        ASSERT_ALWAYS("unhandled parameter set");
        return 0;
    }
    ASSERT_ALWAYS("unhandled parameter set");
    return 0;
}
