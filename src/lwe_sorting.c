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

#include "lwe_sorting.h"
#include "assert_utils.h"
#include <string.h>
#include <stddef.h>

static const char *sorting_label[numSortingMethods] =
{
    "unordered", /* unordered */
    "Plain BKW", /* plainBKW */
    "Lazy Modulus Switching (LMS)", /* LMS */
    "Smooth LMS", /* smooth LMS */
    "Coded BKW" /* coded BKW */
};

const char *sortingAsString(sortingMethod sorting)
{
    return sorting < numSortingMethods ? sorting_label[sorting] : NULL;
}

sortingMethod sortingFromString(const char *str)
{
    for (int i=0; i<numSortingMethods; i++)
    {
        const char *sortingName = sortingAsString(i);
        int len = strlen(sortingName);
        if (!strncmp(sortingName, str, len))
        {
            return i;
        }
    }
    ASSERT_ALWAYS("sorting could not be determined");
    return 0;
}
