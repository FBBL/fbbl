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

#ifndef LWE_SORTING_H
#define LWE_SORTING_H

typedef enum
{
    unordered,
    plainBKW,
    LMS, /* Lazy Modulus Switching */
    smoothLMS, /* Smooth Lazy Modulus Switching */
    codedBKW,
    numSortingMethods
} sortingMethod;

const char *sortingAsString(sortingMethod);
sortingMethod sortingFromString(const char *str);

#endif
