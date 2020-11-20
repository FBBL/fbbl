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

#ifndef LWE_SELECTION_METHOD_H
#define LWE_SELECTION_METHOD_H

typedef enum
{
    LF1, // linear column combination
    LF2, // all pairs columns combination
#if 0
    LF2_unnaturalSelection, // all pairs columns combination
    LSH_LF2 // all pairs in LSH wedges
#endif
} selectionMethod;

#endif
