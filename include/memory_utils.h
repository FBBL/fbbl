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

#ifndef MEMORY_UTILS_H
#define MEMORY_UTILS_H

#include <stdlib.h> /* size_t, malloc, free, realloc */
#include <string.h> /* memcpy, memset */

#define MALLOC  malloc
#define FREE    free
#define CALLOC  calloc
#define REALLOC realloc

#define MEMCPY memcpy
#define MEMSET memset
#define MEMXOR memxor
#define MEMOR memor
#define MEMAND memand
#define MEMRND memrnd

#define MEMBITCPY membitcpy
#define MEMBITSET membitset
#define MEMBITXOR membitxor
#define MEMBITOR  membitor
#define MEMBITAND membitand
#define MEMBITRND membitrnd

#endif

