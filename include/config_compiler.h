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

#ifndef CONFIG_H
#define CONFIG_H

#ifdef __MINGW32__
#define MINGW
#elif defined(__GNUC__)
#define GCC
#else
#error "compiler not determined"
#endif


#if defined(GCC)

#include <time.h>
#define fseeko64 fseek
#define fseeko fseek
#define ftello64 ftell
#define ftello ftell

#define _off64_t __off64_t
// #define mkdir(A, B) mkdir(A)

#define S_IFDIR 16384

#endif

#if defined(MINGW)
#include <stdio.h>
#endif

#endif
