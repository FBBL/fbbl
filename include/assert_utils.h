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

#ifndef ASSERT_UTILS_H
#define ASSERT_UTILS_H
#if defined(_DEBUG) || defined(DEBUG)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CODE_BODY(message) \
  const char *file = strrchr((const char *)__FILE__, '\\'); \
  printf("Assert triggered on line %d in file %s: %s\n", __LINE__, file ? (file + 1) : __FILE__, message); \
  exit(EXIT_FAILURE);

#define ASSERT(condition, message) \
  do { \
    if (!(condition)) { \
      CODE_BODY(message) \
    } \
  } while(0)

#define ASSERT_ALWAYS(message) \
  do { \
      CODE_BODY(message) \
  } while(0)

#else

#define ASSERT(condition, message)

#define ASSERT_ALWAYS(message)

#endif

#endif

