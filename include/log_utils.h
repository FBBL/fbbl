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

#ifndef LOG_UTILS_H
#define LOG_UTILS_H
#include <stdio.h>
#include <time.h>
#include "platform_types.h"

#define LOGSCREEN  1
#define LOGFILE    2
#define LOGFLUSH   4

#define LOGNOFLUSH 3
#define LOGALL     7

void logger(FILE *logFile, int flags, const char *formatString, ...);

// void logBuf(FILE *logFile, int flags, const u8 *buf, int len, int gap, int compactifyLeadingZeros);
void logTimeRaw(FILE *logFile, int flags, time_t start);
void logDateRaw(FILE *logFile, int flags);

void timeStamp(time_t start);

#endif
