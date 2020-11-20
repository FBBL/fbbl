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

#include "log_utils.h"
#include <stdio.h>
#include <stdarg.h>

void logToScreenWithArgList(int flags, const char *formatString, va_list args)
{

    (void)flags;

    vprintf(formatString, args);

    /* Lunarc machines buffer output to stdout and do not flush
     * if a job exceeds the given time limit.
     * Flush all output to stdout immediately to circumvent
     * this problem. */
    fflush(stdout);
}

void logToFileWithArgList(FILE *logFile, int flags, const char *formatString, va_list args)
{

    (void)flags;

    if (vfprintf(logFile, formatString, args) == EOF)
    {
        fprintf(stderr, "Could not write to log file!\n");
    }

    if (fflush(logFile))
    {
        fprintf(stderr, "Could not flush log file!\n");
    }
}

void logger(FILE *logFile, int flags, const char *formatString, ...)
{

    if (logFile && (flags & LOGFILE))
    {
        va_list args;
        va_start(args, formatString);
        logToFileWithArgList(logFile, flags, formatString, args);
        va_end(args);
    }

    if (flags & LOGSCREEN)
    {
        va_list args;
        va_start(args, formatString);
        logToScreenWithArgList(flags, formatString, args);
        va_end(args);
    }
}

/*void logBuf(FILE *logFile, int flags, const u8 *buf, int len, int gap, int compactifyLeadingZeros)
{
    if (gap == 0)
        gap = 1000000;
    if (compactifyLeadingZeros)
    {
        int i;
        int allZero = 0;
        for (i=0; i<len; i++)
        {
            if (allZero == -1)
                logger(logFile, flags, "%s%02X", (i % gap) == 0 ? " " : "", buf[i]);
            else if (buf[i] == 0)
                allZero++;
            else if (allZero == 1)
            {
                logger(logFile, flags, " 00%02X", buf[i]);
                allZero = -1;
            }
            else if (allZero > 1)
            {
                logger(logFile, flags, " [%d zero bytes] 00%02X", allZero, buf[i]);
                allZero = -1;
            }
        }
        if (allZero > 0)
            logger(logFile, flags, " [%d zero bytes]", allZero);
    }
    else
    {
        int i;
        for (i=0; i<len; i++)
        {
            logger(logFile, flags, "%s%02X", (i % gap) == 0 ? " " : "", buf[i]);
        }
    }
}*/

void logTimeRaw(FILE *logFile, int flags, time_t start)
{
    int t = (int)difftime(time(NULL), start);
    int d = t / 3600 / 24;
    int h = (t % (3600 * 24)) / 3600;
    int min = (t % 3600) / 60;
    int sec = t % 60;
    logger(logFile, flags, "%d day%s %2d h %2d min %2d sec", d, d==1 ? " " : "s", h, min, sec);
}

void logDateRaw(FILE *logFile, int flags)
{
    time_t t = time(NULL);
    struct tm *tblock;
    const char *weekday[] = {"Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"};
    tblock = localtime(&t);
    logger(logFile, flags, "%s %2d:%02d:%02d", weekday[tblock->tm_wday], tblock->tm_hour, tblock->tm_min, tblock->tm_sec);
}

void timeStamp(time_t start)
{
    printf("[");
    logDateRaw(NULL, LOGALL);
    printf(",");
    logTimeRaw(NULL, LOGALL, start);
    printf("] ");
}
