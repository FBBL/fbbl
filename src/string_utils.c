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

#include "string_utils.h"
#include "assert_utils.h"
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

char *sprintf_u64_delim(char *s, u64 v)
{
    ASSERT(s, "unexpected parameter usage");
    char vv[256];
    char *p = s;
    char *q = vv;
    sprintf(vv, "%" PRIu64 "", v);
    int len = strlen(vv);
    *p++ = *q++;
    for (int i=1; i<len; i++)
    {
        if ((len-i) % 3 == 0)
        {
            *p++ = ',';
        }
        *p++ = *q++;
    }
    *p = 0;
    return s;
}
