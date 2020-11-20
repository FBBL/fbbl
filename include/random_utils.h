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

#ifndef RANDOM_UTILS_H
#define RANDOM_UTILS_H
#include "platform_types.h"

typedef struct
{
    u64 A1, A2, B1, B2, C1, C2;
    int initialized;
} rand_ctx;

void randomUtilRandomize(void);
void randomUtilInit(rand_ctx *ctx);
u64 randomUtil64(rand_ctx *ctx);
int randomUtilInt(rand_ctx *ctx, int n);
// void randomUtilIntBufUnique(rand_ctx *ctx, int n, int *buf, int len);
double randomUtilDouble(rand_ctx *ctx);
long double randomUtilLongDouble(rand_ctx *ctx);
void randomUtilAppendRandomness(rand_ctx *ctx, u8 *randBuf, int len);

#endif

