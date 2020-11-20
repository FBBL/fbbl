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

#include <stdlib.h>
#include <time.h>
#include "random_utils.h"
#include "assert_utils.h"
#include "rdtsc.h"

static rand_ctx state;

#define SHIFT(REG, steps) ((REG##1 << (64 - (steps))) | (REG##2 >> (steps)))
#define MUL(REG) (SHIFT(REG, 1) & SHIFT(REG, 2))

#define S66(ctx) SHIFT(ctx->A, 93 - 66)
#define S69(ctx) SHIFT(ctx->A, 93 - 69)
#define S93(ctx) ctx->A2

#define S162(ctx) SHIFT(ctx->B, 177 - 162)
#define S171(ctx) SHIFT(ctx->B, 177 - 171)
#define S177(ctx) ctx->B2

#define S243(ctx) SHIFT(ctx->C, 288 - 243)
#define S264(ctx) SHIFT(ctx->C, 288 - 264)
#define S288(ctx) ctx->C2

#define PRE_R(ctx, t1, t2, t3) \
    t1  = S93(ctx); t2  = S177(ctx); t3  = S243(ctx); \
    t1 ^= S66(ctx); t2 ^= S162(ctx); t3 ^= S288(ctx);

#define POST_R(ctx, t1, t2, t3) \
    t1 ^= MUL(ctx->A); t2 ^= MUL(ctx->B); t3 ^= MUL(ctx->C); \
    t1 ^= S171(ctx);   t2 ^= S264(ctx);   t3 ^= S69(ctx); \
    ctx->A2 = (t3 << 29) | ctx->A1; ctx->A1 = t3 >> (64 - 29); \
    ctx->B2 = (t1 << 20) | ctx->B1; ctx->B1 = t1 >> (64 - 20); \
    ctx->C2 = (t2 << 47) | ctx->C1; ctx->C1 = t2 >> (64 - 47);

#define ROUND_WITH_OUTPUT(ctx, t1, t2, t3, u) \
    PRE_R(ctx, t1, t2, t3) u = t1 ^ t2 ^ t3; POST_R(ctx, t1, t2, t3)

#define ROUND_WITHOUT_OUTPUT(ctx, t1, t2, t3) \
    PRE_R(ctx, t1, t2, t3) POST_R(ctx, t1, t2, t3)

void randomUtilRandomize(void)
{
    //randomize();
    u64 t = (u64)time(NULL);
    u64 c = (u64)clock();
    u64 r = rdtsc();
//  printf("t = %Ld\n", t);
//  printf("c = %Ld\n", c);
//  printf("r = %Ld\n", r);
    u64 seed = t ^ c ^ r;
    srand(seed);
}

static u64 random64(void)
{
    u64 t;
    u8 *p = (u8*)&t;
    *p++ = rand() & 0xFF;
    *p++ = rand() & 0xFF;
    *p++ = rand() & 0xFF;
    *p++ = rand() & 0xFF;
    *p++ = rand() & 0xFF;
    *p++ = rand() & 0xFF;
    *p++ = rand() & 0xFF;
    *p = rand() & 0xFF;
    t ^= rdtsc();
    return t;
}

void randomUtilAppendBadRandom(rand_ctx *ctx)
{
    if (!ctx)
    {
        ctx = &state;
    }
    ctx->A1 ^= random64(); // retain register content in case appendRandomness was used
    ctx->A2 ^= random64();
    ctx->B1 ^= random64();
    ctx->B2 ^= random64();
    ctx->C1 ^= random64();
    ctx->C2 ^= random64();
}

void randomUtilInit(rand_ctx *ctx)
{
    u64 t1, t2, t3;
    int i;
    randomUtilRandomize();
    if (!ctx)
    {
        ctx = &state;
    }
    if (ctx != &state && state.initialized)
    {
        ctx->A1 ^= randomUtil64(NULL);
        ctx->A2 ^= randomUtil64(NULL);
        ctx->B1 ^= randomUtil64(NULL);
        ctx->B2 ^= randomUtil64(NULL);
        ctx->C1 ^= randomUtil64(NULL);
        ctx->C2 ^= randomUtil64(NULL);
    }
    else
    {
        randomUtilAppendBadRandom(ctx);
    }
    for (i=0; i<18; i++)
    {
        ROUND_WITHOUT_OUTPUT(ctx, t1, t2, t3)
    }
    ctx->initialized = 1;
}

u64 randomUtil64(rand_ctx *ctx)
{
    u64 t1, t2, t3, u;
    if (!ctx)
    {
        ctx = &state;
    }
    if (ctx == &state && !state.initialized)
    {
        randomUtilInit(&state);
    }
    ASSERT(ctx->initialized, "Not initialized!\n");
    ROUND_WITH_OUTPUT(ctx, t1, t2, t3, u)
    return u;
}

int randomUtilInt(rand_ctx *ctx, int n)
{
    int u = randomUtil64(ctx) % n;
    ASSERT(u >= 0, "Unexpectedly low value!\n");
    ASSERT(u < n, "Unexpectedly high value!\n");
    return u;
}

// static int intInList(int *buf, int len, int num)
// {
//     for (int i=0; i<len; i++)
//     {
//         if (buf[i] == num) return 1;
//     }
//     return 0;
// }

// /* fill buf with unique numbers in [0,n-1] */
// void randomUtilIntBufUnique(rand_ctx *ctx, int n, int *buf, int len)
// {
//     ASSERT(len <= n, "Cannot create so many unique numbers!\n");
//     for (int i=0; i<len; i++)
//     {
//         int num;
//         do
//         {
//             num = randomUtilInt(ctx, n);
//         }
//         while (intInList(buf, i-1, num));
//         buf[i] = num;
//     }
// }

double randomUtilDouble(rand_ctx *ctx)
{
    double d, maskValue = (double)0x00FFFFFFFFFFFF;
    u64 r = randomUtil64(ctx) & 0x00FFFFFFFFFFFF;
    d = r / maskValue;
    ASSERT(d >= 0.0, "Value too small!\n");
    ASSERT(d <= 1.0, "Value too large!\n");
    return d;
}

long double randomUtilLongDouble(rand_ctx *ctx)
{
    long double d, maskValue = (long double)0x00FFFFFFFFFFFF;
    u64 r = randomUtil64(ctx) & 0x00FFFFFFFFFFFF;
    d = r / maskValue;
    ASSERT(d >= 0.0, "Value too small!\n");
    ASSERT(d <= 1.0, "Value too large!\n");
    return d;
}

static void appendRandom(u64 *reg, u8 *rand, int len)
{
    u8 *p = (u8*)reg;
    if (len <= 0) return;
    if (len > 8) len = 8;
    while (len-- > 0) *p++ ^= *rand++;
}

void randomUtilAppendRandomness(rand_ctx *ctx, u8 *randBuf, int len)
{
    if (!ctx)
    {
        ctx = &state;
    }
    if (randBuf)
    {
        appendRandom(&ctx->A1, randBuf, len);
        appendRandom(&ctx->A2, randBuf, len - 8);
        appendRandom(&ctx->B1, randBuf, len - 16);
        appendRandom(&ctx->B2, randBuf, len - 24);
        appendRandom(&ctx->C1, randBuf, len - 32);
        appendRandom(&ctx->C2, randBuf, len - 40);
    }
    else
    {
        randomUtilAppendBadRandom(ctx);
    }
}


