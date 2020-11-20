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

#include "rdtsc.h"

//u64(*const rdtsc)(void)=(u64(*const)(void))"\x0f\x31\xc3";
u64 rdtsc(void)
{
    u32 eax, edx;
    __asm__ volatile ( "rdtsc\n\t" : "=a" (eax), "=d" (edx) );
    return (u64)eax | (u64)edx << 32;
}

#ifdef EXAMPLE_OF_USAGE
#include <stdio.h>
#define N 1000000
#define WARMUP_ROUNDS 100
void testRdtsc(void)
{
    u64 c1, c2;
    long i;
    u64 min = 1000000;
    u64 tot = 0;
    u64 overflow = 0;

    c1 = rdtsc();

    for (i=0; i<WARMUP_ROUNDS+N; i++)
    {
        u64 c3, c4, diff;
        c3 = rdtsc();
        c4 = rdtsc();
        if (i >= WARMUP_ROUNDS)
        {
            if (c3 < c4)   /* the expected case */
            {
                diff = c4 - c3;
                tot += diff;
                if (diff < min) min = diff;
            }
            else /* overflow */
                overflow++;
        }
    }

    c2 = rdtsc();

    printf("RDTSC test\n\n");

    printf("Start           : %Ld clock cycles\n", c1);
    printf("End             : %Ld clock cycles\n\n", c2);

    printf("External total  : %Ld clock cycles\n", c2 - c1);
    printf("External average: %f clock cycles\n\n", (c2 - c1)/(double)N);

    printf("Internal total  : %Ld clock cycles\n", tot);
    printf("Internal average: %f clock cycles\n", tot/(double)(N - overflow));
    printf("Minimum         : %Ld clock cycles\n", min);
    printf("Overflow        : %d times\n", overflow);
}
#endif

