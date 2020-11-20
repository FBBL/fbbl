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

#ifndef PLATFORM_TYPES_H
#define PLATFORM_TYPES_H

#ifdef __BORLANDC__

#ifndef u8
#define u8 unsigned char
#endif
#ifndef s16
#define s16 short
#endif
#ifndef u16
#define u16 unsigned short
#endif
#ifndef s32
#define s32 int
#endif
#ifndef u32
#define u32 unsigned int
#define U32C(x) x##i64
#endif
#ifndef s64
#define s64 __int64
#endif
#ifndef u64
#define u64 unsigned __int64
#define U64C(x) x##ui64
#endif

#else

#define u8 unsigned char
#define s16 short
#define u16 unsigned short
#define s32 int
#define u32 unsigned int
#define U32C(x) x##UL
#define s64 long long
#define u64 unsigned long long
#define U64C(x) x##ULL

#endif

#ifndef BYTE
#define BYTE u8
#endif

#ifndef UINT16
#define UINT16 u16
#endif

#ifndef UINT32
#define UINT32 u32
#endif

#ifndef UINT64
#define UINT64 u64
#endif

#define ASSERT_CONCAT_(a, b) a##b
#define ASSERT_CONCAT(a, b) ASSERT_CONCAT_(a, b)
#define size_check(e) enum { ASSERT_CONCAT(assert_line_, __LINE__) = 1/(!!(e)) }

size_check(sizeof(int)==4);

size_check(sizeof(u8)==1);
size_check(sizeof(u16)==2);
size_check(sizeof(s16)==2);
size_check(sizeof(u32)==4);
size_check(sizeof(s32)==4);
size_check(sizeof(u64)==8);

size_check(sizeof(BYTE)==1);
size_check(sizeof(UINT16)==2);
size_check(sizeof(UINT32)==4);
size_check(sizeof(UINT64)==8);

#endif

