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

#ifndef BITSLICE_UTILS_H
#define BITSLICE_UTILS_H
#include "memory_utils_conf.h"
#include "platform_types.h"
#include "log_utils.h"
#include "fill_utils.h"

void toBitslicedKeyIv(u64 *key, int keySizeInBits, const u8 *key8,
                      u64 *iv, int ivSizeInBits, const u8 *iv8,
                      int *keyBit, int numKeyBits,
                      int *ivBit, int numIvBits,
                      int *numBitsclicedKeyBits,
                      int *numBitsclicedIvBits);

void initializeBitslicedKeyValues(u64 *key, int keyLengthInBits, int keyFill, const u8 *explicitKeyFill); /* obsolete, use toBitslicedKeyIv instead */
void initializeBitslicedIvValues(u64 *iv, int ivLengthInBits, const int *ivBit, int ivFill, const u8 *explicitIvFill, const u64 *key); /* obsolete, use toBitslicedKeyIv instead */

#ifndef THREADING_OFF
void fromBitslicedBufByWeight(u8 *dst, const u64 *src, int numBits);
void fromBitslicedBufByLSB(u8 *dst, const u64 *src, int numBits);
void toBitslicedBuf(u64 *dst, const u8 *src, int numBits);

void logBitslicedBufByWeight(FILE *logFile, const u64 *buf, int len, int gap);
void logBitslicedBufByLSB(FILE *logFile, const u64 *buf, int len, int gap);
#endif

void setKeyBits(u64 *key, u64 j, int numKeyBits, const int *keyBit, int numBitsliceBits);
void setIvBits(u64 *iv, u64 i, int numIvBits, const int *ivBit, int numBitsliceBits);

#endif

