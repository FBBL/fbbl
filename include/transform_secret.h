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

#ifndef TRANSFORM_SECRET_H
#define TRANSFORM_SECRET_H
#include "lwe_instance.h"

/* transform (known) secret vector s */
/* s -> b - A^{-1}s */
/* \vector{s} -> b - A^{-1}\vector{s} */
void transformSecret(lweInstance *lwe, short *s);

/* inverse transform secret */
/* s -> b - A^{-1}s */
/* \vector{s} -> b - A^{-1}\vector{s} */
void inverseTransformSecret(lweInstance *lwe, short *s);

#endif
