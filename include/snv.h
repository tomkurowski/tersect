/*  snv.h

    Copyright (C) 2019 Cranfield University

    Author: Tomasz Kurowski <t.j.kurowski@cranfield.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. */

#ifndef SNV_H
#define SNV_H

#include <stdint.h>

/**
 * Variant codes.
 */
#define V_INDEL 0
#define SNV_A_C 1
#define SNV_A_G 2
#define SNV_A_T 3
#define SNV_C_A 4
#define SNV_C_G 5
#define SNV_C_T 6
#define SNV_G_A 7
#define SNV_G_C 8
#define SNV_G_T 9
#define SNV_T_A 10
#define SNV_T_C 11
#define SNV_T_G 12

extern const char * const variant_format[];

/**
 * @param ref Reference allele
 * @param alt Alt allele
 * @return int single nucleotide variant code
 */
uint8_t snv_type(char ref, char alt);

#endif
