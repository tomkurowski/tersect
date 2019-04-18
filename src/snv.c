/*  snv.c

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

#include "snv.h"

uint8_t snv_type(char ref, char alt) {
    if (ref == 'A' || ref == 'a') {
        if (alt == 'C' || alt == 'c') {
            return SNV_A_C;
        } else if (alt == 'G' || alt == 'g') {
            return SNV_A_G;
        } else if (alt == 'T' || alt == 't') {
            return SNV_A_T;
        }
    } else if (ref == 'C' || ref == 'c') {
        if (alt == 'A' || alt == 'a') {
            return SNV_C_A;
        } else if (alt == 'G' || alt == 'g') {
            return SNV_C_G;
        } else if (alt == 'T' || alt == 't') {
            return SNV_C_T;
        }
    } else if (ref == 'G' || ref == 'g') {
        if (alt == 'A' || alt == 'a') {
            return SNV_G_A;
        } else if (alt == 'C' || alt == 'c') {
            return SNV_G_C;
        } else if (alt == 'T' || alt == 't') {
            return SNV_G_T;
        }
    } else if (ref == 'T' || ref == 't') {
        if (alt == 'A' || alt == 'a') {
            return SNV_T_A;
        } else if (alt == 'C' || alt == 'c') {
            return SNV_T_C;
        } else if (alt == 'G' || alt == 'g') {
            return SNV_T_G;
        }
    }
    return V_INDEL;
}
