/*  alleles.h

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

#ifndef ALLELES_H
#define ALLELES_H

#include <stdint.h>
#include <string.h>

struct allele {
    uint32_t position;
    char *ref;
    char *alt;
};

/*
** Compare SNVs, first by position then (alphabetically) by alternate allele
** letter. Used for sorting SNVs.
*/
inline int allele_cmp(const struct allele *a, const struct allele *b)
{
    int64_t position_diff;
    if ((position_diff = a->position - b->position)) {
        return position_diff;
    }
    int ref_diff;
    if ((ref_diff = strcmp(a->ref, b->ref))) {
        return ref_diff;
    }
    return strcmp(a->alt, b->alt);
}

#endif
