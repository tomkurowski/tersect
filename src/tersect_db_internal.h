/*  tersect_db_internal.h

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

#ifndef TERSECT_DB_INTERNAL
#define TERSECT_DB_INTERNAL

#include "hashmap.h"

#include <stdint.h>

typedef uint64_t tdb_offset;

struct tersect_db {
    char *filename;
    HashMap *sequences;
    uintptr_t mapping;
    struct tersect_db_hdr *hdr;
};

struct chrom_hdr {
    tdb_offset name;
    tdb_offset variants;
    tdb_offset bitarrays;
    uint32_t variant_count;
    uint32_t length;
    tdb_offset next;
};

struct genome_hdr {
    tdb_offset name;
    tdb_offset next;
};

struct tersect_db_hdr {
    char format[14]; // TERSECT_FORMAT_VERSION
    uint64_t db_size;
    uint16_t word_size;
    tdb_offset chromosomes;
    uint32_t chromosome_count;
    tdb_offset genomes;
    uint32_t genome_count;
    tdb_offset free_head;
};

struct bitarray_hdr {
    tdb_offset genome_offset;
    size_t size;
    tdb_offset array;
    bitarray_word start_mask;
    bitarray_word end_mask;
    tdb_offset next;
};

struct variant {
    uint32_t position;
    uint8_t type;
    tdb_offset allele;
};

struct genomic_variant {
    char *chromosome;
    uint32_t position;
    char *ref;
    char *alt;
    int type;
};

#endif
