/*  bitarray.h

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

#ifndef BITARRAY_H
#define BITARRAY_H

#include <stdio.h>
#include <stdint.h>

#ifndef WAH_BITARRAY_WORD
#define WAH_BITARRAY_WORD uint64_t
#endif

typedef WAH_BITARRAY_WORD bitarray_word; // type used internally to store the array
extern const uint16_t bitarray_word_capacity; // number of booleans in word

/*
 * The BitArray is a structure meant for compact storage and fast set
 * theoretical operations on sets of boolean values. The structure stores two
 * "sizes": the size element corresponds to the actual number of boolean values
 * stored, while the word_capacity element is an internal variable which tracks
 * the size of allocated storage (*array) in terms of multiples of an unsigned
 * integer type.
 *
 * The start and end masks are used to delimit the valid bits in a sub-bit array
 * extracted from a larger one. Note that the internal *array of such an
 * extracted bit array points to the original array, not a copy.
 */
struct bitarray {
    size_t size; // Size in terms of bitarray_word variables
    size_t last_word; // Position of previously set bit
    size_t ncompressed; // Number of compressed words
    bitarray_word *array;
    bitarray_word start_mask;
    bitarray_word end_mask;
};

/**
 * Used to store start and end position of bit array intervals.
 */
struct bitarray_interval {
    uint64_t start_index;
    uint64_t end_index;
};

/*
 * Initialisation/allocation, zeroing and deallocation routines.
 */
struct bitarray *init_bitarray(uint64_t bit_size);
struct bitarray *copy_bitarray(const struct bitarray *ba);
void clear_bitarray(struct bitarray *ba);
void free_bitarray(struct bitarray *ba);
void bitarray_shrinkwrap(struct bitarray *ba);
void bitarray_resize(struct bitarray *ba, uint64_t new_size);

/*
 * Routines meant for randomising, printing, and writing the contents of bit
 * arrays.
 */
void print_bitarray(const struct bitarray *ba);
void print_set_indices(const struct bitarray *ba);

/*
 * Routines for set theoretical operations.
 */
void bitarray_intersection(const struct bitarray *a,
                           const struct bitarray *b,
                           struct bitarray **out);
void bitarray_difference(const struct bitarray *a,
                         const struct bitarray *b,
                         struct bitarray **out);
void bitarray_symmetric_difference(const struct bitarray *a,
                                   const struct bitarray *b,
                                   struct bitarray **out);
void bitarray_union(const struct bitarray *a,
                    const struct bitarray *b,
                    struct bitarray **out);
uint64_t bitarray_distance(const struct bitarray *a, const struct bitarray *b);
uint64_t bitarray_weight(const struct bitarray *ba);
void bitarray_extract_region(struct bitarray *dest_ba,
                             const struct bitarray *src_ba,
                             const struct bitarray_interval *region);
void bitarray_extract_bins(struct bitarray *dest_bas,
                           const struct bitarray *src_ba,
                           size_t nbins,
                           const struct bitarray_interval *bins);

/*
 * Routines for manipulating individual bits.
 */
int bitarray_set_bit(struct bitarray *bitset, size_t pos);
int bitarray_get_bit(const struct bitarray *ba, size_t pos);
int bitarray_get_set_indices(const struct bitarray *ba,
                             size_t *nset_indices, size_t **set_indices);

#endif
