/*  bitarray.c

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

#include "bitarray.h"

#include <inttypes.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#define GROWTH_FACTOR 1.5

static const bitarray_word WORD_MAX = (bitarray_word)~0;
static const bitarray_word MSB = (bitarray_word)1
                                 << (CHAR_BIT * sizeof(bitarray_word) - 1);

// Number of boolean bits in the internal bit array storage type.
const uint16_t bitarray_word_capacity = CHAR_BIT * sizeof(bitarray_word) - 1;

static inline size_t bit_to_word_size(size_t bit_size)
{
    return (bit_size + bitarray_word_capacity - 1) / bitarray_word_capacity;
}

/*
 * Allocate and initialise bit array. All the bits are unset (i.e. set to 0).
 */
struct bitarray *init_bitarray(uint64_t bit_size)
{
    // TODO: set masks to match how region extraction works
    struct bitarray *ba = malloc(sizeof *ba);
    // Rounding up
    ba->size = bit_to_word_size(bit_size);
    ba->last_word = 0;
    ba->ncompressed = 0;
    ba->array = calloc(ba->size, sizeof *(ba->array));
    ba->array[0] = ba->size - 1;
    ba->start_mask = WORD_MAX;
    ba->end_mask = WORD_MAX;
    return ba;
}

/*
 * Duplicate bit array. Needs to be freed manually.
 */
struct bitarray *copy_bitarray(const struct bitarray *ba)
{
    struct bitarray *copied_ba = malloc(sizeof *ba);
    size_t array_size = ba->size * sizeof *(ba->array);
    *copied_ba = (struct bitarray) {
        .size = ba->size,
        .last_word = ba->last_word,
        .ncompressed = ba->ncompressed,
        .array = malloc(array_size),
        .start_mask = ba->start_mask,
        .end_mask = ba->end_mask
    };
    memcpy(copied_ba->array, ba->array, array_size);
    return copied_ba;
}

/*
 * Unsets (zeroes) entire bit array.
 */
void clear_bitarray(struct bitarray *ba)
{
    ba->last_word = 0;
    ba->ncompressed = 0;
    ba->start_mask = WORD_MAX;
    ba->end_mask = WORD_MAX;
    memset(ba->array, 0, ba->size * sizeof *(ba->array));
    ba->array[0] = ba->size - 1;
}

/*
 * Free bit array.
 */
void free_bitarray(struct bitarray *ba)
{
    free(ba->array);
    free(ba);
}

void print_bitarray_word(bitarray_word w)
{
    for (uint64_t i = 0; i < CHAR_BIT * sizeof(bitarray_word); ++i) {
        printf("%"PRIu64, w >> i & 1);
    }
    printf(" (%"PRIu64")\n", w);
}

/*
 * Print bit array (as a series of 0s and 1s) in standard output.
 */
void print_bitarray(const struct bitarray *ba)
{
    for (size_t i = 0; i < ba->size; ++i) {
        printf("%zu:\t", i);
        print_bitarray_word(ba->array[i]);
    }
    printf("SM:\t");
    print_bitarray_word(ba->start_mask);
    printf("EM:\t");
    print_bitarray_word(ba->end_mask);
}

/*
 * Print indices of set (i.e. true) bits in array.
 */
void print_set_indices(const struct bitarray *ba)
{
    size_t nset;
    size_t *set_indices;
    bitarray_get_set_indices(ba, &nset, &set_indices);
    for (size_t i = 0; i < nset; ++i) {
        printf("%lu,", set_indices[i]);
    }
    free(set_indices);
    printf("\n");
}

/**
 * Returns the number of words compressed in a zero fill word at a given
 * position in a bit array.
 */
static inline size_t load_zerofill(const struct bitarray *ba, size_t pos)
{
    if (pos == 0) {
        return ba->start_mask + 1;
    } else if (pos + 1 == ba->size) {
        return ba->end_mask + 1;
    } else {
        return ba->array[pos] + 1;
    }
}

static inline void load_masks(const struct bitarray *a,
                              const struct bitarray *b,
                              struct bitarray *out)
{
    if (a->array[0] & MSB) {
        out->start_mask = a->start_mask;
    } else if (b->array[0] & MSB) {
        out->start_mask = b->start_mask;
    }
    if (a->array[a->size - 1] & MSB) {
        out->end_mask = a->end_mask;
    } else if (b->array[b->size - 1] & MSB) {
        out->end_mask = b->end_mask;
    }
}

/**
 * Either adds a zero fill word corresponding to zf_num compressed zero words
 * at the specified position and increments the position by one, or increases
 * the preceding zero fill word if one exists.
 */
static inline void append_zerofills(struct bitarray *ba, size_t *pos,
                                    size_t zf_num)
{
    if (*pos && !(ba->array[(*pos) - 1] & MSB)) {
        ba->array[(*pos) - 1] += zf_num;
    } else {
        ba->array[(*pos)++] = zf_num - 1;
    }
}

void bitarray_union(const struct bitarray *a, const struct bitarray *b,
                    struct bitarray **out)
{
    *out = init_bitarray((a->size + a->ncompressed + b->size + b->ncompressed)
                         * bitarray_word_capacity);
    load_masks(a, b, *out);

    size_t a_pos = 0;
    size_t a_ncomp = 0; // Number of compressed empty words in a
    size_t b_pos = 0;
    size_t b_ncomp = 0; // Number of compressed empty words in b
    size_t out_pos = 0;

    while (a_pos < a->size || b_pos < b->size) {
        if (!(a->array[a_pos] & MSB)) {
            a_ncomp += load_zerofill(a, a_pos++);
        }
        if (!(b->array[b_pos] & MSB)) {
            b_ncomp += load_zerofill(b, b_pos++);
        }
        if (a_ncomp) {
            if (b_ncomp) {
                // Skipping over the smaller set of compressed words
                size_t to_skip = (a_ncomp > b_ncomp) ? b_ncomp : a_ncomp;
                append_zerofills(*out, &out_pos, to_skip);
                a_ncomp -= to_skip;
                b_ncomp -= to_skip;
            } else {
                // Only b word set
                (*out)->array[out_pos++] = b->array[b_pos++];
                --a_ncomp;
            }
            continue;
        } else if (b_ncomp) {
            // Only a word set
            (*out)->array[out_pos++] = a->array[a_pos++];
            --b_ncomp;
            continue;
        }
        (*out)->array[out_pos++] = a->array[a_pos++] | b->array[b_pos++];
    }
    (*out)->last_word = out_pos - 1;
    bitarray_shrinkwrap(*out);
}

void bitarray_intersection(const struct bitarray *a,
                           const struct bitarray *b,
                           struct bitarray **out)
{
    *out = init_bitarray((a->size + a->ncompressed + b->size + b->ncompressed)
                         * bitarray_word_capacity);
    load_masks(a, b, *out);

    size_t a_pos = 0;
    size_t a_ncomp = 0; // Number of compressed empty words in a
    size_t b_pos = 0;
    size_t b_ncomp = 0; // Number of compressed empty words in b
    size_t out_pos = 0;

    while (a_pos < a->size || b_pos < b->size) {
        if (a_pos < a->size) {
            if (!(a->array[a_pos] & MSB)) {
                a_ncomp += load_zerofill(a, a_pos++);
            }
        }
        if (b_pos < b->size) {
            if (!(b->array[b_pos] & MSB)) {
                b_ncomp += load_zerofill(b, b_pos++);
            }
        }
        if (a_ncomp) {
            if (b_ncomp) {
                // Skipping over the smaller set of compressed words
                size_t to_skip = (a_ncomp > b_ncomp) ? b_ncomp : a_ncomp;
                append_zerofills(*out, &out_pos, to_skip);
                a_ncomp -= to_skip;
                b_ncomp -= to_skip;
            } else {
                // Only b word set
                append_zerofills(*out, &out_pos, 1);
                --a_ncomp;
                ++b_pos;
            }
            continue;
        } else if (b_ncomp) {
            // Only a word set
            append_zerofills(*out, &out_pos, 1);
            --b_ncomp;
            ++a_pos;
            continue;
        }
        bitarray_word res = a->array[a_pos++] & b->array[b_pos++];
        if (res == MSB) {
            // Only MSB set, no intersection
            append_zerofills(*out, &out_pos, 1);
        } else {
            (*out)->array[out_pos++] = res;
        }
    }
    (*out)->last_word = out_pos - 1;
    bitarray_shrinkwrap(*out);
}

void bitarray_difference(const struct bitarray *a, const struct bitarray *b,
                         struct bitarray **out)
{
    *out = init_bitarray((a->size + a->ncompressed + b->size + b->ncompressed)
                         * bitarray_word_capacity);
    load_masks(a, b, *out);

    size_t a_pos = 0;
    size_t a_ncomp = 0; // Number of compressed empty words in a
    size_t b_pos = 0;
    size_t b_ncomp = 0; // Number of compressed empty words in b
    size_t out_pos = 0;

    while (a_pos < a->size || b_pos < b->size) {
        if (!(a->array[a_pos] & MSB)) {
            a_ncomp += load_zerofill(a, a_pos++);
        }
        if (!(b->array[b_pos] & MSB)) {
            b_ncomp += load_zerofill(b, b_pos++);
        }
        if (a_ncomp) {
            if (b_ncomp) {
                // Skipping over the smaller set of compressed words
                size_t to_skip = (a_ncomp > b_ncomp) ? b_ncomp : a_ncomp;
                append_zerofills(*out, &out_pos, to_skip);
                a_ncomp -= to_skip;
                b_ncomp -= to_skip;
            } else {
                // Only b word set
                append_zerofills(*out, &out_pos, 1);
                --a_ncomp;
                ++b_pos;
            }
            continue;
        } else if (b_ncomp) {
            // Only a word set
            (*out)->array[out_pos++] = a->array[a_pos++];
            --b_ncomp;
            continue;
        }
        bitarray_word res = a->array[a_pos++] & ~b->array[b_pos++];
        if (res == 0) {
            // MSB is not set due to ~
            append_zerofills(*out, &out_pos, 1);
        } else {
            (*out)->array[out_pos++] = res | MSB;
        }
    }
    (*out)->last_word = out_pos - 1;
    bitarray_shrinkwrap(*out);
}

void bitarray_symmetric_difference(const struct bitarray *a,
                                   const struct bitarray *b,
                                   struct bitarray **out)
{
    *out = init_bitarray((a->size + a->ncompressed + b->size + b->ncompressed)
                         * bitarray_word_capacity);
    load_masks(a, b, *out);

    size_t a_pos = 0;
    size_t a_ncomp = 0; // Number of compressed empty words in a
    size_t b_pos = 0;
    size_t b_ncomp = 0; // Number of compressed empty words in b
    size_t out_pos = 0;

    while (a_pos < a->size || b_pos < b->size) {
        if (!(a->array[a_pos] & MSB)) {
            a_ncomp += load_zerofill(a, a_pos++);
        }
        if (!(b->array[b_pos] & MSB)) {
            b_ncomp += load_zerofill(b, b_pos++);
        }
        if (a_ncomp) {
            if (b_ncomp) {
                // Skipping over the smaller set of compressed words
                size_t to_skip = (a_ncomp > b_ncomp) ? b_ncomp : a_ncomp;
                append_zerofills(*out, &out_pos, to_skip);
                a_ncomp -= to_skip;
                b_ncomp -= to_skip;
            } else {
                // Only b word set
                (*out)->array[out_pos++] = b->array[b_pos++];
                --a_ncomp;
            }
            continue;
        } else if (b_ncomp) {
            // Only a word set
            (*out)->array[out_pos++] = a->array[a_pos++];
            --b_ncomp;
            continue;
        }
        bitarray_word res = a->array[a_pos++] ^ b->array[b_pos++];
        if (res == 0) {
            append_zerofills(*out, &out_pos, 1);
        } else {
            (*out)->array[out_pos++] = res | MSB;
        }
    }
    (*out)->last_word = out_pos - 1;
    bitarray_shrinkwrap(*out);
}

/**
 * Calculate the Hamming distance (number of different bits) between two bit
 * arrays.
 */
uint64_t bitarray_distance(const struct bitarray *a, const struct bitarray *b)
{
    uint64_t distance = 0;

    size_t a_pos = 0;
    size_t a_ncomp = 0; // Number of compressed empty words in a
    size_t b_pos = 0;
    size_t b_ncomp = 0; // Number of compressed empty words in b

    while (a_pos < a->size || b_pos < b->size) {
        if (!(a->array[a_pos] & MSB)) {
            a_ncomp += load_zerofill(a, a_pos++);
        }
        if (!(b->array[b_pos] & MSB)) {
            b_ncomp += load_zerofill(b, b_pos++);
        }
        if (a_ncomp) {
            if (b_ncomp) {
                // Skipping over the smaller set of compressed words
                size_t to_skip = (a_ncomp > b_ncomp) ? b_ncomp : a_ncomp;
                a_ncomp -= to_skip;
                b_ncomp -= to_skip;
            } else {
                // Only b word set
                distance += __builtin_popcountll(b->array[b_pos++]) - 1;
                --a_ncomp;
            }
            continue;
        } else if (b_ncomp) {
            // Only a word set
            distance += __builtin_popcountll(a->array[a_pos++]) - 1;
            --b_ncomp;
            continue;
        }
        distance += __builtin_popcountll(a->array[a_pos++]
                                            ^ b->array[b_pos++]);
    }

    // Removing bits masked by start
    if ((a->array[0] & MSB) && (b->array[0] & MSB)) {
        distance -= __builtin_popcountll((a->array[0] ^ b->array[0])
                                         & ~a->start_mask);
    } else if (a->array[0] & MSB) {
        distance -= __builtin_popcountll(a->array[0] & ~a->start_mask);
    } else if (b->array[0] & MSB) {
        distance -= __builtin_popcountll(b->array[0] & ~b->start_mask);
    }

    // Removing bits masked by end
    size_t a_last = a->size - 1;
    size_t b_last = b->size - 1;
    if ((a->array[a_last] & MSB) && (b->array[b_last] & MSB)) {
        distance -= __builtin_popcountll((a->array[a_last] ^ b->array[b_last])
                                         & ~a->end_mask);
    } else if (a->array[a_last] & MSB) {
        distance -= __builtin_popcountll(a->array[a_last] & ~a->end_mask);
    } else if (b->array[b_last] & MSB) {
        distance -= __builtin_popcountll(b->array[b_last] & ~b->end_mask);
    }

    return distance;
}

/*
 * Get the number of non-zero bits (Hamming weight) in the bit array.
 */
uint64_t bitarray_weight(const struct bitarray *ba)
{
    uint64_t weight = 0;
    if (ba->size == 1 && (ba->array[0] & MSB)) {
        return __builtin_popcountll(ba->array[0]
                                    & ba->start_mask
                                    & ba->end_mask) - 1;
    }
    if (ba->array[0] & MSB) {
        weight += __builtin_popcountll(ba->array[0] & ba->start_mask) - 1;
    }
    if (ba->array[ba->size - 1] & MSB) {
        weight += __builtin_popcountll(ba->array[ba->size - 1]
                                       & ba->end_mask) - 1;
    }
    for (size_t i = 1; i < ba->size - 1; ++i) {
        if (ba->array[i] & MSB) {
            weight += __builtin_popcountll(ba->array[i]) - 1;
        }
    }
    return weight;
}

static inline void bitarray_resize_internal(struct bitarray *ba,
                                            size_t new_size_words)
{
    size_t old_size = ba->size;
    ba->size = new_size_words;
    ba->array = realloc(ba->array, ba->size * sizeof *(ba->array));
    if (new_size_words > old_size) {
        memset(&ba->array[old_size], 0,
            (ba->size - old_size) * sizeof *(ba->array));
        ba->array[ba->last_word + 1] = ba->size - ba->last_word - 2;
    } else {
        // TODO: Handle shrinking (primarily adjusting last word and end mask)
    }
}

static inline void bitarray_grow(struct bitarray *ba)
{
    bitarray_resize_internal(ba, ba->size * GROWTH_FACTOR);
}

void bitarray_resize(struct bitarray *ba, uint64_t new_size)
{
    size_t new_size_words = bit_to_word_size(new_size);
    if (new_size_words != ba->size) {
        bitarray_resize_internal(ba, new_size_words);
    } else {
        // TODO: Adjust end mask
    }
}

/*
 * Set the bit at the specified position to 1.
 * Returns 0 on success, -1 on failure.
 */
int bitarray_set_bit(struct bitarray *ba, size_t pos)
{
    // TODO: simplify this
    size_t word_pos = pos / bitarray_word_capacity;
    if (ba->last_word + ba->ncompressed > word_pos) {
        // Cannot set bits in words prior to the last_word
        return -1;
    }
    // TODO: may need to handle word_diff being > WORD_MAX/2 (unlikely)
    size_t word_diff = word_pos - ba->last_word - ba->ncompressed;
    if (word_diff == 1) {
        // New bit is in the next word
        if (ba->last_word + 1 >= ba->size) {
            bitarray_grow(ba);
        }
        if (ba->last_word == 0 && !(ba->array[0] & MSB)) {
            ba->array[0] = 0;
        }
        ba->last_word += 1;
    } else if (word_diff > 1) {
        // New bit is further than one word away from the previous bit
        if (ba->last_word == 0 && !(ba->array[0] & MSB)) {
            // Adding compressed gap in first word (edge case) only if empty
            if (ba->last_word + 1 >= ba->size) {
                bitarray_grow(ba);
            }
            ba->array[0] = word_diff - 1;
            ba->ncompressed += word_diff - 1;
            ba->last_word += 1;
        } else {
            // Adding compressed gap
            if (ba->last_word + 2 >= ba->size) {
                bitarray_grow(ba);
            }
            ba->array[ba->last_word + 1] = word_diff - 2;
            ba->ncompressed += word_diff - 2;
            ba->last_word += 2;
        }
    }
    if (!(ba->array[ba->last_word] & MSB)) {
        ba->array[ba->last_word] = MSB;
        if (ba->last_word + 1 < ba->size) {
            /* Except for the last word, store the number of remaining
            empty words in the following word */
            ba->array[ba->last_word + 1] = ba->size
                                                   - ba->last_word - 2;
        }
    }
    ba->array[ba->last_word] |= (bitarray_word)1
                                        << pos % bitarray_word_capacity;
    return 0;
}

/*
 * Get the value (0 or 1) of the bit at a particular position.
 */
int bitarray_get_bit(const struct bitarray *ba, size_t pos)
{
    size_t word_index = pos / bitarray_word_capacity;
    bitarray_word mask = ((bitarray_word)1 << pos % bitarray_word_capacity);
    if (word_index == 0) {
        if (!(ba->array[0] & MSB) || !(ba->array[0] & mask & ba->start_mask)) {
            // Bit in first word and not set or not covered by start mask
            return 0;
        } else {
            return 1;
        }
    }
    if (word_index == ba->size + ba->ncompressed - 1) {
        if (!(ba->array[ba->size - 1] & MSB)
            || !(ba->array[ba->size - 1] & mask & ba->end_mask)) {
            // Bit in last word and not set or not covered by end mask
            return 0;
        } else {
            return 1;
        }
    }
    size_t ncompressed = 0;
    if (!(ba->array[0] & MSB)) {
        // Include the compression at the start
        // (or only part of it, depending on the start mask)
        ncompressed += (ba->array[0] > ba->start_mask) ? ba->start_mask
                                                       : ba->array[0];
    }
    for (size_t i = 1; i < ba->size - 1; ++i) {
        if (i + ncompressed > word_index) {
            // Index was in the compressed interval
            return 0;
        }
        if (i + ncompressed == word_index) {
            return (ba->array[i] & MSB) && (ba->array[i] & mask);
        }
        if (!(ba->array[i] & MSB)) {
            ncompressed += ba->array[i];
        }
    }
    return 0;
}

/*
 * Get array of the indices of set (i.e. true) bits in the bit array.
 * Allocates memory to the first argument; this needs to be freed independently.
 */
int bitarray_get_set_indices(const struct bitarray *ba,
                             size_t *nset_indices, size_t **set_indices)
{
    *set_indices = malloc(bitarray_word_capacity * ba->size
                          * sizeof **set_indices);
    *nset_indices = 0;
    uint64_t ncompressed = 0;
    bitarray_word current_word;
    for (uint64_t i = 0; i < ba->size; ++i) {
        if (!(ba->array[i] & MSB))  {
            // MSB is 0, fill word (run-length of zeroes)
            ncompressed += ba->array[i];
            continue;
        }
        current_word = ba->array[i];
        if (!i) {
            current_word &= ba->start_mask;
        }
        if (i + 1 == ba->size) {
            current_word &= ba->end_mask;
        }
        for (uint16_t j = 0; j < bitarray_word_capacity; ++j) {
            if (current_word & ((bitarray_word)1 << j)) {
                (*set_indices)[(*nset_indices)++] = bitarray_word_capacity
                                                    * (i + ncompressed) + j;
            }
        }
    }
    // Truncate to final size
    *set_indices = realloc(*set_indices,
                           (*nset_indices) * sizeof **set_indices);
    return 0;
}

void bitarray_shrinkwrap(struct bitarray *ba)
{
    ba->size = ba->last_word + 1;
    ba->array = realloc(ba->array, ba->size * sizeof *(ba->array));
    if (!(ba->array[0] & MSB)) {
        ba->start_mask = ba->array[0];
        if (ba->size == 1) {
            ba->end_mask = ba->start_mask;
        }
    }
    if (!(ba->array[ba->size - 1] & MSB)) {
        ba->end_mask = ba->array[ba->size - 1];
    }
}

/**
 * Extracts region from bit array starting from the (*index) word and assuming
 * (*ncompressed) words were compressed (in fill words) before (*index).
 * This is helpful for speeding up the extraction of successive regions from the
 * same bit array, e.g. in binning, as we keep track of the number of traversed
 * words instead of traversing from the beginning for each bin.
 */
static inline void extract_region(struct bitarray *dest_ba,
                                  const struct bitarray *src_ba,
                                  const struct bitarray_interval *region,
                                  size_t *index, size_t *ncompressed)
{
    // The 'internal' indices are in terms of the storage word type
    uint64_t internal_start_index = region->start_index
                                    / bitarray_word_capacity;
    uint64_t internal_end_index = region->end_index
                                  / bitarray_word_capacity;
    // Shifting the internal indices by the number of preceding fill words
    size_t i;
    for (i = *index; i <= internal_end_index; ++i) {
        if (src_ba->array[i] & MSB) continue; // no compression
        if (internal_start_index >= i) {
            if (internal_start_index <= i + src_ba->array[i]) {
                dest_ba->start_mask = src_ba->array[i]
                                      - (internal_start_index - i);
                internal_start_index = i;
            } else {
                internal_start_index -= src_ba->array[i];
            }
        }
        if (internal_end_index <= i + src_ba->array[i]) {
            dest_ba->end_mask = internal_end_index - i;
            internal_end_index = i;
        } else {
            internal_end_index -= src_ba->array[i];
        }
        *ncompressed += src_ba->array[i];
    }
    *index = i;
    dest_ba->size = 1 + internal_end_index - internal_start_index;
    dest_ba->last_word = 0;
    dest_ba->ncompressed = *ncompressed;
    dest_ba->array = &(src_ba->array[internal_start_index]);
    if (src_ba->array[internal_start_index] & MSB) {
        dest_ba->start_mask = WORD_MAX << region->start_index
                                          % bitarray_word_capacity;
    }
    if (src_ba->array[internal_end_index] & MSB) {
        dest_ba->end_mask = WORD_MAX >> (bitarray_word_capacity
                                         - region->end_index
                                           % bitarray_word_capacity)
                            | MSB;
    }
}

/**
 * Extract successive bit array regions (bins) from a single bit array.
 * The start index of each region should be 1 higher than the end index of the
 * preceding region.
 * Does not allocate memory. The caller should ensure the dest_bas output array
 * is at least nbins elements long.
 */
void bitarray_extract_bins(struct bitarray *dest_bas,
                           const struct bitarray *src_ba,
                           size_t nbins,
                           const struct bitarray_interval *bins)
{
    size_t index = 0;
    size_t ncompressed = 0;
    for (size_t i = 0; i < nbins; ++i) {
        extract_region(&dest_bas[i], src_ba, &bins[i], &index, &ncompressed);
    }
}

/**
 * Extract a bit array representing a region of a larger bit array.
 *
 * If an extracted region starts/ends on a literal word, the start/end mask is
 * a normal mask on that word.
 * If an extracted region starts/ends on a fill woed, the start/end mask is the
 * number of words in that fill included in the region minus one.
 * e.g. with a 64-bit word, if the first word is a zero-fill of ten words
 * (i.e. 630 consecutive 0 bits) and the region includes all 630, the start_mask
 * will be 9 (10 - 1).
 */
void bitarray_extract_region(struct bitarray *dest_ba,
                             const struct bitarray *src_ba,
                             const struct bitarray_interval *region)
{
    bitarray_extract_bins(dest_ba, src_ba, 1, region);
}
