/*  tersect_db.h

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

#ifndef TERSECT_DB
#define TERSECT_DB

#include "alleles.h"
#include "bitarray.h"
#include "errorc.h"

#include <stdbool.h>
#include <stdint.h>

/* Database operation flags */
#define TDB_FORCE       2
#define TDB_VERBOSE     4

/* Opaque header handles */
typedef struct chrom_hdr chrom_hdr;
typedef struct genome_hdr genome_hdr;
struct variant;
struct tersect_db_interval;

struct chromosome {
    char *name;
    uint32_t length;
    uint32_t variant_count;
    struct variant *variants;
    chrom_hdr *hdr;
};

struct genome {
    char *name;
    genome_hdr *hdr;
};

/**
 * Stores a genomic interval as stored in the database - by the chromosome
 * object and bit array indices it represents.
 */
struct tersect_db_interval {
    struct chromosome chromosome;
    size_t nvariants;
    struct variant *variants;
    struct bitarray_interval interval;
};

/**
 * Stores a genomic interval in terms of the name of the chromosome and the
 * start and end positions in bases (inclusive on both sides).
 */
struct genomic_interval {
    char *chromosome;
    uint32_t start_base;
    uint32_t end_base;
};

typedef struct tersect_db tersect_db;

error_t tersect_db_create(const char *filename, int flags, tersect_db **tdb);
tersect_db *tersect_db_open(const char *filename);
void tersect_db_close(tersect_db *tdb);
error_t tersect_db_insert_allele(tersect_db *tdb, const struct allele *allele,
                                 struct variant *out);
void tersect_db_add_chromosome(tersect_db *tdb,
                               const char *chr_name,
                               const struct variant *variants,
                               uint32_t variant_num,
                               uint32_t length);
void tersect_db_add_genome(tersect_db *tdb,
                           const char *genome_name);
uint32_t tersect_db_get_genome_count(const tersect_db *tdb);
uint32_t tersect_db_get_chromosome_count(const tersect_db *tdb);
error_t tersect_db_get_genomes(const tersect_db *tdb,
                               size_t nmatch, char *const *matches,
                               size_t ncont, char *const *contains,
                               size_t *ngenomes, struct genome **genomes);
void tersect_db_add_bitarray(tersect_db *tdb, const char *genome,
                             const char *chromosome, const struct bitarray *ba);
void tersect_db_get_bitarray(const tersect_db *tdb,
                             const struct genome *gen,
                             const struct chromosome *chr,
                             struct bitarray *output);
void tersect_db_get_chromosomes(const tersect_db *tdb,
                                size_t *nchroms, struct chromosome **chroms);
void tersect_db_get_chromosome(const tersect_db *tdb, const char *name,
                               struct chromosome *chrom);
bool tersect_db_contains_chromosome(const tersect_db *tdb, const char *name);
void tersect_db_get_interval(const tersect_db *tdb,
                             const struct genomic_interval *gi,
                             struct tersect_db_interval *ti);
void tersect_db_get_bin_intervals(const tersect_db *tdb,
                                  const struct genomic_interval *gi,
                                  uint32_t bin_size,
                                  size_t *nbins,
                                  struct tersect_db_interval **bins);
error_t tersect_db_rename_genome(tersect_db *tdb, const char *old_name,
                                 const char *new_name);

/**
 * Extracts a genomic interval structure from a region string.
 * Valid region strings take the following forms:
 *
 *  chromosome:start-end        (e.g. 'ch1:1-10000' for the interval from base 1
 *                               to base 10000 on chromosome 'ch1')
 *
 *  chromosome                  (e.g. 'ch1' for an interval covering the
 *                               entirety of chromosome 'ch1')
 *
 * The position bounds are inclusive on both sides.
 */
error_t tersect_db_parse_regions(const tersect_db *tdb, size_t nregions,
                                 char **region_strings,
                                 struct genomic_interval **output);

error_t tersect_db_get_regions(const tersect_db *tdb,
                               size_t *nregions,
                               struct genomic_interval **output);

#endif
