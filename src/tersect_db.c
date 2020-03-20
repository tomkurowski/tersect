/*  tersect_db.c

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

#include "tersect_db.h"
#include "tersect_db_internal.h"

#include "snv.h"
#include "version.h"
#include "vcf_writer.h"

#include <fcntl.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#define PAGE_SIZE 4096
#define INITIAL_DB_SIZE 4096

/* Types used for queries */
#define QUERY_SINGLE_MATCH    0   // Selecting only a single element (exact match)
#define QUERY_WILDCARD_MATCH  1   // Selecting multiple elements
#define QUERY_ALL_MATCH       2   // Selecting all elements

// Capacity of sequence hashmap used during database build
#define SEQUENCE_MAP_CAPACITY 50000000

/**
 * Initialize header in new file. Note that the file size was already set by
 * tersect_db_resize_file.
 *
 * Returns 0 success, 1 on failure.
 */
static error_t tersect_db_init_header(tersect_db *tdb)
{
    strcpy(tdb->hdr->format, TERSECT_FORMAT_VERSION);
    tdb->hdr->free_head = sizeof *tdb->hdr;
    tdb->hdr->chromosome_count = 0;
    tdb->hdr->chromosomes = 0;
    tdb->hdr->genome_count = 0;
    tdb->hdr->genomes = 0;
    tdb->hdr->word_size = CHAR_BIT * sizeof(bitarray_word);
    return SUCCESS;
}

error_t tersect_db_resize_file(tersect_db *tdb, off_t new_size)
{
    if ((void *)tdb->mapping != NULL) {
        // Removing previous mapping
        munmap((void *)tdb->mapping, tdb->hdr->db_size);
    }
    int fd;
    if ((fd = open(tdb->filename, O_RDWR | O_CREAT, 0664)) == -1) return FAILURE;
    if (ftruncate(fd, new_size)) return FAILURE;
    tdb->mapping = (uintptr_t)mmap(NULL, new_size, PROT_READ | PROT_WRITE,
                                   MAP_SHARED, fd, 0);
    if ((void *)tdb->mapping == MAP_FAILED) return FAILURE;
    if (close(fd) == -1) {
        munmap((void *)tdb->mapping, new_size);
        return FAILURE;
    }
    tdb->hdr = (struct tersect_db_hdr *)tdb->mapping;
    tdb->hdr->db_size = new_size;
    return SUCCESS;
}

static tdb_offset tersect_db_malloc(tersect_db *tdb, size_t size)
{
    if (tdb->hdr->free_head + size >= tdb->hdr->db_size) {
        // Rounded up
        size_t page_num = (tdb->hdr->free_head + size + PAGE_SIZE - 1)
                          / PAGE_SIZE;
        if (tersect_db_resize_file(tdb, page_num * PAGE_SIZE)) return 0;
    }
    tdb_offset offset = tdb->hdr->free_head;
    tdb->hdr->free_head += size;
    return offset;
}

/**
 * Verifies file existence and write permissions. Adds .tsi extension if not
 * present. Allocates memory for the output.
 */
static inline error_t validate_filename(const char *filename, int flags,
                                        char **output_filename)
{
    error_t rc;
    size_t original_length = strlen(filename);
    if (!strcmp(&filename[original_length - 4], ".tsi")) {
        // Name already has the extension
        *output_filename = malloc(original_length + 1);
        if (*output_filename == NULL) return E_ALLOC;
        strcpy(*output_filename, filename);
    } else {
        // Adding extension
        *output_filename = malloc(original_length + 5);
        if (*output_filename == NULL) return E_ALLOC;
        sprintf(*output_filename, "%s.tsi", filename);
    }
    if (access(*output_filename, F_OK) == 0) {
        if (flags & TDB_FORCE) {
            if (access(*output_filename, W_OK) == -1) {
                rc = E_BUILD_NO_WRITE;
                goto cleanup;
            }
            // TODO: Warning: overwriting existing file
        } else {
            rc = E_BUILD_DB_EXISTS;
            goto cleanup;
        }
    }
    return SUCCESS;
cleanup:
    free(*output_filename);
    return rc;
}

error_t tersect_db_create(const char *filename, int flags, tersect_db **tdb)
{
    error_t rc = SUCCESS;
    size_t size = INITIAL_DB_SIZE;
    if (size <= sizeof(struct tersect_db_hdr)) {
        size = sizeof(struct tersect_db_hdr) + 1;
    }
    *tdb = malloc(sizeof **tdb);
    if (!(*tdb)) return E_BUILD_CREATE;
    **tdb = (tersect_db) {
        .mapping = 0,
        .sequences = init_hashmap(SEQUENCE_MAP_CAPACITY)
    };
    rc = validate_filename(filename, flags, &(*tdb)->filename);
    if (rc != SUCCESS) {
        goto cleanup_1;
    }
    if (tersect_db_resize_file(*tdb, size)) goto cleanup_2;
    if (tersect_db_init_header(*tdb)) goto cleanup_2;
    return rc;
cleanup_2:
    free((*tdb)->filename);
cleanup_1:
    free_hashmap((*tdb)->sequences);
    free(*tdb);
    return rc;
}

tersect_db *tersect_db_open(const char *filename)
{
    tersect_db *tdb = malloc(sizeof *tdb);
    if (!tdb) return NULL;
    *tdb = (tersect_db) {
        .mapping = 0,
        .sequences = NULL
    };
    if (validate_filename(filename, TDB_FORCE, &tdb->filename) != SUCCESS) {
        goto cleanup_1;
    }
    int fd;
    if ((fd = open(tdb->filename, O_RDWR, 0664)) == -1) goto cleanup_2;
    struct stat st;
    if (fstat(fd, &st) == -1) goto cleanup_3;
    tdb->mapping = (uintptr_t)mmap(NULL, st.st_size, PROT_READ | PROT_WRITE,
                              MAP_SHARED, fd, 0);
    if ((void *)tdb->mapping == MAP_FAILED) goto cleanup_3;
    if (close(fd) == -1) goto cleanup_4;
    tdb->hdr = (struct tersect_db_hdr *)tdb->mapping;
    return tdb;
cleanup_4:
    munmap((void *)tdb->mapping, st.st_size);
cleanup_3:
    close(fd);
cleanup_2:
    free(tdb->filename);
cleanup_1:
    free(tdb);
    return NULL;
}

void tersect_db_close(tersect_db *tdb)
{
    munmap((void *)tdb->mapping, tdb->hdr->db_size);
    if (tdb->sequences != NULL) {
        // Free allelic sequences
        HashIterator it = hashmap_iterator(tdb->sequences);
        hashmap_iterator_next(&it);
        while (it.value) {
            free(it.value);
            hashmap_iterator_next(&it);
        }
        free_hashmap(tdb->sequences);
    }
    free(tdb->filename);
    free(tdb);
}

static tdb_offset tersect_db_add_string(tersect_db *tdb,
                                        const char *string)
{
    size_t len = strlen(string) + 1;
    tdb_offset offset = tersect_db_malloc(tdb, len);
    strcpy((char *)(tdb->mapping + offset), string);
    return offset;
}

static tdb_offset tersect_db_add_variants(tersect_db *tdb,
                                          uint32_t variant_count,
                                          const struct variant *variants)
{
    size_t size = variant_count * sizeof *variants;
    tdb_offset offset = tersect_db_malloc(tdb, size);
    memcpy((void *)(tdb->mapping + offset), variants, size);
    return offset;
}

/**
 * Add the raw content of a bit array to the database.
 */
static tdb_offset tersect_db_add_raw_bitarray(tersect_db *tdb,
                                              const struct bitarray *ba)
{
    size_t size = ba->size * sizeof *ba->array;
    tdb_offset offset = tersect_db_malloc(tdb, size);
    memcpy((void *)(tdb->mapping + offset), ba->array, size);
    return offset;
}

/**
 * Finds genome header by name. Returns NULL if not found.
 */
static struct genome_hdr *tersect_db_find_genome(const tersect_db *tdb,
                                                 const char *genome)
{
    tdb_offset offset = tdb->hdr->genomes;
    while (offset) {
        struct genome_hdr *gen_hdr = (struct genome_hdr *)(tdb->mapping
                                                           + offset);
        if (!strcmp((char *)(tdb->mapping + gen_hdr->name), genome)) {
            return (struct genome_hdr *)(tdb->mapping + offset);
        }
        offset = gen_hdr->next;
    }
    return NULL;
}

/**
 * Finds chromosome header by name. Returns NULL if not found.
 */
static struct chrom_hdr *tersect_db_find_chromosome(const tersect_db *tdb,
                                                    const char *chromosome)
{
    tdb_offset offset = tdb->hdr->chromosomes;
    while (offset) {
        struct chrom_hdr *chr_hdr = (struct chrom_hdr *)(tdb->mapping + offset);
        if (!strcmp((char *)(tdb->mapping + chr_hdr->name), chromosome)) {
            return (struct chrom_hdr *)(tdb->mapping + offset);
        }
        offset = chr_hdr->next;
    }
    return NULL;
}

/**
 * Finds bitarray by genome name and chromosome. Returns NULL if not found.
 */
static struct bitarray_hdr *tersect_db_find_bitarray(const tersect_db *tdb,
                                                     const struct genome *gen,
                                                     const struct chromosome *chr)
{
    tdb_offset genome_offset = (uintptr_t)tersect_db_find_genome(tdb, gen->name)
                               - tdb->mapping;
    tdb_offset offset = chr->hdr->bitarrays;
    while (offset) {
        struct bitarray_hdr *ba_hdr = (struct bitarray_hdr *)(tdb->mapping
                                                              + offset);
        if (ba_hdr->genome_offset == genome_offset) {
            return ba_hdr;
        }
        offset = ba_hdr->next;
    }
    return NULL;
}

error_t tersect_db_insert_allele(tersect_db *tdb, const struct allele *allele,
                                 struct variant *out)
{
    size_t ref_len = strlen(allele->ref);
    size_t alt_len = strlen(allele->alt);
    if (ref_len > 1 || alt_len > 1) {
        // Indel
        char *allele_string = malloc(ref_len + alt_len + 2);
        if (allele_string == NULL) return E_ALLOC;
        sprintf(allele_string, "%s\t%s", allele->ref, allele->alt);
        tdb_offset *allele_offset = (tdb_offset *)hashmap_get(tdb->sequences,
                                                              allele_string);
        if (allele_offset == NULL) {
            // New allelic sequence
            allele_offset = malloc(sizeof *allele_offset);
            *allele_offset = tersect_db_add_string(tdb, allele_string);
            hashmap_insert(tdb->sequences, allele_string, allele_offset);
        }
        *out = (struct variant) {
            .position = allele->position,
            .type = 0,
            .allele = *allele_offset
        };
        free(allele_string);
        return SUCCESS;
    }
    // SNV
    *out = (struct variant) {
        .position = allele->position,
        .type = snv_type(allele->ref[0], allele->alt[0]),
        .allele = 0
    };
    return SUCCESS;
}

void tersect_db_add_bitarray(tersect_db *tdb, const char *genome,
                             const char *chromosome,
                             const struct bitarray *ba)
{
    tdb_offset array_offset = tersect_db_add_raw_bitarray(tdb, ba);
    tdb_offset ba_offset = tersect_db_malloc(tdb, sizeof(struct bitarray_hdr));
    tdb_offset genome_offset = (uintptr_t)tersect_db_find_genome(tdb, genome)
                               - tdb->mapping;
    struct bitarray_hdr *ba_hdr = (struct bitarray_hdr *)(tdb->mapping
                                                          + ba_offset);
    struct chrom_hdr *chr_hdr = tersect_db_find_chromosome(tdb, chromosome);
    *ba_hdr = (struct bitarray_hdr) {
        .genome_offset = genome_offset,
        .size = ba->size,
        .array = array_offset,
        .start_mask = ba->start_mask,
        .end_mask = ba->end_mask,
        .next = chr_hdr->bitarrays
    };
    chr_hdr->bitarrays = ba_offset;
}

void tersect_db_get_bitarray(const tersect_db *tdb,
                             const struct genome *gen,
                             const struct chromosome *chr,
                             struct bitarray *output)
{
    struct bitarray_hdr *ba_hdr = tersect_db_find_bitarray(tdb, gen, chr);
    *output = (struct bitarray) {
        .size = ba_hdr->size,
        .array = (bitarray_word *)(tdb->mapping + ba_hdr->array),
        .start_mask = ba_hdr->start_mask,
        .end_mask = ba_hdr->end_mask
    };
}

void tersect_db_add_chromosome(tersect_db *tdb,
                               const char *chr_name,
                               const struct variant *variants,
                               uint32_t variant_count,
                               uint32_t length)
{
    tdb_offset name_offset = tersect_db_add_string(tdb, chr_name);
    tdb_offset var_offset = tersect_db_add_variants(tdb, variant_count,
                                                    variants);
    tdb_offset chr_offset = tersect_db_malloc(tdb, sizeof(struct chrom_hdr));
    struct chrom_hdr *chr_hdr = (struct chrom_hdr *)(tdb->mapping + chr_offset);
    *chr_hdr = (struct chrom_hdr) {
        .name = name_offset,
        .variants = var_offset,
        .variant_count = variant_count,
        .length = length ? length : variants[variant_count - 1].position,
        .next = tdb->hdr->chromosomes
    };
    tdb->hdr->chromosomes = chr_offset;
    ++tdb->hdr->chromosome_count;
}

void tersect_db_add_genome(tersect_db *tdb, const char *genome_name)
{
    tdb_offset name_offset = tersect_db_add_string(tdb, genome_name);
    tdb_offset genome_offset = tersect_db_malloc(tdb, sizeof(struct genome_hdr));
    struct genome_hdr *gen_hdr = (struct genome_hdr *)(tdb->mapping
                                                       + genome_offset);
    *gen_hdr = (struct genome_hdr) {
        .name = name_offset,
        .next = tdb->hdr->genomes
    };
    tdb->hdr->genomes = genome_offset;
    ++tdb->hdr->genome_count;
}

uint32_t tersect_db_get_genome_count(const tersect_db *tdb)
{
    return tdb->hdr->genome_count;
}

uint32_t tersect_db_get_chromosome_count(const tersect_db *tdb)
{
    return tdb->hdr->chromosome_count;
}

static inline bool wildcard_match(const char *query, const char *pattern)
{
    char *pattern_copy = malloc(strlen(pattern) + 1);
    strcpy(pattern_copy, pattern);
    char *context;
    char *token = strtok_r(pattern_copy, "*", &context);
    while (token != NULL) {
        query = strstr(query, token);
        if (query == NULL) {
            free(pattern_copy);
            return false;
        }
        query += strlen(token);
        token = strtok_r(NULL, "*", &context);
    }
    free(pattern_copy);
    // False if there is trailing text and the last pattern character is not *
    if (strlen(query) && pattern[strlen(pattern) - 1] != '*') return false;
    return true;
}

/**
 * Finds variant by string (e.g. ch02:204:A:G)
 */
error_t parse_variant(struct genomic_variant *output,
                      const tersect_db *tdb,
                      const char *variant)
{
    error_t rc = SUCCESS;
    char *context;
    char *variant_temp = malloc(strlen(variant) + 1);
    if (variant_temp == NULL) {
        rc = E_ALLOC;
        goto cleanup;
    }
    strcpy(variant_temp, variant);
    char *chr_name = strtok_r(variant_temp, ":", &context);
    char *position = strtok_r(NULL, ":", &context);
    char *ref = strtok_r(NULL, ":", &context);
    char *alt = strtok_r(NULL, ":", &context);
    if (chr_name == NULL || position == NULL || ref == NULL || alt == NULL) {
        rc = E_PARSE_ALLELE;
        goto cleanup;
    }
    struct chrom_hdr *chr_hdr = tersect_db_find_chromosome(tdb, chr_name);
    if (chr_hdr == NULL) {
        rc = E_PARSE_ALLELE_NO_CHROMOSOME;
        goto cleanup;
    }
    output->chromosome = (char *)(tdb->mapping + chr_hdr->name);
    char *endptr;
    output->position = strtol(position, &endptr, 10);
    if (endptr == position || *endptr != '\0') {
        rc = E_PARSE_ALLELE;
        goto cleanup;
    }
    output->ref = ref;
    output->alt = alt;
    output->type = snv_type(ref[0], alt[0]);
cleanup:
    if (variant_temp != NULL) {
        free(variant_temp);
    }
    return rc;
}

/**
 * Converts an array of variant strings (each of which may contain several
 * comma-separated variants, e.g. "ch02:100:A:G,ch05:4031:C:T") into an
 * array of strings representing each individual variant.
 */
void flatten_contains_queries(char *const *contains, size_t ncont,
                              char ***out, size_t *nout)
{
    *nout = 0;
    *out = malloc(ncont * sizeof **out);
    size_t out_capacity = ncont;
    for (size_t i = 0; i < ncont; ++i) {
        if (contains[i] == NULL) continue;
        char *variant;
        char *context;
        char *copy = malloc(strlen(contains[i]) + 1);
        strcpy(copy, contains[i]);
        variant = strtok_r(copy, ",", &context);
        while (variant != NULL) {
            (*out)[*nout] = malloc(strlen(variant) + 1);
            strcpy((*out)[*nout], variant);
            ++(*nout);
            if (*nout >= out_capacity) {
                // Ran out of space, expanding
                out_capacity *= 2;
                *out = realloc(*out, out_capacity * sizeof **out);
            }
            variant = strtok_r(NULL, ",", &context);
        }
        free(copy);
    }
    // Clean-up / shrinkwrapping
    if (!(*nout)) {
        free(*out);
        *out = NULL;
    } else {
        *out = realloc(*out, *nout * sizeof **out);
    }
}

/**
 * variant_index and intervals are parallel arrays of size nvars
 * the former represents the position of each respective variant in the latter
 */
error_t parse_contains_queries(const tersect_db *tdb,
                               char *const *contains, size_t ncont,
                               uint64_t **variant_index,
                               struct tersect_db_interval **out_intervals,
                               size_t *nvars)
{
    // TODO: error handling
    error_t rc = SUCCESS;
    char **required_var_strings;
    size_t n_required_vars;
    flatten_contains_queries(contains, ncont, &required_var_strings,
                             &n_required_vars);
    *nvars = 0;
    *variant_index = malloc(n_required_vars * sizeof **variant_index);
    *out_intervals = malloc(n_required_vars * sizeof **out_intervals);
    for (size_t i = 0; i < n_required_vars; ++i) {
        struct genomic_variant variant;
        rc = parse_variant(&variant, tdb, required_var_strings[i]);
        if (rc != SUCCESS) {
            goto cleanup;
        }
        struct genomic_interval gi = {
            .chromosome = variant.chromosome,
            .start_base = variant.position,
            .end_base = variant.position
        };
        struct tersect_db_interval ti;
        tersect_db_get_interval(tdb, &gi, &ti);
        if (ti.interval.start_index > ti.interval.end_index) {
            // site not found
            rc = E_PARSE_ALLELE_UNKNOWN;
            goto cleanup;
        }
        // Checking all variants on the same position
        for (size_t j = ti.interval.start_index;
             j <= ti.interval.end_index; ++j) {
            if (variant.type == ti.chromosome.variants[j].type) {
                // Found variant
                (*variant_index)[*nvars] = j;
                (*out_intervals)[*nvars] = ti;
                ++(*nvars);
                break;
            }
            if (j == ti.interval.end_index) {
                // Variant not in database, so no sample can contain it
                rc = E_PARSE_ALLELE_UNKNOWN;
                goto cleanup;
            }
        }
    }
    // Clean-up / shrinkwrapping
cleanup:
    if (!(*nvars)) {
        free(*variant_index);
        free(*out_intervals);
        *variant_index = NULL;
        *out_intervals = NULL;
    } else {
        *variant_index = realloc(*variant_index,
                                 *nvars * sizeof **variant_index);
        *out_intervals = realloc(*out_intervals,
                                 *nvars * sizeof **out_intervals);
    }
    return rc;
}

/**
 * Returns true if str matches at least one pattern, no patterns are provided,
 * or all patterns are NULL.
 */
bool matches_any_pattern(const char *str,
                         char *const *patterns,
                         size_t npatterns)
{
    if (!npatterns) return true;
    size_t null_patterns = 0;
    for (size_t i = 0; i < npatterns; ++i) {
        if (patterns[i] == NULL) {
            ++null_patterns;
            continue;
        }
        if (wildcard_match(str, patterns[i])) {
            return true;
        }
    }
    if (null_patterns == npatterns) return true;
    return false;
}

/**
 * variant_index and intervals are parallel arrays of size nvars
 * the former represents the position of each respective variant in the latter
 */
bool contains_all_variants(const tersect_db *tdb, const struct genome *gen,
                           const uint64_t *variant_index,
                           const struct tersect_db_interval *intervals,
                           size_t nvars)
{
    struct bitarray ba;
    for (size_t i = 0; i < nvars; ++i) {
        tersect_db_get_bitarray(tdb, gen, &intervals[i].chromosome, &ba);
        if (!bitarray_get_bit(&ba, variant_index[i])) {
            return false;
        }
    }
    return true;
}

error_t tersect_db_get_genomes(const tersect_db *tdb,
                               size_t nmatch, char *const *matches,
                               size_t ncont, char *const *contains,
                               size_t *ngenomes, struct genome **genomes)
{
    error_t rc = SUCCESS;
    struct genome_hdr *gen_hdr = (struct genome_hdr *)(tdb->mapping
                                                       + tdb->hdr->genomes);
    uint64_t *vars_index;
    struct tersect_db_interval *vars_intervals;
    size_t n_contains_vars;
    rc = parse_contains_queries(tdb, contains, ncont, &vars_index,
                                &vars_intervals, &n_contains_vars);
    if (rc != SUCCESS) {
        if (rc == E_PARSE_ALLELE_UNKNOWN) {
            // Variant not in database, so no sample can contain it
            // TODO: perhaps print a warning?
            *ngenomes = 0;
            *genomes = NULL;
            return SUCCESS;
        } else {
            return rc;
        }
    }
    *ngenomes = 0;
    *genomes = malloc(tdb->hdr->genome_count * sizeof **genomes);
    for (uint32_t i = 0; i < tdb->hdr->genome_count; ++i) {
        char *genome_name = (char *)(tdb->mapping + gen_hdr->name);
        (*genomes)[*ngenomes] = (struct genome) {
            .name = genome_name,
            .hdr = gen_hdr
        };
        if ((nmatch == 0 || matches_any_pattern(genome_name, matches, nmatch))
            && contains_all_variants(tdb, &(*genomes)[*ngenomes], vars_index,
                                     vars_intervals, n_contains_vars)) {
            ++(*ngenomes);
        }
        if (gen_hdr->next) {
            gen_hdr = (struct genome_hdr *)(tdb->mapping + gen_hdr->next);
        }
    }
    if (!(*ngenomes)) {
        // No genomes found
        free(*genomes);
        *genomes = NULL;
    } else {
        *genomes = realloc(*genomes, *ngenomes * sizeof **genomes);
    }
    return rc;
}

/**
 * Loads chromosome structure based on a chromosome header.
 */
static inline void load_chromosome(const tersect_db *tdb,
                                   struct chrom_hdr *chr_hdr,
                                   struct chromosome *chrom)
{
    *chrom = (struct chromosome) {
        .name = (char *)(tdb->mapping + chr_hdr->name),
        .length = chr_hdr->length,
        .variant_count = chr_hdr->variant_count,
        .variants = (struct variant *)(tdb->mapping + chr_hdr->variants),
        .hdr = chr_hdr
    };
}

void tersect_db_get_chromosomes(const tersect_db *tdb,
                                size_t *nchroms, struct chromosome **chroms)
{
    struct chrom_hdr *chr_hdr = (struct chrom_hdr *)(tdb->mapping
                                                       + tdb->hdr->chromosomes);
    *nchroms = tdb->hdr->chromosome_count;
    *chroms = malloc(tdb->hdr->chromosome_count * sizeof **chroms);
    for (uint32_t i = tdb->hdr->chromosome_count; i; --i) {
        load_chromosome(tdb, chr_hdr, &(*chroms)[i - 1]);
        if (chr_hdr->next) {
            chr_hdr = (struct chrom_hdr *)(tdb->mapping + chr_hdr->next);
        }
    }
}

void tersect_db_get_chromosome(const tersect_db *tdb, const char *name,
                               struct chromosome *chrom)
{
    struct chrom_hdr *chr_hdr = tersect_db_find_chromosome(tdb, name);
    load_chromosome(tdb, chr_hdr, chrom);
}

bool tersect_db_contains_chromosome(const tersect_db *tdb, const char *name)
{
    return tersect_db_find_chromosome(tdb, name) != NULL;
}

void tersect_db_get_interval(const tersect_db *tdb,
                             const struct genomic_interval *gi,
                             struct tersect_db_interval *ti)
{
    tersect_db_get_chromosome(tdb, gi->chromosome, &ti->chromosome);
    // TODO: binary search instead of linear search for interval boundaries
    for (uint32_t i = 0; i < ti->chromosome.variant_count; ++i) {
        if (ti->chromosome.variants[i].position >= gi->start_base) {
            ti->interval.start_index = i;
            break;
        }
    }
    for (uint32_t i = ti->chromosome.variant_count; i > 0; --i) {
        if (ti->chromosome.variants[i - 1].position <= gi->end_base) {
            ti->interval.end_index = i - 1;
            break;
        }
    }
    ti->nvariants = 1 + ti->interval.end_index - ti->interval.start_index;
    // Rounding down to word index
    ti->variants = &ti->chromosome.variants[(ti->interval.start_index
                                             / bitarray_word_capacity)
                                            * bitarray_word_capacity];
}

void tersect_db_get_bin_intervals(const tersect_db *tdb,
                                  const struct genomic_interval *gi,
                                  uint32_t bin_size,
                                  size_t *nbins,
                                  struct tersect_db_interval **bins)
{
    struct chromosome chrom;
    tersect_db_get_chromosome(tdb, gi->chromosome, &chrom);

    uint32_t region_size = gi->end_base - gi->start_base + 1;
    *nbins = (region_size + bin_size - 1) / bin_size;
    *bins = calloc(*nbins, sizeof **bins);
    for (size_t i = 0; i < *nbins; ++i) {
        (*bins)[i].chromosome = chrom;
    }

    // Finding the first bin
    for (uint32_t i = 0; i < chrom.variant_count; ++i) {
        if (chrom.variants[i].position < gi->start_base) continue;
        (*bins)[0].interval.start_index = i;
        (*bins)[0].variants = &chrom.variants[(i / bitarray_word_capacity)
                                              * bitarray_word_capacity];
        (*bins)[0].nvariants = 1;
        break;
    }
    size_t prev_bin = 0;

    // Went through all variants without finding bin
    if (!(*bins)[0].nvariants) return;

    // Finding further bins
    for (uint32_t i = (*bins)[0].interval.start_index + 1;
         i < chrom.variant_count; ++i) {
        if (chrom.variants[i].position > gi->end_base) {
            // Close last bin
            (*bins)[prev_bin].interval.end_index = i - 1;
            break;
        }
        size_t bin = (chrom.variants[i].position - gi->start_base) / bin_size;
        if (bin != prev_bin) {
            // Close previous bin
            (*bins)[prev_bin].interval.end_index = i - 1;
            if (prev_bin + 1 == *nbins) {
                // Closed last bin, we're done here
                break;
            }
            // Start next bin
            (*bins)[bin].interval.start_index = i;
            (*bins)[bin].variants = &chrom.variants[(i / bitarray_word_capacity)
                                                    * bitarray_word_capacity];
            prev_bin = bin;
        }
        ++((*bins)[bin].nvariants);
    }
    if ((*bins)[prev_bin].interval.end_index == 0) {
        // Last bin was not closed since the last variant was still within it
        (*bins)[prev_bin].interval.end_index = chrom.variant_count - 1;
    }
}

error_t tersect_db_rename_genome(tersect_db *tdb, const char *old_name,
                                 const char *new_name)
{
    struct genome_hdr *gen_hdr = tersect_db_find_genome(tdb, old_name);
    if (gen_hdr == NULL) {
        return E_NO_GENOME;
    }
    tdb_offset new_name_offset = tersect_db_add_string(tdb, new_name);
    // Have to find genome header again as adding the string can update mapping
    gen_hdr = tersect_db_find_genome(tdb, old_name);
    gen_hdr->name = new_name_offset;
    return SUCCESS;
}

error_t parse_region(struct genomic_interval *output, const tersect_db *tdb,
                     const char *region)
{
    error_t rc = SUCCESS;
    char *context;
    char *region_temp = malloc(strlen(region) + 1);
    if (region_temp == NULL) {
        rc = E_ALLOC;
        goto cleanup;
    }
    strcpy(region_temp, region);
    char *chr_name = strtok_r(region_temp, ":", &context);
    char *bounds = strtok_r(NULL, ":", &context);
    struct chrom_hdr *chr_hdr = tersect_db_find_chromosome(tdb, chr_name);
    if (chr_hdr == NULL) {
        rc = E_PARSE_REGION_NO_CHROMOSOME;
        goto cleanup;
    }
    output->chromosome = (char *)(tdb->mapping + chr_hdr->name);
    if (bounds == NULL) {
        // No bounds, interval covers entire chromosome
        output->start_base = 1;
        output->end_base = chr_hdr->length;
    } else {
        char *start_base = strtok_r(bounds, "-", &context);
        char *end_base = strtok_r(NULL, "-", &context);
        if (start_base == NULL || end_base == NULL) {
            rc = E_PARSE_REGION_BAD_BOUNDS;
            goto cleanup;
        }
        char *endptr;
        output->start_base = strtol(start_base, &endptr, 10);
        if (endptr == start_base || *endptr != '\0') {
            rc = E_PARSE_REGION_BAD_BOUNDS;
            goto cleanup;
        }
        output->end_base = strtol(end_base, &endptr, 10);
        if (endptr == end_base || *endptr != '\0') {
            rc = E_PARSE_REGION_BAD_BOUNDS;
            goto cleanup;
        }
    }
cleanup:
    if (region_temp != NULL) {
        free(region_temp);
    }
    return rc;
}

error_t tersect_db_parse_regions(const tersect_db *tdb, size_t nregions,
                                 char **region_strings,
                                 struct genomic_interval **output)
{
    error_t rc = SUCCESS;
    *output = malloc(nregions * sizeof **output);
    for (size_t i = 0; i < nregions; ++i) {
        rc = parse_region(&(*output)[i], tdb, region_strings[i]);
        if (rc != SUCCESS) {
            free(*output);
            return rc;
        }
    }
    return rc;
}

/**
 * Returns all regions (covering entire chromosomes) in the database.
 */
error_t tersect_db_get_regions(const tersect_db *tdb,
                               size_t *nregions,
                               struct genomic_interval **output)
{
    error_t rc = SUCCESS;
    struct chromosome *chroms;
    tersect_db_get_chromosomes(tdb, nregions, &chroms);
    *output = malloc(*nregions * sizeof **output);
    for (size_t i = 0; i < *nregions; ++i) {
        (*output)[i] = (struct genomic_interval) {
            .chromosome = chroms[i].name,
            .start_base = 1,
            .end_base = chroms[i].length
        };
    }
    free(chroms);
    return rc;
}
