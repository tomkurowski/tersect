/*  distance.c

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

#include "distance.h"

#include "bitarray.h"
#include "tersect_db.h"

#include <getopt.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local flags for distance */
#define JSON_OUTPUT     2
static int local_flags = 0;

/* Argument options without a short equivalent */
#define A_CONTAINS           1000
#define B_CONTAINS           1001
#define A_MATCHLIST_FILE     1002
#define B_MATCHLIST_FILE     1003
#define MATCHLIST_FILE       1004

struct distance_matrix {
    char **row_samples;
    char **col_samples;
    size_t nrows;
    size_t ncols;
    size_t nmatrices;
    uint64_t ***distance;
    uint32_t bin_size;
    bool symmetric;
};

static void usage(FILE *stream)
{
    fprintf(stream,
            "\n"
            "Usage:    tersect dist [options] <db.tsi> [region]...\n"
            "          tersect dist [options] <db.tsi> [-a <ma>] [-b <mb>] [region]...\n"
            "          tersect dist [options] <db.tsi> [--ac <va>] [--bc <vb>] [region]...\n\n"
            "Options:\n"
            "    -a, --a-match STR       name pattern to be matched by samples in set A\n"
            "    -b, --b-match STR       name pattern to be matched by samples in set B\n"
            "    --ac STR                variants required for sample inclusion in set A\n"
            "    --ac STR                variants required for sample inclusion in set B\n"
            "    --a-list-file STR       file containing list of names to include in set A\n"
            "    --b-list-file STR       file containing list of names to include in set B\n"
            "    --list-file STR         file containing list of names to include in any set\n"
            "    -c, --contains STR      variants required for sample inclusion in any set\n"
            "    -m, --match STR         name pattern to be matched by samples in any set\n"
            "    -B, --bin-size INT      size of bins into which the region is split\n"
            "    -h, --help              print this help message\n"
            "    -j, --json              output JSON; implied if match/contains settings for\n"
            "                            set A and set B differ\n"
            "\n");
}

static inline void init_distance_matrix(size_t nmatrices,
                                        uint32_t bin_size,
                                        size_t nrows,
                                        const struct genome *row_samples,
                                        size_t ncols,
                                        const struct genome *col_samples,
                                        struct distance_matrix *matrix)
{
    matrix->bin_size = bin_size;
    // TODO: optimize for partial symmetry, i.e. some but not all shared genomes
    matrix->symmetric = row_samples == col_samples;
    matrix->nmatrices = nmatrices;
    matrix->nrows = nrows;
    matrix->ncols = ncols;

    matrix->row_samples = malloc(nrows * sizeof *matrix->row_samples);
    for (size_t i = 0; i < nrows; ++i) {
        matrix->row_samples[i] = strdup(row_samples[i].name);
    }
    if (matrix->symmetric) {
        matrix->col_samples = matrix->row_samples;
    } else {
        matrix->col_samples = malloc(ncols * sizeof *matrix->col_samples);
        for (size_t i = 0; i < ncols; ++i) {
            matrix->col_samples[i] = strdup(col_samples[i].name);
        }
    }

    matrix->distance = malloc(nmatrices * sizeof *matrix->distance);
    for (size_t i = 0; i < nmatrices; ++i) {
        matrix->distance[i] = malloc(nrows * sizeof **matrix->distance);
        for (size_t j = 0; j < nrows; ++j) {
            matrix->distance[i][j] = calloc(ncols, sizeof ***matrix->distance);
        }
    }
}

static inline void calculate_distance_matrix(size_t nrows,
                                             const struct genome *row_samples,
                                             const struct bitarray *row_bas,
                                             size_t ncols,
                                             const struct genome *col_samples,
                                             const struct bitarray *col_bas,
                                             bool symmetric,
                                             uint64_t * const *output_dist)
{
    for (size_t j = 0; j < nrows; ++j) {
        for (size_t k = symmetric ? j : 0; k < ncols; ++k) {
            uint64_t dist = 0;
            if (row_samples[j].hdr != col_samples[k].hdr) {
                // Only calculating distance for distinct samples
                dist = bitarray_distance(&row_bas[j], &col_bas[k]);
            }
            output_dist[j][k] += dist;
            if (symmetric) {
                output_dist[k][j] += dist;
            }
        }
    }
}

static inline void build_bin_distance_matrix(const tersect_db *tdb,
                                             size_t nrows,
                                             const struct genome *row_samples,
                                             size_t ncols,
                                             const struct genome *col_samples,
                                             uint32_t bin_size,
                                             const struct genomic_interval *region,
                                             struct distance_matrix *matrix)
{
    size_t nbins;
    struct tersect_db_interval *bins;

    tersect_db_get_bin_intervals(tdb, region, bin_size,
                                 &nbins, &bins);

    // The chromosome is the same for all bins
    struct chromosome *chrom = &bins[0].chromosome;

    // Extract bit array intevals
    struct bitarray_interval *ba_intervals = malloc(nbins
                                                    * sizeof *ba_intervals);
    for (size_t i = 0; i < nbins; ++i) {
        ba_intervals[i] = bins[i].interval;
    }

    // Set up bin iterators for each row/column sample
    ba_bin_it **row_its = malloc(nrows * sizeof *row_its);
    for (size_t i = 0; i < nrows; ++i) {
        struct bitarray tmp;
        tersect_db_get_bitarray(tdb, &row_samples[i], chrom, &tmp);
        row_its[i] = init_bitarray_bin_iterator(&tmp, nbins, ba_intervals);
    }
    ba_bin_it **col_its = malloc(ncols * sizeof *col_its);
    for (size_t i = 0; i < ncols; ++i) {
        struct bitarray tmp;
        tersect_db_get_bitarray(tdb, &col_samples[i], chrom, &tmp);
        col_its[i] = init_bitarray_bin_iterator(&tmp, nbins, ba_intervals);
    }

    init_distance_matrix(nbins, bin_size, nrows, row_samples,
                         ncols, col_samples, matrix);

    struct bitarray *row_bas = malloc(nrows * sizeof *row_bas);
    struct bitarray *col_bas = malloc(ncols * sizeof *col_bas);

    for (size_t i = 0; i < nbins; ++i) {
        // Extracting region bitarrays for rows and cols
        for (size_t j = 0; j < nrows; ++j) {
            bitarray_bin_iterator_next(row_its[j], &row_bas[j]);
        }
        for (size_t j = 0; j < ncols; ++j) {
            bitarray_bin_iterator_next(col_its[j], &col_bas[j]);
        }
        // Calculate distances
        calculate_distance_matrix(nrows, row_samples, row_bas,
                                  ncols, col_samples, col_bas,
                                  matrix->symmetric,
                                  matrix->distance[i]);
    }
    free(bins);
    free(row_bas);
    free(col_bas);
    for (size_t i = 0; i < nrows; ++i) {
        free_bitarray_bin_iterator(row_its[i]);
    }
    free(row_its);
    for (size_t i = 0; i < ncols; ++i) {
        free_bitarray_bin_iterator(col_its[i]);
    }
    free(col_its);
    free(ba_intervals);
}

static inline void build_distance_matrix(const tersect_db *tdb,
                                         size_t nrows,
                                         const struct genome *row_samples,
                                         size_t ncols,
                                         const struct genome *col_samples,
                                         size_t nregions,
                                         const struct genomic_interval *regions,
                                         struct distance_matrix *matrix)
{
    struct tersect_db_interval *intervals;

    intervals = malloc(nregions * sizeof *intervals);
    for (size_t i = 0; i < nregions; ++i) {
        tersect_db_get_interval(tdb, &regions[i], &intervals[i]);
    }
    init_distance_matrix(1, 0, nrows, row_samples, ncols, col_samples, matrix);

    struct bitarray *row_bas = malloc(nrows * sizeof *row_bas);
    struct bitarray *col_bas = malloc(ncols * sizeof *col_bas);

    for (size_t i = 0; i < nregions; ++i) {
        // Extracting region bitarrays for rows and cols
        for (size_t j = 0; j < nrows; ++j) {
            struct bitarray tmp;
            tersect_db_get_bitarray(tdb, &row_samples[j],
                                    &intervals[i].chromosome, &tmp);
            bitarray_extract_region(&row_bas[j], &tmp, &intervals[i].interval);
        }
        for (size_t j = 0; j < ncols; ++j) {
            struct bitarray tmp;
            tersect_db_get_bitarray(tdb, &col_samples[j],
                                    &intervals[i].chromosome, &tmp);
            bitarray_extract_region(&col_bas[j], &tmp, &intervals[i].interval);
        }
        // Calculate distances
        calculate_distance_matrix(nrows, row_samples, row_bas,
                                  ncols, col_samples, col_bas,
                                  matrix->symmetric,
                                  matrix->distance[0]);
    }
    free(intervals);
    free(row_bas);
    free(col_bas);
}

static inline void print_distance_matrix_phylip(const struct distance_matrix *matrix)
{
    // TODO: truncate/pad sample name to ten characters
    // Note: the format is only valid for symmetric matrices
    // http://evolution.genetics.washington.edu/phylip/doc/distance.html
    printf("%lu\n", matrix->nrows);
    for (size_t i = 0; i < matrix->nrows; ++i) {
        printf("%s ", matrix->row_samples[i]);
        for (size_t j = 0; j < matrix->ncols; ++j) {
            printf("%"PRIu64, matrix->distance[0][i][j]);
            if (j + 1 < matrix->ncols) {
                printf(" ");
            }
        }
        printf("\n");
    }
}

static inline void print_single_matrix_json(size_t nrows,
                                            size_t ncols,
                                            uint64_t * const *dist)
{
    printf("\t[\n");
    for (size_t i = 0; i < nrows; ++i) {
        printf("\t\t[");
        for (size_t j = 0; j < ncols; ++j) {
            if (j + 1 < ncols) {
                printf("%"PRIu64", ", dist[i][j]);
            } else {
                printf("%"PRIu64, dist[i][j]);
            }
        }
        if (i + 1 < nrows) {
            printf("],\n");
        } else {
            printf("]\n");
        }
    }
    printf("\t]\n");
}

static inline void print_distance_matrix_json(const struct distance_matrix *matrix)
{
    printf("{\n");
    printf("\t\"rows\": [\n");
    for (size_t i = 0; i < matrix->nrows; ++i) {
        if (i + 1 < matrix->nrows) {
            printf("\t\t\"%s\",\n", matrix->row_samples[i]);
        } else {
            printf("\t\t\"%s\"\n", matrix->row_samples[i]);
        }
    }
    printf("\t],\n");
    printf("\t\"columns\": [\n");
    for (size_t i = 0; i < matrix->ncols; ++i) {
        if (i + 1 < matrix->ncols) {
            printf("\t\t\"%s\",\n", matrix->col_samples[i]);
        } else {
            printf("\t\t\"%s\"\n", matrix->col_samples[i]);
        }
    }
    printf("\t],\n");
    printf("\t\"matrix\":\n");
    if (matrix->nmatrices > 1) {
        // Binning
        printf("\t[\n");
        for (size_t i = 0; i < matrix->nmatrices; ++i) {
            print_single_matrix_json(matrix->nrows,
                                     matrix->ncols,
                                     matrix->distance[i]);
            if (i + 1 < matrix->nmatrices) {
                printf(",\n");
            }
        }
        printf("\t]\n");
    } else {
        print_single_matrix_json(matrix->nrows,
                                 matrix->ncols,
                                 matrix->distance[0]);
    }
    printf("}\n");
}

static inline void dealloc_distance_matrix(struct distance_matrix *matrix)
{
    for (size_t i = 0; i < matrix->nrows; ++i) {
        free(matrix->row_samples[i]);
    }
    free(matrix->row_samples);
    if (!matrix->symmetric) {
        for (size_t i = 0; i < matrix->ncols; ++i) {
            free(matrix->col_samples[i]);
        }
        free(matrix->col_samples);
    }
    for (size_t i = 0; i < matrix->nmatrices; ++i) {
        for (size_t j = 0; j < matrix->nrows; ++j) {
            free(matrix->distance[i][j]);
        }
        free(matrix->distance[i]);
    }
    free(matrix->distance);
}

/**
 * Merge the provided sample name match string array and match strings
 * contained in the specified files into a single match string array.
 * Allocated memory to the output
 */
static inline error_t merge_match_queries(size_t nmatch,
                                          const char *const match[],
                                          size_t nfiles,
                                          const char *const filenames[],
                                          size_t *nmerged, char ***merged)
{
    *nmerged = nmatch;
    *merged = calloc(nmatch, sizeof **merged);
    if (*merged == NULL) return E_ALLOC;
    for (size_t i = 0; i < nmatch; ++i) {
        if (match[i] != NULL) {
            (*merged)[i] = strdup(match[i]);
            if ((*merged)[i] == NULL) {
                free(*merged);
                return E_ALLOC;
            }
        } else {
            (*merged)[i] = NULL;
        }
    }

    for (size_t i = 0; i < nfiles; ++i) {
        if (filenames[i] == NULL) continue;
        FILE *fh = fopen(filenames[i], "r");
        if (fh == NULL) return E_DIST_LIST_NOPEN;
        char *line_buffer = NULL;
        size_t buffer_size = 0;
        while (getline(&line_buffer, &buffer_size, fh) != -1) {
            ++(*nmerged);
            *merged = realloc(*merged, *nmerged * sizeof **merged);
            if (*merged == NULL) return E_ALLOC;
            (*merged)[*nmerged - 1] = strndup(line_buffer,
                                              strlen(line_buffer) - 1);
        }
        fclose(fh);
        if (line_buffer != NULL) {
            free(line_buffer);
        }
    }
    return SUCCESS;
}

static inline void free_match_queries(size_t nmerged, char **merged)
{
    for (size_t i = 0; i < nmerged; ++i) {
        if (merged[i] != NULL) {
            free(merged[i]);
        }
    }
    free(merged);
}

error_t tersect_distance(int argc, char **argv)
{
    error_t rc = SUCCESS;
    char *db_filename = NULL;
    char **region_strings = NULL;
    size_t nregions = 0;
    char *a_match = NULL;
    char *a_contains = NULL;
    char *a_matchlist_filename = NULL;
    char *b_match = NULL;
    char *b_contains = NULL;
    char *b_matchlist_filename = NULL;
    char *contains = NULL;
    char *match = NULL;
    char *matchlist_filename = NULL;
    bool binning = false;
    uint32_t bin_size = 0;
    bool symmetric = true;
    static struct option loptions[] = {
        {"a-match", required_argument, NULL, 'a'},
        {"b-match", required_argument, NULL, 'b'},
        {"ac", required_argument, NULL, A_CONTAINS},
        {"bc", required_argument, NULL, B_CONTAINS},
        {"a-list-file", required_argument, NULL, A_MATCHLIST_FILE},
        {"b-list-file", required_argument, NULL, B_MATCHLIST_FILE},
        {"list-file", required_argument, NULL, MATCHLIST_FILE},
        {"contains", required_argument, NULL, 'c'},
        {"match", required_argument, NULL, 'm'},
        {"help", no_argument, NULL, 'h'},
        {"json", no_argument, NULL, 'j'},
        {"bin-size", required_argument, NULL, 'B'},
        {NULL, 0, NULL, 0}
    };
    int c;
    while ((c = getopt_long(argc, argv, ":a:b:B:c:m:hj", loptions, NULL)) != -1) {
        switch(c) {
        case 'a':
            a_match = optarg;
            break;
        case 'b':
            b_match = optarg;
            break;
        case 'B':
            binning = true;
            bin_size = strtol(optarg, NULL, 10);
            local_flags |= JSON_OUTPUT;
            break;
        case A_CONTAINS:
            a_contains = optarg;
            break;
        case B_CONTAINS:
            b_contains = optarg;
            break;
        case A_MATCHLIST_FILE:
            a_matchlist_filename = optarg;
            break;
        case B_MATCHLIST_FILE:
            b_matchlist_filename = optarg;
            break;
        case MATCHLIST_FILE:
            matchlist_filename = optarg;
            break;
        case 'c':
            contains = optarg;
            break;
        case 'm':
            match = optarg;
            break;
        case 'h':
            usage(stdout);
            return SUCCESS;
        case 'j':
            local_flags |= JSON_OUTPUT;
            break;
        default:
            usage(stderr);
            return SUCCESS;
        }
    }
    argc -= optind;
    argv += optind;
    if (!argc) {
        // Missing Tersect index file
        usage(stderr);
        return E_NO_TSI_FILE;
    }
    db_filename = argv[0];
    argc -= 1;
    argv += 1;
    if (argc) {
        region_strings = argv;
        nregions = argc;
    }
    if ((a_match != b_match) || (a_contains != b_contains)
        || (a_matchlist_filename != b_matchlist_filename)) {
        // Force JSON for non-symmetric distance matrices
        local_flags |= JSON_OUTPUT;
        symmetric = false;
    }
    // End parsing options

    tersect_db *tdb = tersect_db_open(db_filename);
    if (tdb == NULL) return E_TSI_NOPEN;

    struct genomic_interval *regions;
    if (nregions) {
        rc = tersect_db_parse_regions(tdb, nregions, region_strings, &regions);
    } else {
        rc = tersect_db_get_regions(tdb, &nregions, &regions);
    }
    if (rc != SUCCESS) goto cleanup_1;

    if (binning && nregions > 1) {
        rc = E_DIST_BIN_REGIONS;
        goto cleanup_2;
    }

    char *contains_queries_a[2] = { contains, a_contains };
    char *contains_queries_b[2] = { contains, b_contains };

    size_t nmatches_a;
    char **merged_matches_a;
    rc = merge_match_queries(2, (const char *const[]){ match, a_match },
                             2, (const char *const[]){ matchlist_filename,
                                                       a_matchlist_filename },
                             &nmatches_a, &merged_matches_a);
    if (rc != SUCCESS) goto cleanup_2;

    size_t nmatches_b;
    char **merged_matches_b;
    rc = merge_match_queries(2, (const char *const[]){ match, b_match },
                             2, (const char *const[]){ matchlist_filename,
                                                       b_matchlist_filename },
                             &nmatches_b, &merged_matches_b);
    if (rc != SUCCESS) goto cleanup_3;

    size_t count_a;
    size_t count_b;
    struct genome *samples_a;
    struct genome *samples_b;

    rc = tersect_db_get_genomes(tdb, nmatches_a, merged_matches_a,
                                2, contains_queries_a,
                                &count_a, &samples_a);
    if (rc != SUCCESS) goto cleanup_4;

    if (symmetric) {
        samples_b = samples_a;
        count_b = count_a;
    } else {
        rc = tersect_db_get_genomes(tdb, nmatches_b, merged_matches_b,
                                    2, contains_queries_b,
                                    &count_b, &samples_b);
        if (rc != SUCCESS) goto cleanup_5;
    }

    struct distance_matrix matrix;

    if (!binning) {
        build_distance_matrix(tdb, count_a, samples_a, count_b, samples_b,
                              nregions, regions, &matrix);
    } else {
        build_bin_distance_matrix(tdb, count_a, samples_a, count_b, samples_b,
                                  bin_size, regions, &matrix);
    }

    if (local_flags & JSON_OUTPUT) {
        print_distance_matrix_json(&matrix);
    } else {
        print_distance_matrix_phylip(&matrix);
    }

    dealloc_distance_matrix(&matrix);

    if (samples_b != samples_a) {
        free(samples_b);
    }
cleanup_5:
    free(samples_a);
cleanup_4:
    free(merged_matches_b);
cleanup_3:
    free(merged_matches_a);
cleanup_2:
    free(regions);
cleanup_1:
    tersect_db_close(tdb);
    return rc;
}
