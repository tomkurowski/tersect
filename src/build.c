/*  build.c

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

#include "build.h"

#include "hashmap.h"
#include "heap.h"
#include "rename.h"
#include "snv.h"
#include "tersect.h"
#include "tersect_db.h"
#include "tersect_db_internal.h"
#include "vcf_parser.h"

#include <getopt.h>
#include <stdlib.h>
#include <string.h>

/**
 * Maximum number of alleles per chromosome
 * TODO: match this with INITIAL_ALLELE_NUM and expand var_container dynamically
 */
#define MAX_ALLELES 50000000

/**
 * Maximum allele size in base pairs. Used to set up a parsing buffer.
 */
#define MAX_ALLELE_SIZE 20000

#define SAMPLES_PER_FILE 10000

/**
 * Initial size of allele bitarray.
 */
#define INITIAL_ALLELE_NUM 10000

static int tdb_flags = 0;
static int parser_flags = 0;

/**
 * Wrapper for a parser and an associated bit array to record variants present
 * in a specific genome file.
 */
struct parser_wrapper {
    VCF_PARSER parser;
    struct bitarray **ba;
};

static void usage(FILE *stream)
{
    fprintf(stream,
            "\n"
            "Usage:    tersect build [options] <out.tsi> <in1.vcf>...\n\n"
            "Options:\n"
            "    -f, --force             overwrite database file if necessary\n"
            "    -H, --homozygous        include only homozygous variants\n"
            "    -h, --help              print this help message\n"
            "    -n, --name-file         tsv file containing sample names\n"
            "    -t, --types             include snps, indels, or both (default)\n"
            "    -v, --verbose           run in verbose mode\n"
            "\n");
}

/**
 * Comparator for wrapped parsers.
 */
int wrapped_parser_cmp(const void *a, const void *b)
{
    return parser_allele_cmp(&((struct parser_wrapper *)a)->parser,
                             &((struct parser_wrapper *)b)->parser);
}

static error_t import_files(tersect_db *tdb,
                            const int file_num,
                            char **filenames,
                            int parser_flags);
static inline int load_chromosome_queue(const char *chromosome,
                                        int parser_count,
                                        struct parser_wrapper *parsers,
                                        Heap *queue);
static inline int load_next_chromosome_queue(Heap *queue,
                                             char *chromosome,
                                             struct parser_wrapper *parsers,
                                             int parser_num);
static inline uint32_t process_chromosome_queue(tersect_db *tdb, Heap *queue,
                                                struct variant *var_container);

error_t tersect_build_database(int argc, char **argv)
{
    error_t rc = SUCCESS;
    char *db_filename = NULL;
    char *name_filename = NULL;
    static struct option loptions[] = {
        {"force", no_argument, NULL, 'f'},
        {"help", no_argument, NULL, 'h'},
        {"homozygous", no_argument, NULL, 'H'},
        {"name-file", required_argument, NULL, 'n'},
        {"types", required_argument, NULL, 't'},
        {"verbose", no_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}
    };
    int c;
    while ((c = getopt_long(argc, argv, ":fHhn:t:v", loptions, NULL)) != -1) {
        switch(c) {
        case 'f':
            tdb_flags |= TDB_FORCE;
            break;
        case 'h':
            usage(stdout);
            return SUCCESS;
        case 'H':
            parser_flags |= VCF_ONLY_HOMOZYGOUS;
            break;
        case 'n':
            name_filename = optarg;
            break;
        case 't':
            if (!strcmp(optarg, "snps")) {
                parser_flags |= VCF_ONLY_SNPS;
            } else if (!strcmp(optarg, "indels")) {
                parser_flags |= VCF_ONLY_INDELS;
            } else if (!strcmp(optarg, "both")) {
                // Default, no flag needed
            } else {
                usage(stderr);
                return SUCCESS;
            }
            break;
        case 'v':
            tdb_flags |= TDB_VERBOSE;
            break;
        default:
            usage(stderr);
            return SUCCESS;
        }
    }
    argc -= optind;
    argv += optind;
    if (argc) {
        db_filename = argv[0];
        --argc;
        ++argv;
    } else {
        // Missing output filename
        usage(stderr);
        return E_BUILD_NO_OUTNAME;
    }
    if (!argc) {
        return E_BUILD_NO_FILES;
    }
    tersect_db *tdb;
    rc = tersect_db_create(db_filename, tdb_flags, &tdb);
    if (rc != SUCCESS) return rc;
    rc = import_files(tdb, argc, argv, parser_flags);
    tersect_db_close(tdb);
    if (name_filename != NULL) {
        tdb = tersect_db_open(db_filename);
        if (tdb == NULL) return E_TSI_NOPEN;
        rc = tersect_load_name_file(tdb, name_filename);
        tersect_db_close(tdb);
    }
    return rc;
}

/**
 * Builds a database out of k files via a queue-based k-way merge.
 */
static error_t import_files(tersect_db *tdb,
                            const int file_num,
                            char **filenames,
                            int parser_flags)
{
    error_t rc = SUCCESS;
    if (!file_num) return E_BUILD_NO_FILES;
    struct variant *var_container = malloc(MAX_ALLELES * sizeof *var_container);
    if (!var_container) return E_ALLOC;
    struct parser_wrapper *parsers = malloc(file_num * sizeof *parsers);
    if (!parsers) {
        rc = E_ALLOC;
        goto cleanup_1;
    }
    char current_chromosome[MAX_CHROMOSOME_NAME_LENGTH] = "";
    Heap *queue = init_heap(file_num, wrapped_parser_cmp);
    if (!queue) {
        rc = E_ALLOC;
        goto cleanup_2;
    }
    // Open parsers & prepare bit arrays
    HashMap *sample_names = init_hashmap(SAMPLES_PER_FILE * file_num);
    for (int i = 0; i < file_num; ++i) {
        if (init_parser(filenames[i], parser_flags, &parsers[i].parser)
            != VCF_PARSER_INIT_SUCCESS) {
            rc = E_VCF_PARSE_FILE;
            goto cleanup_3;

        }
        parsers[i].ba = malloc(parsers[i].parser.sample_num
                               * sizeof *parsers[i].ba);
        for (size_t j = 0; j < parsers[i].parser.sample_num; ++j) {
            if (hashmap_get(sample_names,
                            parsers[i].parser.samples[j]) != NULL) {
                rc = E_BUILD_DUPSAMPLE;
                goto cleanup_3;
            } else {
                hashmap_insert(sample_names, parsers[i].parser.samples[j],
                               parsers[i].parser.samples[j]);
            }
            parsers[i].ba[j] = init_bitarray(INITIAL_ALLELE_NUM);
            tersect_db_add_genome(tdb, parsers[i].parser.samples[j]);
        }
        goto_next_chromosome(&parsers[i].parser);
    }
    while (load_next_chromosome_queue(queue, current_chromosome,
                                      parsers, file_num)) {
        uint32_t var_count = process_chromosome_queue(tdb, queue,
                                                      var_container);
        // Taking the position of the last variant in the chromosome as proxy
        // for the chromosome size
        tersect_db_add_chromosome(tdb, current_chromosome,
                                  var_container, var_count,
                                  var_container[var_count - 1].position);
        // The bit arrays stored in parsers are larger and get re-used for each
        // chromosome. The interval is used to extract chromosome-specific
        // bit arrays.
        struct bitarray_interval chr_interval = {
            .start_index = 0,
            .end_index = var_count - 1
        };
        for (int i = 0; i < file_num; ++i) {
            for (size_t j = 0; j < parsers[i].parser.sample_num; ++j) {
                struct bitarray ba;
                bitarray_resize(parsers[i].ba[j], var_count);
                bitarray_extract_region(&ba, parsers[i].ba[j], &chr_interval);
                tersect_db_add_bitarray(tdb, parsers[i].parser.samples[j],
                                        current_chromosome, &ba);
                clear_bitarray(parsers[i].ba[j]);
            }
        }
    }
    // Close parsers
    for (int i = 0; i < file_num; ++i) {
        for (size_t j = 0; j < parsers[i].parser.sample_num; ++j) {
            free_bitarray(parsers[i].ba[j]);
        }
        close_parser(&parsers[i].parser);
        free(parsers[i].ba);
    }
    free_heap(queue);
cleanup_3:
    free_hashmap(sample_names);
cleanup_2:
    free(parsers);
cleanup_1:
    free(var_container);
    return rc;
}

static inline int load_next_chromosome_queue(Heap *queue,
                                             char *chromosome,
                                             struct parser_wrapper *parsers,
                                             int parser_num)
{
    // Save currently set chromosome in order to skip it
    char previous_chromosome[MAX_CHROMOSOME_NAME_LENGTH];
    strcpy(previous_chromosome, chromosome);
    strcpy(chromosome, ""); // Clear chromosome
    clear_heap(queue);
    for (int i = 0; i < parser_num; ++i) {
        if (!strcmp(parsers[i].parser.current_chromosome, previous_chromosome)) {
            goto_next_chromosome(&parsers[i].parser);
        }
        if (parsers[i].parser.current_result == ALLELE_NOT_FETCHED) {
            continue; // End of file reached, skip parser
        }
        if (!strlen(chromosome)) {
            // Set next chromosome
            strcpy(chromosome, parsers[i].parser.current_chromosome);
            heap_push(queue, &parsers[i]);
            continue;
        }
        // Go to set chromosome and put parser on heap if found
        while (strcmp(parsers[i].parser.current_chromosome, chromosome)
               && parsers[i].parser.current_result != ALLELE_NOT_FETCHED) {
            goto_next_chromosome(&parsers[i].parser);
        }
        if (!strcmp(parsers[i].parser.current_chromosome, chromosome)) {
            heap_push(queue, &parsers[i]);
        }
    }
    return queue->size;
}

static inline int load_chromosome_queue(const char *chromosome,
                                        int parser_count,
                                        struct parser_wrapper *parsers,
                                        Heap *queue)
{
    clear_heap(queue);
    for (int i = 0; i < parser_count; ++i) {
        if (goto_chromosome(&parsers[i].parser, chromosome) != NULL) {
            heap_push(queue, &parsers[i]);
        }
    }
    return queue->size;
}

static inline uint32_t process_chromosome_queue(tersect_db *tdb, Heap *queue,
                                                struct variant *var_container)
{
    if (!(queue->size)) return 0;
    struct allele previous_allele = {
        .position = 0,
        .ref = calloc(MAX_ALLELE_SIZE + 1, 1),
        .alt = calloc(MAX_ALLELE_SIZE + 1, 1)
    };

    char chromosome[MAX_CHROMOSOME_NAME_LENGTH];
    strcpy(chromosome,
           ((struct parser_wrapper *)heap_peek(queue))->parser.current_chromosome);
    uint32_t var_count = 0;
    while (queue->size) {
        // Alternatively could pop & push, but it is slower
        struct parser_wrapper *pwr = heap_peek(queue);
        if (allele_cmp(&previous_allele, &pwr->parser.current_allele)) {
            if (tersect_db_insert_allele(tdb, &pwr->parser.current_allele,
                                         &var_container[var_count]) == SUCCESS) {
                previous_allele.position = pwr->parser.current_allele.position;
                strcpy(previous_allele.ref, pwr->parser.current_allele.ref);
                strcpy(previous_allele.alt, pwr->parser.current_allele.alt);
                for (size_t i = 0; i < pwr->parser.sample_num; ++i) {
                    if (parser_flags & VCF_ONLY_HOMOZYGOUS) {
                        if (pwr->parser.genotypes[i] != GENOTYPE_HOM_ALT) {
                            continue;
                        }
                    }
                    if (pwr->parser.genotypes[i] != GENOTYPE_HOM_REF) {
                        bitarray_set_bit(pwr->ba[i], var_count);
                    }
                }
                var_count++;
            }
        } else {
            // The same allele as previously
            for (size_t i = 0; i < pwr->parser.sample_num; ++i) {
                if (parser_flags & VCF_ONLY_HOMOZYGOUS) {
                    if (pwr->parser.genotypes[i] != GENOTYPE_HOM_ALT) continue;
                }
                if (pwr->parser.genotypes[i] != GENOTYPE_HOM_REF) {
                    bitarray_set_bit(pwr->ba[i], var_count - 1);
                }
            }
        }
        if ((fetch_next_allele(&pwr->parser) == ALLELE_NOT_FETCHED)
            || strcmp(chromosome, pwr->parser.current_chromosome)) {
            // end of file or end of chromosome
            heap_pop(queue); // Remove parser from queue
        } else {
            sift_down(queue); // Sift new allele down
        }
    }
    free(previous_allele.ref);
    free(previous_allele.alt);
    return var_count;
}
