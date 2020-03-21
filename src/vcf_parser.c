/*  vcf_parser.c

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

#include "vcf_parser.h"

#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define VCF_NUM_COLUMNS     9

#define CHROM_COLUMN        0
#define POS_COLUMN          1
#define ID_COLUMN           2
#define REF_COLUMN          3
#define ALT_COLUMN          4
#define QUAL_COLUMN         5
#define FILTER_COLUMN       6
#define INFO_COLUMN         7
#define FORMAT_COLUMN       8

/**
 * Size of the header line prior to the genotype list (46 characters).
 * Sample names start past this position.
 *
 * #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t
 *
 */
#define HEADER_LINE_SIZE    46

/**
 * Parsers metadata lines until the one starting with #CHROM
 */
static inline void load_metadata(VCF_PARSER *parser)
{
    parser->samples = malloc(sizeof *parser->samples);
    char *line_context;
    while (getline(&parser->line_buffer, &parser->buffer_size,
           parser->file_handle) != -1) {
        if (!strncmp(parser->line_buffer, "#CHROM", 6)) {
            char *sample_columns = &parser->line_buffer[HEADER_LINE_SIZE];
            char *sample_name = strtok_r(sample_columns, "\t\n", &line_context);
            parser->samples[0] = strdup(sample_name);
            parser->sample_num = 1;
            while ((sample_name = strtok_r(NULL, "\t\n", &line_context))
                   != NULL) {
                ++parser->sample_num;
                parser->samples = realloc(parser->samples,
                                          parser->sample_num
                                          * sizeof *parser->samples);
                parser->samples[parser->sample_num - 1] = strdup(sample_name);
            }
            break;
        }
    }
}

static inline int get_genotype_code(char *gt)
{
    if (gt[0] == '0' && gt[2] == '0') {
        return GENOTYPE_HOM_REF;
    } else if (gt[0] == gt[2]) {
        return GENOTYPE_HOM_ALT;
    } else {
        return GENOTYPE_HET;
    }
}

int is_gzipped(const char *filename)
{
    FILE *fh = fopen(filename, "rb");
    char byte1 = getc(fh);
    char byte2 = getc(fh);
    fclose(fh);
    return byte1 == (char)0x1f && byte2 == (char)0x8b;
}

static inline void reset_vcf_position(VCF_PARSER *parser)
{
    strcpy(parser->current_chromosome, "");
    parser->current_allele_index = 0;
    if (parser->line_buffer != NULL) {
        free(parser->line_buffer);
    }
    parser->line_buffer = NULL;
    parser->buffer_size = 0;
    parser->n_alts = 0;
    parser->current_result = ALLELE_NOT_FETCHED;
}

static inline int open_vcf_file(VCF_PARSER *parser, const char *filename)
{
    if (access(filename, F_OK) != 0) {
        return VCF_PARSER_INIT_FAILURE;
    }
    if (is_gzipped(filename)) {
        char *cmd = malloc(strlen(filename) + 8); // strlen of "zcat \'%s\'"
        sprintf(cmd, "zcat \'%s\'", filename);
        parser->file_handle = popen(cmd, "r");
        free(cmd);
    } else {
        parser->file_handle = fopen(filename, "r");
    }
    if (parser->file_handle == NULL) {
        return VCF_PARSER_INIT_FAILURE;
    }
    strcpy(parser->filename, filename);
    parser->line_buffer = NULL;
    reset_vcf_position(parser);
    return VCF_PARSER_INIT_SUCCESS;
}

static inline void close_vcf_file(VCF_PARSER *parser)
{
    if (parser->line_buffer != NULL) {
        free(parser->line_buffer);
    }
    fclose(parser->file_handle);
}

static inline void rewind_vcf_file(VCF_PARSER *parser)
{
    if (is_gzipped(parser->filename)) {
        close_vcf_file(parser);
        open_vcf_file(parser, parser->filename);
    } else {
        reset_vcf_position(parser);
        rewind(parser->file_handle);
    }
}

int init_parser(const char *filename, int flags, VCF_PARSER *parser)
{
    if (open_vcf_file(parser, filename) != VCF_PARSER_INIT_SUCCESS) {
        return VCF_PARSER_INIT_FAILURE;
    }
    parser->flags = flags;
    load_metadata(parser);
    parser->genotypes = malloc(parser->sample_num * sizeof *parser->genotypes);
    if (parser->genotypes == NULL) {
        return VCF_PARSER_INIT_FAILURE;
    }
    parser->chromosome_names = init_stringset();
    return VCF_PARSER_INIT_SUCCESS;
}

/* Compare most recent alleles from two parsers */
int parser_allele_cmp(const void *a, const void *b)
{
    VCF_PARSER *parser_a = (VCF_PARSER *)a;
    VCF_PARSER *parser_b = (VCF_PARSER *)b;
    return allele_cmp(&parser_a->current_allele, &parser_b->current_allele);
}

/* Compare alleles (for sorting them alphabetically) */
static inline int alt_comp(const void *a, const void *b)
{
    return strcmp(*(char *const *)b, *(char *const *)a);
}

int fetch_next_allele(VCF_PARSER *parser)
{
    // TODO: add error handling in case of incorrect file contents
    // Fetching successive ALT alleles at the same position
    while (parser->n_alts) {
        // Pop ALT allele (LIFO)
        parser->current_allele.alt = parser->alt_alleles[--(parser->n_alts)];
        if (parser->flags & VCF_ONLY_SNPS) {
            if (parser->current_allele.alt[1] != '\0') continue;
        } else if (parser->flags & VCF_ONLY_INDELS) {
            if (parser->current_allele.alt[1] == '\0'
                && parser->current_allele.ref[1] == '\0') continue;
        }
        // Increase allele index (ordinal position of allee in chromosome)
        ++(parser->current_allele_index);
        return parser->current_result = ALLELE_FETCHED;
    }
    char *columns[VCF_NUM_COLUMNS];
    char *line_context;
    while (getline(&parser->line_buffer,
                   &parser->buffer_size,
                   parser->file_handle) != -1) {
        if (parser->line_buffer[0] != '#') {
            columns[CHROM_COLUMN] = strtok_r(parser->line_buffer, "\t",
                                             &line_context);
            for (int i = 1; i < VCF_NUM_COLUMNS; ++i) {
                columns[i] = strtok_r(NULL, "\t", &line_context);
            }
            if (strcmp(parser->current_chromosome, columns[CHROM_COLUMN])) {
                // New chromosome
                strcpy(parser->current_chromosome, columns[CHROM_COLUMN]);
                stringset_add(parser->chromosome_names, columns[CHROM_COLUMN]);
                // Reset allele index
                parser->current_allele_index = 0;
            }
            parser->current_allele.position = atoi(columns[POS_COLUMN]);
            parser->current_allele.ref = columns[REF_COLUMN];

            if (parser->flags & VCF_ONLY_SNPS) {
                if (parser->current_allele.ref[1] != '\0') continue;
            }

            for (size_t i = 0; i < parser->sample_num; ++i) {
                char *gt = strtok_r(NULL, "\t", &line_context);
                parser->genotypes[i] = get_genotype_code(gt);
            }
            parser->alt_alleles[0] = strtok_r(columns[ALT_COLUMN], ",",
                                              &parser->allele_context);
            parser->n_alts = 1;
            while ((parser->alt_alleles[parser->n_alts]
                    = strtok_r(NULL, ",", &(parser->allele_context)))
                   != NULL) {
                ++(parser->n_alts);
            }

            qsort(parser->alt_alleles, parser->n_alts,
                  sizeof(char *), alt_comp);

            return fetch_next_allele(parser);
        }
    }
    return parser->current_result = ALLELE_NOT_FETCHED;
}

const char *goto_next_chromosome(VCF_PARSER *parser)
{
    char previous_chromosome[MAX_CHROMOSOME_NAME_LENGTH];
    strcpy(previous_chromosome, parser->current_chromosome);
    while (fetch_next_allele(parser) != ALLELE_NOT_FETCHED) {
        if (strcmp(parser->current_chromosome, previous_chromosome)) {
            // Reached next chromosome
            return parser->current_chromosome;
        }
    }
    return NULL;
}

const char *goto_chromosome(VCF_PARSER *parser, const char *chromosome)
{
    // TODO: Keep track of previously seen chromosomes to make this more
    // sensible and return to their positions instead of re-reading the file.
    if (!strcmp(parser->current_chromosome, chromosome)
        && parser->current_allele_index <= 1) {
        // TODO: may need to rewind to the start
        return parser->current_chromosome;
    }
    while (goto_next_chromosome(parser) != NULL) {
        if (!strcmp(parser->current_chromosome, chromosome)) {
            return parser->current_chromosome;
        }
    }
    // Rewinding and retrying once if the file end was reached
    rewind_vcf_file(parser);
    while (goto_next_chromosome(parser) != NULL) {
        if (!strcmp(parser->current_chromosome, chromosome)) {
            return parser->current_chromosome;
        }
    }
    // Chromosome could not be found
    return NULL;
}

void close_parser(VCF_PARSER *parser)
{
    free(parser->genotypes);
    for (size_t i = 0; i < parser->sample_num; ++i) {
        free(parser->samples[i]);
    }
    free_stringset(parser->chromosome_names);
    free(parser->samples);
    close_vcf_file(parser);
}
