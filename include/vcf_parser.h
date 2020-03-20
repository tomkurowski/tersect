/*  vcf_parser.h

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

#ifndef VCF_PARSER_H
#define VCF_PARSER_H

#include "alleles.h"

#include <stdio.h>

#define VCF_PARSER_INIT_SUCCESS     0
#define VCF_PARSER_INIT_FAILURE     1

#define VCF_PARSER_CHROM_NOT_FOUND  0
#define VCF_PARSER_CHROM_FOUND      1

#define ALLELE_FETCHED              1
#define ALLELE_NOT_FETCHED          0

#define MAX_ALT_ALLELES             100
#define MAX_CHROMOSOME_NAME_LENGTH  100
#define MAX_FILENAME_LENGTH         500
#define MAX_SAMPLE_NAME_LENGTH      250

/* Genotype codes */
#define GENOTYPE_HOM_REF            0
#define GENOTYPE_HOM_ALT            1
#define GENOTYPE_HET                2

/* Parser flags */
#define VCF_ONLY_HOMOZYGOUS         2
#define VCF_ONLY_SNPS               4
#define VCF_ONLY_INDELS             8

typedef struct ParserHandle_t {
    char filename[MAX_FILENAME_LENGTH];
    int *genotypes;
    char **samples;
    size_t sample_num;
    int flags;
    FILE *file_handle;
    char current_chromosome[MAX_CHROMOSOME_NAME_LENGTH];
    char *alt_alleles[MAX_ALT_ALLELES];
    int n_alts;
    char *line_buffer;
    size_t buffer_size;
    struct allele current_allele;
    char *allele_context;
    int current_result;
} VCF_PARSER;

int init_parser(const char *filename, int flags, VCF_PARSER *parser);
int fetch_next_allele(VCF_PARSER *parser);
const char *goto_next_chromosome(VCF_PARSER *parser);
const char *goto_chromosome(VCF_PARSER *parser, const char *chromosome);
int parser_allele_cmp(const void *a, const void *b);
void close_parser(VCF_PARSER *parser);

#endif
