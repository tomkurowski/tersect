/*  view.c

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

#include "view.h"

#include "ast.h"
#include "query.h"
#include "tersect_db.h"
#include "vcf_writer.h"

#include <getopt.h>
#include <stdlib.h>

/* Local flags for view */
#define NO_HEADERS      2
static int local_flags = 0;

static void usage(FILE *stream)
{
    fprintf(stream,
            "\n"
            "Usage:    tersect view [options] <db.tsi> <query> [region]...\n\n"
            "Options:\n"
            "    -h, --help              print this help message\n"
            "    -n, --no-header         skip VCF header\n"
            "\n");
}

error_t tersect_view_set(int argc, char **argv)
{
    error_t rc = SUCCESS;
    char *db_filename = NULL;
    char *query = NULL;
    char **region_strings = NULL;
    size_t nregions = 0;
    static struct option loptions[] = {
        {"help", no_argument, NULL, 'h'},
        {"no-headers", no_argument, NULL, 'n'},
        {NULL, 0, NULL, 0}
    };
    int c;
    while ((c = getopt_long(argc, argv, ":hn", loptions, NULL)) != -1) {
        switch(c) {
        case 'h':
            usage(stdout);
            return SUCCESS;
        case 'n':
            local_flags |= NO_HEADERS;
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
    if (argc == 1) {
        // Missing set query
        return E_VIEW_NO_QUERY;
    }
    query = argv[1];
    argc -= 2;
    argv += 2;
    if (argc) {
        region_strings = argv;
        nregions = argc;
    }
    tersect_db *tdb = tersect_db_open(db_filename);
    if (tdb == NULL) return E_TSI_NOPEN;
    struct genomic_interval *regions;
    if (nregions) {
        rc = tersect_db_parse_regions(tdb, nregions, region_strings, &regions);
    } else {
        rc = tersect_db_get_regions(tdb, &nregions, &regions);
    }
    if (rc != SUCCESS) goto cleanup_1;
    struct ast_node *command = run_set_parser(query, tdb);
    if (command == NULL) goto cleanup_2;
    if (!(local_flags & NO_HEADERS)) {
        vcf_print_header(query, nregions, region_strings);
    }
    for (size_t i = 0; i < nregions; ++i) {
        struct tersect_db_interval ti;
        tersect_db_get_interval(tdb, &regions[i], &ti);
        struct bitarray *result = eval_ast(command, tdb, &ti);
        if (result == NULL) goto cleanup_3;
        vcf_print_bitarray(tdb, result, &ti);
        free_bitarray(result);
    }
cleanup_3:
    free_ast(command);
cleanup_2:
    free(regions);
cleanup_1:
    tersect_db_close(tdb);
    return rc;
}
