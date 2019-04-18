/*  chroms.c

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

#include "chroms.h"

#include "tersect_db.h"

#include <getopt.h>
#include <stdint.h>
#include <stdlib.h>

/* Local flags for chroms */
#define NO_HEADERS      2
static int local_flags = 0;

static void usage(FILE *stream)
{
    fprintf(stream,
            "\n"
            "Usage:    tersect chroms [options] <db.tsi>\n\n"
            "Options:\n"
            "    -h, --help              print this help message\n"
            "    -n, --no-headers        skip column headers\n"
            "\n");
}

error_t tersect_print_chromosomes(int argc, char **argv)
{
    char *db_filename = NULL;
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
    } else if (argc > 1) {
        // Too many arguments
        usage(stderr);
        return SUCCESS;
    }
    db_filename = argv[0];
    tersect_db *tdb = tersect_db_open(db_filename);
    if (tdb == NULL) return E_TSI_NOPEN;
    size_t count;
    struct chromosome *chroms;
    tersect_db_get_chromosomes(tdb, &count, &chroms);
    if (!(local_flags & NO_HEADERS)) {
        printf("Chromosome\tLength\tVariants\n");
    }
    for (size_t i = 0; i < count; ++i) {
        printf("%s\t%u\t%u\n", chroms[i].name,
               chroms[i].length, chroms[i].variant_count);
    }
    tersect_db_close(tdb);
    free(chroms);
    return SUCCESS;
}
