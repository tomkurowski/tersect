/*  samples.c

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

#include "samples.h"

#include "tersect_db.h"

#include <getopt.h>
#include <stdint.h>
#include <stdlib.h>

/* Local flags for samples */
#define NO_HEADERS      2
static int local_flags = 0;

static void usage(FILE *stream)
{
    fprintf(stream,
            "\n"
            "Usage:    tersect samples [options] <db.tsi>\n\n"
            "Options:\n"
            "    -c, --contains STR      print only samples containing each variant from an\n"
            "                            input list (e.g. \"ch02:100:A:G,ch05:4031:C:T\")\n"
            "    -h, --help              print this help message\n"
            "    -m, --match STR         print only samples matching a wildcard pattern\n"
            "                            (e.g. \"S.chi*\" to match all samples beginning\n"
            "                             with \"S.chi\")\n"
            "    -n, --no-headers        skip column headers\n"
            "\n");
}

error_t tersect_print_samples(int argc, char **argv)
{
    error_t rc = SUCCESS;
    char *db_filename = NULL;
    char *contains = NULL;
    char *pattern = NULL;
    static struct option loptions[] = {
        {"contains", required_argument, NULL, 'c'},
        {"help", no_argument, NULL, 'h'},
        {"match", required_argument, NULL, 'm'},
        {"no-headers", no_argument, NULL, 'n'},
        {NULL, 0, NULL, 0}
    };
    int c;
    while ((c = getopt_long(argc, argv, ":c:hm:n", loptions, NULL)) != -1) {
        switch(c) {
        case 'c':
            contains = optarg;
            break;
        case 'h':
            usage(stdout);
            return SUCCESS;
        case 'm':
            pattern = optarg;
            break;
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
    struct genome *samples = NULL;
    if (pattern == NULL && contains == NULL) {
        rc = tersect_db_get_genomes(tdb, 0, NULL, 0, NULL,
                                    &count, &samples);
    } else if (pattern == NULL) {
        rc = tersect_db_get_genomes(tdb, 0, NULL, 1, &contains,
                                    &count, &samples);
    } else if (contains == NULL) {
        rc = tersect_db_get_genomes(tdb, 1, &pattern, 0, NULL,
                                    &count, &samples);
    } else {
        rc = tersect_db_get_genomes(tdb, 1, &pattern, 1, &contains,
                                    &count, &samples);
    }
    if (rc != SUCCESS) goto cleanup;

    if (!(local_flags & NO_HEADERS)) {
        printf("Sample\n");
    }
    for (uint32_t i = 0; i < count; ++i) {
        printf("%s\n", samples[i].name);
    }
    if (count) {
        free(samples);
    }
cleanup:
    tersect_db_close(tdb);
    return rc;
}
