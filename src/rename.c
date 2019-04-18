/*  rename.c

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

#include "rename.h"

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void usage(FILE *stream)
{
    fprintf(stream,
            "\n"
            "Usage:    tersect rename [options] <db.tsi> <oldname> <newname>\n"
            "          tersect rename [options] <db.tsi> -n <names.tsv>\n\n"
            "Options:\n"
            "    -h, --help              print this help message\n"
            "    -n, --name-file         tsv file containing sample names\n"
            "\n");
}

error_t tersect_rename_sample(int argc, char **argv)
{
    error_t rc = SUCCESS;
    char *db_filename = NULL;
    char *oldname = NULL;
    char *newname = NULL;
    char *name_filename = NULL;
    static struct option loptions[] = {
        {"help", no_argument, NULL, 'h'},
        {"name-file", required_argument, NULL, 'n'},
        {NULL, 0, NULL, 0}
    };
    int c;
    while ((c = getopt_long(argc, argv, ":hn:", loptions, NULL)) != -1) {
        switch(c) {
        case 'h':
            usage(stdout);
            return SUCCESS;
        case 'n':
            name_filename = optarg;
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
    tersect_db *tdb = tersect_db_open(db_filename);
    if (tdb == NULL) return E_TSI_NOPEN;

    if (name_filename != NULL) {
        rc = tersect_load_name_file(tdb, name_filename);
    } else {
        if (argc != 3) {
            // Missing old name/new name/too many arguments
            usage(stderr);
            return SUCCESS;
        }
        oldname = argv[1];
        newname = argv[2];
        rc = tersect_db_rename_genome(tdb, oldname, newname);
    }

    tersect_db_close(tdb);
    return rc;
}

error_t tersect_load_name_file(tersect_db *tdb, const char *filename)
{
    error_t rc = SUCCESS;
    FILE *fh = fopen(filename, "r");
    if (fh == NULL) return E_RENAME_NOPEN;
    char *line_buffer = NULL;
    size_t buffer_size = 0;
    char *line_context = NULL;
    while (getline(&line_buffer, &buffer_size, fh) != -1) {
        char *old_name = strtok_r(line_buffer, "\t\n", &line_context);
        if (old_name == NULL) {
            rc = E_RENAME_PARSE;
            goto cleanup;
        }
        char *new_name = strtok_r(NULL, "\t\n", &line_context);
        if (new_name == NULL) {
            rc = E_RENAME_PARSE;
            goto cleanup;
        }
        rc = tersect_db_rename_genome(tdb, old_name, new_name);
        if (rc != SUCCESS && rc != E_NO_GENOME) {
            goto cleanup;
        } else {
            rc = SUCCESS; // considering E_NO_GENOME a success as well
        }
    }
cleanup:
    if (line_buffer != NULL) {
        free(line_buffer);
    }
    fclose(fh);
    return rc;
}
