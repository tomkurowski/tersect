/*  tersect.c

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

#include "tersect.h"

#include "build.h"
#include "errorc.h"
#include "view.h"
#include "chroms.h"
#include "samples.h"
#include "distance.h"
#include "rename.h"
#include "version.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "alleles.h"

static void usage(FILE *stream)
{
    fprintf(stream,
            "\n"
            "Version:  "TERSECT_VERSION"\n"
            "Usage:    tersect <command> [options]\n\n"
            "Commands:\n"
            "    build       build new VCF database\n"
            "    chroms      list chromosomes in the database\n"
            "    dist        calculate distance matrix for samples\n"
            "    help        print this help message\n"
            "    rename      rename sample\n"
            "    samples     list samples in the database\n"
            "    view        display variants belonging to a sample\n"
            "\n");
}

int main(int argc, char *argv[])
{
    error_t rc = SUCCESS;
    char *command = "help";
    if (argc > 1) {
        command = argv[1];
        --argc;
        ++argv;
    }
    if (!strcmp(command, "build")) {
        rc = tersect_build_database(argc, argv);
    } else if (!strcmp(command, "view")) {
        rc = tersect_view_set(argc, argv);
        if (rc != SUCCESS) goto cleanup;
    } else if (!strcmp(command, "chroms")) {
        rc = tersect_print_chromosomes(argc, argv);
    } else if (!strcmp(command, "rename")) {
        rc = tersect_rename_sample(argc, argv);
    } else if (!strcmp(command, "samples")) {
        rc = tersect_print_samples(argc, argv);
    } else if (!strcmp(command, "dist")) {
        rc = tersect_distance(argc, argv);
    } else if (!strcmp(command, "help")) {
        usage(stdout);
    } else {
        usage(stderr);
    }
cleanup:
    if (rc != SUCCESS) {
        report_error(rc);
    }
    return rc;
}
