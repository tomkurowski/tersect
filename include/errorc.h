/*  errorc.h

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

#ifndef TERSECT_ERROR_H
#define TERSECT_ERROR_H

typedef enum error_t {
    SUCCESS = 0,
    FAILURE = 1,
    E_ALLOC = 500,
    E_NO_GENOME = 600,
    E_NO_TSI_FILE = 700,
    E_TSI_NOPEN = 701,
    E_BUILD_NO_OUTNAME = 5000,
    E_BUILD_NO_FILES = 5001,
    E_BUILD_CREATE = 5002,
    E_BUILD_DB_EXISTS = 5003,
    E_BUILD_NO_WRITE = 5004,
    E_BUILD_DUPSAMPLE = 5005,
    E_PARSE_REGION = 6000,
    E_PARSE_REGION_NO_CHROMOSOME = 6001,
    E_PARSE_REGION_BAD_BOUNDS = 6002,
    E_PARSE_ALLELE = 7000,
    E_PARSE_ALLELE_NO_CHROMOSOME = 7001,
    E_PARSE_ALLELE_BAD_POSITION = 7002,
    E_PARSE_ALLELE_UNKNOWN = 7003,
    E_VCF_PARSE_FILE = 7100,
    E_VIEW_NO_QUERY = 8000,
    E_RENAME_NOPEN = 9000,
    E_RENAME_PARSE = 9001,
    E_DIST_BIN_REGIONS = 10000,
    E_DIST_LIST_NOPEN = 10001
} error_t;

extern struct error_desc {
    error_t code;
    char *description;
} error_desc[];

void report_error(error_t code);

#endif
