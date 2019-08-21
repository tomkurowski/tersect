/*  errorc.c

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

#include "errorc.h"

#include <stdio.h>

struct error_desc error_desc[] = {
    { SUCCESS, "No error" },
    { FAILURE, "General failure" },
    { E_ALLOC, "Memory allocation error"},
    { E_NO_GENOME, "Sample not found"},
    { E_NO_TSI_FILE, "No Tersect index (.tsi) file specified"},
    { E_TSI_NOPEN, "Could not open specified Tersect index (.tsi) file"},
    { E_BUILD_NO_OUTNAME, "Output filename missing"},
    { E_BUILD_NO_FILES, "No input files specified"},
    { E_BUILD_CREATE, "Tersect database file could not be created"},
    { E_BUILD_DB_EXISTS, "Output file already exists (use -f to overwrite)"},
    { E_BUILD_NO_WRITE, "No write permissions on specified output file"},
    { E_BUILD_DUPSAMPLE, "Duplicate sample in input data"},
    { E_PARSE_REGION, "Region could not be parsed"},
    { E_PARSE_REGION_NO_CHROMOSOME, "Requested chromosome is not in the database"},
    { E_PARSE_REGION_BAD_BOUNDS, "Incorrect region bounds specified"},
    { E_PARSE_ALLELE, "Allele could not be parsed"},
    { E_PARSE_ALLELE_NO_CHROMOSOME, "Requested chromosome is not in the database"},
    { E_PARSE_ALLELE_BAD_POSITION, "Incorrect position specified"},
    { E_VCF_PARSE_FILE, "Failed to parse VCF/VCF.GZ file"},
    { E_PARSE_ALLELE_UNKNOWN, "Allele not in database"},
    { E_VIEW_NO_QUERY, "No set query specified"},
    { E_RENAME_NOPEN, "Coult not open specified name file"},
    { E_RENAME_PARSE, "Name file could not be parsed"},
    { E_DIST_BIN_REGIONS, "Only one region allowed if binning is enabled"},
    { E_DIST_LIST_NOPEN, "Match string list file could not be opened"}
};

void report_error(error_t code) {
    for (size_t i = 0; i < sizeof error_desc / sizeof error_desc[0]; ++i) {
        if (error_desc[i].code == code) {
            fprintf(stderr, "Error %d: %s\n", code, error_desc[i].description);
            return;
        }
    }
    fprintf(stderr, "Error %d: Unknown error code\n", code);
}
