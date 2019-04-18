/*  vcf_writer.c

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

#include "vcf_writer.h"

#include "tersect_db_internal.h"
#include "version.h"

#include <stdlib.h>

/**
 * Format strings for printing out variants as VCF lines.
 * The positions in the array correspond to codes defined in snv.h.
 *
 * These are the 8 fixed fields required by the VCF format:
 * CHROM    (string)    chromosome name
 * POS      (integer)   reference position
 * ID       (string)    identifier, not used by tersect (".")
 * REF      (string)    reference base, only one character for SNVs
 * ALT      (string)    alternate base, only one character for SNVs
 * QUAL     (float)     quality, not used by tersect (".")
 * FILTER   (string)    filter status, not used by tersect (".")
 * INFO     (string)    additional information, not used by tersect (".")
 */
const char * const variant_format[] = {
    "%s\t%u\t.\t%s\t.\t.\t.\n",     // 0    non-SNV variant
    "%s\t%u\t.\tA\tC\t.\t.\t.\n",   // 1    SNV_A_C
    "%s\t%u\t.\tA\tG\t.\t.\t.\n",   // 2    SNV_A_G
    "%s\t%u\t.\tA\tT\t.\t.\t.\n",   // 3    SNV_A_T
    "%s\t%u\t.\tC\tA\t.\t.\t.\n",   // 4    SNV_C_A
    "%s\t%u\t.\tC\tG\t.\t.\t.\n",   // 5    SNV_C_G
    "%s\t%u\t.\tC\tT\t.\t.\t.\n",   // 6    SNV_C_T
    "%s\t%u\t.\tG\tA\t.\t.\t.\n",   // 7    SNV_G_A
    "%s\t%u\t.\tG\tC\t.\t.\t.\n",   // 8    SNV_G_C
    "%s\t%u\t.\tG\tT\t.\t.\t.\n",   // 9    SNV_G_T
    "%s\t%u\t.\tT\tA\t.\t.\t.\n",   // 10   SNV_T_A
    "%s\t%u\t.\tT\tC\t.\t.\t.\n",   // 11   SNV_T_C
    "%s\t%u\t.\tT\tG\t.\t.\t.\n",   // 12   SNV_T_G
};

void vcf_print_header(const char *command, size_t nregions,
                      char **region_strings)
{
    printf("##fileformat="VCF_FORMAT"\n");
    printf("##tersectVersion="TERSECT_VERSION"\n");
    printf("##tersectCommand=%s\n", command);
    if (region_strings != NULL) {
        printf("##tersectRegion=");
        for (size_t i = 0; i < nregions; ++i) {
            printf("%s", region_strings[i]);
            if (i + 1 < nregions) printf(" ");
        }
        printf("\n");
    }
    printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
}

static inline void print_snv(const tersect_db *tdb, const struct variant v,
                             const char *chr_name)
{
    if (v.type) {
        // SNV
        printf(variant_format[v.type], chr_name, v.position);
    } else {
        // InDel
        printf(variant_format[0], chr_name, v.position,
               (char *)(tdb->mapping + v.allele));
    }
}

void vcf_print_bitarray(const tersect_db *tdb, const struct bitarray *ba,
                        const struct tersect_db_interval *ti)
{
    size_t allele_num;
    uint64_t *allele_indices;
    bitarray_get_set_indices(ba, &allele_num, &allele_indices);
    for (size_t i = 0; i < allele_num; ++i) {
        print_snv(tdb, ti->variants[allele_indices[i]], ti->chromosome.name);
    }
    free(allele_indices); // Allocated by bitarray_get_set_indices
}
