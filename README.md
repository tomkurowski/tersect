# Tersect

Tersect is a command-line utility for conducting fast set theoretical operations and genetic distance estimation on biological sequence variant data. The tool [generates index files](#building-a-tersect-index) based on provided variant data (VCF files) which can then be used to rapidly [execute flexible set theoretical queries](#set-operations) and output the resulting lists of variants in selected regions.

Tersect is intended to allow for highly responsive, exploratory interaction with variant data as well as for integration with larger tools and pipelines. It follows the Samtools/tabix convention for specifying [genomic regions](#regions) which allows for much faster operations and more manageable output sizes.

Tersect can also be used to provide estimates of genetic distance between sets of samples, using the number of differing sites as a proxy for distance measures.

## Table of Contents

- [Tersect](#tersect)
  - [Table of Contents](#table-of-contents)
  - [Installation](#installation)
    - [Pre-compiled releases](#pre-compiled-releases)
      - [Linux](#linux)
      - [macOS](#macos)
    - [Building Tersect from source](#building-tersect-from-source)
      - [1. Cloning the repository](#1-cloning-the-repository)
      - [2. Building](#2-building)
      - [3. Installing](#3-installing)
  - [Example data](#example-data)
  - [Building a Tersect index](#building-a-tersect-index)
  - [Inspecting a Tersect index](#inspecting-a-tersect-index)
  - [Set operations](#set-operations)
    - [Overview](#overview)
    - [Queries](#queries)
      - [Genomes](#genomes)
      - [Binary operators](#binary-operators)
      - [Genome list](#genome-list)
      - [Functional operators](#functional-operators)
    - [Regions](#regions)

## Installation

### Pre-compiled releases

Tersect packages and binaries are available for download below:

#### Linux

- 64-bit binaries: [tersect-0.12.0-Linux.tar.gz](https://github.com/tomkurowski/tersect/releases/download/v0.12.0/tersect-0.12.0-Linux.tar.gz)
- 64-bit .deb package (Debian, Ubuntu): [tersect-0.12.0-Linux.deb](https://github.com/tomkurowski/tersect/releases/download/v0.12.0/tersect-0.12.0-Linux.deb)
- 64-bit .rpm package (Fedora, openSUSE): [tersect-0.12.0-Linux.rpm](https://github.com/tomkurowski/tersect/releases/download/v0.12.0/tersect-0.12.0-Linux.rpm)

#### macOS

- 64-bit binaries: [tersect-0.12.0-macOS.tar.gz](https://github.com/tomkurowski/tersect/releases/download/v0.12.0/tersect-0.12.0-macOS.tar.gz)

### Building Tersect from source

Building Tersect from source requires CMake version 3.1+ as well as Flex (lexical analyzer) version 2.5+ and Bison (parser generator) version 2.6+.

#### 1. Cloning the repository

```bash
git clone https://github.com/tomkurowski/tersect.git
```

#### 2. Building

For an out-of-source build after cloning the repository use the following commands:

```bash
cd tersect
mkdir build
cd build
cmake ..
make
```

#### 3. Installing

This step may require elevated permissions (e.g. prefacing the command with ``sudo``). The default installation location for Tersect is `/usr/local/bin`.

```bash
make install
```

## Example data

Two archives containing example Tersect index files (.tsi) are available for download below to allow you to try out the utility without needing to create an index file yourself.

The first index contains human genomic variant data for 2504 individuals from the 1000 Genomes Project. While Tersect is capable of handling the entire human genome, the index below is limited to chromosome 1 to make the example archive smaller and quicker to download.

The second index contains tomato genomic variant data for 360 tomato accessions from the AGIS Tomato 360 Resequencing Project and 84 accessions from the Wageningen UR 150 Tomato Genome Resequencing Project for a combined data set of 444 accessions. Samples have been renamed according to a provided key ([accession_names.tsv](https://github.com/tomkurowski/tersect/releases/download/v0.12.0/accession_names.tsv)) to make them more informative and consistent between the two source data sets.

**Note:** the index files provided below are compressed using gzip and need to be uncompressed before use.

- [2504 human genomes, chromosome 1](https://github.com/tomkurowski/tersect/releases/download/v0.12.0/human_chr1.tsi.gz)
- [444 tomato genomes](https://github.com/tomkurowski/tersect/releases/download/v0.12.0/tomato.tsi.gz) [ sample names: [accession_names.tsv](https://github.com/tomkurowski/tersect/releases/download/v0.12.0/accession_names.tsv) ]

## Building a Tersect index

You can build your own Tersect index based on a set of VCF files using the `tersect build` command. You need to provide a name for your index file (a .tsi extension will be added if you omit it) as the first argument, followed by any number of input VCF files (which may be compressed using gzip) to be included in the index. 

Please note that although from a technical point of you, Tersect would work even if your VCF files were called against different reference genomes or versions of the same reference, the biological context of your theoretical operations won't be accurate (depending on how different the reference genomes used). Therefore, we strongly recommend using VCF files called against the same reference version.

**Example:**

The command below builds a Tersect index file named *index.tsi* which includes variants from all *.vcf.gz* files in the *data* directory. Depending on the input size this can take several minutes.

```console
foo@bar:~$ tersect build tomato.tsi ./data/*.vcf.gz
```

Optionally, you can also provide a ``--name-file`` input file containing custom sample names to be used by Tersect. These names will replace the default sample IDs defined in the input VCF header lines. The ``--name-file`` should be a tab-delimited file containing two columns, the first with the sample IDs to be replaced and the second with the names to be used by Tersect. An example is shown below:

```console
TS-1	S.lyc B TS-1
TS-10	S.lyc B TS-10
TS-100	S.lyc B TS-100
TS-101	S.lyc B TS-101
TS-102	S.lyc B TS-102
TS-103	S.lyc B TS-103
TS-104	S.lyc B TS-104
TS-108	S.lyc B TS-108
TS-11	S.lyc B TS-11
TS-110	S.lyc B TS-110
```

You can also modify sample names in an existing Tersect index file by using the `tersect rename` command.

It is worth noting that the descriptive fields of the VCF files are not stored within the Tersect database. The reason for that is once an operation is performed on two of more VCF files, these fields will be discarded anyway as they are genotype-specific. However, you should be able to retrieve it back by intesecting Tersect's output with any VCF files from this list.

## Inspecting a Tersect index

The data contained in a Tersect index file can be inspected using several commands. The `tersect chroms` command prints information on the number of variants present in each of the reference chromosomes as well as the chromosome names and size. **Note:** In the absence of a reference file, the *length* of each chromosome is represented by the position of the last variant, which will always be an underestimate.

**Example:**

The command below prints the per-chromosome variant content of the example Tersect index file named *tomato.tsi* (you can download it [here](#example-data)).

```console
foo@bar:~$ tersect chroms tomato.tsi
Chromosome	Length	Variants
SL2.50ch00	21805702	1343815
SL2.50ch01	98543411	9965680
SL2.50ch02	55340384	5189338
SL2.50ch03	70787603	6741448
SL2.50ch04	66470926	7257520
SL2.50ch05	65875078	6830857
SL2.50ch06	49751619	4870941
SL2.50ch07	68044764	6868152
SL2.50ch08	65866627	6504025
SL2.50ch09	72481975	7102356
SL2.50ch10	65527500	6712293
SL2.50ch11	56302478	5367032
SL2.50ch12	67145147	6719621
```

The `tersect samples` command prints the names of samples present in a Tersect index file. These can be either all samples or a subset based on a naming pattern (the `--match` parameter) and/or on the presence of specific variants (the `--contains` parameter).

Sample name patterns can include wildcard symbols (\*) which match zero or more characters. For example, a pattern like "S.lyc\*" will match all samples whose names begin with "S.lyc". A lone wildcard character matches all samples stored in the Tersect index file.

If you specify a list of variants via the `--contains` paramter, only samples which contain each of those variants will be printed. The variant format should look as follows: `chromosome:position:ref:alt` where `ref` and `alt` are reference and alternate alleles. Multiple variants can be included, separated by commas (e.g. `chr1:1245:A:G,chr8:5300:T:A`).

**Examples:**

The command below prints the names of samples matching the *"S.gal\*"* wildcard pattern contained in the example Tersect index file *tomato.tsi*.

```console
foo@bar:~$ tersect samples tomato.tsi -m "S.gal*"
Sample
S.gal W TS-208
S.gal LA1044
S.gal LA1401
S.gal LA0483
```

The command below prints the names of all samples containing both a T/G SNV at position 100642 on chromosome 3 and an A/G SNV at position 5001015 on chromosome 6 contained in the example Tersect index file *tomato.tsi*.

```console
foo@bar:~$ tersect samples tomato.tsi -c "SL2.50ch03:100642:T:G,SL2.50ch06:5001015:A:G"
Sample
S.lyc LA1479
S.pen LA0716
S.hab LYC4
S.hab LA0407
S.hab LA1777
S.hab LA1718
S.hab CGN15792
S.hab PI134418
S.hab CGN15791
S.chm LA2695
```

## Set operations

### Overview

Tersect can interpret and display the results of set theoretical commands using the `tersect view` command. This is the primary and most flexible functionality of the application and allows the user to construct arbitrarily complex queries. The expected format of a `tersect view` query is as follows:

```console
tersect view [options] index.tsi QUERY [REGION1...]
```

### Queries

A query is a command interpreted and evaluated by Tersect which (if successful) prints either a list of variants (if the result is a single genome or virtual genome) or a list of genome sample names (if the result is a list of genomes). The simplest query consists of a genome name and prints out the variants belonging to that genome. More advanced queries can contain complex combinations of operations described in the sections below.

**Note:** The term *virtual genome* refers to a collection of variants not representing a specific genome - for example, the symmetric difference of two genomes (the collection of variants which appear in one but not both of the genomes). Tersect treats these *virtual genomes* the same way it treats "real" genomes so they can be used as operands in nested queries.

#### Genomes

Genomes can be referred to by their sample name, which is either taken from the header line of the source VCF file or set by the user either manually (see `tersect rename`) or through a tab-delimited name file (see `--name-file` in `tersect build` and `tersect rename`). A sample name can be of any length and can include any characters (including whitespace) except for single quotes ('). However, if a sample name includes whitespace, parentheses, or characters used as Tersect operators (-^&|>,\\), it has to be surrounded by single quotes (').

If the query is (or results in) a single genome or virtual genome, the variants contained by that one genome are printed out.

**Example:**

Print out all the variants belonging to the "S.hab LYC4" genome in the *tomato.tsi* Tersect index file:

```console
foo@bar:~$ tersect view tomato.tsi "'S.hab LYC4'"
##fileformat=VCFv4.3
##tersectVersion=0.11.0
##tersectCommand='S.hab LYC4'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
SL2.50ch00	391	.	C	T	.	.	.
SL2.50ch00	416	.	T	A	.	.	.
SL2.50ch00	734	.	T	G	.	.	.
SL2.50ch00	759	.	C	T	.	.	.
SL2.50ch00	771	.	A	G	.	.	.
SL2.50ch00	778	.	T	A	.	.	.
...
```

**Note:** The sample name had to be surrounded by single quotes because it contains a whitespace character.

#### Binary operators

Tersect supports four basic binary operators. Each operand has to be a **single** genome. All four operators have the same precedence and are left-associative. You can use parentheses to enforce precedence other than simple left-to-right.

| Operator  | Name                 | Usage              | Result                                                       |
|:---------:|:--------------------:|:------------------:|:------------------------------------------------------------:|
| &         | intersection         | GENOME1 & GENOME2  | Virtual genome containing variants found in both GENOME1 and GENOME2 |
| \|        | union                | GENOME1 \| GENOME2 | Virtual genome containing variants found in GENOME1, GENOME2, or both |
| \         | difference           | GENOME1 \ GENOME2  | Virtual genome containing variants found in GENOME1 but not in GENOME2 |
| ^         | symmetric difference | GENOME1 ^ GENOME2  | Virtual genome containing variants found in GENOME1 or GENOME2 but not in both |

The result of a binary operation is treated as a single genome (though it does not have a sample name) called a *virtual genome*, which can be used in further operations.

**Examples:**

Print out the variants shared by 'S.hua LA1983' and 'S.pim LYC2798':

```console
foo@bar:~$ tersect view tomato.tsi "'S.hua LA1983' & 'S.pim LYC2798'"
##fileformat=VCFv4.3
##tersectVersion=0.11.0
##tersectCommand='S.hua LA1983' & 'S.pim LYC2798'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
SL2.50ch00	3235	.	A	G	.	.	.
SL2.50ch00	3277	.	A	G	.	.	.
SL2.50ch00	3873	.	C	T	.	.	.
SL2.50ch00	4083	.	A	G	.	.	.
SL2.50ch00	4112	.	T	G	.	.	.
SL2.50ch00	4314	.	A	C	.	.	.
...
```

Print out the variants which appear in only one of 'S.gal LA1044' or 'S.gal W TS-208':

```console
foo@bar:~$ tersect view tomato.tsi "'S.gal LA1044' ^ 'S.gal W TS-208'"
##fileformat=VCFv4.3
##tersectVersion=0.11.0
##tersectCommand='S.gal LA1044' ^ 'S.gal W TS-208'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
SL2.50ch00	362	.	G	T	.	.	.
SL2.50ch00	867	.	G	T	.	.	.
SL2.50ch00	1198	.	G	A	.	.	.
SL2.50ch00	3235	.	A	G	.	.	.
SL2.50ch00	3567	.	T	G	.	.	.
SL2.50ch00	3873	.	C	T	.	.	.
...
```

Print out the variants which appear in 'S.chi CGN15532' but not 'S.chi CGN15530' or 'S.chi W TS-408':

```console
foo@bar:~$ tersect view tomato.tsi "'S.chi CGN15532' \ 'S.chi CGN15530' \ 'S.chi W TS-408'"
##fileformat=VCFv4.3
##tersectVersion=0.11.0
##tersectCommand='S.chi CGN15532' \ 'S.chi CGN15530' \ 'S.chi W TS-408'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
SL2.50ch00	1163	.	C	G	.	.	.
SL2.50ch00	1811	.	C	G	.	.	.
SL2.50ch00	1818	.	C	A	.	.	.
SL2.50ch00	1818	.	C	G	.	.	.
SL2.50ch00	2432	.	G	A	.	.	.
SL2.50ch00	2544	.	A	T	.	.	.
...
```

**Note:** A more convenient way to conduct the same operation on many genomes is by using [functional operators](#functional-operators).

#### Genome list

Instead of individual genomes, Tersect can also operate on lists of genomes. These can be selected using wildcard patterns matching genome sample names, with the most general pattern of a lone wildcard operator (`*`) matching *all* the genomes in the Tersect index file. Individual genomes can also be appended to lists using commas (`,`) or removed from lists using minus signs (`-`).

Genome lists can also be filtered (using the `>` operator) by whether they contain a specified list of variants. The variant format should look as follows: `chromosome:position:ref:alt` where `ref` and `alt` are reference and alternate alleles. Multiple variants can be included, separated by commas (e.g. `chr1:1245:A:G,chr8:5300:T:A`).

| Operator | Name     | Usage     | Result       |
|:--------:|:--------:|:---------:|:------------:|
| *        | wildcard | PATTERN | Genome list containing all genomes whose sample names match the provided wildcard pattern |
| ,        | append   | GENOMELIST, GENOME | Genome list containing all genomes in GENOMELIST and GENOME |
| -        | remove   | GENOMELIST - GENOME | Genome list containing all genomes in GENOMELIST except GENOME |
| >        | superset | GENOMELIST > VARIANTLIST | Genome list containing all genomes in GENOMELIST which contain all variants in VARIANTLIST

**Note:** Tersect does not distinguish between a genome list which contains only one genome and a single genome. The former can be used in binary operations and the latter can be used in functional operations or in constructing genome lists.

If the query is (or results in) a genome list, the list of their genome sample names is printed out.

**Examples:**

Print out all the names of genomes which begin with "S.pim":

```console
foo@bar:~$ tersect view tomato.tsi "S.pim*"
S.pim P TS-92
S.pim P TS-79
S.pim P TS-77
S.pim P TS-50
S.pim P TS-441
S.pim P TS-440
S.pim P TS-439
S.pim P TS-438
...
```

Print out all the names of genomes which contain an A/G single nucleotide polymorphism at position 828587 in chromosome 7:

```console
foo@bar:~$ tersect view tomato.tsi "* > SL2.50ch07:828587:A:G"
S.lyc C TS-97
S.lyc C TS-94
S.pim P TS-79
S.pim P TS-77
S.lyc B TS-68
S.lyc C TS-53
S.pim P TS-50
S.pim P TS-441
...
```

Print out all the names of genomes which contain a G/A SNP at position 1590608 in chromosome 5 and a T/C SNP at position 5230 in chromosome 12, except for 'S.gal LA1401' and those whose names begin with "S.pim":

```console
foo@bar:~$ tersect view tomato.tsi "* > SL2.50ch05:1590608:G:A,SL2.50ch12:5230:T:C - ('S.pim*','S.gal LA1401')"
S.lyc C TS-431
S.lyc C TS-430
S.lyc LA1314
```

#### Functional operators

Functional operators are used to conduct operations on genome lists instead of individual genomes.

| Operator | Name | Usage | Result |
|:--------:|:----:|:-----:|:------:|
| union() <br> u() | arbitrary union | union(GENOMELIST) <br> u(GENOMELIST) | Virtual genome containing all variants contained in any of the genomes in GENOMELIST |
| intersect() <br> i() | arbitrary intersection | intersect(GENOMELIST) <br> i(GENOMELIST) | Virtual genome containing all variants which appear in every genome in GENOMELIST |

The result of a functional operation is treated as a single genome (though it does not have a sample name).

**Examples:**

Union of all genomes, containing every variant recorded in the *tomato.tsi* Tersect index file:

```console
foo@bar:~$ tersect view tomato.tsi "u(*)"
##fileformat=VCFv4.3
##tersectVersion=0.11.0
##tersectCommand=u(*)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
SL2.50ch00	280	.	A	C	.	.	.
SL2.50ch00	284	.	A	G	.	.	.
SL2.50ch00	316	.	C	T	.	.	.
SL2.50ch00	323	.	C	T	.	.	.
SL2.50ch00	332	.	A	T	.	.	.
SL2.50ch00	362	.	G	T	.	.	.
...
```

Intersection of all genomes which contain a T/A single nucleotide polymorphism at position 12547 in chromosome 12, containing all variants that are shared by each of those genomes:

```console
foo@bar:~$ tersect view tomato.tsi "i(* > SL2.50ch12:12547:T:A)"
##fileformat=VCFv4.3
##tersectVersion=0.11.0
##tersectCommand=i(* > SL2.50ch12:12547:T:A)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
SL2.50ch00	16576	.	T	C	.	.	.
SL2.50ch00	26171	.	G	T	.	.	.
SL2.50ch00	29880	.	A	G	.	.	.
SL2.50ch00	37486	.	T	G	.	.	.
SL2.50ch00	40476	.	G	T	.	.	.
SL2.50ch00	436850	.	A	G	.	.	.
...
```

Print all the variants which appear only in genome S.hab CGN15792. This is achieved by finding the difference of that genome and the union of all genomes except S.hab CGN15792:

```console
foo@bar:~$ tersect view tomato.tsi "'S.hab CGN15792' \ u(* - 'S.hab CGN15792')"
##fileformat=VCFv4.3
##tersectVersion=0.11.0
##tersectCommand='S.hab CGN15792' \ u(* - 'S.hab CGN15792')
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
SL2.50ch00	1163	.	C	T	.	.	.
SL2.50ch00	1596	.	G	A	.	.	.
SL2.50ch00	2048	.	G	A	.	.	.
SL2.50ch00	2933	.	G	A	.	.	.
SL2.50ch00	2987	.	A	T	.	.	.
SL2.50ch00	4349	.	C	T	.	.	.
...
```

### Regions

By default, queries are executed and results are returned for the entire genome. However, it is possible to selectively execute a query only on a specified region. The familiar tabix/samtools format `chromosome:beginPos-endPos` is used to specify those regions. The coordinates are one-based and inclusive.

Limiting queries to regions allows for much faster execution since far fewer positions need to be processed and printed, capturing only intervals of interest. This feature makes it possible to use Tersect's flexible queries as a high-performance part of a larger pipeline or the back-end of a highly responsive, interactive application.

**Example:**

Print a union, that is all the variants appearing either in genome 'S.lyc SG16', 'S.lyc LA1421', or both, from the first 90 kbp of chromosome 2 in the *tomato.tsi* index file:

```console
foo@bar:~$ tersect view tomato.tsi "'S.lyc SG16' | 'S.lyc LA1421'" SL2.50ch02:1-90000
##fileformat=VCFv4.3
##tersectVersion=0.11.0
##tersectCommand='S.lyc SG16' | 'S.lyc LA1421'
##tersectRegion=SL2.50ch02:1-90000
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
SL2.50ch02      204     .       A       G       .       .       .
SL2.50ch02      255     .       TCC     TCCC    .       .       .
SL2.50ch02      255     .       TCC     TCCCC   .       .       .
SL2.50ch02      2382    .       G       A       .       .       .
SL2.50ch02      13383   .       G       A       .       .       .
SL2.50ch02      21752   .       C       T       .       .       .
SL2.50ch02      24538   .       T       C       .       .       .
SL2.50ch02      29276   .       G       T       .       .       .
SL2.50ch02      71245   .       A       C       .       .       .
SL2.50ch02      73326   .       C       T       .       .       .
SL2.50ch02      86236   .       C       A       .       .       .
SL2.50ch02      86601   .       A       G       .       .       .
SL2.50ch02      86635   .       T       A       .       .       .
SL2.50ch02      86695   .       T       C       .       .       .
SL2.50ch02      86769   .       G       A       .       .       .
SL2.50ch02      87079   .       T       A       .       .       .
```
