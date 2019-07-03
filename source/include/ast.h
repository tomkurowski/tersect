/*  ast.h

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

#ifndef AST_H
#define AST_H

#include "bitarray.h"
#include "tersect_db.h"

/**
 * AST node types (set theoretical operations / bit arrays / genomes).
 */
#define AST_INTERSECTION 1
#define AST_UNION 2
#define AST_DIFFERENCE 3
#define AST_SYMMETRIC_DIFFERENCE 4
#define AST_GENOME 10

struct ast_node {
    int type;
    struct ast_node *l;
    struct ast_node *r;
    struct genome *genome;
};

/**
 * Create AST subtree describing an intersection/union/symmetric difference of a
 * list of genomes and return its root.
 */
struct ast_node *create_subtree(int operation_type, size_t ngenomes,
                                struct genome *genomes);
struct ast_node *create_ast_node(int operation_type, struct ast_node *l,
                                 struct ast_node *r);
struct ast_node *create_genome_node(struct genome *genome);
struct bitarray *eval_ast(struct ast_node *root, const tersect_db *tdb,
                          const struct tersect_db_interval *ti);
void free_ast(struct ast_node *root);

#endif
