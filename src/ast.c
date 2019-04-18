/*  ast.c

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

#include "ast.h"

#include <stdbool.h>
#include <stdlib.h>

static struct bitarray *eval_node(struct ast_node *node,
                                  const tersect_db *tdb,
                                  const struct tersect_db_interval *ti);

/**
 * Allocate and initialise abstract syntax tree node for a binary operation.
 */
struct ast_node *create_ast_node(int operation_type, struct ast_node *l,
                                 struct ast_node *r)
{
    struct ast_node *node = malloc(sizeof *node);
    node->type = operation_type;
    node->l = l;
    node->r = r;
    return node;
}

struct ast_node *create_genome_node(struct genome *genome)
{
    struct ast_node *node = malloc(sizeof *node);
    node->type = AST_GENOME;
    node->genome = malloc(sizeof *node->genome);
    *(node->genome) = *genome;
    return node;
}

struct ast_node *create_subtree(int operation_type, size_t ngenomes,
                                struct genome *genomes)
{
    struct ast_node *root = create_genome_node(&genomes[0]);
    for (size_t i = 1; i < ngenomes; ++i) {
        root = create_ast_node(operation_type, root,
                               create_genome_node(&genomes[i]));
    }
    return root;
}

/**
 * Allocates memory for the container struct.
 */
static struct bitarray *load_bitarray(const tersect_db *tdb, struct genome *genome,
                                      const struct tersect_db_interval *ti)
{
    struct bitarray ba;
    tersect_db_get_bitarray(tdb, genome, &ti->chromosome, &ba);
    struct bitarray *region_ba = malloc(sizeof *region_ba);
    bitarray_extract_region(region_ba, &ba, &ti->interval);
    return region_ba;
}

static inline struct bitarray *ast_node_operation(struct ast_node *node,
                                                  void (*op)(const struct bitarray*,
                                                             const struct bitarray*,
                                                             struct bitarray**),
                                                  const tersect_db *tdb,
                                                  const struct tersect_db_interval *ti)
{
    struct bitarray *ba = eval_node(node->l, tdb, ti);
    struct bitarray *bb = eval_node(node->r, tdb, ti);
    struct bitarray *out;
    op(ba, bb, &out);
    if (node->l->type == AST_GENOME) {
        free(ba);
    } else {
        free_bitarray(ba);
    }
    if (node->r->type == AST_GENOME) {
        free(bb);
    } else {
        free_bitarray(bb);
    }
    return out;
}

static struct bitarray *eval_node(struct ast_node *node, const tersect_db *tdb,
                                  const struct tersect_db_interval *ti)
{
    switch (node->type) {
    case AST_INTERSECTION:
        return ast_node_operation(node, &bitarray_intersection, tdb, ti);
    case AST_UNION:
        return ast_node_operation(node, &bitarray_union, tdb, ti);
    case AST_DIFFERENCE:
        return ast_node_operation(node, &bitarray_difference, tdb, ti);
    case AST_SYMMETRIC_DIFFERENCE:
        return ast_node_operation(node, &bitarray_symmetric_difference, tdb, ti);
    case AST_GENOME:
        return load_bitarray(tdb, node->genome, ti);
    }
    return NULL;
}

struct bitarray *eval_ast(struct ast_node *root, const tersect_db *tdb,
                   const struct tersect_db_interval *ti)
{
    if (root->type == AST_GENOME) {
        struct bitarray *ba = load_bitarray(tdb, root->genome, ti);
        struct bitarray *out = copy_bitarray(ba);
        free(ba);
        return out;
    } else {
        return eval_node(root, tdb, ti);
    }
}

/**
 * Free the entire abstract syntax tree starting from the root.
 */
void free_ast(struct ast_node *root)
{
    if (root->type != AST_GENOME) {
        free_ast(root->l);
        free_ast(root->r);
    } else {
        free(root->genome);
    }
    free(root);
}
