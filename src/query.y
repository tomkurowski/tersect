%{
    #include "query.h"

    #include "ast.h"

    #include <stdarg.h>
    #include <stdlib.h>
    #include <string.h>

    struct id_list {
        char **ids;
        size_t count;
    };

    struct gen_list {
        struct genome *genomes;
        size_t count;
    };

    int yylex(void);
    void yy_scan_string(const char *);
    int yylex_destroy(void);
    struct id_list *merge_id_lists(struct id_list *list_a,
                                   struct id_list *list_b);
    void load_genomes(const struct id_list *gen_id_list,
                  const struct id_list *var_id_list,
                  struct gen_list *gen_list);
    struct gen_list *diff_gen_lists(struct gen_list *list_a,
                                    struct gen_list *list_b);
    void free_id_list(struct id_list *id_list);
    void free_gen_list(struct gen_list *gen_list);
    const tersect_db *PARSE_TERSECT_DB;
    struct ast_node *PARSE_OUTPUT;
%}

%code provides {
    void yyerror(const char *, ...);
}

%union {
    char *name;
    struct ast_node *ast;
    struct id_list *id_list;
    struct gen_list *gen_list;
}

%token <name> IDENT
%token UNION
%token INTER
%token SYMDIFF
%type <ast> expr
%type <id_list> list
%type <gen_list> genlist;
%left '&' '|' '-' '^' '\\'
%right '>'
%left ','

// Needed to handle parantheses for list / genlist / expr
%expect 2

%%

program:
        program expr            {
                                    PARSE_OUTPUT = $2;
                                }
        | %empty
        ;

list:
        IDENT                   {
                                    struct id_list *list = malloc(sizeof *list);
                                    list->ids = malloc(sizeof *list->ids);
                                    list->ids[0] = $1;
                                    list->count = 1;
                                    $$ = list;
                                }
        | list ',' list         {
                                    $$ = merge_id_lists($1, $3);
                                }
        | '(' list ')'          {
                                    $$ = $2;
                                }
        ;

genlist:
        list                    {
                                    struct id_list *gen_id_list = $1;
                                    struct gen_list *gen_list = malloc(sizeof *gen_list);
                                    load_genomes(gen_id_list, NULL, gen_list);
                                    free_id_list(gen_id_list);
                                    $$ = gen_list;
                                }
        | list '>' list         {
                                    struct id_list *gen_id_list = $1;
                                    struct id_list *var_id_list = $3;
                                    struct gen_list *gen_list = malloc(sizeof *gen_list);
                                    load_genomes(gen_id_list,
                                                 var_id_list,
                                                 gen_list);
                                    free_id_list(gen_id_list);
                                    free_id_list(var_id_list);
                                    $$ = gen_list;
                                }
        | genlist '-' genlist   {
                                    $$ = diff_gen_lists($1, $3);
                                }
        | '(' genlist ')'       {
                                    $$ = $2;
                                }
        ;

expr:
        genlist                 {
                                    struct gen_list *gen_list = $1;
                                    if (gen_list->count > 1) {
                                        // yyerror("Genome lists are only allowed within functions");
                                        for (size_t i = 0; i < gen_list->count; ++i) {
                                            printf("%s\n", gen_list->genomes[i].name);
                                        }
                                        $$ = NULL;
                                    } else if (gen_list->count == 0) {
                                        yyerror("Empty genome list");
                                        $$ = NULL;
                                    } else {
                                        $$ = create_genome_node(&gen_list->genomes[0]);
                                    }
                                    free_gen_list(gen_list);
                                }
        | expr '&' expr         {
                                    if ($1 == NULL || $3 == NULL) {
                                        yyerror("Invalid operand in intersection");
                                    }
                                    $$ = create_ast_node(AST_INTERSECTION,
                                                         $1, $3);
                                }
        | expr '|' expr         {
                                    if ($1 == NULL || $3 == NULL) {
                                        yyerror("Invalid operand in union");
                                    }
                                    $$ = create_ast_node(AST_UNION, $1, $3);
                                }
        | expr '\\' expr        {
                                    if ($1 == NULL || $3 == NULL) {
                                        yyerror("Invalid operand in difference");
                                    }
                                    $$ = create_ast_node(AST_DIFFERENCE, $1, $3);
                                }
        | expr '^' expr         {
                                    if ($1 == NULL || $3 == NULL) {
                                        yyerror("Invalid operand in symmetric difference");
                                    }
                                    $$ = create_ast_node(AST_SYMMETRIC_DIFFERENCE,
                                                         $1, $3);
                                }
        | '(' expr ')'          {
                                    $$ = $2;
                                }
        | UNION '(' genlist ')' {
                                    struct gen_list *gen_list = $3;
                                    $$ = create_subtree(AST_UNION,
                                                        gen_list->count,
                                                        gen_list->genomes);
                                    free_gen_list(gen_list);
                                }
        | INTER '(' genlist ')' {
                                    struct gen_list *gen_list = $3;
                                    $$ = create_subtree(AST_INTERSECTION,
                                                        gen_list->count,
                                                        gen_list->genomes);
                                    free_gen_list(gen_list);
                                }
        | SYMDIFF '(' genlist ')'  {
                                    struct gen_list *gen_list = $3;
                                    $$ = create_subtree(AST_SYMMETRIC_DIFFERENCE,
                                                        gen_list->count,
                                                        gen_list->genomes);
                                    free_gen_list(gen_list);
                                }
        ;

%%

void yyerror(const char *s, ...)
{
    va_list ap;
    va_start(ap, s);
    fprintf(stderr, "Error: ");
    vfprintf(stderr, s, ap);
    fprintf(stderr, "\n");
    exit(1);
}

/**
 * Naive implementation,
 */
struct gen_list *diff_gen_lists(struct gen_list *list_a,
                                struct gen_list *list_b)
{
    for (size_t i = 0; i < list_a->count; ++i) {
        for (size_t j = 0; j < list_b->count; ++j) {
            if (list_a->genomes[i].hdr == list_b->genomes[j].hdr) {
                --list_a->count;
                j = 0;
                if (i < list_a->count) {
                    memmove(&list_a->genomes[i], &list_a->genomes[i + 1],
                            (list_a->count - i) * sizeof *list_a->genomes);
                }
                --i;
            }
        }
    }
    list_a->genomes = realloc(list_a->genomes,
                              list_a->count * sizeof *list_a->genomes);
    free(list_b->genomes);
    free(list_b);
    return list_a;
}

/**
 * Merges two id lists. Saves the result in list_a and frees list_b.
 * Returns list_a.
 */
struct id_list *merge_id_lists(struct id_list *list_a, struct id_list *list_b)
{
    size_t old_count = list_a->count;
    list_a->count += list_b->count;
    list_a->ids = realloc(list_a->ids, list_a->count * sizeof *list_a->ids);
    memcpy(&list_a->ids[old_count], list_b->ids,
           list_b->count * sizeof *list_b->ids);
    free(list_b->ids);
    free(list_b);
    return list_a;
}

void load_genomes(const struct id_list *gen_id_list,
                  const struct id_list *var_id_list,
                  struct gen_list *gen_list)
{
    error_t rc;
    if (var_id_list == NULL) {
        rc = tersect_db_get_genomes(PARSE_TERSECT_DB,
                                    gen_id_list->count, gen_id_list->ids,
                                    0, NULL,
                                    &gen_list->count, &gen_list->genomes);
    } else {
        rc = tersect_db_get_genomes(PARSE_TERSECT_DB,
                                    gen_id_list->count, gen_id_list->ids,
                                    var_id_list->count, var_id_list->ids,
                                    &gen_list->count, &gen_list->genomes);
    }
    if (rc != SUCCESS || gen_list->count == 0) {
        /* TODO: may need to make "tersect_db_get_genomes" verbose here to
        print specific missing name(s) */
        yyerror("Could not find specified genome(s)");
    }
}

void free_id_list(struct id_list *id_list)
{
    for (size_t i = 0; i < id_list->count; ++i) {
        free(id_list->ids[i]);
    }
    free(id_list->ids);
    free(id_list);
}

void free_gen_list(struct gen_list *gen_list)
{
    free(gen_list->genomes);
    free(gen_list);
}

struct ast_node *run_set_parser(const char *query, const tersect_db *tdb)
{
    PARSE_TERSECT_DB = tdb;
    yy_scan_string(query);
    yyparse();
    yylex_destroy();
    return PARSE_OUTPUT;
}
