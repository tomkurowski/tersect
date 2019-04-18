%option noyywrap nodefault noinput nounput
%{
    #include "query.tab.h"

    #include "errorc.h"
    #include "tersect_db.h"

    #include <stdlib.h>

    const tersect_db *PARSE_TERSECT_DB;
    struct tersect_db_interval PARSE_TERSECT_REGION;

    static char *strip_single_quotes(char *str);
%}

%%

"union"|"u"     {
                    return UNION;
                }

"intersect"|"i" {
                    return INTER;
                }

"symdiff"|"s"   {
                    return SYMDIFF;
                }

([^-^&|()>,\\ \t\n']+)|('[^']+') {
                    yylval.name = strdup(strip_single_quotes(yytext));
                    return IDENT;
                }

[-^&|()>,\\]      return *yytext;

[ \t\n]         ; /* skip whitespace */

.               {
                    yyerror("Unknown character in query string: %c", *yytext);
                }

%%

/**
 * Used to strip single quotations marks from the name of a GENOME token.
 * NOTE: This modifies the original string.
 */
static char *strip_single_quotes(char *str)
{
    size_t last = strlen(str) - 1;
    if (str[0] == '\'' && str[last] == '\'') {
        str[last] = '\0';
        return str + 1;
    }
    return str;
}
