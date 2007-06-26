/*!
\file
\brief A top down parser for function evaluation.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

*/

#ifndef PSS_PARSER_H
#define PSS_PARSER_H

#include "../headers/standardtypes.h"

struct _ST_NODE;

struct _ST_NODE* pss_parse(CHAR* funct);
void pss_parse_cleanup(struct _ST_NODE* head);

/* global interpreter call for space functions */
DOUBLE pss_evaluate_funct(struct _ST_NODE* funct, DOUBLE x, DOUBLE y, DOUBLE z);

/* global interpreter call for time curves */
DOUBLE pss_evaluate_curve(struct _ST_NODE* funct, DOUBLE t);

/* print the parsed function */
void pss_parser_print(FILE* out, struct _ST_NODE* node);

#endif
