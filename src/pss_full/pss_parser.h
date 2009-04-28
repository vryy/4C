/*!
\file
\brief A top down parser for function evaluation.

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich
              
Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed, 
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions and of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de) 
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de                   

-------------------------------------------------------------------------
<\pre>

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

/* a Heaviside function */
double heaviside();

#endif
