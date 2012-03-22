/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#ifndef STANDARDTYPES_H
#define STANDARDTYPES_H

/*----------------------------------------------------------------------*
 | includes of Ansi C standard headers                   m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
/*----------------------------------------------------------------------*
 | definitions file                                       m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "definitions.h"
/*----------------------------------------------------------------------*
 | structure types used by the array management           m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "am.h"
/*----------------------------------------------------------------------*
 | various types of enums                                 m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "enums.h"

/*!
\addtogroup FRSYSTEM
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*
 | general problem-variables                              m.gee 4/01    |
 | General information is held here                                     |
 *----------------------------------------------------------------------*/
typedef struct _GENPROB
{
  enum _PROBLEM_TYP probtyp;       /* type of problem, see enum.h */

  INT               ndim;          /* dimension of problem (2 or 3) */

  INT               restart;       /* is restart or not */

  INT               numsf;         /* actual number of struct-field */
  INT               numff;         /* actual number of fluid field */
  INT               numaf;         /* actual number of ale field */
  INT               numtf;         /* actual number of the thermal field */
  INT               numscatra;     /* actual number of the scalar transport field */
  INT               numartf;       /* actual number of the 1D_Artery field */
  INT               numawf;        /* actual number of the reduced dimensional airways field */
  INT               numof;         /* actual number of optimization field */

} GENPROB;




#endif
