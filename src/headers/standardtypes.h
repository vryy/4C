/*---------------------------------------------------------------------*/
/*! \file

\brief Definition of PI

\level 0

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#ifndef STANDARDTYPES_H
#define STANDARDTYPES_H

/*----------------------------------------------------------------------*
 | includes of Ansi C standard headers                   m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include <math.h> /* include this before(!) M_PI is tried to be accessed below */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

/*----------------------------------------------------------------------*
 | value for pi (if not provided by math.h)                             |
 *----------------------------------------------------------------------*/
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
/* since a lot of routines still use the old PI macro instead of M_PI: */
#define PI M_PI
/*#define PI               (3.1415926535897932)  used till 18.07.2012 */
/* further options being unused since a long time: */
/*#define PI               (asin(1.0)*2.0) */
/*#define PI  (3.141592653589793238462643383279502884197169399375)*/

/*----------------------------------------------------------------------*
 | definitions file                                       m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "definitions.h"


#endif
