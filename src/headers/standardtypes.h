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
/*----------------------------------------------------------------------*
 | structures concerning domain decomposition             m.gee 8/00    |
 | and intra-communicators                                              |
 *----------------------------------------------------------------------*/
#include "partition.h"

/*!
\addtogroup FRSYSTEM
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>

*----------------------------------------------------------------------*/
typedef struct _FILES
{
/*------------------------------------------------------------- file I/O */
char             *inputfile_name;         /* input file name             */

char             *outputfile_kenner;      /* output file kenner          */
char              outputfile_name[100];   /* output file name            */

FILE             *out_err;                /* file-pointer .err  file     */
} FILES;
/*! @} (documentation module close)*/




/*----------------------------------------------------------------------*
 | general problem-variables                              m.gee 4/01    |
 | General information is held here                                     |
 *----------------------------------------------------------------------*/
typedef struct _GENPROB
{
  enum _PROBLEM_TYP probtyp;       /* type of problem, see enum.h */
  enum _TIME_TYP    timetyp;       /* type of time, see enum.h */

  INT               ndim;          /* dimension of problem (2 or 3) */

  INT               restart;       /* is restart or not */

  INT               adaptive;      /* adaptive mesh algorithms */
  
  INT               numsf;         /* actual number of struct-field */
  INT               numff;         /* actual number of fluid field */
  INT               numaf;         /* actual number of ale field */
  INT               numls;         /* actual number of ls field */
  INT               numtf;         /* actual number of the thermal field */
  INT               numscatra;     /* actual number of the scalar transport field */
  INT               numartf;       /* actual number of the 1D_Artery field */
  INT               numawf;        /* actual number of the reduced dimensional airways field */

} GENPROB;




#endif
