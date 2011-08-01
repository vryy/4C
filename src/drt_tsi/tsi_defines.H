/*!----------------------------------------------------------------------
\file tsi_defines.H

\brief A set of preprocessor defines for tsi methods

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef TSI_DEFINES_H
#define TSI_DEFINES_H

/************************************************************************/
/* TSI algorithm options                                                */
/************************************************************************/

//#define COUPLEINITTEMPERATURE       /* flag for THR: constant temperature T_0 in couple term  */

/************************************************************************/
/* Debugging options                                                    */
/************************************************************************/

// GENERAL DEBUGGING OPTIONS
//#define TSIASOUTPUT          /* flag for detailed active set output */
//#define TSIMONOLITHASOUTPUT  /* flag for detailed output in monolithic TSI */
//#define TSIBLOCKMATRIXMERGE  /* flag to activate merge of systemmatrix in monolithic TSI */
//#define TSIPARALLEL          /* flag for parallel TSI */
//#define THRASOUTPUT          /* flag for output in the thermal field */
//#define DEBUGMATERIAL        /* flag for material debug output */
//#define MonTSIwithoutTHR     /* flag to comment out all coupling terms in STR field */
//#define MonTSIwithoutSTR     /* flag to comment out all coupling terms in THR field */

#endif  // #ifndef TSI_DEFINES_H
#endif  // #ifdef CCADISCRET
