/*----------------------------------------------------------------------*/
/*! \file

\brief A set of preprocessor defines for tsi methods

\level 1

\maintainer Christoph Meier

*----------------------------------------------------------------------*/
#ifndef TSI_DEFINES_H
#define TSI_DEFINES_H

/************************************************************************/
/* Debugging options for TSI algorithm                                  */
/************************************************************************/
//#define COUPLEINITTEMPERATURE  /* flag for THR: constant temperature T_0 in couple term  */
//#define TSIPARALLEL            /* flag for parallel TSI */
//#define MonTSIwithoutTHR       /* flag to comment out all coupling terms in STR field */
//#define MonTSIwithoutSTR       /* flag to comment out all coupling terms in THR field */
//#define TFSI                   /* flag to reduce the output to screen */
//#define TSI_DEBUG              /* output to screen information */


/************************************************************************/
/* Debugging options dependent on problem type                          */
/************************************************************************/

// GENERAL DEBUGGING OPTIONS
//#define TSIASOUTPUT            /* flag for detailed active set output */
//#define TSIMONOLITHASOUTPUT    /* flag for detailed output in monolithic TSI */
//#define TSIPARTITIONEDASOUTPUT /* flag for detailed output in partitioned TSI */


/************************************************************************/
/* Output options for thermal problems                                  */
/************************************************************************/

//#define THRASOUTPUT            /* flag for output in the thermal field */
//#define CALCSTABILOFREACTTERM  /* flag for output of stabilisation parameter in THR equation */


/************************************************************************/
/* Debugging options for plastic materials                              */
/************************************************************************/
//#define DEBUGMATERIAL          /* flag for material debug output */

/************************************************************************/
/* Debugging options for pseudo-SLM in TSI                              */
/************************************************************************/
//#define TSISLMFDCHECK         /* flag for enabling FDChecks*/
//#define TSISLMFDCHECKDEBUG    /* prints element matrices in FDcheck to a file "FDCheck_capa.txt"*/

#endif  // #ifndef TSI_DEFINES_H
