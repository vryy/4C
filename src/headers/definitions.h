/*!---------------------------------------------------------------------
\file definitions.h
\brief C-style definitions of some frequently used code idioms
\level 3
\maintainer Martin Kronbichler
---------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | NOTE:                                                                |
 | - if changes or additions are made to this file, a complete recompile|
 |   of the whole code is recommended                                   |
 | - always use strict upper case letters                               |
 |                                                                      |
 *----------------------------------------------------------------------*/

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

/*----------------------------------------------------------------------*
 | sign of an integer (used in drt_geometry/tetrahedradecomposition.cpp)|
 *----------------------------------------------------------------------*/
#define SIGN(x)    ((x) <  0  ? (-1) : (1))

/*----------------------------------------------------------------------*
 | square of a double (widely used)                                     |
 *----------------------------------------------------------------------*/
#define DSQR(a) ((a)*(a))

/*----------------------------------------------------------------------*
 | a set of different tolerances                                        |
 *----------------------------------------------------------------------*/
#define EPS1             (1.0E-01)
#define EPS2             (1.0E-02)
#define EPS3             (1.0E-03)
#define EPS4             (1.0E-04)
#define EPS5             (1.0E-05)
#define EPS6             (1.0E-06)
#define EPS7             (1.0E-07)
#define EPS8             (1.0E-08)
#define EPS9             (1.0E-09)
#define EPS10            (1.0E-10)
#define EPS11            (1.0E-11)
#define EPS12            (1.0E-12)
#define EPS13            (1.0E-13)
#define EPS14            (1.0E-14)
#define EPS15            (1.0E-15)

/*----------------------------------------------------------------------*
 | a set of numbers (still used in shell8)                              |
 *----------------------------------------------------------------------*/
#define VERYLARGEINT     (1000000000)
#define VERYLARGEREAL    (1000000000.0)

#endif /* DEFINITIONS_H */
