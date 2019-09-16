/*---------------------------------------------------------------------*/
/*! \file

\brief C-style definitions of data types

\level 3

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | NOTE:                                                                |
 | - if changes or additions are made to this file, a complete recompile|
 |   of the whole code is recommended                                   |
 | - always use strict upper case letters                               |
 |                                                                      |
 *----------------------------------------------------------------------*/

#ifndef CCARAT_DATATYPES_H
#define CCARAT_DATATYPES_H

/*----------------------------------------------------------------------*
 | basic data types, do not use INT, DOUBLE or CHAR in baci !           |
 | only still here to allow the compilation of shell8 and some C stuff  |
 *----------------------------------------------------------------------*/
#ifdef INT
#undef INT
#endif
typedef int INT;
#ifdef DOUBLE
#undef DOUBLE
#endif
typedef double DOUBLE;
#ifdef CHAR
#undef CHAR
#endif
typedef char CHAR;


#endif /* CCARAT_DATATYPES_H */
