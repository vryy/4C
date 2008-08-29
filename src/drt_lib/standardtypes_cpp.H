/*!
 * \file standardtypes_cpp.H
 *
 * Use this header when you need access to standardtypes.h 
 * 
 * (this is to ensure that the C++ mpi.h header is included before the c version)
 * <pre>
 * Maintainer: Axel Gerstenberger
 *             gerstenberger@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15236
 * </pre>
 */
#ifndef STANDARDTYPES_CPP_H
#define STANDARDTYPES_CPP_H

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include <mpi.h>
#endif
extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}

#endif
