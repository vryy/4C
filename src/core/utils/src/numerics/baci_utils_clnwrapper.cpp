/*---------------------------------------------------------------------*/
/*! \file

\brief Custom memory allocator user for CLN data type

\level 3


*----------------------------------------------------------------------*/
#include "baci_utils_clnwrapper.H"

#include "baci_cut_tolerance.H"

#include <cmath>
#include <iomanip>
#ifdef CUT_CLN_CALC
#include <cln/malloc.h>
#endif

BACI_NAMESPACE_OPEN

// initial value of precision_
unsigned int CORE::CLN::ClnWrapper::precision_ = CLN_START_PRECISION;

BACI_NAMESPACE_CLOSE
