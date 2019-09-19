/*---------------------------------------------------------------------*/
/*! \file

\brief Custom memory allocator user for CLN data type

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de

*----------------------------------------------------------------------*/
#include "cut_clnwrapper.H"
#include <iomanip>
#include "cut_tolerance.H"
#include <boost/unordered_map.hpp>
#include <cmath>
#ifdef CLN_CALC
#include <cln/malloc.h>
#endif

// initial value of precision_
unsigned int GEO::CUT::ClnWrapper::precision_ = CLN_START_PRECISION;


#ifdef CUSTOM_MEMORY_ALLOCATOR


// wrappers to call the method (function pointers and class method pointers are not convertible, so
// we need to use singleton here)
void* mallocwrap(size_t size) { return (GEO::CUT::MemorySingleton::getInstance().Allocate(size)); }

void deallocwrap(void* ptr) { return (GEO::CUT::MemorySingleton::getInstance().Free(ptr)); }
// just function pointers to overwrite the ones provided by cln ( malloc and free)
namespace cln
{
  void* (*malloc_hook)(size_t size) = mallocwrap;
  void (*free_hook)(void* ptr) = deallocwrap;
}  // namespace cln

#endif

#ifndef CLN_CALC

// If this is true, floating point underflow returns zero instead of throwing an exception.
cln::cl_inhibit_floating_point_underflow = true;

// to satisfy linalg fixed matrix
cln::float_format_t::float_format_t() { dserror("Should not be called"); }

#endif
