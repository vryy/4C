/*---------------------------------------------------------------------*/
/*! \file

\brief Custom memory allocator user for CLN data type

\level 3


*----------------------------------------------------------------------*/
#include "baci_cut_tolerance.H"
#include "baci_utils_clnwrapper.H"

#include <cmath>
#include <iomanip>
#ifdef CUT_CLN_CALC
#include <cln/malloc.h>
#endif

BACI_NAMESPACE_OPEN

#ifdef CUSTOM_MEMORY_ALLOCATOR


// wrappers to call the method (function pointers and class method pointers are not convertible, so
// we need to use singleton here)
void* mallocwrap(size_t size)
{
  return (CORE::GEO::CUT::MemorySingleton::getInstance().Allocate(size));
}

void deallocwrap(void* ptr) { return (CORE::GEO::CUT::MemorySingleton::getInstance().Free(ptr)); }
// just function pointers to overwrite the ones provided by cln ( malloc and free)
namespace cln
{
  void* (*malloc_hook)(size_t size) = mallocwrap;
  void (*free_hook)(void* ptr) = deallocwrap;
}  // namespace cln

#endif

BACI_NAMESPACE_CLOSE

#ifndef CUT_CLN_CALC

// If this is true, floating point underflow returns zero instead of throwing an exception.
cln::cl_inhibit_floating_point_underflow = true;

// to satisfy linalg fixed matrix
cln::float_format_t::float_format_t() { dserror("Should not be called"); }

#endif
