/*----------------------------------------------------------------------*/
/*! \file
 \brief helper for accessing FunctionOfAnything

\level 3

*----------------------------------------------------------------------*/
#ifndef LIBGETFUNCTIONOFANYTHING_H
#define LIBGETFUNCTIONOFANYTHING_H

#include "lib_function.H"
#include "lib_globalproblem.H"

namespace DRT::UTILS
{
  inline DRT::UTILS::FunctionOfAnything& GetFunctionOfAnything(int functnum)
  {
    return DRT::Problem::Instance()->FunctionById<DRT::UTILS::FunctionOfAnything>(functnum);
  }
}  // namespace DRT::UTILS


#endif