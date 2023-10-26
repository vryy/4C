/*----------------------------------------------------------------------*/
/*! \file
 \brief helper for accessing FunctionOfAnything

\level 3

*----------------------------------------------------------------------*/
#ifndef BACI_UTILS_GET_FUNCTIONOFANYTHING_H
#define BACI_UTILS_GET_FUNCTIONOFANYTHING_H

#include "baci_lib_globalproblem.H"
#include "baci_utils_function.H"

namespace DRT::UTILS
{
  inline const CORE::UTILS::FunctionOfAnything& GetFunctionOfAnything(int functnum)
  {
    return DRT::Problem::Instance()->FunctionById<CORE::UTILS::FunctionOfAnything>(functnum);
  }
}  // namespace DRT::UTILS


#endif
