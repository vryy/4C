/*----------------------------------------------------------------------*/
/*! \file
 \brief helper for accessing FunctionOfAnything

\level 3

*----------------------------------------------------------------------*/
#ifndef BACI_LIB_GET_FUNCTIONOFANYTHING_H
#define BACI_LIB_GET_FUNCTIONOFANYTHING_H

#include "baci_lib_function.H"
#include "baci_lib_globalproblem.H"

namespace DRT::UTILS
{
  inline const DRT::UTILS::FunctionOfAnything& GetFunctionOfAnything(int functnum)
  {
    return DRT::Problem::Instance()->FunctionById<DRT::UTILS::FunctionOfAnything>(functnum);
  }
}  // namespace DRT::UTILS


#endif
