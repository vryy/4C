/*-----------------------------------------------------------*/
/*! \file

\brief Setting of general fluid parameter for element evaluation


\level 1

*/
/*-----------------------------------------------------------*/

#include "baci_fluid_ele_parameter_std.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::FluidEleParameterStd* DRT::ELEMENTS::FluidEleParameterStd::Instance(
    CORE::UTILS::SingletonAction action)
{
  static auto singleton_owner = CORE::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<DRT::ELEMENTS::FluidEleParameterStd>(
            new DRT::ELEMENTS::FluidEleParameterStd());
      });

  return singleton_owner.Instance(action);
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleParameterStd::FluidEleParameterStd()
    : DRT::ELEMENTS::FluidEleParameter::FluidEleParameter()
{
}

FOUR_C_NAMESPACE_CLOSE
