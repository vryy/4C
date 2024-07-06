/*-----------------------------------------------------------*/
/*! \file

\brief Setting of general fluid parameter for element evaluation


\level 1

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_ele_parameter_std.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::FluidEleParameterStd* Discret::ELEMENTS::FluidEleParameterStd::instance(
    Core::UTILS::SingletonAction action)
{
  static auto singleton_owner = Core::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<Discret::ELEMENTS::FluidEleParameterStd>(
            new Discret::ELEMENTS::FluidEleParameterStd());
      });

  return singleton_owner.instance(action);
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
Discret::ELEMENTS::FluidEleParameterStd::FluidEleParameterStd()
    : Discret::ELEMENTS::FluidEleParameter::FluidEleParameter()
{
}

FOUR_C_NAMESPACE_CLOSE
