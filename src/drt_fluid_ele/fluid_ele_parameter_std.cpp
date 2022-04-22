/*-----------------------------------------------------------*/
/*! \file

\brief Setting of general fluid parameter for element evaluation


\level 1

*/
/*-----------------------------------------------------------*/

#include "fluid_ele_parameter_std.H"

DRT::ELEMENTS::FluidEleParameterStd* DRT::ELEMENTS::FluidEleParameterStd::Instance(
    ::UTILS::SingletonAction action)
{
  static ::UTILS::SingletonOwner<DRT::ELEMENTS::FluidEleParameterStd> singleton_owner([]() {
    return std::unique_ptr<DRT::ELEMENTS::FluidEleParameterStd>(
        new DRT::ELEMENTS::FluidEleParameterStd());
  });

  return singleton_owner.Instance(action);
}

//----------------------------------------------------------------------*/
//    destruction method
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameterStd::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(::UTILS::SingletonAction::destruct);
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleParameterStd::FluidEleParameterStd()
    : DRT::ELEMENTS::FluidEleParameter::FluidEleParameter()
{
}
