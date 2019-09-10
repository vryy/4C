/*-----------------------------------------------------------*/
/*! \file

\brief Setting of general fluid parameter for element evaluation

\maintainer Martin Kronbichler

\level 1

*/
/*-----------------------------------------------------------*/

#include "fluid_ele_parameter_std.H"

//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleParameterStd* DRT::ELEMENTS::FluidEleParameterStd::Instance(bool create)
{
  static FluidEleParameterStd* instance;
  if (create)
  {
    if (instance == NULL)
    {
      instance = new FluidEleParameterStd();
    }
  }
  else
  {
    if (instance != NULL) delete instance;
    instance = NULL;
  }
  return instance;
}

//----------------------------------------------------------------------*/
//    destruction method
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameterStd::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleParameterStd::FluidEleParameterStd()
    : DRT::ELEMENTS::FluidEleParameter::FluidEleParameter()
{
}
