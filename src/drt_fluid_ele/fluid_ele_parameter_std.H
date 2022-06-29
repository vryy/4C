/*-----------------------------------------------------------*/
/*! \file

\brief Setting of general fluid parameter for element evaluation


\level 1

*/
/*-----------------------------------------------------------*/

#ifndef FLUID_ELE_PARAMETER_STD_H
#define FLUID_ELE_PARAMETER_STD_H

#include "fluid_ele_parameter.H"

#include "../headers/singleton_owner.H"

namespace DRT
{
  namespace ELEMENTS
  {
    class FluidEleParameterStd : public FluidEleParameter
    {
     public:
      /// Singleton access method
      static FluidEleParameterStd* Instance(
          ::UTILS::SingletonAction action = ::UTILS::SingletonAction::create);

     private:
     protected:
      /// protected Constructor since we are a Singleton.
      FluidEleParameterStd();
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
