/*-----------------------------------------------------------*/
/*! \file

\brief Setting of general fluid parameter for element evaluation


\level 1

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_PARAMETER_STD_HPP
#define FOUR_C_FLUID_ELE_PARAMETER_STD_HPP

#include "4C_config.hpp"

#include "4C_fluid_ele_parameter.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    class FluidEleParameterStd : public FluidEleParameter
    {
     public:
      /// Singleton access method
      static FluidEleParameterStd* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

     private:
     protected:
      /// protected Constructor since we are a Singleton.
      FluidEleParameterStd();
    };

  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
