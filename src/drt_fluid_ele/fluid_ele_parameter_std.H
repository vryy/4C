/*-----------------------------------------------------------*/
/*!

\brief Setting of general fluid parameter for element evaluation

\maintainer Martin Kronbichler

\level 1

*/
/*-----------------------------------------------------------*/

#ifndef FLUID_ELE_PARAMETER_STD_H
#define FLUID_ELE_PARAMETER_STD_H

#include "fluid_ele_parameter.H"


namespace DRT
{
  namespace ELEMENTS
  {
    class FluidEleParameterStd : public FluidEleParameter
    {
     public:
      virtual ~FluidEleParameterStd() {}

      /// Singleton access method
      static FluidEleParameterStd* Instance(bool create = true);

      /// called upon destruction
      virtual void Done();

     private:
     protected:
      /// protected Constructor since we are a Singleton.
      FluidEleParameterStd();
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
