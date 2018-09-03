/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_parameter_std.H

\brief Setting of general fluid parameter for element evaluation

<pre>
Maintainers: Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/

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
