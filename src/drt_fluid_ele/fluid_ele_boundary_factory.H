/*!----------------------------------------------------------------------
\file fluid_ele_boundary_factory.H

\brief factory for fluid boundary evaluation

<pre>
Maintainers: Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>
*----------------------------------------------------------------------*/


#ifndef FLUID_ELE_BOUNDARY_FACTORY_H
#define FLUID_ELE_BOUNDARY_FACTORY_H

#include "../drt_lib/drt_element.H"

namespace DRT
{

namespace ELEMENTS
{
  class FluidBoundaryImplInterface;

class FluidBoundaryFactory
{
public:
  //! ctor
  FluidBoundaryFactory(){ return; };

  //! dtor
  virtual ~FluidBoundaryFactory(){ return; };

  //! ProvideImpl
  static FluidBoundaryImplInterface* ProvideImpl(DRT::Element::DiscretizationType distype, std::string problem);

private:
  //! define FluidEleBoundaryCalc instances dependent on problemtype
  template<DRT::Element::DiscretizationType distype>
  static FluidBoundaryImplInterface* DefineProblemType(std::string problem);


};

}
}

#endif
