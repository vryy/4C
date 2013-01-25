/*!
\file fluid_ele_factory.H

<pre>
Maintainer: Ursula Rasthofer & Volker Gravemeier
            {rasthofer,vgravem}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_FACTORY_H
#define FLUID_ELE_FACTORY_H

#include "fluid_ele_interface.H"

namespace DRT
{

namespace ELEMENTS
{
class FluidFactory
{
public:
  FluidFactory();

  virtual ~FluidFactory(){}

  //!
  static FluidEleInterface* ProvideImpl(DRT::Element::DiscretizationType distype, string problem);

  //! special implementation of ProvideImpl for XFEM problems to reduce created template combinations
  static FluidEleInterface* ProvideImplXFEM(DRT::Element::DiscretizationType distype, string problem);

private:
  //! define FluidEle instances dependent on problem
  template<DRT::Element::DiscretizationType distype>
  static FluidEleInterface* DefineProblemType(string problem);

  //! create FluidEleCalcXFEM instance for problem type xfem
  template<DRT::Element::DiscretizationType distype>
  static FluidEleInterface* DefineProblemTypeXFEM(string problem);

};

}
}

#endif
