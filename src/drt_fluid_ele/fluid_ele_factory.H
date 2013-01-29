/*--------------------------------------------------------------------------*/
/*!
\file fluid_ele_factory.H

\brief Factory of fluid elements

<pre>
Maintainer: Ursula Rasthofer & Volker Gravemeier
            {rasthofer,vgravem}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15236/245
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef FLUID_ELE_FACTORY_H
#define FLUID_ELE_FACTORY_H

#include "fluid_ele_interface.H"

namespace DRT
{

namespace ELEMENTS
{

/*--------------------------------------------------------------------------*/
/*!
 * \brief to do
 *
 * \author Ursula Rasthofer (rasthofer@lnm.mw.tum.de)
 *
 * \date March, 2012
 */
/*--------------------------------------------------------------------------*/
class FluidFactory
{
public:
  //! ctor
  FluidFactory(){ return; }

  //! dtor
  virtual ~FluidFactory(){ return; }

  //! ProvideImpl
  static FluidEleInterface* ProvideImpl(DRT::Element::DiscretizationType distype, string problem);

  //! special ProvideImpl for XFEM problems to reduce created template combinations
  static FluidEleInterface* ProvideImplXFEM(DRT::Element::DiscretizationType distype, string problem);

  //! special ProvideImpl for meshfree problems to reduce created template combinations
  static FluidEleInterface* ProvideImplMeshfree(DRT::Element::DiscretizationType distype, string problem);

private:
  //! define FluidEle instances dependent on problem
  template<DRT::Element::DiscretizationType distype>
  static FluidEleInterface* DefineProblemType(string problem);

}; // end class FluidFactory

} // end ELEMENTS

} // end DRT

#endif
