/*--------------------------------------------------------------------------*/
/*!
\file fluid_ele_factory.H

\brief Factory of fluid elements

<pre>
\maintainer Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15236/245
</pre>

\level 1
*/
/*--------------------------------------------------------------------------*/

#ifndef FLUID_ELE_FACTORY_H
#define FLUID_ELE_FACTORY_H

#include "../drt_lib/drt_element.H"

namespace DRT
{

namespace ELEMENTS
{
  class FluidEleInterface;

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
  static FluidEleInterface* ProvideImpl(DRT::Element::DiscretizationType distype, std::string problem);

  //! special ProvideImpl for XFEM problems to reduce created template combinations
  static FluidEleInterface* ProvideImplXFEM(DRT::Element::DiscretizationType distype, std::string problem);

private:
  //! define FluidEle instances dependent on problem
  template<DRT::Element::DiscretizationType distype>
  static FluidEleInterface* DefineProblemType(std::string problem);

  //! special DefineProblemTypeXFEM for XFEM problems
  template<DRT::Element::DiscretizationType distype>
  static FluidEleInterface* DefineProblemTypeXFEM(std::string problem);

}; // end class FluidFactory

} // end ELEMENTS

} // end DRT

#endif
