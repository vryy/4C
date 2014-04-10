/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_factory.H

\brief factory for fluid scatra boundary evaluation

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 */
/*----------------------------------------------------------------------*/


#ifndef SCATRA_ELE_BOUNDARY_FACTORY_H
#define SCATRA_ELE_BOUNDARY_FACTORY_H

#include "../drt_lib/drt_element.H"
#include "../drt_inpar/inpar_scatra.H"

namespace DRT
{

namespace ELEMENTS
{

class ScaTraBoundaryInterface;

class ScaTraBoundaryFactory
{

public:
  //! ctor
  ScaTraBoundaryFactory(){ return; };

  //! dtor
  virtual ~ScaTraBoundaryFactory(){ return; };

  //! ProvideImpl
  static ScaTraBoundaryInterface* ProvideImpl(const DRT::Element* ele,
                                              const enum INPAR::SCATRA::ScaTraType scatratype,
                                              const int numdofpernode,
                                              const int numscal);


private:

  //! define FluidEleBoundaryCalc instances dependent on problemtype
  template<DRT::Element::DiscretizationType distype>
  static ScaTraBoundaryInterface* DefineProblemType(const enum INPAR::SCATRA::ScaTraType scatratype,
                                                    const int numdofpernode,
                                                    const int numscal);

}; // class ScaTraBoundaryFactory

} // namespace ELEMENTS
} // namespace DRT

#endif
