/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_factory.H

\brief Factory of scatra elements

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_FACTORY_H
#define SCATRA_ELE_FACTORY_H

#include "scatra_ele_interface.H"

#include "../drt_meshfree_discret/meshfree_scatra_cell_interface.H"

#include "../drt_lib/drt_element.H"

namespace DRT
{

namespace ELEMENTS
{

class ScaTraFactory
{
public:
  //! ctor
  ScaTraFactory(){ return; }

  //! dtor
  virtual ~ScaTraFactory(){ return; }

  //! ProvideImpl
  static ScaTraEleInterface* ProvideImpl(
    DRT::Element::DiscretizationType distype,
    std::string problem,
    const int numdofpernode,
    const int numscal);

  //! ProvideMeshfreeImpl
  static MeshfreeScaTraCellInterface* ProvideMeshfreeImpl(
    DRT::Element::DiscretizationType distype,
    std::string problem,
    const int numdofpernode,
    const int numscal);

private:
  //! define ScatraEle instances dependent on problem
  template<DRT::Element::DiscretizationType distype>
  static ScaTraEleInterface* DefineProblemType(
    std::string problem,
    const int numdofpernode,
    const int numscal);

  //! define MeshfreeScatraCell instances dependent on problem
  template<DRT::Element::DiscretizationType distype>
  static MeshfreeScaTraCellInterface* DefineMeshfreeProblemType(
    std::string problem,
    const int numdofpernode,
    const int numscal);


}; // end class ScaTraFactory

} // end ELEMENTS

} // end DRT

#endif
