/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_factory.H

\brief Factory of scatra elements

<pre>
\level 1

\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_FACTORY_H
#define SCATRA_ELE_FACTORY_H

#include "scatra_ele_interface.H"

#include "../drt_meshfree_discret/meshfree_scatra_cell_interface.H"
#include "../drt_inpar/inpar_scatra.H"

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
    INPAR::SCATRA::ImplType problem,
    const int numdofpernode,
    const int numscal,
    const std::string& disname);

  //! ProvideImplHDG
  static ScaTraEleInterface* ProvideImplHDG(
    DRT::Element::DiscretizationType distype,
    INPAR::SCATRA::ImplType problem,
    const int numdofpernode,
    const int numscal,
    const std::string& disname);

  //! ProvideMeshfreeImpl
  static MeshfreeScaTraCellInterface* ProvideMeshfreeImpl(
    DRT::Element::DiscretizationType distype,
    INPAR::SCATRA::ImplType problem,
    const int numdofpernode,
    const int numscal,
    const std::string& disname);

private:
  //! define ScatraEle instances dependent on problem
  template<DRT::Element::DiscretizationType distype,int probdim>
  static ScaTraEleInterface* DefineProblemType(
    INPAR::SCATRA::ImplType problem,
    const int numdofpernode,
    const int numscal,
    const std::string& disname);

  //! define ScatraEle instances dependent on problem
  template<DRT::Element::DiscretizationType distype,int probdim>
  static ScaTraEleInterface* DefineProblemTypeHDG(
    INPAR::SCATRA::ImplType problem,
    const int numdofpernode,
    const int numscal,
    const std::string& disname);

  //! define MeshfreeScatraCell instances dependent on problem
  template<DRT::Element::DiscretizationType distype>
  static MeshfreeScaTraCellInterface* DefineMeshfreeProblemType(
    INPAR::SCATRA::ImplType impltype,
    const int numdofpernode,
    const int numscal,
    const std::string& disname);


}; // end class ScaTraFactory

} // end ELEMENTS

} // end DRT

#endif
