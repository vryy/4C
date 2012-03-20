#ifndef XFEM_FLUIDWIZARD_H
#define XFEM_FLUIDWIZARD_H

#include <map>
#include <set>


#include "../drt_geometry/geo_utils.H"
#include "../drt_geometry/geo_intersection.H"

#include "physics.H"


class Epetra_Vector;
class Epetra_IntVector;

namespace DRT
{
  class Discretization;
}

namespace Teuchos
{
  class ParameterList;
}

namespace XFEM
{
  class FieldEnr;
  class ElementAnsatz;
  class FluidDofSet;

class FluidWizard
{
public:

  /*!
  \brief The interface that bridges the cut libraries and fluid part of XFEM
  */
  FluidWizard( DRT::Discretization & backdis, DRT::Discretization & cutterdis )
    : backdis_( backdis ),
      cutterdis_( cutterdis )
  {
  }

  /*!
  \brief Cut routine for the new XFEM framework (XFSI and XFLUIDFLUID)
  */
  void Cut(  bool include_inner, const Epetra_Vector & idispcol, bool parallel, std::string VCellgausstype,
      std::string BCellgausstype, bool positions=true  );

  Teuchos::RCP<FluidDofSet> DofSet(int maxNumMyReservedDofs);

  GEO::CutWizard & CutWizard()
  {
    return *cut_;
  }

  GEO::CUT::ElementHandle * GetElement( DRT::Element * ele ) { return cut_->GetElement( ele ); }

  GEO::CUT::Node * GetNode( int nid ) { return cut_->GetNode( nid ); }

private:

  /*!
  \brief Writes GMSH output of integrationcells for old cut algorithm
   */
  void DumpGmshIntegrationCells( std::map< int, GEO::DomainIntCells > & domainintcells,
                                 std::map< int, GEO::BoundaryIntCells > & boundaryintcells );

  Teuchos::RCP<GEO::CutWizard> cut_;

  /// background discretization
  DRT::Discretization & backdis_;

  /// cutter discretization
  DRT::Discretization & cutterdis_;

  Teuchos::RCP<Epetra_IntVector> position_;

#if 0
  /// domainintegrationcells for each intersected element
  std::map< int, GEO::DomainIntCells > domainintcells_;

  /// boundaryintegrationcells for each intersected element
  std::map< int, GEO::BoundaryIntCells > boundaryintcells_;

  /// which physical fields are used in this problem
  std::set<XFEM::PHYSICS::Field> fieldset_;

  /// element dofs for this problem
  Teuchos::RCP<XFEM::ElementAnsatz> elementAnsatz_;

  Teuchos::RCP<XFEM::DofManager> dofmanager_;
#endif
};

}

#endif
