#ifndef XFEM_FLUIDWIZARD_H
#define XFEM_FLUIDWIZARD_H

#include <map>
#include <set>

#include "../drt_geometry/integrationcell.H"
#include "../drt_geometry/geo_intersection.H"

#include "physics.H"

namespace DRT
{
  class Discretization;
}

namespace XFEM
{
  class FieldEnr;
  class ElementAnsatz;
  class FluidDofSet;

class FluidWizard
{
public:

  FluidWizard( DRT::Discretization & backdis, DRT::Discretization & cutterdis )
    : backdis_( backdis ),
      cutterdis_( cutterdis )
  {
  }

  void Cut(  bool include_inner, const Epetra_Vector & idispcol, std::string gausstype, bool positions=true  );

  void Cut( const Epetra_Vector & idispcol,
            std::map< int, GEO::DomainIntCells > & domainintcells,
            std::map< int, GEO::BoundaryIntCells > & boundaryintcells,
            const std::map<int,int>& labelPerElementId,
            const std::vector<int>& MovingFluideleGIDs,
            std::string gausstype,
            bool positions=true );

  void CreateDofMap( std::map<int, const std::set<XFEM::FieldEnr> >&   nodalDofSetFinal,    ///< enriched fields per node
                     std::map<int, const std::set<XFEM::FieldEnr> >&   elementalDofsFinal,  ///< enriched fields per element
                     const std::set<XFEM::PHYSICS::Field>&             fieldset,
                     const XFEM::ElementAnsatz&                        elementAnsatz,
                     const Teuchos::ParameterList&                     params );

  Teuchos::RCP<FluidDofSet> DofSet();

  GEO::CutWizard & CutWizard()
  {
    return *cut_;
  }

  GEO::CUT::ElementHandle * GetElement( DRT::Element * ele ) { return cut_->GetElement( ele ); }

  GEO::CUT::Node * GetNode( int nid ) { return cut_->GetNode( nid ); }

private:

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
