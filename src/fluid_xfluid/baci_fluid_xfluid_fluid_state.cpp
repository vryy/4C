/*----------------------------------------------------------------------*/
/*! \file

\brief State class for (in)stationary XFEM fluid problems involving embedded
fluid meshes

\level 2

 */
/*----------------------------------------------------------------------*/

#include "baci_fluid_xfluid_fluid_state.hpp"

#include "baci_fluid_utils.hpp"
#include "baci_fluid_utils_mapextractor.hpp"
#include "baci_fluid_xfluid_state.hpp"
#include "baci_global_data.hpp"
#include "baci_linalg_mapextractor.hpp"
#include "baci_linalg_utils_sparse_algebra_create.hpp"
#include "baci_linalg_utils_sparse_algebra_manipulation.hpp"
#include "baci_xfem_condition_manager.hpp"
#include "baci_xfem_xfield_state_utils.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Constructor for XFluidFluidState                         kruse 01/15 |
 *----------------------------------------------------------------------*/
FLD::XFluidFluidState::XFluidFluidState(
    const Teuchos::RCP<XFEM::ConditionManager>& condition_manager,
    const Teuchos::RCP<CORE::GEO::CutWizard>& wizard, const Teuchos::RCP<XFEM::XFEMDofSet>& dofset,
    const Teuchos::RCP<const Epetra_Map>& xfluiddofrowmap,
    const Teuchos::RCP<const Epetra_Map>& xfluiddofcolmap,
    const Teuchos::RCP<const Epetra_Map>& embfluiddofrowmap)
    : XFluidState(condition_manager, wizard, dofset, xfluiddofrowmap, xfluiddofcolmap),
      xffluiddofrowmap_(CORE::LINALG::MergeMap(xfluiddofrowmap, embfluiddofrowmap, false)),
      xffluidsplitter_(Teuchos::rcp(new FLD::UTILS::XFluidFluidMapExtractor())),
      xffluidvelpressplitter_(Teuchos::rcp(new CORE::LINALG::MapExtractor())),
      embfluiddofrowmap_(embfluiddofrowmap)
{
  xffluidsplitter_->Setup(*xffluiddofrowmap_, xfluiddofrowmap, embfluiddofrowmap);
  InitSystemMatrix();
  InitStateVectors();
}

/*----------------------------------------------------------------------*
 |  Initialize (coupling) matrices & rhs-vectors
 |                                                         kruse 01/15
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidState::InitSystemMatrix()
{
  // the combined fluid system matrix is not of FECrs-type - it is solely composed out of
  // fully assembled submatrices
  xffluidsysmat_ =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*xffluiddofrowmap_, 108, false, true));
}

/*----------------------------------------------------------------------*
 |  Initialize state vectors                                kruse 01/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidState::InitStateVectors()
{
  // matrices & vectors for merged background & embedded fluid
  xffluidresidual_ = CORE::LINALG::CreateVector(*xffluiddofrowmap_, true);
  xffluidincvel_ = CORE::LINALG::CreateVector(*xffluiddofrowmap_, true);
  xffluidvelnp_ = CORE::LINALG::CreateVector(*xffluiddofrowmap_, true);
  xffluidveln_ = CORE::LINALG::CreateVector(*xffluiddofrowmap_, true);
  xffluidzeros_ = CORE::LINALG::CreateVector(*xffluiddofrowmap_, true);
}

/*----------------------------------------------------------------------*
 |  Access system matrix                                    kruse 01/15 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> FLD::XFluidFluidState::SystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(xffluidsysmat_);
}

/*----------------------------------------------------------------------*
 | Complete coupling matrices and rhs vectors              schott 12/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidState::CompleteCouplingMatricesAndRhs()
{
  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat_block =
      Teuchos::rcp_dynamic_cast<CORE::LINALG::BlockSparseMatrixBase>(xffluidsysmat_, false);

  // in case that fluid-fluid sysmat is merged (no block matrix), we have to complete the coupling
  // blocks (e.g. fluid-structure) w.r.t. xff sysmat instead of just the xfluid block

  if (sysmat_block == Teuchos::null)  // merged matrix
  {
    XFluidState::CompleteCouplingMatricesAndRhs(xffluiddofrowmap_);
  }
  else  // block matrix
  {
    XFluidState::CompleteCouplingMatricesAndRhs(xfluiddofrowmap_);
  }
}

/*----------------------------------------------------------------------*
 |  Create merged DBC map extractor                         kruse 01/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidState::CreateMergedDBCMapExtractor(
    Teuchos::RCP<const CORE::LINALG::MapExtractor> embfluiddbcmaps)
{
  // create merged dbc map from both fluids
  std::vector<Teuchos::RCP<const Epetra_Map>> dbcmaps;
  dbcmaps.push_back(XFluidState::dbcmaps_->CondMap());
  dbcmaps.push_back(embfluiddbcmaps->CondMap());

  Teuchos::RCP<const Epetra_Map> xffluiddbcmap =
      CORE::LINALG::MultiMapExtractor::MergeMaps(dbcmaps);

  std::vector<Teuchos::RCP<const Epetra_Map>> othermaps;
  othermaps.push_back(XFluidState::dbcmaps_->OtherMap());
  othermaps.push_back(embfluiddbcmaps->OtherMap());
  Teuchos::RCP<const Epetra_Map> xffluidothermap =
      CORE::LINALG::MultiMapExtractor::MergeMaps(othermaps);

  xffluiddbcmaps_ = Teuchos::rcp(
      new CORE::LINALG::MapExtractor(*xffluiddofrowmap_, xffluiddbcmap, xffluidothermap));
}

/*----------------------------------------------------------------------*
 |  Set dirichlet- and velocity/pressure-map extractor      kruse 01/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidState::SetupMapExtractors(
    const Teuchos::RCP<DRT::Discretization>& xfluiddiscret,
    const Teuchos::RCP<DRT::Discretization>& embfluiddiscret, const double& time)
{
  // create merged dirichlet map extractor
  XFluidState::SetupMapExtractors(xfluiddiscret, time);
  xffluidsplitter_->Setup(*xffluiddofrowmap_, embfluiddofrowmap_, XFluidState::xfluiddofrowmap_);

  FLD::UTILS::SetupFluidFluidVelPresSplit(*xfluiddiscret, GLOBAL::Problem::Instance()->NDim(),
      *embfluiddiscret, *xffluidvelpressplitter_, xffluiddofrowmap_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool FLD::XFluidFluidState::Destroy()
{
  // destroy system matrix
  std::cout << "Destroying the xffluidsysmat_ is not possible at the moment. Internally more "
               "strong RCPs point to the EpetraMatrix. This has to be checked!!!"
            << std::endl;

  XFEM::DestroyRCPObject(xffluidvelnp_);
  XFEM::DestroyRCPObject(xffluidveln_);

  XFEM::DestroyRCPObject(xffluidresidual_);

  XFEM::DestroyRCPObject(xffluidzeros_);
  XFEM::DestroyRCPObject(xffluidincvel_);


  XFEM::DestroyRCPObject(xffluidsplitter_);
  XFEM::DestroyRCPObject(xffluidvelpressplitter_);
  XFEM::DestroyRCPObject(xffluiddbcmaps_);

  // destroy dofrowmap
  XFEM::DestroyRCPObject(embfluiddofrowmap_);

  // TODO: actually it should be possible to delete the dofrowmap, however this causes problems in
  // xffsi applications! (CHECK THIS!!!)
  // DofRowMap() in Xfluidfluid currently returns a strong RCP
  if (xffluiddofrowmap_.strong_count() == 1)
    xffluiddofrowmap_ = Teuchos::null;
  else  // FOUR_C_THROW("could not destroy object: %i!=1 pointers",
        // xffluiddofrowmap_.strong_count());
    std::cout << "could not destroy xffluiddofrowmap_: number of pointers is "
              << xffluiddofrowmap_.strong_count() << "!=1";

  XFluidState::Destroy();

  return true;
}

FOUR_C_NAMESPACE_CLOSE
