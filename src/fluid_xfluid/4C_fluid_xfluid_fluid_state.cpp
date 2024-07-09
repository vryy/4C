/*----------------------------------------------------------------------*/
/*! \file

\brief State class for (in)stationary XFEM fluid problems involving embedded
fluid meshes

\level 2

 */
/*----------------------------------------------------------------------*/

#include "4C_fluid_xfluid_fluid_state.hpp"

#include "4C_fluid_utils.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fluid_xfluid_state.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_xfem_condition_manager.hpp"
#include "4C_xfem_xfield_state_utils.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Constructor for XFluidFluidState                         kruse 01/15 |
 *----------------------------------------------------------------------*/
FLD::XFluidFluidState::XFluidFluidState(
    const Teuchos::RCP<XFEM::ConditionManager>& condition_manager,
    const Teuchos::RCP<Core::Geo::CutWizard>& wizard, const Teuchos::RCP<XFEM::XFEMDofSet>& dofset,
    const Teuchos::RCP<const Epetra_Map>& xfluiddofrowmap,
    const Teuchos::RCP<const Epetra_Map>& xfluiddofcolmap,
    const Teuchos::RCP<const Epetra_Map>& embfluiddofrowmap)
    : XFluidState(condition_manager, wizard, dofset, xfluiddofrowmap, xfluiddofcolmap),
      xffluiddofrowmap_(Core::LinAlg::MergeMap(xfluiddofrowmap, embfluiddofrowmap, false)),
      xffluidsplitter_(Teuchos::rcp(new FLD::UTILS::XFluidFluidMapExtractor())),
      xffluidvelpressplitter_(Teuchos::rcp(new Core::LinAlg::MapExtractor())),
      embfluiddofrowmap_(embfluiddofrowmap)
{
  xffluidsplitter_->setup(*xffluiddofrowmap_, xfluiddofrowmap, embfluiddofrowmap);
  init_system_matrix();
  init_state_vectors();
}

/*----------------------------------------------------------------------*
 |  Initialize (coupling) matrices & rhs-vectors
 |                                                         kruse 01/15
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidState::init_system_matrix()
{
  // the combined fluid system matrix is not of FECrs-type - it is solely composed out of
  // fully assembled submatrices
  xffluidsysmat_ =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*xffluiddofrowmap_, 108, false, true));
}

/*----------------------------------------------------------------------*
 |  Initialize state vectors                                kruse 01/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidState::init_state_vectors()
{
  // matrices & vectors for merged background & embedded fluid
  xffluidresidual_ = Core::LinAlg::CreateVector(*xffluiddofrowmap_, true);
  xffluidincvel_ = Core::LinAlg::CreateVector(*xffluiddofrowmap_, true);
  xffluidvelnp_ = Core::LinAlg::CreateVector(*xffluiddofrowmap_, true);
  xffluidveln_ = Core::LinAlg::CreateVector(*xffluiddofrowmap_, true);
  xffluidzeros_ = Core::LinAlg::CreateVector(*xffluiddofrowmap_, true);
}

/*----------------------------------------------------------------------*
 |  Access system matrix                                    kruse 01/15 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix> FLD::XFluidFluidState::system_matrix()
{
  return Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(xffluidsysmat_);
}

/*----------------------------------------------------------------------*
 | Complete coupling matrices and rhs vectors              schott 12/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidState::complete_coupling_matrices_and_rhs()
{
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmat_block =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(xffluidsysmat_, false);

  // in case that fluid-fluid sysmat is merged (no block matrix), we have to complete the coupling
  // blocks (e.g. fluid-structure) w.r.t. xff sysmat instead of just the xfluid block

  if (sysmat_block == Teuchos::null)  // merged matrix
  {
    XFluidState::complete_coupling_matrices_and_rhs(xffluiddofrowmap_);
  }
  else  // block matrix
  {
    XFluidState::complete_coupling_matrices_and_rhs(xfluiddofrowmap_);
  }
}

/*----------------------------------------------------------------------*
 |  Create merged DBC map extractor                         kruse 01/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidState::create_merged_dbc_map_extractor(
    Teuchos::RCP<const Core::LinAlg::MapExtractor> embfluiddbcmaps)
{
  // create merged dbc map from both fluids
  std::vector<Teuchos::RCP<const Epetra_Map>> dbcmaps;
  dbcmaps.push_back(XFluidState::dbcmaps_->cond_map());
  dbcmaps.push_back(embfluiddbcmaps->cond_map());

  Teuchos::RCP<const Epetra_Map> xffluiddbcmap =
      Core::LinAlg::MultiMapExtractor::merge_maps(dbcmaps);

  std::vector<Teuchos::RCP<const Epetra_Map>> othermaps;
  othermaps.push_back(XFluidState::dbcmaps_->other_map());
  othermaps.push_back(embfluiddbcmaps->other_map());
  Teuchos::RCP<const Epetra_Map> xffluidothermap =
      Core::LinAlg::MultiMapExtractor::merge_maps(othermaps);

  xffluiddbcmaps_ = Teuchos::rcp(
      new Core::LinAlg::MapExtractor(*xffluiddofrowmap_, xffluiddbcmap, xffluidothermap));
}

/*----------------------------------------------------------------------*
 |  Set dirichlet- and velocity/pressure-map extractor      kruse 01/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidState::setup_map_extractors(
    const Teuchos::RCP<Core::FE::Discretization>& xfluiddiscret,
    const Teuchos::RCP<Core::FE::Discretization>& embfluiddiscret, const double& time)
{
  // create merged dirichlet map extractor
  XFluidState::setup_map_extractors(xfluiddiscret, time);
  xffluidsplitter_->setup(*xffluiddofrowmap_, embfluiddofrowmap_, XFluidState::xfluiddofrowmap_);

  FLD::UTILS::SetupFluidFluidVelPresSplit(*xfluiddiscret, Global::Problem::instance()->n_dim(),
      *embfluiddiscret, *xffluidvelpressplitter_, xffluiddofrowmap_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool FLD::XFluidFluidState::destroy()
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
  // dof_row_map() in Xfluidfluid currently returns a strong RCP
  if (xffluiddofrowmap_.strong_count() == 1)
    xffluiddofrowmap_ = Teuchos::null;
  else  // FOUR_C_THROW("could not destroy object: %i!=1 pointers",
        // xffluiddofrowmap_.strong_count());
    std::cout << "could not destroy xffluiddofrowmap_: number of pointers is "
              << xffluiddofrowmap_.strong_count() << "!=1";

  XFluidState::destroy();

  return true;
}

FOUR_C_NAMESPACE_CLOSE
