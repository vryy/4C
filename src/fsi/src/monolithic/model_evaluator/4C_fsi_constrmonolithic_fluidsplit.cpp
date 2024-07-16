/*----------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problem with constraints

\level 2

*/
/*----------------------------------------------------------------------*/
#include "4C_fsi_constrmonolithic_fluidsplit.hpp"

#include "4C_adapter_ale_fsi.hpp"
#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_adapter_str_structure.hpp"
#include "4C_ale_utils_mapextractor.hpp"
#include "4C_constraint_manager.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_structure_aux.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::ConstrMonolithicFluidSplit::ConstrMonolithicFluidSplit(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : ConstrMonolithic(comm, timeparams)
{
  // ---------------------------------------------------------------------------
  // FSI specific check of Dirichlet boundary conditions
  // ---------------------------------------------------------------------------
  // Create intersection of slave DOFs that hold a Dirichlet boundary condition
  // and are located at the FSI interface
  std::vector<Teuchos::RCP<const Epetra_Map>> intersectionmaps;
  intersectionmaps.push_back(fluid_field()->get_dbc_map_extractor()->cond_map());
  intersectionmaps.push_back(fluid_field()->interface()->fsi_cond_map());
  Teuchos::RCP<Epetra_Map> intersectionmap =
      Core::LinAlg::MultiMapExtractor::intersect_maps(intersectionmaps);

  // Check whether the intersection is empty
  if (intersectionmap->NumGlobalElements() != 0)
  {
    // It is not allowed, that slave DOFs at the interface hold a Dirichlet
    // boundary condition. Thus --> ToDO: Error message

    // We do not have to care whether ALE interface DOFs carry DBCs in the
    // input file since they do not occur in the monolithic system and, hence,
    // do not cause a conflict.

    std::stringstream errormsg;
    errormsg << "  "
                "+---------------------------------------------------------------------------------"
                "------------+"
             << std::endl
             << "  |                DIRICHLET BOUNDARY CONDITIONS ON SLAVE SIDE OF FSI INTERFACE   "
                "              |"
             << std::endl
             << "  "
                "+---------------------------------------------------------------------------------"
                "------------+"
             << std::endl
             << "  | NOTE: The slave side of the interface is not allowed to carry Dirichlet "
                "boundary conditions.|"
             << std::endl
             << "  |                                                                               "
                "              |"
             << std::endl
             << "  | This is a fluid split scheme. Hence, master and slave field are chosen as "
                "follows:          |"
             << std::endl
             << "  |     MASTER  = STRUCTURE                                                       "
                "              |"
             << std::endl
             << "  |     SLAVE   = FLUID                                                           "
                "              |"
             << std::endl
             << "  |                                                                               "
                "              |"
             << std::endl
             << "  | Dirichlet boundary conditions were detected on slave interface degrees of "
                "freedom. Please   |"
             << std::endl
             << "  | remove Dirichlet boundary conditions from the slave side of the FSI "
                "interface.              |"
             << std::endl
             << "  | Only the master side of the FSI interface is allowed to carry Dirichlet "
                "boundary conditions.|"
             << std::endl
             << "  "
                "+---------------------------------------------------------------------------------"
                "------------+"
             << std::endl;

    std::cout << errormsg.str();
  }
  // ---------------------------------------------------------------------------

  scon_t_ =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*conman_->get_constraint_map(), 81, false, true));

  fggtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);
  fgitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  figtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  fmiitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  fmgitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);
  fmigtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  fmggtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);
  aigtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicFluidSplit::setup_system()
{
  general_setup();

  // create combined map
  create_combined_dof_row_map();

  fluid_field()->use_block_matrix(true);

  // build ale system matrix in splitted system
  ale_field()->create_system_matrix(ale_field()->interface());

  create_system_matrix(false);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicFluidSplit::create_combined_dof_row_map()
{
  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;
  vecSpaces.push_back(structure_field()->dof_row_map());
  vecSpaces.push_back(fluid_field()->dof_row_map());
  vecSpaces.push_back(ale_field()->interface()->other_map());
  vecSpaces.push_back(conman_->get_constraint_map());

  if (vecSpaces[0]->NumGlobalElements() == 0)
    FOUR_C_THROW("No inner structural equations. Splitting not possible. Panic.");

  set_dof_row_maps(vecSpaces);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicFluidSplit::setup_rhs_residual(Epetra_Vector& f)
{
  const double scale = fluid_field()->residual_scaling();

  setup_vector(f, structure_field()->rhs(), fluid_field()->rhs(), ale_field()->rhs(),
      conman_->get_error(), scale);

  // add additional ale residual
  extractor().add_vector(*aleresidual_, 2, f);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicFluidSplit::setup_rhs_lambda(Epetra_Vector& f)
{
  // ToDo: We still need to implement this.

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicFluidSplit::setup_rhs_firstiter(Epetra_Vector& f)
{
  Teuchos::RCP<const Core::LinAlg::BlockSparseMatrixBase> blockf =
      fluid_field()->block_system_matrix();

  const Core::LinAlg::SparseMatrix& fig = blockf->matrix(0, 1);
  const Core::LinAlg::SparseMatrix& fgg = blockf->matrix(1, 1);

  Teuchos::RCP<const Epetra_Vector> fveln = fluid_field()->extract_interface_veln();
  const double timescale = fluid_field()->time_scaling();
  const double scale = fluid_field()->residual_scaling();

  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(fig.row_map()));

  fig.Apply(*fveln, *rhs);
  rhs->Scale(timescale * dt());

  rhs = fluid_field()->interface()->insert_other_vector(rhs);
  extractor().add_vector(*rhs, 1, f);

  rhs = Teuchos::rcp(new Epetra_Vector(fgg.row_map()));

  fgg.Apply(*fveln, *rhs);
  rhs->Scale(scale * timescale * dt());

  rhs = fluid_to_struct(rhs);
  rhs = structure_field()->interface()->insert_fsi_cond_vector(rhs);
  extractor().add_vector(*rhs, 0, f);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicFluidSplit::setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::setup_system_matrix");

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W

  const Core::Adapter::Coupling& coupsf = structure_fluid_coupling();
  const Core::Adapter::Coupling& coupsa = structure_ale_coupling();
  const Core::Adapter::Coupling& coupfa = fluid_ale_coupling();

  /*----------------------------------------------------------------------*/

  double scale = fluid_field()->residual_scaling();
  double timescale = fluid_field()->time_scaling();

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockf = fluid_field()->block_system_matrix();


  Core::LinAlg::SparseMatrix& fii = blockf->matrix(0, 0);
  Core::LinAlg::SparseMatrix& fig = blockf->matrix(0, 1);
  Core::LinAlg::SparseMatrix& fgi = blockf->matrix(1, 0);
  Core::LinAlg::SparseMatrix& fgg = blockf->matrix(1, 1);
  /*----------------------------------------------------------------------*/

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> a = ale_field()->block_system_matrix();

  if (a == Teuchos::null) FOUR_C_THROW("expect ale block matrix");

  Core::LinAlg::SparseMatrix& aii = a->matrix(0, 0);
  Core::LinAlg::SparseMatrix& aig = a->matrix(0, 1);
  /*----------------------------------------------------------------------*/
  // structure part

  Teuchos::RCP<Core::LinAlg::SparseMatrix> s = structure_field()->system_matrix();
  // const std::string fname = "cfsstructmatrix.mtl";
  // Core::LinAlg::PrintMatrixInMatlabFormat(fname,*(s->EpetraMatrix()));


  /*----------------------------------------------------------------------*/
  // structure constraint part
  Teuchos::RCP<Core::LinAlg::SparseOperator> tmp = conman_->get_constr_matrix();
  Core::LinAlg::SparseMatrix scon = *(Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(tmp));

  scon.complete(*conman_->get_constraint_map(), s->range_map());

  mat.assign(0, 3, Core::LinAlg::View, scon);

  /*----------------------------------------------------------------------*/
  // fluid part

  //    // uncomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  //   this just once...
  s->un_complete();


  (*fggtransform_)(fgg, scale * timescale, Core::Adapter::CouplingSlaveConverter(coupsf),
      Core::Adapter::CouplingSlaveConverter(coupsf), *s, true, true);

  mat.assign(0, 0, Core::LinAlg::View, *s);

  (*fgitransform_)(fgi, scale, Core::Adapter::CouplingSlaveConverter(coupsf), mat.matrix(0, 1));

  (*figtransform_)(blockf->full_row_map(), blockf->full_col_map(), fig, timescale,
      Core::Adapter::CouplingSlaveConverter(coupsf), mat.matrix(1, 0));

  mat.matrix(1, 1).add(fii, false, 1., 0.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> eye =
      Core::LinAlg::CreateIdentityMatrix(*fluid_field()->interface()->fsi_cond_map());
  mat.matrix(1, 1).add(*eye, false, 1., 1.0);

  (*aigtransform_)(a->full_row_map(), a->full_col_map(), aig, 1.,
      Core::Adapter::CouplingSlaveConverter(coupsa), mat.matrix(2, 0));
  mat.assign(2, 2, Core::LinAlg::View, aii);

  /*----------------------------------------------------------------------*/
  // add optional fluid linearization with respect to mesh motion block

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> mmm = fluid_field()->shape_derivatives();
  if (mmm != Teuchos::null)
  {
    Core::LinAlg::SparseMatrix& fmii = mmm->matrix(0, 0);
    Core::LinAlg::SparseMatrix& fmgi = mmm->matrix(1, 0);

    Core::LinAlg::SparseMatrix& fmig = mmm->matrix(0, 1);
    Core::LinAlg::SparseMatrix& fmgg = mmm->matrix(1, 1);

    // reuse transform objects to add shape derivative matrices to structural blocks

    (*figtransform_)(blockf->full_row_map(), blockf->full_col_map(), fmig, 1.,
        Core::Adapter::CouplingSlaveConverter(coupsf), mat.matrix(1, 0), false, true);

    (*fggtransform_)(fmgg, scale, Core::Adapter::CouplingSlaveConverter(coupsf),
        Core::Adapter::CouplingSlaveConverter(coupsf), mat.matrix(0, 0), false, true);

    // We cannot copy the pressure value. It is not used anyway. So no exact
    // match here.
    (*fmiitransform_)(mmm->full_row_map(), mmm->full_col_map(), fmii, 1.,
        Core::Adapter::CouplingMasterConverter(coupfa), mat.matrix(1, 2), false);

    (*fmgitransform_)(fmgi, scale, Core::Adapter::CouplingSlaveConverter(coupsf),
        Core::Adapter::CouplingMasterConverter(coupfa), mat.matrix(0, 2), false, false);
  }
  /*----------------------------------------------------------------------*/
  // constraint part -> structure

  scon = *(Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(tmp));

  scon_t_->add(scon, true, 1.0, 0.0);
  scon_t_->complete(scon.range_map(), scon.domain_map());
  mat.assign(3, 0, Core::LinAlg::View, *scon_t_);

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.complete();

  // Finally, take care of Dirichlet boundary conditions
  mat.apply_dirichlet(*(dbcmaps_->cond_map()), true);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicFluidSplit::initial_guess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::initial_guess");

  Teuchos::RCP<Epetra_Vector> ConstraintInitialGuess =
      Teuchos::rcp(new Epetra_Vector(*(conman_->get_constraint_map()), true));

  setup_vector(*ig, structure_field()->initial_guess(), fluid_field()->initial_guess(),
      ale_field()->initial_guess(), ConstraintInitialGuess, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicFluidSplit::setup_vector(Epetra_Vector& f,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv,
    Teuchos::RCP<const Epetra_Vector> av, Teuchos::RCP<const Epetra_Vector> cv, double fluidscale)
{
  // extract the inner and boundary dofs of all three fields

  Teuchos::RCP<Epetra_Vector> fov = fluid_field()->interface()->extract_other_vector(fv);
  fov = fluid_field()->interface()->insert_other_vector(fov);
  Teuchos::RCP<Epetra_Vector> aov = ale_field()->interface()->extract_other_vector(av);

  if (fabs(fluidscale) >= 1.0E-10)
  {
    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> fcv = fluid_field()->interface()->extract_fsi_cond_vector(fv);
    Teuchos::RCP<Epetra_Vector> modsv =
        structure_field()->interface()->insert_fsi_cond_vector(fluid_to_struct(fcv));
    modsv->Update(1.0, *sv, fluidscale);

    extractor().insert_vector(*modsv, 0, f);
  }
  else
  {
    extractor().insert_vector(*sv, 0, f);
  }

  extractor().insert_vector(*fov, 1, f);
  extractor().insert_vector(*aov, 2, f);
  Epetra_Vector modcv = *cv;
  modcv.Scale(-1.0);
  extractor().insert_vector(modcv, 3, f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicFluidSplit::extract_field_vectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
    Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::ConstrMonolithicFluidSplit::extract_field_vectors");

  // We have overlap at the interface. Thus we need the interface part of the
  // structure vector and append it to the fluid and ale vector. (With the
  // right translation.)

  sx = extractor().extract_vector(x, 0);
  Teuchos::RCP<const Epetra_Vector> scx =
      structure_field()->interface()->extract_fsi_cond_vector(sx);

  // process fluid unknowns

  Teuchos::RCP<const Epetra_Vector> fox = extractor().extract_vector(x, 1);
  fox = fluid_field()->interface()->extract_other_vector(fox);
  Teuchos::RCP<Epetra_Vector> fcx = struct_to_fluid(scx);

  fluid_field()->displacement_to_velocity(fcx);

  Teuchos::RCP<Epetra_Vector> f = fluid_field()->interface()->insert_other_vector(fox);
  fluid_field()->interface()->insert_fsi_cond_vector(fcx, f);
  fx = f;

  // process ale unknowns

  Teuchos::RCP<const Epetra_Vector> aox = extractor().extract_vector(x, 2);
  Teuchos::RCP<Epetra_Vector> acx = struct_to_ale(scx);
  Teuchos::RCP<Epetra_Vector> a = ale_field()->interface()->insert_other_vector(aox);
  ale_field()->interface()->insert_vector(acx, 1, a);

  ax = a;
}

FOUR_C_NAMESPACE_CLOSE
