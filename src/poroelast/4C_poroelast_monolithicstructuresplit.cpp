/*----------------------------------------------------------------------*/
/*! \file

 \brief  monolithic structure split poroelasticity algorithm

\level 2

 *------------------------------------------------------------------------------------------------*/

#include "4C_poroelast_monolithicstructuresplit.hpp"

#include "4C_adapter_fld_poro.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_structure_aux.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


PoroElast::MonolithicStructureSplit::MonolithicStructureSplit(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams,
    Teuchos::RCP<Core::LinAlg::MapExtractor> porosity_splitter)
    : MonolithicSplit(comm, timeparams, porosity_splitter)
{
  sggtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);
  sgitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  sigtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  csggtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  cfggtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  csgitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  cfigtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);

  // Recovering of Lagrange multiplier happens on structure field
  lambda_ = Teuchos::rcp(new Epetra_Vector(*structure_field()->Interface()->FSICondMap()));
}

void PoroElast::MonolithicStructureSplit::SetupSystem()
{
  {
    // create combined map
    std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;

    vecSpaces.push_back(structure_field()->Interface()->OtherMap());
    vecSpaces.push_back(fluid_field()->dof_row_map());

    if (vecSpaces[0]->NumGlobalElements() == 0) FOUR_C_THROW("No structure equation. Panic.");
    if (vecSpaces[1]->NumGlobalElements() == 0) FOUR_C_THROW("No fluid equation. Panic.");

    // full Poroelasticity-map
    fullmap_ = Core::LinAlg::MultiMapExtractor::MergeMaps(vecSpaces);
    // full Poroelasticity-blockmap
    blockrowdofmap_->setup(*fullmap_, vecSpaces);
  }

  // Use splitted structure matrix
  structure_field()->use_block_matrix();

  setup_coupling_and_matrices();

  build_combined_dbc_map();

  SetupEquilibration();
}

void PoroElast::MonolithicStructureSplit::setup_rhs(bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("PoroElast::MonolithicStructureSplit::setup_rhs");

  // create full monolithic rhs vector
  rhs_ = Teuchos::rcp(new Epetra_Vector(*dof_row_map(), true));

  setup_vector(
      *rhs_, structure_field()->RHS(), fluid_field()->RHS(), fluid_field()->residual_scaling());

  // store interface force onto the structure to know it in the next time step as previous force
  // in order to recover the Lagrange multiplier
  // fgpre_ = fgcur_;
  fgcur_ = structure_field()->Interface()->ExtractFSICondVector(structure_field()->RHS());
}

void PoroElast::MonolithicStructureSplit::setup_system_matrix(
    Core::LinAlg::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("PoroElast::MonolithicStructureSplit::setup_system_matrix");

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> s = structure_field()->block_system_matrix();
  if (s == Teuchos::null) FOUR_C_THROW("expect structure block matrix");
  Teuchos::RCP<Core::LinAlg::SparseMatrix> f = fluid_field()->SystemMatrix();
  if (f == Teuchos::null) FOUR_C_THROW("expect fluid matrix");

  // just to play it safe ...
  mat.Zero();

  /*----------------------------------------------------------------------*/

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  /*----------------------------------------------------------------------*/
  // structural part k_sf (3nxn)
  // build mechanical-fluid block

  // create empty matrix
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> k_sf = struct_fluid_coupling_block_matrix();
  if (k_sf == Teuchos::null) FOUR_C_THROW("expect coupling block matrix");

  // call the element and calculate the matrix block
  apply_str_coupl_matrix(k_sf);

  /*----------------------------------------------------------------------*/
  // fluid part k_fs ( (3n+1)x3n )
  // build fluid-mechanical block

  // create empty matrix
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> k_fs = fluid_struct_coupling_block_matrix();
  if (k_fs == Teuchos::null) FOUR_C_THROW("expect coupling block matrix");

  // call the element and calculate the matrix block
  apply_fluid_coupl_matrix(k_fs);

  /*----------------------------------------------------------------------*/

  k_fs->Complete();
  k_sf->Complete();

  f->UnComplete();

  /*----------------------------------------------------------------------*/
  if (evaluateinterface_)
  {
    double scale = fluid_field()->residual_scaling();
    double timescale = fluid_field()->TimeScaling();

    // get time integration parameters of structure an fluid time integrators
    // to enable consistent time integration among the fields
    double stiparam = structure_field()->TimIntParam();
    double ftiparam = fluid_field()->TimIntParam();

    (*sigtransform_)(s->FullRowMap(), s->FullColMap(), s->Matrix(0, 1), 1. / timescale,
        Core::Adapter::CouplingMasterConverter(*icoupfs_),
        k_sf->Matrix(0, 1),  // k_sf->Matrix(0,1),mat.Matrix(0,1)
        true, true);

    (*sggtransform_)(s->Matrix(1, 1), (1.0 - ftiparam) / ((1.0 - stiparam) * scale * timescale),
        Core::Adapter::CouplingMasterConverter(*icoupfs_),
        Core::Adapter::CouplingMasterConverter(*icoupfs_), *f, true, true);

    (*sgitransform_)(s->Matrix(1, 0), (1.0 - ftiparam) / ((1.0 - stiparam) * scale),
        Core::Adapter::CouplingMasterConverter(*icoupfs_),
        k_fs->Matrix(1, 0),  // k_fs->Matrix(1,0), mat.Matrix(1,0)
        true);

    (*cfggtransform_)(s->FullRowMap(),  // k_fs->FullRowMap(),
        s->FullColMap(),                // k_fs->FullColMap(),
        k_fs->Matrix(1, 1), 1. / timescale, Core::Adapter::CouplingMasterConverter(*icoupfs_), *f,
        true, true);

    (*cfigtransform_)(s->FullRowMap(),  // k_fs->FullRowMap(),
        s->FullColMap(),                // k_fs->FullColMap(),
        k_fs->Matrix(0, 1), 1. / timescale, Core::Adapter::CouplingMasterConverter(*icoupfs_), *f,
        true, true);

    (*csggtransform_)(k_sf->Matrix(1, 1), (1.0 - ftiparam) / ((1.0 - stiparam) * scale),
        Core::Adapter::CouplingMasterConverter(*icoupfs_), *f, true);

    (*csgitransform_)(k_sf->Matrix(1, 0), (1.0 - ftiparam) / ((1.0 - stiparam) * scale),
        Core::Adapter::CouplingMasterConverter(*icoupfs_), *f, true);
  }

  /*----------------------------------------------------------------------*/
  // pure structural part
  mat.Matrix(0, 0).Add(s->Matrix(0, 0), false, 1., 0.0);

  // structure coupling part
  mat.Matrix(0, 1).Add(k_sf->Matrix(0, 0), false, 1.0, 0.0);
  mat.Matrix(0, 1).Add(k_sf->Matrix(0, 1), false, 1.0, 1.0);

  // pure fluid part
  mat.Assign(1, 1, Core::LinAlg::View, *f);

  // fluid coupling part
  mat.Matrix(1, 0).Add(k_fs->Matrix(0, 0), false, 1.0, 0.0);
  mat.Matrix(1, 0).Add(k_fs->Matrix(1, 0), false, 1.0, 1.0);
  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();

  sgicur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(s->Matrix(1, 0)));
  sggcur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(s->Matrix(1, 1)));
  cgicur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(k_sf->Matrix(1, 0)));
  cggcur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(k_sf->Matrix(1, 1)));
}

void PoroElast::MonolithicStructureSplit::setup_vector(Epetra_Vector& f,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv, double fluidscale)
{
  // extract the inner and boundary dofs of all three fields

  Teuchos::RCP<Epetra_Vector> sov = structure_field()->Interface()->ExtractOtherVector(sv);

  if (fluidscale != 0.0)
  {
    // get time integration parameters of structure an fluid time integrators
    // to enable consistent time integration among the fields
    double stiparam = structure_field()->TimIntParam();
    double ftiparam = fluid_field()->TimIntParam();

    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> scv = structure_field()->Interface()->ExtractFSICondVector(sv);
    Teuchos::RCP<Epetra_Vector> modfv =
        fluid_field()->Interface()->InsertFSICondVector(structure_to_fluid_at_interface(scv));
    modfv->Update(1.0, *fv, (1.0 - ftiparam) / ((1.0 - stiparam) * fluidscale));

    // add contribution of Lagrange multiplier from previous time step
    if (lambda_ != Teuchos::null)
      modfv->Update(-ftiparam + stiparam * (1.0 - ftiparam) / (1.0 - stiparam),
          *structure_to_fluid_at_interface(lambda_), 1.0);

    Extractor()->InsertVector(*modfv, 1, f);
  }
  else
  {
    Extractor()->InsertVector(*fv, 1, f);
  }

  Extractor()->InsertVector(*sov, 0, f);
}

void PoroElast::MonolithicStructureSplit::extract_field_vectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("PoroElast::MonolithicStructureSplit::extract_field_vectors");

  // process fluid unknowns of the second field
  fx = Extractor()->ExtractVector(x, 1);

  // process structure unknowns
  if (evaluateinterface_)
  {
    Teuchos::RCP<Epetra_Vector> fcx = fluid_field()->Interface()->ExtractFSICondVector(fx);

    {
      double timescale = 1. / fluid_field()->TimeScaling();
      fcx->Scale(timescale);
    }

    Teuchos::RCP<Epetra_Vector> scx = fluid_to_structure_at_interface(fcx);
    Teuchos::RCP<const Epetra_Vector> sox = Extractor()->ExtractVector(x, 0);

    Teuchos::RCP<Epetra_Vector> s = structure_field()->Interface()->InsertOtherVector(sox);
    structure_field()->Interface()->InsertFSICondVector(scx, s);

    auto zeros = Teuchos::rcp(new const Epetra_Vector(s->Map(), true));
    Core::LinAlg::apply_dirichlet_to_system(
        *s, *zeros, *(structure_field()->GetDBCMapExtractor()->CondMap()));

    sx = s;

    Teuchos::RCP<Epetra_Vector> fox = fluid_field()->Interface()->ExtractOtherVector(fx);

    // Store field vectors to know them later on as previous quantities
    if (solipre_ != Teuchos::null)
      ddiinc_->Update(1.0, *sox, -1.0, *solipre_, 0.0);  // compute current iteration increment
    else
      ddiinc_ = Teuchos::rcp(new Epetra_Vector(*sox));  // first iteration increment

    solipre_ = sox;  // store current step increment

    if (solgpre_ != Teuchos::null)
      ddginc_->Update(1.0, *scx, -1.0, *solgpre_, 0.0);  // compute current iteration increment
    else
      ddginc_ = Teuchos::rcp(new Epetra_Vector(*scx));  // first iteration increment

    solgpre_ = scx;  // store current step increment

    if (solivelpre_ != Teuchos::null)
      duiinc_->Update(1.0, *fox, -1.0, *solivelpre_, 0.0);  // compute current iteration increment
    else
      duiinc_ = Teuchos::rcp(new Epetra_Vector(*fox));  // first iteration increment

    solivelpre_ = fox;  // store current step increment
  }
  else
    sx = Extractor()->ExtractVector(x, 0);
}

void PoroElast::MonolithicStructureSplit::recover_lagrange_multiplier_after_time_step()
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PoroElast::MonolithicStructureSplit::recover_lagrange_multiplier_after_time_step");
  if (evaluateinterface_)
  {
    // get time integration parameter of structural time integrator
    // to enable consistent time integration among the fields
    double stiparam = structure_field()->TimIntParam();
    double timescale = fluid_field()->TimeScaling();

    // store the product S_{\GammaI} \Delta d_I^{n+1} in here
    Teuchos::RCP<Epetra_Vector> sgiddi =
        Core::LinAlg::CreateVector(*structure_field()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    sgicur_->Multiply(false, *ddiinc_, *sgiddi);

    // store the product C_{\GammaI} \Delta u_I^{n+1} in here
    Teuchos::RCP<Epetra_Vector> fgiddi =
        Core::LinAlg::CreateVector(*structure_field()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    cgicur_->Multiply(false, *duiinc_, *fgiddi);

    // store the product S_{\Gamma\Gamma} \Delta d_\Gamma^{n+1} in here
    Teuchos::RCP<Epetra_Vector> sggddg =
        Core::LinAlg::CreateVector(*structure_field()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    sggcur_->Multiply(false, *ddginc_, *sggddg);

    // store the prodcut C_{\Gamma\Gamma} \Delta u_\Gamma^{n+1} in here
    Teuchos::RCP<Epetra_Vector> cggddg =
        Core::LinAlg::CreateVector(*structure_field()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    cggcur_->Multiply(false, *ddginc_, *cggddg);
    cggddg->Scale(timescale);

    // Update the Lagrange multiplier:
    /* \lambda^{n+1} =  1/b * [ - a*\lambda^n - f_\Gamma^S
     *                          - S_{\Gamma I} \Delta d_I
     *                          - C_{\Gamma I} \Delta u_I
     *                          - S_{\Gamma\Gamma} \Delta d_\Gamma]
     *                          - C_{\Gamma\Gamma} * 2 / \Delta t * \Delta d_\Gamma]
     */
    // lambda_->Update(1.0, *fgpre_, -stiparam);
    lambda_->Update(1.0, *fgcur_, -stiparam);
    lambda_->Update(-1.0, *sgiddi, -1.0, *sggddg, 1.0);
    lambda_->Update(-1.0, *fgiddi, -1.0, *cggddg, 1.0);
    lambda_->Scale(1 / (1.0 - stiparam));  // entire Lagrange multiplier is divided by (1.-stiparam)
  }
}

FOUR_C_NAMESPACE_CLOSE
