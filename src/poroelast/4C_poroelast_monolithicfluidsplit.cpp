/*----------------------------------------------------------------------*/
/*! \file

 \brief  monolithic fluid split poroelasticity algorithms

\level 2

 *------------------------------------------------------------------------------------------------*/

#include "4C_poroelast_monolithicfluidsplit.hpp"

#include "4C_adapter_fld_poro.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fsi_overlapprec_fsiamg.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_structure_aux.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

#define FLUIDSPLITAMG


PoroElast::MonolithicFluidSplit::MonolithicFluidSplit(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams,
    Teuchos::RCP<Core::LinAlg::MapExtractor> porosity_splitter)
    : MonolithicSplit(comm, timeparams, porosity_splitter)
{
  fggtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);
  fgitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  figtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  cfggtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  csggtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  cfgitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  csigtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);

  // Recovering of Lagrange multiplier happens on structure field
  lambda_ = Teuchos::rcp(new Epetra_Vector(*fluid_field()->Interface()->FSICondMap()));
}

void PoroElast::MonolithicFluidSplit::SetupSystem()
{
  {
    // create combined map
    std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;

    vecSpaces.push_back(structure_field()->dof_row_map());
#ifdef FLUIDSPLITAMG
    vecSpaces.push_back(fluid_field()->dof_row_map());
#else
    vecSpaces.push_back(fluid_field()->Interface()->OtherMap());
#endif

    if (vecSpaces[0]->NumGlobalElements() == 0) FOUR_C_THROW("No structure equation. Panic.");
    if (vecSpaces[1]->NumGlobalElements() == 0) FOUR_C_THROW("No fluid equation. Panic.");

    // full Poroelasticity-map
    fullmap_ = Core::LinAlg::MultiMapExtractor::MergeMaps(vecSpaces);
    // full Poroelasticity-blockmap
    blockrowdofmap_->setup(*fullmap_, vecSpaces);
  }

  // Switch fluid to interface split block matrix
  fluid_field()->use_block_matrix(true);

  setup_coupling_and_matrices();

  build_combined_dbc_map();

  SetupEquilibration();
}

void PoroElast::MonolithicFluidSplit::setup_rhs(bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("PoroElast::MonolithicFluidSplit::setup_rhs");

  // create full monolithic rhs vector
  rhs_ = Teuchos::rcp(new Epetra_Vector(*dof_row_map(), true));

  setup_vector(
      *rhs_, structure_field()->RHS(), fluid_field()->RHS(), fluid_field()->residual_scaling());

  if (firstcall and evaluateinterface_)
  {
    // add additional rhs-terms depending on the solution on the interface
    // of the previous time step

    // get time integration parameters of structure an fluid time integrators
    // to enable consistent time integration among the fields
    double stiparam = structure_field()->TimIntParam();
    double ftiparam = fluid_field()->TimIntParam();

    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockf = fluid_field()->BlockSystemMatrix();
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> k_sf = struct_fluid_coupling_block_matrix();
    if (k_sf == Teuchos::null) FOUR_C_THROW("expect coupling block matrix");

    Core::LinAlg::SparseMatrix& fig = blockf->Matrix(0, 1);
    Core::LinAlg::SparseMatrix& fgg = blockf->Matrix(1, 1);
    Core::LinAlg::SparseMatrix& kig = k_sf->Matrix(0, 1);
    Core::LinAlg::SparseMatrix& kgg = k_sf->Matrix(1, 1);

    Teuchos::RCP<Epetra_Vector> fveln = fluid_field()->extract_interface_veln();
    double timescale = fluid_field()->TimeScaling();
    double scale = fluid_field()->residual_scaling();

    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(fig.RowMap()));

    fig.Apply(*fveln, *rhs);
    rhs->Scale(timescale * Dt());

#ifdef FLUIDSPLITAMG
    rhs = fluid_field()->Interface()->InsertOtherVector(rhs);
#endif

    Extractor()->AddVector(*rhs, 1, *rhs_);  // add fluid contributions to 'f'

    rhs = Teuchos::rcp(new Epetra_Vector(fgg.RowMap()));

    fgg.Apply(*fveln, *rhs);
    rhs->Scale(scale * timescale * Dt());
    rhs->Scale(
        (1.0 - stiparam) / (1.0 - ftiparam));  // scale 'rhs' due to consistent time integration

    rhs = fluid_to_structure_at_interface(rhs);
    rhs = structure_field()->Interface()->InsertFSICondVector(rhs);

    Extractor()->AddVector(*rhs, 0, *rhs_);  // add structure contributions to 'f'

    rhs = Teuchos::rcp(new Epetra_Vector(kig.RowMap()));

    kig.Apply(*fveln, *rhs);
    rhs->Scale(timescale * Dt());

    rhs = structure_field()->Interface()->InsertOtherVector(rhs);

    Extractor()->AddVector(*rhs, 0, *rhs_);  // add structure contributions to 'f'

    rhs = Teuchos::rcp(new Epetra_Vector(kgg.RowMap()));

    kgg.Apply(*fveln, *rhs);
    rhs->Scale(timescale * Dt());

    rhs = structure_field()->Interface()->InsertFSICondVector(rhs);

    Extractor()->AddVector(*rhs, 0, *rhs_);  // add structure contributions to 'f'
  }

  // store interface force onto the structure to know it in the next time step as previous force
  // in order to recover the Lagrange multiplier
  fgcur_ = fluid_field()->Interface()->ExtractFSICondVector(fluid_field()->RHS());
}

void PoroElast::MonolithicFluidSplit::setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("PoroElast::MonolithicFluidSplit::setup_system_matrix");

  Teuchos::RCP<Core::LinAlg::SparseMatrix> s = structure_field()->system_matrix();
  if (s == Teuchos::null) FOUR_C_THROW("expect structure matrix");
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> f = fluid_field()->BlockSystemMatrix();
  if (f == Teuchos::null) FOUR_C_THROW("expect fluid block matrix");

  mat.Matrix(0, 1).Zero();
  mat.Matrix(1, 0).Zero();
#ifdef FLUIDSPLITAMG
  mat.Matrix(1, 1).Zero();
#endif

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

  s->UnComplete();

  /*----------------------------------------------------------------------*/

  if (evaluateinterface_)
  {
    double timescale = fluid_field()->TimeScaling();

    (*figtransform_)(f->FullRowMap(), f->FullColMap(), f->Matrix(0, 1), timescale,
        Core::Adapter::CouplingSlaveConverter(*icoupfs_), k_fs->Matrix(0, 1), true, true);

    (*csggtransform_)(f->FullRowMap(), f->FullColMap(), k_sf->Matrix(1, 1), timescale,
        Core::Adapter::CouplingSlaveConverter(*icoupfs_), *s, true, true);

    (*csigtransform_)(f->FullRowMap(), f->FullColMap(), k_sf->Matrix(0, 1), timescale,
        Core::Adapter::CouplingSlaveConverter(*icoupfs_), *s, true, true);
  }

  /*----------------------------------------------------------------------*/
  // pure fluid part
  // uncomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  // this just once...
#ifdef FLUIDSPLITAMG
  mat.Matrix(1, 1).Add(f->Matrix(0, 0), false, 1., 0.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> eye =
      Core::LinAlg::Eye(*fluid_field()->Interface()->FSICondMap());
  mat.Matrix(1, 1).Add(*eye, false, 1., 1.0);
#else
  f->Matrix(0, 0).UnComplete();
  mat.Assign(1, 1, View, f->Matrix(0, 0));
#endif

  // fluid coupling part
  mat.Matrix(1, 0).Add(k_fs->Matrix(0, 0), false, 1.0, 0.0);
  mat.Matrix(1, 0).Add(k_fs->Matrix(0, 1), false, 1.0, 1.0);

  // pure structure part
  mat.Assign(0, 0, Core::LinAlg::View, *s);

  // structure coupling part
  mat.Matrix(0, 1).Add(k_sf->Matrix(0, 0), false, 1.0, 0.0);
  mat.Matrix(0, 1).Add(k_sf->Matrix(1, 0), false, 1.0, 1.0);
  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();

  fgicur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(f->Matrix(1, 0)));
  fggcur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(f->Matrix(1, 1)));
  cgicur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(k_fs->Matrix(1, 0)));
  cggcur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(k_fs->Matrix(1, 1)));
}

void PoroElast::MonolithicFluidSplit::setup_vector(Epetra_Vector& f,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv, double fluidscale)
{
  // extract the inner and boundary dofs of all fields

  Teuchos::RCP<Epetra_Vector> fov = fluid_field()->Interface()->ExtractOtherVector(fv);
#ifdef FLUIDSPLITAMG
  fov = fluid_field()->Interface()->InsertOtherVector(fov);
#endif

  Extractor()->InsertVector(*sv, 0, f);

  Extractor()->InsertVector(*fov, 1, f);  // add fluid contributions to 'f'
}

void PoroElast::MonolithicFluidSplit::extract_field_vectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("PoroElast::MonolithicFluidSplit::extract_field_vectors");

  // process structure unknowns
  sx = Extractor()->ExtractVector(x, 0);

  // process fluid unknowns
  if (evaluateinterface_)
  {
    Teuchos::RCP<const Epetra_Vector> scx =
        structure_field()->Interface()->ExtractFSICondVector(sx);

    Teuchos::RCP<Epetra_Vector> fcx = structure_to_fluid_at_interface(scx);
    Teuchos::RCP<const Epetra_Vector> fox = Extractor()->ExtractVector(x, 1);
#ifdef FLUIDSPLITAMG
    fox = fluid_field()->Interface()->ExtractOtherVector(fox);
#endif

    {
      double timescale = fluid_field()->TimeScaling();
      fcx->Scale(timescale);
    }

    Teuchos::RCP<Epetra_Vector> f = fluid_field()->Interface()->InsertOtherVector(fox);
    fluid_field()->Interface()->InsertFSICondVector(fcx, f);

    auto zeros = Teuchos::rcp(new const Epetra_Vector(f->Map(), true));
    Core::LinAlg::apply_dirichlet_to_system(
        *f, *zeros, *(fluid_field()->GetDBCMapExtractor()->CondMap()));

    fx = f;

    // Store field vectors to know them later on as previous quantities
    Teuchos::RCP<Epetra_Vector> sox = structure_field()->Interface()->ExtractOtherVector(sx);
    if (solipre_ != Teuchos::null)
      ddiinc_->Update(1.0, *sox, -1.0, *solipre_, 0.0);  // compute current iteration increment
    else
      ddiinc_ = Teuchos::rcp(new Epetra_Vector(*sox));  // first iteration increment

    solipre_ = sox;  // store current step increment

    if (solgvelpre_ != Teuchos::null)
      duginc_->Update(1.0, *fcx, -1.0, *solgvelpre_, 0.0);  // compute current iteration increment
    else
      duginc_ = Teuchos::rcp(new Epetra_Vector(*fcx));  // first iteration increment

    solgvelpre_ = fcx;  // store current step increment

    if (solivelpre_ != Teuchos::null)
      duiinc_->Update(1.0, *fox, -1.0, *solivelpre_, 0.0);  // compute current iteration increment
    else
      duiinc_ = Teuchos::rcp(new Epetra_Vector(*fox));  // first iteration increment

    solivelpre_ = fox;  // store current step increment
  }
  else
    fx = Extractor()->ExtractVector(x, 1);
}

void PoroElast::MonolithicFluidSplit::recover_lagrange_multiplier_after_time_step()
{
  TEUCHOS_FUNC_TIME_MONITOR("PoroElast::MonolithicFluidSplit::recover_lagrange_multiplier");

  if (evaluateinterface_)
  {
    // get time integration parameter of structural time integrator
    // to enable consistent time integration among the fields
    double ftiparam = fluid_field()->TimIntParam();
    double timescale = fluid_field()->TimeScaling();

    // store the product F_{\GammaI} \Delta u_I^{n+1} in here
    Teuchos::RCP<Epetra_Vector> fgiddi =
        Core::LinAlg::CreateVector(*fluid_field()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    fgicur_->Multiply(false, *duiinc_, *fgiddi);

    // store the product C_{\GammaI} \Delta d_I^{n+1} in here
    Teuchos::RCP<Epetra_Vector> sgiddi =
        Core::LinAlg::CreateVector(*fluid_field()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    cgicur_->Multiply(false, *ddiinc_, *sgiddi);

    // store the product F_{\Gamma\Gamma} \Delta u_\Gamma^{n+1} in here
    Teuchos::RCP<Epetra_Vector> sggddg =
        Core::LinAlg::CreateVector(*fluid_field()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    fggcur_->Multiply(false, *duginc_, *sggddg);

    // store the prodcut C_{\Gamma\Gamma} \Delta d_\Gamma^{n+1} in here
    Teuchos::RCP<Epetra_Vector> cggddg =
        Core::LinAlg::CreateVector(*fluid_field()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    cggcur_->Multiply(false, *duginc_, *cggddg);
    cggddg->Scale(1.0 / timescale);

    // Update the Lagrange multiplier:
    /* \lambda^{n+1} =  1/b * [ - a*\lambda^n - f_\Gamma^S
     *                          - F_{\Gamma I} \Delta u_I
     *                          - C_{\Gamma I} \Delta d_I
     *                          - F_{\Gamma\Gamma} \Delta u_\Gamma]
     *                          - C_{\Gamma\Gamma} * Delta t / 2 * \Delta u_\Gamma]
     */
    lambda_->Update(1.0, *fgcur_, -ftiparam);
    lambda_->Update(-1.0, *fgiddi, -1.0, *sggddg, 1.0);
    lambda_->Update(-1.0, *fgiddi, -1.0, *cggddg, 1.0);
    lambda_->Scale(1 / (1.0 - ftiparam));  // entire Lagrange multiplier is divided by (1.-ftiparam)
  }
}

FOUR_C_NAMESPACE_CLOSE
