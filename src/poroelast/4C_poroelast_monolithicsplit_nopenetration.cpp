/*----------------------------------------------------------------------*/
/*! \file

 \brief porous medium algorithm with matrix split for condensation of
      no-penetration constraint

\level 2

 *----------------------------------------------------------------------*/

#include "4C_poroelast_monolithicsplit_nopenetration.hpp"

#include "4C_adapter_coupling_nonlin_mortar.hpp"
#include "4C_adapter_fld_poro.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_contact_interface.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_assemblestrategy.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_structure_aux.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


PoroElast::MonolithicSplitNoPenetration::MonolithicSplitNoPenetration(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams,
    Teuchos::RCP<Core::LinAlg::MapExtractor> porosity_splitter)
    : MonolithicSplit(comm, timeparams, porosity_splitter), normrhs_nopenetration_(-1.0)
{
  // Initialize Transformation Objects
  k_d_transform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  k_inv_d_transform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);

  k_d_lin_transform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);

  // Recovering of Lagrange multiplier happens on fluid field
  lambda_ = Teuchos::rcp(new Epetra_Vector(*structure_field()->interface()->fsi_cond_map()));
  lambdanp_ = Teuchos::rcp(new Epetra_Vector(*structure_field()->interface()->fsi_cond_map()));

  k_dn_ = Teuchos::null;

  mortar_adapter_ = Teuchos::rcp(new Adapter::CouplingNonLinMortar(
      Global::Problem::instance()->n_dim(), Global::Problem::instance()->mortar_coupling_params(),
      Global::Problem::instance()->contact_dynamic_params(),
      Global::Problem::instance()->spatial_approximation_type()));
}

void PoroElast::MonolithicSplitNoPenetration::setup_system()
{
  {
    const int ndim = Global::Problem::instance()->n_dim();
    std::vector<int> coupleddof(ndim + 1, 1);
    coupleddof[ndim] = 0;

    mortar_adapter_->setup(structure_field()->discretization(), fluid_field()->discretization(),
        coupleddof, "FSICoupling");
  }

  // use full maps of both fields. Only Lagrange multipliers are condensed
  {
    // create combined map
    std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;

    vecSpaces.push_back(structure_field()->dof_row_map());
    vecSpaces.push_back(fluid_field()->dof_row_map());

    if (vecSpaces[0]->NumGlobalElements() == 0) FOUR_C_THROW("No structure equation. Panic.");
    if (vecSpaces[1]->NumGlobalElements() == 0) FOUR_C_THROW("No fluid equation. Panic.");

    // full Poroelasticity-map
    fullmap_ = Core::LinAlg::MultiMapExtractor::merge_maps(vecSpaces);
    // full Poroelasticity-blockmap
    blockrowdofmap_->setup(*fullmap_, vecSpaces);
  }

  // Switch fluid to interface split block matrix
  fluid_field()->use_block_matrix(true);

  // setup coupling objects, system and coupling matrices
  setup_coupling_and_matrices();

  // build map of dofs subjected to a DBC of whole problem
  build_combined_dbc_map();

  setup_equilibration();
}

void PoroElast::MonolithicSplitNoPenetration::setup_rhs(bool firstcall)
{
  // only Lagrange multipliers are condensed -> use unchanged maps from single fields
  TEUCHOS_FUNC_TIME_MONITOR("PoroElast::MonolithicSplitNoPenetration::setup_rhs");

  // create full monolithic rhs vector
  if (rhs_ == Teuchos::null) rhs_ = Teuchos::rcp(new Epetra_Vector(*dof_row_map(), true));

  setup_vector(*rhs_, structure_field()->rhs(), fluid_field()->rhs());
}

void PoroElast::MonolithicSplitNoPenetration::setup_vector(
    Epetra_Vector& f, Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv)
{
  // extract dofs of the two fields
  // and put the structural/fluid field vector into the global vector f
  // noticing the block number

  extractor()->insert_vector(*sv, 0, f);

  Teuchos::RCP<Epetra_Vector> fov = fluid_field()->interface()->extract_other_vector(fv);
  Teuchos::RCP<Epetra_Vector> fcv = fluid_field()->interface()->extract_fsi_cond_vector(fv);

  Teuchos::RCP<Epetra_Vector> Dlam =
      Teuchos::rcp(new Epetra_Vector(*fluid_field()->interface()->fsi_cond_map(), true));
  Teuchos::RCP<Epetra_Vector> couprhs =
      Teuchos::rcp(new Epetra_Vector(*fluid_field()->interface()->fsi_cond_map(), true));
  if (k_dn_ != Teuchos::null)
  {
    double stiparam = structure_field()->tim_int_param();

    k_dn_->multiply(false, *lambda_, *Dlam);  // D(n)*lambda(n)

    Dlam->Scale(stiparam);  //*(1-b)
  }
  Dlam->Update(-1.0, *fcv, 1.0);
  k_lambdainv_d_->multiply(false, *Dlam, *couprhs);

  couprhs->Update(1.0, *nopenetration_rhs_, 1.0);

  // std::cout << "nopenetration_rhs_: " << *nopenetration_rhs_ << std::endl;

  Teuchos::RCP<Epetra_Vector> fullcouprhs =
      Teuchos::rcp(new Epetra_Vector(*fluid_field()->dof_row_map(), true));
  Core::LinAlg::Export(*couprhs, *fullcouprhs);
  extractor()->insert_vector(*fullcouprhs, 1, f);

  Teuchos::RCP<Epetra_Vector> fullfov =
      Teuchos::rcp(new Epetra_Vector(*fluid_field()->dof_row_map(), true));
  Core::LinAlg::Export(*fov, *fullfov);
  extractor()->add_vector(*fullfov, 1, f, 1.0);

  rhs_fgcur_ = fcv;  // Store interface rhs for recovering of lagrange multiplier
}

void PoroElast::MonolithicSplitNoPenetration::recover_lagrange_multiplier_after_newton_step(
    Teuchos::RCP<const Epetra_Vector> x)
{
  // call base class
  Monolithic::recover_lagrange_multiplier_after_newton_step(x);


  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  extract_field_vectors(x, sx, fx);

  Teuchos::RCP<Epetra_Vector> sox = structure_field()->interface()->extract_other_vector(sx);
  Teuchos::RCP<Epetra_Vector> scx = structure_field()->interface()->extract_fsi_cond_vector(sx);
  Teuchos::RCP<Epetra_Vector> fox = fluid_field()->interface()->extract_other_vector(fx);
  Teuchos::RCP<Epetra_Vector> fcx = fluid_field()->interface()->extract_fsi_cond_vector(fx);

  ddiinc_ = Teuchos::rcp(new Epetra_Vector(*sox));  // first iteration increment

  ddginc_ = Teuchos::rcp(new Epetra_Vector(*scx));  // first iteration increment

  duiinc_ = Teuchos::rcp(new Epetra_Vector(*fox));  // first iteration increment

  duginc_ = Teuchos::rcp(new Epetra_Vector(*fcx));  // first iteration increment

  double stiparam = structure_field()->tim_int_param();

  // store the product Cfs_{\GammaI} \Delta d_I^{n+1} in here
  Teuchos::RCP<Epetra_Vector> cfsgiddi =
      Core::LinAlg::CreateVector(*fluid_field()->interface()->fsi_cond_map(), true);
  // compute the above mentioned product
  cfsgicur_->multiply(false, *ddiinc_, *cfsgiddi);

  // store the product F_{\GammaI} \Delta u_I^{n+1} in here
  Teuchos::RCP<Epetra_Vector> fgiddi =
      Core::LinAlg::CreateVector(*fluid_field()->interface()->fsi_cond_map(), true);
  // compute the above mentioned product
  fgicur_->multiply(false, *duiinc_, *fgiddi);

  // store the product Cfs_{\Gamma\Gamma} \Delta d_\Gamma^{n+1} in here
  Teuchos::RCP<Epetra_Vector> cfsggddg =
      Core::LinAlg::CreateVector(*fluid_field()->interface()->fsi_cond_map(), true);
  // compute the above mentioned product
  cfsggcur_->multiply(false, *ddginc_, *cfsggddg);

  // store the prodcut F_{\Gamma\Gamma} \Delta u_\Gamma^{n+1} in here
  Teuchos::RCP<Epetra_Vector> fggddg =
      Core::LinAlg::CreateVector(*fluid_field()->interface()->fsi_cond_map(), true);
  // compute the above mentioned product
  fggcur_->multiply(false, *duginc_, *fggddg);

  // Update the Lagrange multiplier:
  /* \lambda^{n+1}_{i} =  -1/b * invD^{n+1} * [
   *                          + CFS_{\Gamma I} \Delta d_I
   *                          + CFS_{\Gamma \Gamma} \Delta d_\Gamma
   *                          + F_{\Gamma I} \Delta u_I
   *                          + F_{\Gamma\Gamma} \Delta u_\Gamma
   *                          - f_{\Gamma}^f]
   *                          - (1-b)/b * invD^{n+1} * D^n * \lambda^n
   */

  Teuchos::RCP<Epetra_Vector> tmplambda =
      Teuchos::rcp(new Epetra_Vector(*fluid_field()->interface()->fsi_cond_map(), true));

  tmplambda->Update(1.0, *cfsgiddi, 0.0);
  tmplambda->Update(1.0, *fgiddi, 1.0);
  tmplambda->Update(1.0, *cfsggddg, 1.0);
  tmplambda->Update(1.0, *fggddg, 1.0);
  tmplambda->Update(-1.0, *rhs_fgcur_, 1.0);

  if (k_dn_ != Teuchos::null)  // for first timestep lambda = 0 !
  {
    Teuchos::RCP<Epetra_Vector> Dlam =
        Teuchos::rcp(new Epetra_Vector(*fluid_field()->interface()->fsi_cond_map(), true));
    k_dn_->Apply(*lambda_, *Dlam);  // D(n)*lambda(n)
    Dlam->Scale(stiparam);          //*(1-b)
    tmplambda->Update(1.0, *Dlam, 1.0);
  }

  k_inv_d_->Apply(*tmplambda, *lambdanp_);
  lambdanp_->Scale(-1 / (1.0 - stiparam));  //*-1/b
}

void PoroElast::MonolithicSplitNoPenetration::setup_system_matrix(
    Core::LinAlg::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("PoroElast::MonolithicSplitNoPenetration::setup_system_matrix");

  Teuchos::RCP<Core::LinAlg::SparseMatrix> s = structure_field()->system_matrix();
  if (s == Teuchos::null) FOUR_C_THROW("expect structure matrix");
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> f = fluid_field()->block_system_matrix();
  if (f == Teuchos::null) FOUR_C_THROW("expect fluid block matrix");

  // Get Idx of fluid and structure field map extractor
  const int& fidx_other = FLD::UTILS::MapExtractor::cond_other;
  const int& fidx_nopen = FLD::UTILS::MapExtractor::cond_fsi;

  const int& sidx_other = Solid::MapExtractor::cond_other;
  const int& sidx_nopen = Solid::MapExtractor::cond_fsi;

  /*----------------------------------------------------------------------*/

  // just to play it safe ...
  mat.reset();

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  /*----------------------------------------------------------------------*/
  // structural part k_sf (3nxn)
  // build mechanical-fluid block

  // create empty matrix
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> k_sf = struct_fluid_coupling_block_matrix();

  // call the element and calculate the matrix block
  apply_str_coupl_matrix(k_sf);

  /*----------------------------------------------------------------------*/
  // fluid part k_fs ( (3n+1)x3n )
  // build fluid-mechanical block

  // create empty matrix
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> k_fs = fluid_struct_coupling_block_matrix();

  // call the element and calculate the matrix block
  apply_fluid_coupl_matrix(k_fs);

  /*----------------------------------------------------------------------*/

  k_fs->complete();
  k_sf->complete();

  /*----------------------------------------------------------------------*/
  // pure structural part
  mat.assign(0, 0, Core::LinAlg::View, *s);

  // structure coupling part
  mat.matrix(0, 1).add(k_sf->matrix(sidx_other, fidx_other), false, 1.0, 0.0);
  mat.matrix(0, 1).add(k_sf->matrix(sidx_other, fidx_nopen), false, 1.0, 1.0);
  mat.matrix(0, 1).add(k_sf->matrix(sidx_nopen, fidx_other), false, 1.0, 1.0);
  mat.matrix(0, 1).add(k_sf->matrix(sidx_nopen, fidx_nopen), false, 1.0, 1.0);
  /*----------------------------------------------------------------------*/
  // pure fluid part
  // incomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  // this just once...
  // f->UnComplete();

  mat.matrix(1, 1).add(f->matrix(fidx_other, fidx_other), false, 1.0, 0.0);
  mat.matrix(1, 1).add(f->matrix(fidx_other, fidx_nopen), false, 1.0, 1.0);

  // fluid coupling part
  mat.matrix(1, 0).add(k_fs->matrix(fidx_other, fidx_other), false, 1.0, 0.0);
  mat.matrix(1, 0).add(k_fs->matrix(fidx_other, fidx_nopen), false, 1.0, 1.0);

  /*----------------------------------------------------------------------*/
  /*Add lines for poro nopenetration condition*/

  fgicur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(f->matrix(fidx_nopen, fidx_other)));
  fggcur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(f->matrix(fidx_nopen, fidx_nopen)));
  cfsgicur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(k_fs->matrix(fidx_nopen, sidx_other)));
  cfsggcur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(k_fs->matrix(fidx_nopen, sidx_nopen)));

  Teuchos::RCP<Core::LinAlg::SparseMatrix> tanginvDkfsgi =
      Core::LinAlg::MLMultiply(*k_lambdainv_d_, *cfsgicur_, true);  // T*D^-1*K^FS_gi;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tanginvDfgi =
      Core::LinAlg::MLMultiply(*k_lambdainv_d_, *fgicur_, true);  // T*D^-1*Fgi;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tanginvDfgg =
      Core::LinAlg::MLMultiply(*k_lambdainv_d_, *fggcur_, true);  // T*D^-1*Fgg;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tanginvDkfsgg =
      Core::LinAlg::MLMultiply(*k_lambdainv_d_, *cfsggcur_, true);  // T*D^-1*K^FS_gg;

  mat.matrix(1, 0).add(*tanginvDkfsgi, false, -1.0, 1.0);
  mat.matrix(1, 0).add(*tanginvDkfsgg, false, -1.0, 1.0);
  mat.matrix(1, 0).add(*k_struct_, false, 1.0, 1.0);
  mat.matrix(1, 0).add(k_porodisp_->matrix(1, 0), false, 1.0, 1.0);
  mat.matrix(1, 0).add(k_porodisp_->matrix(1, 1), false, 1.0, 1.0);
  mat.matrix(1, 1).add(*tanginvDfgi, false, -1.0, 1.0);
  mat.matrix(1, 1).add(*k_fluid_, false, 1.0, 1.0);
  mat.matrix(1, 1).add(*tanginvDfgg, false, -1.0, 1.0);
  mat.matrix(1, 1).add(*k_porofluid_, false, 1.0, 1.0);

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.complete();
}

void PoroElast::MonolithicSplitNoPenetration::apply_fluid_coupl_matrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator> k_fs)
{
  // call base class
  Monolithic::apply_fluid_coupl_matrix(k_fs);

  // reset
  k_fluid_->zero();
  k_d_->zero();
  k_inv_d_->zero();
  k_struct_->zero();
  k_lambda_->zero();
  k_porodisp_->zero();
  k_porofluid_->zero();
  nopenetration_rhs_->PutScalar(0.0);

  Teuchos::RCP<Core::LinAlg::SparseMatrix> tmp_k_D = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *(fluid_field()->interface()->fsi_cond_map()), 81, false, false));

  // fill diagonal blocks
  {
    // create the parameters for the discretization
    Teuchos::ParameterList params;
    // action for elements
    params.set<int>("action", FLD::poro_splitnopenetration);
    params.set("total time", time());
    params.set("delta time", dt());
    params.set("timescale", fluid_field()->residual_scaling());
    params.set<int>("Physical Type", fluid_field()->physical_type());

    fluid_field()->discretization()->clear_state();
    fluid_field()->discretization()->set_state(0, "dispnp", fluid_field()->dispnp());
    fluid_field()->discretization()->set_state(0, "gridv", fluid_field()->grid_vel());
    fluid_field()->discretization()->set_state(0, "velnp", fluid_field()->velnp());
    fluid_field()->discretization()->set_state(0, "scaaf", fluid_field()->scaaf());

    // fluid_field()->discretization()->set_state(0,"lambda",
    //    fluid_field()->Interface()->insert_fsi_cond_vector(structure_to_fluid_at_interface(lambdanp_)));

    // build specific assemble strategy for the fluid-mechanical system matrix
    // from the point of view of fluid_field:
    // fluiddofset = 0, structdofset = 1

    Core::FE::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
        0,                                       // fluiddofset for column
        k_fluid_, Teuchos::null, nopenetration_rhs_, Teuchos::null, Teuchos::null);
    fluid_field()->discretization()->evaluate_condition(params, fluidstrategy, "FSICoupling");

    fluid_field()->discretization()->clear_state();
  }

  Teuchos::RCP<Epetra_Vector> disp_interface =
      fluid_field()->interface()->extract_fsi_cond_vector(fluid_field()->dispnp());
  mortar_adapter_->integrate_lin_d(
      "displacement", disp_interface, structure_to_fluid_at_interface(lambdanp_));
  tmp_k_D = mortar_adapter_->get_mortar_matrix_d();

  // fill off diagonal blocks
  {
    // create the parameters for the discretization
    Teuchos::ParameterList params;
    // action for elements
    params.set<int>("action", FLD::poro_splitnopenetration_OD);
    params.set("total time", time());
    params.set("delta time", dt());
    params.set("timescale", fluid_field()->residual_scaling());
    params.set<int>("Physical Type", fluid_field()->physical_type());

    fluid_field()->discretization()->clear_state();
    fluid_field()->discretization()->set_state(0, "dispnp", fluid_field()->dispnp());
    fluid_field()->discretization()->set_state(0, "gridv", fluid_field()->grid_vel());
    fluid_field()->discretization()->set_state(0, "velnp", fluid_field()->velnp());
    fluid_field()->discretization()->set_state(0, "scaaf", fluid_field()->scaaf());

    fluid_field()->discretization()->set_state(0, "lambda",
        fluid_field()->interface()->insert_fsi_cond_vector(
            structure_to_fluid_at_interface(lambdanp_)));

    // build specific assemble strategy for the fluid-mechanical system matrix
    // from the point of view of fluid_field:
    // fluiddofset = 0, structdofset = 1
    Core::FE::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
        1,                                       // structdofset for column
        k_struct_,                               // fluid-mechanical matrix
        k_lambda_, Teuchos::null, Teuchos::null, Teuchos::null);
    fluid_field()->discretization()->evaluate_condition(params, fluidstrategy, "FSICoupling");

    fluid_field()->discretization()->clear_state();
  }

  // fill off diagonal blocks
  {
    // create the parameters for the discretization
    Teuchos::ParameterList params;
    // action for elements
    params.set<int>("action", FLD::poro_splitnopenetration_ODdisp);
    params.set("total time", time());
    params.set("delta time", dt());
    params.set("timescale", fluid_field()->residual_scaling());
    params.set<int>("Physical Type", fluid_field()->physical_type());

    fluid_field()->discretization()->clear_state();
    fluid_field()->discretization()->set_state(0, "dispnp", fluid_field()->dispnp());
    fluid_field()->discretization()->set_state(0, "gridv", fluid_field()->grid_vel());
    fluid_field()->discretization()->set_state(0, "velnp", fluid_field()->velnp());
    //  fluid_field()->discretization()->set_state(0,"scaaf",fluid_field()->Scaaf());

    // build specific assemble strategy for the fluid-mechanical system matrix
    // from the point of view of fluid_field:
    // fluiddofset = 0, structdofset = 1
    Core::FE::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
        1,                                       // structdofset for column
        k_porodisp_,                             // fluid-mechanical matrix
        Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
    fluid_field()->discretization()->evaluate_condition(params, fluidstrategy, "FSICoupling");

    fluid_field()->discretization()->clear_state();
  }

  // fill off diagonal blocks
  {
    // create the parameters for the discretization
    Teuchos::ParameterList params;
    // action for elements
    params.set<int>("action", FLD::poro_splitnopenetration_ODpres);
    params.set("total time", time());
    params.set("delta time", dt());
    params.set("timescale", fluid_field()->residual_scaling());
    params.set<int>("Physical Type", fluid_field()->physical_type());

    fluid_field()->discretization()->clear_state();
    fluid_field()->discretization()->set_state(0, "dispnp", fluid_field()->dispnp());
    fluid_field()->discretization()->set_state(0, "gridv", fluid_field()->grid_vel());
    fluid_field()->discretization()->set_state(0, "velnp", fluid_field()->velnp());
    //  fluid_field()->discretization()->set_state(0,"scaaf",fluid_field()->Scaaf());

    // build specific assemble strategy for the fluid-mechanical system matrix
    // from the point of view of fluid_field:
    // fluiddofset = 0, structdofset = 1
    Core::FE::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
        0,                                       // fluiddofset for column
        k_porofluid_,                            // fluid-mechanical matrix
        Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
    fluid_field()->discretization()->evaluate_condition(params, fluidstrategy, "FSICoupling");

    fluid_field()->discretization()->clear_state();
  }

  // Complete Coupling matrices which should be *.Add later!
  k_struct_->complete(
      *structure_field()->interface()->fsi_cond_map(), *fluid_field()->interface()->fsi_cond_map());
  k_fluid_->complete();
  k_porofluid_->complete();
  k_porodisp_->complete();

  //------------------------------invert D Matrix!-----------------------------------------------
  tmp_k_D->complete();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> invd =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*tmp_k_D, Core::LinAlg::Copy));
  // invd->Complete();

  Teuchos::RCP<Epetra_Vector> diag =
      Core::LinAlg::CreateVector(*fluid_field()->interface()->fsi_cond_map(), true);

  int err = 0;

  // extract diagonal of invd into diag
  invd->extract_diagonal_copy(
      *diag);  // That the Reason, why tmp_k_D has to have Fluid Maps for Rows & Columns!!!

  // set zero diagonal values to dummy 1.0 ??
  for (int i = 0; i < diag->MyLength(); ++i)
  {
    if ((*diag)[i] == 0.0)
    {
      (*diag)[i] = 1.0;
      std::cout << "--- --- --- WARNING: D-Matrix Diagonal Element " << i
                << " is zero!!! --- --- ---" << std::endl;
    }
  }

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err > 0) FOUR_C_THROW("ERROR: Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd->replace_diagonal_values(*diag);
  invd->complete();
  //------------------------------End of invert D
  // Matrix!-----------------------------------------------

  // Transform also colum map of D-Matrix
  (*k_d_transform_)(*fluid_field()->interface()->fsi_cond_map(),
      fluid_field()->block_system_matrix()->matrix(1, 1).col_map(), *tmp_k_D, 1.0,
      Core::Adapter::CouplingSlaveConverter(*icoupfs_), *k_d_);

  (*k_inv_d_transform_)(
      *invd, 1.0, Core::Adapter::CouplingSlaveConverter(*icoupfs_), *k_inv_d_, false);

  double stiparam = structure_field()->tim_int_param();

  Teuchos::RCP<Core::LinAlg::SparseMatrix> tmp_k_DLin = mortar_adapter_->d_lin_matrix();
  tmp_k_DLin->complete();

  // Transform also column map of D-Matrix
  (*k_d_lin_transform_)(*fluid_field()->interface()->fsi_cond_map(),
      fluid_field()->block_system_matrix()->matrix(1, 1).col_map(), *tmp_k_DLin,
      1.0 - stiparam,  // *= b
      Core::Adapter::CouplingSlaveConverter(*icoupfs_),
      (Teuchos::rcp_static_cast<Core::LinAlg::BlockSparseMatrixBase>(k_fs))->matrix(1, 1), true,
      true);

  k_lambda_->complete(
      *structure_field()->interface()->fsi_cond_map(), *fluid_field()->interface()->fsi_cond_map());
  k_inv_d_->complete(
      *fluid_field()->interface()->fsi_cond_map(), *structure_field()->interface()->fsi_cond_map());

  // Calculate 1/b*Tangent*invD
  k_lambdainv_d_ = Core::LinAlg::MLMultiply(*k_lambda_, *k_inv_d_, true);
  k_lambdainv_d_->scale(1.0 / (1.0 - stiparam));  // *= 1/b
}

void PoroElast::MonolithicSplitNoPenetration::apply_str_coupl_matrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator> k_sf  //!< off-diagonal tangent matrix term
)
{
  // call base class
  Monolithic::apply_str_coupl_matrix(k_sf);
}

void PoroElast::MonolithicSplitNoPenetration::recover_lagrange_multiplier_after_time_step()
{
  // we do not need to recover after a time step, it is done after every newton step
}

void PoroElast::MonolithicSplitNoPenetration::update()
{
  // call base class
  MonolithicSplit::update();

  // update lagrangean multiplier
  lambda_->Update(1.0, *lambdanp_, 0.0);

  // copy D matrix from current time step to old D matrix
  k_dn_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *k_d_, Core::LinAlg::Copy));  // store D-Matrix from last timestep
}

void PoroElast::MonolithicSplitNoPenetration::output(bool forced_writerestart)
{
  // call base class
  MonolithicSplit::output(forced_writerestart);

  // for now, we always write the lagrange multiplier
  Teuchos::RCP<Epetra_Vector> fulllambda =
      Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(*structure_field()->dof_row_map()));
  Core::LinAlg::Export(*lambdanp_, *fulllambda);
  structure_field()->disc_writer()->write_vector("poronopencond_lambda", fulllambda);
}

void PoroElast::MonolithicSplitNoPenetration::setup_coupling_and_matrices()
{
  const int ndim = Global::Problem::instance()->n_dim();
  icoupfs_->setup_condition_coupling(*structure_field()->discretization(),
      structure_field()->interface()->fsi_cond_map(), *fluid_field()->discretization(),
      fluid_field()->interface()->fsi_cond_map(), "FSICoupling", ndim);

  evaluateinterface_ = false;

  // initialize Poroelasticity-systemmatrix_
  systemmatrix_ =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *extractor(), *extractor(), 81, false, true));

  // initialize coupling matrices
  k_fs_ =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *(structure_field()->interface()), *(fluid_field()->interface()), 81, false, true));

  k_sf_ =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *(fluid_field()->interface()), *(structure_field()->interface()), 81, false, true));

  // initialize no penetration coupling matrices
  k_struct_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *(fluid_field()->interface()->fsi_cond_map()), 81, true, true));

  k_fluid_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *(fluid_field()->interface()->fsi_cond_map()), 81, false, false));

  k_lambda_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *(fluid_field()->interface()->fsi_cond_map()), 81, true, true));

  k_d_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *(fluid_field()->interface()->fsi_cond_map()), 81, true, true));

  k_inv_d_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *(structure_field()->interface()->fsi_cond_map()), 81, true, true));

  k_porodisp_ =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *(structure_field()->interface()), *(fluid_field()->interface()), 81, true, true));

  k_porofluid_ =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*(fluid_field()->dof_row_map()), 81, true, true));

  nopenetration_rhs_ =
      Teuchos::rcp(new Epetra_Vector(*fluid_field()->interface()->fsi_cond_map(), true));
}

void PoroElast::MonolithicSplitNoPenetration::prepare_time_step()
{
  // call base class
  PoroElast::Monolithic::prepare_time_step();
}

void PoroElast::MonolithicSplitNoPenetration::read_restart(const int step)
{
  // call base class
  PoroElast::PoroBase::read_restart(step);

  // get lagrange multiplier and D matrix
  if (step)
  {
    // get the structure reader (this is where the lagrange multiplier was saved)
    Core::IO::DiscretizationReader reader(structure_field()->discretization(),
        Global::Problem::instance()->input_control_file(), structure_field()->step());
    Teuchos::RCP<Epetra_Vector> fulllambda =
        Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(*structure_field()->dof_row_map()));

    // this is the lagrange multiplier on the whole structure field
    reader.read_vector(fulllambda, "poronopencond_lambda");

    // extract lambda on fsi interface vector
    lambda_ = structure_field()->interface()->extract_fsi_cond_vector(fulllambda);
    lambdanp_->Update(1.0, *lambda_, 0.0);

    // call an additional evaluate to get the old D matrix
    setup_system();
    // call evaluate to recalculate D matrix
    evaluate(zeros_, false);

    // copy D matrix from current time step to old D matrix
    k_dn_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
        *k_d_, Core::LinAlg::Copy));  // store D-Matrix from last timestep
  }
}

void PoroElast::MonolithicSplitNoPenetration::print_newton_iter_header_stream(
    std::ostringstream& oss)
{
  Monolithic::print_newton_iter_header_stream(oss);

  oss << std::setw(20) << "abs-crhs-res";
}

void PoroElast::MonolithicSplitNoPenetration::print_newton_iter_text_stream(std::ostringstream& oss)
{
  Monolithic::print_newton_iter_text_stream(oss);

  oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhs_nopenetration_;
}

void PoroElast::MonolithicSplitNoPenetration::build_convergence_norms()
{
  Monolithic::build_convergence_norms();

  normrhs_nopenetration_ = UTILS::calculate_vector_norm(vectornormfres_, nopenetration_rhs_);
}

FOUR_C_NAMESPACE_CLOSE
