/*----------------------------------------------------------------------*/
/*!

 \brief porous medium algorithm with matrix split for condensation of
      no-penetration constraint

\level 2

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/

#include "poro_monolithicsplit_nopenetration.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/ad_str_fpsiwrapper.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_coupling_nonlin_mortar.H"

#include "../drt_contact/contact_interface.H"

#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fsi/fsi_matrixtransform.H"

#include "../drt_io/io.H"

#include "../drt_structure/stru_aux.H"

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_blocksparsematrix.H"


/*----------------------------------------------------------------------*
 |                                                         vuong 09/14  |
 *----------------------------------------------------------------------*/
POROELAST::MonolithicSplitNoPenetration::MonolithicSplitNoPenetration(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : MonolithicSplit(comm, timeparams), normrhs_nopenetration_(-1.0)
{
  // Initialize Transformation Objects
  k_D_transform_ = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
  k_invD_transform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);

  k_DLin_transform_ = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);

  // Recovering of Lagrange multiplier happens on fluid field
  lambda_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->FSICondMap()));
  lambdanp_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->FSICondMap()));

  k_Dn_ = Teuchos::null;

  mortar_adapter_ = Teuchos::rcp(new ADAPTER::CouplingNonLinMortar);

  return;
}

/*----------------------------------------------------------------------*
 | setup system (called in porolast.cpp)               vuong 09/14     |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::SetupSystem()
{
  {
    const int ndim = DRT::Problem::Instance()->NDim();
    std::vector<int> coupleddof(ndim + 1, 1);
    coupleddof[ndim] = 0;

    mortar_adapter_->Setup(StructureField()->Discretization(), FluidField()->Discretization(),
        coupleddof, "FSICoupling");
  }

  // use full maps of both fields. Only Lagrange multipliers are condensed
  {
    // create combined map
    std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;

    vecSpaces.push_back(StructureField()->DofRowMap());
    vecSpaces.push_back(FluidField()->DofRowMap());

    if (vecSpaces[0]->NumGlobalElements() == 0) dserror("No structure equation. Panic.");
    if (vecSpaces[1]->NumGlobalElements() == 0) dserror("No fluid equation. Panic.");

    // full Poroelasticity-map
    fullmap_ = LINALG::MultiMapExtractor::MergeMaps(vecSpaces);
    // full Poroelasticity-blockmap
    blockrowdofmap_->Setup(*fullmap_, vecSpaces);
  }

  // initialize vectors for row and column sums of global system matrix if necessary
  if (rowequilibration_)
    invrowsums_ = Teuchos::rcp(new Epetra_Vector(*blockrowdofmap_->FullMap(), false));

  // Switch fluid to interface split block matrix
  FluidField()->UseBlockMatrix(true);

  // setup coupling objects, system and coupling matrices
  SetupCouplingAndMatrices();

  // build map of dofs subjected to a DBC of whole problem
  BuildCombinedDBCMap();
}  // SetupSystem()

/*----------------------------------------------------------------------*
 |                                                         vuong 09/14  |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::SetupRHS(bool firstcall)
{
  // only Lagrange multipliers are condensed -> use unchanged maps from single fields
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicSplitNoPenetration::SetupRHS");

  // create full monolithic rhs vector
  if (rhs_ == Teuchos::null) rhs_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  SetupVector(*rhs_, StructureField()->RHS(), FluidField()->RHS());
}

/*----------------------------------------------------------------------*
 | setup vector of the structure and fluid field                        |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::SetupVector(
    Epetra_Vector& f, Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv)
{
  // extract dofs of the two fields
  // and put the structural/fluid field vector into the global vector f
  // noticing the block number

  Extractor()->InsertVector(*sv, 0, f);

  Teuchos::RCP<Epetra_Vector> fov = FluidField()->Interface()->ExtractOtherVector(fv);
  Teuchos::RCP<Epetra_Vector> fcv = FluidField()->Interface()->ExtractFSICondVector(fv);

  Teuchos::RCP<Epetra_Vector> Dlam =
      Teuchos::rcp(new Epetra_Vector(*FluidField()->Interface()->FSICondMap(), true));
  Teuchos::RCP<Epetra_Vector> couprhs =
      Teuchos::rcp(new Epetra_Vector(*FluidField()->Interface()->FSICondMap(), true));
  if (k_Dn_ != Teuchos::null)
  {
    double stiparam = StructureField()->TimIntParam();

    k_Dn_->Multiply(false, *lambda_, *Dlam);  // D(n)*lambda(n)

    Dlam->Scale(stiparam);  //*(1-b)
  }
  Dlam->Update(-1.0, *fcv, 1.0);
  k_lambdainvD_->Multiply(false, *Dlam, *couprhs);

  couprhs->Update(1.0, *nopenetration_rhs_, 1.0);

  // std::cout << "nopenetration_rhs_: " << *nopenetration_rhs_ << std::endl;

  Teuchos::RCP<Epetra_Vector> fullcouprhs =
      Teuchos::rcp(new Epetra_Vector(*FluidField()->DofRowMap(), true));
  LINALG::Export(*couprhs, *fullcouprhs);
  Extractor()->InsertVector(*fullcouprhs, 1, f);

  Teuchos::RCP<Epetra_Vector> fullfov =
      Teuchos::rcp(new Epetra_Vector(*FluidField()->DofRowMap(), true));
  LINALG::Export(*fov, *fullfov);
  Extractor()->AddVector(*fullfov, 1, f, 1.0);

  rhs_fgcur_ = fcv;  // Store interface rhs for recovering of lagrange multiplier
  return;
}

/*----------------------------------------------------------------------*
| Recover the Lagrange multipliers for poro no pen.          ager 09/14 |
*----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::RecoverLagrangeMultiplierAfterNewtonStep(
    Teuchos::RCP<const Epetra_Vector> x)
{
  // call base class
  Monolithic::RecoverLagrangeMultiplierAfterNewtonStep(x);


  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  ExtractFieldVectors(x, sx, fx);

  Teuchos::RCP<Epetra_Vector> sox = StructureField()->Interface()->ExtractOtherVector(sx);
  Teuchos::RCP<Epetra_Vector> scx = StructureField()->Interface()->ExtractFSICondVector(sx);
  Teuchos::RCP<Epetra_Vector> fox = FluidField()->Interface()->ExtractOtherVector(fx);
  Teuchos::RCP<Epetra_Vector> fcx = FluidField()->Interface()->ExtractFSICondVector(fx);

  ddiinc_ = Teuchos::rcp(new Epetra_Vector(*sox));  // first iteration increment

  ddginc_ = Teuchos::rcp(new Epetra_Vector(*scx));  // first iteration increment

  duiinc_ = Teuchos::rcp(new Epetra_Vector(*fox));  // first iteration increment

  duginc_ = Teuchos::rcp(new Epetra_Vector(*fcx));  // first iteration increment

  double stiparam = StructureField()->TimIntParam();

  // store the product Cfs_{\GammaI} \Delta d_I^{n+1} in here
  Teuchos::RCP<Epetra_Vector> cfsgiddi =
      LINALG::CreateVector(*FluidField()->Interface()->FSICondMap(), true);
  // compute the above mentioned product
  cfsgicur_->Multiply(false, *ddiinc_, *cfsgiddi);

  // store the product F_{\GammaI} \Delta u_I^{n+1} in here
  Teuchos::RCP<Epetra_Vector> fgiddi =
      LINALG::CreateVector(*FluidField()->Interface()->FSICondMap(), true);
  // compute the above mentioned product
  fgicur_->Multiply(false, *duiinc_, *fgiddi);

  // store the product Cfs_{\Gamma\Gamma} \Delta d_\Gamma^{n+1} in here
  Teuchos::RCP<Epetra_Vector> cfsggddg =
      LINALG::CreateVector(*FluidField()->Interface()->FSICondMap(), true);
  // compute the above mentioned product
  cfsggcur_->Multiply(false, *ddginc_, *cfsggddg);

  // store the prodcut F_{\Gamma\Gamma} \Delta u_\Gamma^{n+1} in here
  Teuchos::RCP<Epetra_Vector> fggddg =
      LINALG::CreateVector(*FluidField()->Interface()->FSICondMap(), true);
  // compute the above mentioned product
  fggcur_->Multiply(false, *duginc_, *fggddg);

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
      Teuchos::rcp(new Epetra_Vector(*FluidField()->Interface()->FSICondMap(), true));

  tmplambda->Update(1.0, *cfsgiddi, 0.0);
  tmplambda->Update(1.0, *fgiddi, 1.0);
  tmplambda->Update(1.0, *cfsggddg, 1.0);
  tmplambda->Update(1.0, *fggddg, 1.0);
  tmplambda->Update(-1.0, *rhs_fgcur_, 1.0);

  if (k_Dn_ != Teuchos::null)  // for first timestep lambda = 0 !
  {
    Teuchos::RCP<Epetra_Vector> Dlam =
        Teuchos::rcp(new Epetra_Vector(*FluidField()->Interface()->FSICondMap(), true));
    k_Dn_->Apply(*lambda_, *Dlam);  // D(n)*lambda(n)
    Dlam->Scale(stiparam);          //*(1-b)
    tmplambda->Update(1.0, *Dlam, 1.0);
  }

  k_invD_->Apply(*tmplambda, *lambdanp_);
  lambdanp_->Scale(-1 / (1.0 - stiparam));  //*-1/b

  return;
}

/*----------------------------------------------------------------------*
 |                                                         vuong 09/14  |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicSplitNoPenetration::SetupSystemMatrix");

  // Teuchos::RCP<LINALG::BlockSparseMatrixBase> s = StructureField()->BlockSystemMatrix();
  Teuchos::RCP<LINALG::SparseMatrix> s = StructureField()->SystemMatrix();
  if (s == Teuchos::null) dserror("expect structure matrix");
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> f = FluidField()->BlockSystemMatrix();
  if (f == Teuchos::null) dserror("expect fluid block matrix");

  // Get Idx of fluid and structure field map extractor
  const int& fidx_other = FLD::UTILS::MapExtractor::cond_other;
  const int& fidx_nopen = FLD::UTILS::MapExtractor::cond_fsi;

  const int& sidx_other = STR::AUX::MapExtractor::cond_other;
  const int& sidx_nopen = STR::AUX::MapExtractor::cond_fsi;

  /*----------------------------------------------------------------------*/

  // just to play it save ...
  mat.Reset();

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  /*----------------------------------------------------------------------*/
  // structural part k_sf (3nxn)
  // build mechanical-fluid block

  // create empty matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> k_sf = StructFluidCouplingBlockMatrix();

  // call the element and calculate the matrix block
  ApplyStrCouplMatrix(k_sf);

  /*----------------------------------------------------------------------*/
  // fluid part k_fs ( (3n+1)x3n )
  // build fluid-mechanical block

  // create empty matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> k_fs = FluidStructCouplingBlockMatrix();

  // call the element and calculate the matrix block
  ApplyFluidCouplMatrix(k_fs);

  /*----------------------------------------------------------------------*/

  k_fs->Complete();
  k_sf->Complete();

  /*----------------------------------------------------------------------*/
  // pure structural part
  mat.Assign(0, 0, LINALG::View, *s);

  // structure coupling part
  mat.Matrix(0, 1).Add(k_sf->Matrix(sidx_other, fidx_other), false, 1.0, 0.0);
  mat.Matrix(0, 1).Add(k_sf->Matrix(sidx_other, fidx_nopen), false, 1.0, 1.0);
  mat.Matrix(0, 1).Add(k_sf->Matrix(sidx_nopen, fidx_other), false, 1.0, 1.0);
  mat.Matrix(0, 1).Add(k_sf->Matrix(sidx_nopen, fidx_nopen), false, 1.0, 1.0);
  /*----------------------------------------------------------------------*/
  // pure fluid part
  // uncomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  // this just once...
  // f->UnComplete();

  mat.Matrix(1, 1).Add(f->Matrix(fidx_other, fidx_other), false, 1.0, 0.0);
  mat.Matrix(1, 1).Add(f->Matrix(fidx_other, fidx_nopen), false, 1.0, 1.0);

  // fluid coupling part
  mat.Matrix(1, 0).Add(k_fs->Matrix(fidx_other, fidx_other), false, 1.0, 0.0);
  mat.Matrix(1, 0).Add(k_fs->Matrix(fidx_other, fidx_nopen), false, 1.0, 1.0);

  /*----------------------------------------------------------------------*/
  /*Add lines for poro nopenetration condition*/

  fgicur_ = Teuchos::rcp(new LINALG::SparseMatrix(f->Matrix(fidx_nopen, fidx_other)));
  fggcur_ = Teuchos::rcp(new LINALG::SparseMatrix(f->Matrix(fidx_nopen, fidx_nopen)));
  cfsgicur_ = Teuchos::rcp(new LINALG::SparseMatrix(k_fs->Matrix(fidx_nopen, sidx_other)));
  cfsggcur_ = Teuchos::rcp(new LINALG::SparseMatrix(k_fs->Matrix(fidx_nopen, sidx_nopen)));

  Teuchos::RCP<LINALG::SparseMatrix> tanginvDkfsgi =
      LINALG::MLMultiply(*k_lambdainvD_, *cfsgicur_, true);  // T*D^-1*K^FS_gi;
  Teuchos::RCP<LINALG::SparseMatrix> tanginvDfgi =
      LINALG::MLMultiply(*k_lambdainvD_, *fgicur_, true);  // T*D^-1*Fgi;
  Teuchos::RCP<LINALG::SparseMatrix> tanginvDfgg =
      LINALG::MLMultiply(*k_lambdainvD_, *fggcur_, true);  // T*D^-1*Fgg;
  Teuchos::RCP<LINALG::SparseMatrix> tanginvDkfsgg =
      LINALG::MLMultiply(*k_lambdainvD_, *cfsggcur_, true);  // T*D^-1*K^FS_gg;

  mat.Matrix(1, 0).Add(*tanginvDkfsgi, false, -1.0, 1.0);
  mat.Matrix(1, 0).Add(*tanginvDkfsgg, false, -1.0, 1.0);
  mat.Matrix(1, 0).Add(*k_struct_, false, 1.0, 1.0);
  mat.Matrix(1, 0).Add(k_porodisp_->Matrix(1, 0), false, 1.0, 1.0);
  mat.Matrix(1, 0).Add(k_porodisp_->Matrix(1, 1), false, 1.0, 1.0);
  mat.Matrix(1, 1).Add(*tanginvDfgi, false, -1.0, 1.0);
  mat.Matrix(1, 1).Add(*k_fluid_, false, 1.0, 1.0);
  mat.Matrix(1, 1).Add(*tanginvDfgg, false, -1.0, 1.0);
  mat.Matrix(1, 1).Add(*k_porofluid_, false, 1.0, 1.0);

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();

  return;
}

/*----------------------------------------------------------------------*
 |    evaluate fluid-structural system matrix at state      vuong 09/14   |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::ApplyFluidCouplMatrix(
    Teuchos::RCP<LINALG::SparseOperator> k_fs  //!< off-diagonal tangent matrix term
)
{
  // call base class
  Monolithic::ApplyFluidCouplMatrix(k_fs);

  // reset
  k_fluid_->Zero();
  k_D_->Zero();
  k_invD_->Zero();
  k_struct_->Zero();
  k_lambda_->Zero();
  k_porodisp_->Zero();
  k_porofluid_->Zero();
  nopenetration_rhs_->Scale(0.0);

  Teuchos::RCP<LINALG::SparseMatrix> tmp_k_D = Teuchos::rcp(
      new LINALG::SparseMatrix(*(FluidField()->Interface()->FSICondMap()), 81, false, false));

  // fill diagonal blocks
  {
    // create the parameters for the discretization
    Teuchos::ParameterList params;
    // action for elements
    params.set<int>("action", FLD::poro_splitnopenetration);
    params.set("total time", Time());
    params.set("delta time", Dt());
    params.set("timescale", FluidField()->ResidualScaling());
    params.set<int>("Physical Type", FluidField()->PhysicalType());

    FluidField()->Discretization()->ClearState();
    FluidField()->Discretization()->SetState(0, "dispnp", FluidField()->Dispnp());
    FluidField()->Discretization()->SetState(0, "gridv", FluidField()->GridVel());
    FluidField()->Discretization()->SetState(0, "velnp", FluidField()->Velnp());
    FluidField()->Discretization()->SetState(0, "scaaf", FluidField()->Scaaf());

    // FluidField()->Discretization()->SetState(0,"lambda",
    //    FluidField()->Interface()->InsertFSICondVector(StructureToFluidAtInterface(lambdanp_)));

    // build specific assemble strategy for the fluid-mechanical system matrix
    // from the point of view of FluidField:
    // fluiddofset = 0, structdofset = 1

    DRT::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
        0,                                  // fluiddofset for column
        k_fluid_, Teuchos::null, nopenetration_rhs_, Teuchos::null, Teuchos::null);
    FluidField()->Discretization()->EvaluateCondition(params, fluidstrategy, "FSICoupling");

    FluidField()->Discretization()->ClearState();
  }

  Teuchos::RCP<Epetra_Vector> disp_interface =
      FluidField()->Interface()->ExtractFSICondVector(FluidField()->Dispnp());
  mortar_adapter_->IntegrateLinD(
      "displacement", disp_interface, StructureToFluidAtInterface(lambdanp_));
  tmp_k_D = mortar_adapter_->GetDMatrix();

  // fill off diagonal blocks
  {
    // create the parameters for the discretization
    Teuchos::ParameterList params;
    // action for elements
    params.set<int>("action", FLD::poro_splitnopenetration_OD);
    params.set("total time", Time());
    params.set("delta time", Dt());
    params.set("timescale", FluidField()->ResidualScaling());
    params.set<int>("Physical Type", FluidField()->PhysicalType());

    FluidField()->Discretization()->ClearState();
    FluidField()->Discretization()->SetState(0, "dispnp", FluidField()->Dispnp());
    FluidField()->Discretization()->SetState(0, "gridv", FluidField()->GridVel());
    FluidField()->Discretization()->SetState(0, "velnp", FluidField()->Velnp());
    FluidField()->Discretization()->SetState(0, "scaaf", FluidField()->Scaaf());

    FluidField()->Discretization()->SetState(0, "lambda",
        FluidField()->Interface()->InsertFSICondVector(StructureToFluidAtInterface(lambdanp_)));

    // build specific assemble strategy for the fluid-mechanical system matrix
    // from the point of view of FluidField:
    // fluiddofset = 0, structdofset = 1
    DRT::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
        1,                                  // structdofset for column
        k_struct_,                          // fluid-mechanical matrix
        k_lambda_, Teuchos::null, Teuchos::null, Teuchos::null);
    FluidField()->Discretization()->EvaluateCondition(params, fluidstrategy, "FSICoupling");

    FluidField()->Discretization()->ClearState();
  }

  // fill off diagonal blocks
  {
    // create the parameters for the discretization
    Teuchos::ParameterList params;
    // action for elements
    params.set<int>("action", FLD::poro_splitnopenetration_ODdisp);
    params.set("total time", Time());
    params.set("delta time", Dt());
    params.set("timescale", FluidField()->ResidualScaling());
    params.set<int>("Physical Type", FluidField()->PhysicalType());

    FluidField()->Discretization()->ClearState();
    FluidField()->Discretization()->SetState(0, "dispnp", FluidField()->Dispnp());
    FluidField()->Discretization()->SetState(0, "gridv", FluidField()->GridVel());
    FluidField()->Discretization()->SetState(0, "velnp", FluidField()->Velnp());
    //  FluidField()->Discretization()->SetState(0,"scaaf",FluidField()->Scaaf());

    // build specific assemble strategy for the fluid-mechanical system matrix
    // from the point of view of FluidField:
    // fluiddofset = 0, structdofset = 1
    DRT::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
        1,                                  // structdofset for column
        k_porodisp_,                        // fluid-mechanical matrix
        Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
    FluidField()->Discretization()->EvaluateCondition(params, fluidstrategy, "FSICoupling");

    FluidField()->Discretization()->ClearState();
  }

  // fill off diagonal blocks
  {
    // create the parameters for the discretization
    Teuchos::ParameterList params;
    // action for elements
    params.set<int>("action", FLD::poro_splitnopenetration_ODpres);
    params.set("total time", Time());
    params.set("delta time", Dt());
    params.set("timescale", FluidField()->ResidualScaling());
    params.set<int>("Physical Type", FluidField()->PhysicalType());

    FluidField()->Discretization()->ClearState();
    FluidField()->Discretization()->SetState(0, "dispnp", FluidField()->Dispnp());
    FluidField()->Discretization()->SetState(0, "gridv", FluidField()->GridVel());
    FluidField()->Discretization()->SetState(0, "velnp", FluidField()->Velnp());
    //  FluidField()->Discretization()->SetState(0,"scaaf",FluidField()->Scaaf());

    // build specific assemble strategy for the fluid-mechanical system matrix
    // from the point of view of FluidField:
    // fluiddofset = 0, structdofset = 1
    DRT::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
        0,                                  // fluiddofset for column
        k_porofluid_,                       // fluid-mechanical matrix
        Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
    FluidField()->Discretization()->EvaluateCondition(params, fluidstrategy, "FSICoupling");

    FluidField()->Discretization()->ClearState();
  }

  // Complete Coupling matrices which should be *.Add later!
  k_struct_->Complete(
      *StructureField()->Interface()->FSICondMap(), *FluidField()->Interface()->FSICondMap());
  k_fluid_->Complete();
  k_porofluid_->Complete();
  k_porodisp_->Complete();

  //------------------------------Invert D Matrix!-----------------------------------------------
  tmp_k_D->Complete();
  Teuchos::RCP<LINALG::SparseMatrix> invd =
      Teuchos::rcp(new LINALG::SparseMatrix(*tmp_k_D, LINALG::Copy));
  // invd->Complete();

  Teuchos::RCP<Epetra_Vector> diag =
      LINALG::CreateVector(*FluidField()->Interface()->FSICondMap(), true);

  int err = 0;

  // extract diagonal of invd into diag
  invd->ExtractDiagonalCopy(
      *diag);  // That the Reason, why tmp_k_D has to have Fluid Maps for Rows & Columns!!!

  // set zero diagonal values to dummy 1.0 ??
  for (int i = 0; i < diag->MyLength(); ++i)
    if ((*diag)[i] == 0.0)
    {
      (*diag)[i] = 1.0;
      std::cout << "--- --- --- WARNING: D-Matrix Diagonal Element " << i
                << " is zero!!! --- --- ---" << std::endl;
    }

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err > 0) dserror("ERROR: Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd->ReplaceDiagonalValues(*diag);
  invd->Complete();
  //------------------------------End of Invert D
  // Matrix!-----------------------------------------------

  // Transform also colum map of D-Matrix
  (*k_D_transform_)(*FluidField()->Interface()->FSICondMap(),
      FluidField()->BlockSystemMatrix()->Matrix(1, 1).ColMap(), *tmp_k_D, 1.0,
      ADAPTER::CouplingSlaveConverter(*icoupfs_), *k_D_);

  (*k_invD_transform_)(*invd, 1.0, ADAPTER::CouplingSlaveConverter(*icoupfs_), *k_invD_, false);

  double stiparam = StructureField()->TimIntParam();

  Teuchos::RCP<LINALG::SparseMatrix> tmp_k_DLin = mortar_adapter_->DLinMatrix();
  tmp_k_DLin->Complete();

  // Transform also column map of D-Matrix
  (*k_DLin_transform_)(*FluidField()->Interface()->FSICondMap(),
      FluidField()->BlockSystemMatrix()->Matrix(1, 1).ColMap(), *tmp_k_DLin,
      1.0 - stiparam,  // *= b
      ADAPTER::CouplingSlaveConverter(*icoupfs_),
      (Teuchos::rcp_static_cast<LINALG::BlockSparseMatrixBase>(k_fs))->Matrix(1, 1), true, true);

  k_lambda_->Complete(
      *StructureField()->Interface()->FSICondMap(), *FluidField()->Interface()->FSICondMap());
  k_invD_->Complete(
      *FluidField()->Interface()->FSICondMap(), *StructureField()->Interface()->FSICondMap());

  // Calculate 1/b*Tangent*invD
  k_lambdainvD_ = LINALG::MLMultiply(*k_lambda_, *k_invD_, true);
  k_lambdainvD_->Scale(1.0 / (1.0 - stiparam));  // *= 1/b
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate mechanical-fluid system matrix at state        vuong 09/14    |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::ApplyStrCouplMatrix(
    Teuchos::RCP<LINALG::SparseOperator> k_sf  //!< off-diagonal tangent matrix term
)
{
  // call base class
  Monolithic::ApplyStrCouplMatrix(k_sf);

  // done.
  return;
}

/*----------------------------------------------------------------------*
 | RecoverLagrangeMultiplier (protected)              vuong 09/14        |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::RecoverLagrangeMultiplierAfterTimeStep()
{
  // we do not need to recover after a time step, it is done after every newton step
  return;
}

/*----------------------------------------------------------------------*
 | RecoverLagrangeMultiplier (protected)              vuong 09/14        |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::Update()
{
  // call base class
  MonolithicSplit::Update();

  // update lagrangean multiplier
  lambda_->Update(1.0, *lambdanp_, 0.0);

  // copy D matrix from current time step to old D matrix
  k_Dn_ = Teuchos::rcp(
      new LINALG::SparseMatrix(*k_D_, LINALG::Copy));  // store D-Matrix from last timestep
  return;
}

/*----------------------------------------------------------------------*
| RecoverLagrangeMultiplier (protected)               ager 09/14        |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::Output(bool forced_writerestart)
{
  // call base class
  MonolithicSplit::Output(forced_writerestart);

  // for now, we always write the lagrange multiplier
  Teuchos::RCP<Epetra_Vector> fulllambda =
      Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(*StructureField()->DofRowMap()));
  LINALG::Export(*lambdanp_, *fulllambda);
  StructureField()->DiscWriter()->WriteVector("poronopencond_lambda", fulllambda);
}  // MonolithicSplitNoPenetration::Output()

/*----------------------------------------------------------------------*
 | RecoverLagrangeMultiplier (protected)              vuong 09/14        |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::SetupCouplingAndMatrices()
{
  const int ndim = DRT::Problem::Instance()->NDim();
  icoupfs_->SetupConditionCoupling(*StructureField()->Discretization(),
      StructureField()->Interface()->FSICondMap(), *FluidField()->Discretization(),
      FluidField()->Interface()->FSICondMap(), "FSICoupling", ndim);

  evaluateinterface_ = false;

  // initialize Poroelasticity-systemmatrix_
  systemmatrix_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      *Extractor(), *Extractor(), 81, false, true));

  // initialize coupling matrices
  k_fs_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      *(StructureField()->Interface()), *(FluidField()->Interface()), 81, false, true));

  k_sf_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      *(FluidField()->Interface()), *(StructureField()->Interface()), 81, false, true));

  // initialize no penetration coupling matrices
  k_struct_ = Teuchos::rcp(
      new LINALG::SparseMatrix(*(FluidField()->Interface()->FSICondMap()), 81, true, true));

  k_fluid_ = Teuchos::rcp(
      new LINALG::SparseMatrix(*(FluidField()->Interface()->FSICondMap()), 81, false, false));

  k_lambda_ = Teuchos::rcp(
      new LINALG::SparseMatrix(*(FluidField()->Interface()->FSICondMap()), 81, true, true));

  k_D_ = Teuchos::rcp(
      new LINALG::SparseMatrix(*(FluidField()->Interface()->FSICondMap()), 81, true, true));

  k_invD_ = Teuchos::rcp(
      new LINALG::SparseMatrix(*(StructureField()->Interface()->FSICondMap()), 81, true, true));

  k_porodisp_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      *(StructureField()->Interface()), *(FluidField()->Interface()), 81, true, true));

  k_porofluid_ =
      Teuchos::rcp(new LINALG::SparseMatrix(*(FluidField()->DofRowMap()), 81, true, true));

  nopenetration_rhs_ =
      Teuchos::rcp(new Epetra_Vector(*FluidField()->Interface()->FSICondMap(), true));

  return;
}


/*----------------------------------------------------------------------*
 | RecoverLagrangeMultiplier (protected)              vuong 09/14        |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::PrepareTimeStep()
{
  // call base class
  POROELAST::Monolithic::PrepareTimeStep();
}

/*----------------------------------------------------------------------*
 | Read Restart (public)                          vuong 09/14        |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::ReadRestart(const int step)
{
  // call base class
  POROELAST::PoroBase::ReadRestart(step);

  // get lagrange multiplier and D matrix
  if (step)
  {
    // get the structure reader (this is where the lagrange multiplier was saved)
    IO::DiscretizationReader reader(StructureField()->Discretization(), StructureField()->Step());
    Teuchos::RCP<Epetra_Vector> fulllambda =
        Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(*StructureField()->DofRowMap()));

    // this is the lagrange multiplier on the whole structure field
    reader.ReadVector(fulllambda, "poronopencond_lambda");

    // extract lambda on fsi interface vector
    lambda_ = StructureField()->Interface()->ExtractFSICondVector(fulllambda);
    lambdanp_->Update(1.0, *lambda_, 0.0);

    // call an additional evaluate to get the old D matrix
    SetupSystem();
    // call evaluate to recalculate D matrix
    Evaluate(zeros_,false);

    // copy D matrix from current time step to old D matrix
    k_Dn_ = Teuchos::rcp(
        new LINALG::SparseMatrix(*k_D_, LINALG::Copy));  // store D-Matrix from last timestep
  }
}

/*----------------------------------------------------------------------*
 *                                                    vuong 11/14        |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::PrintNewtonIterHeaderStream(std::ostringstream& oss)
{
  Monolithic::PrintNewtonIterHeaderStream(oss);

  oss << std::setw(20) << "abs-crhs-res";

  return;
}

/*----------------------------------------------------------------------*
 *                                                    vuong 11/14        |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::PrintNewtonIterTextStream(std::ostringstream& oss)
{
  Monolithic::PrintNewtonIterTextStream(oss);

  oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhs_nopenetration_;
}

/*----------------------------------------------------------------------*
 *                                                    vuong 11/14        |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::BuildConvergenceNorms()
{
  Monolithic::BuildConvergenceNorms();

  normrhs_nopenetration_ = UTILS::CalculateVectorNorm(vectornormfres_, nopenetration_rhs_);

  return;
}
