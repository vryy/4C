/*----------------------------------------------------------------------*/
/*!

 \brief  monolithic structure split poroelasticity algorithm

\level 2

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
 *------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/

#include "poro_monolithicstructuresplit.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_adapter/ad_str_fpsiwrapper.H"
#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/adapter_coupling.H"

#include "../drt_fsi/fsi_matrixtransform.H"

#include "../drt_fluid/fluid_utils_mapextractor.H"

#include "../drt_structure/stru_aux.H"

#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |                                                         vuong 11/12  |
 *----------------------------------------------------------------------*/
POROELAST::MonolithicStructureSplit::MonolithicStructureSplit(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : MonolithicSplit(comm, timeparams)
{
  sggtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform);
  sgitransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  sigtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
  csggtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  cfggtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
  csgitransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  cfigtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);

  // Recovering of Lagrange multiplier happens on structure field
  lambda_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->FSICondMap()));

  return;
}

/*----------------------------------------------------------------------*
 | setup system (called in porolast.cpp)                 vuong 11/12    |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicStructureSplit::SetupSystem()
{
  {
    // create combined map
    std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;

    vecSpaces.push_back(StructureField()->Interface()->OtherMap());
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

  // Use splitted structure matrix
  StructureField()->UseBlockMatrix();

  SetupCouplingAndMatrices();

  BuildCombinedDBCMap();
}  // SetupSystem()

/*----------------------------------------------------------------------*
 |                                                         vuong 11/12  |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicStructureSplit::SetupRHS(bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicStructureSplit::SetupRHS");

  // create full monolithic rhs vector
  rhs_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  SetupVector(*rhs_, StructureField()->RHS(), FluidField()->RHS(), FluidField()->ResidualScaling());

// i am not sure, if this is needed
#if 0
  if (firstcall and evaluateinterface_)
  {
     // get time integration parameters of structure an fluid time integrators
    // to enable consistent time integration among the fields
    double stiparam = StructureField()->TimIntParam();
    double ftiparam = FluidField()->TimIntParam();

    Teuchos::RCP<Epetra_Vector> fveln = FluidField()->ExtractInterfaceVeln();
    Teuchos::RCP<Epetra_Vector> sveln = FluidToStructureAtInterface(fveln);

    // additional rhs term for structure equations
    Teuchos::RCP<Epetra_Vector> veln = StructureField()->Interface()->InsertFSICondVector(sveln);
    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(veln->Map()));

    Teuchos::RCP<LINALG::BlockSparseMatrixBase> s = StructureField()->BlockSystemMatrix();
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> k_fs = FluidStructCouplingBlockMatrix();
    if (k_fs==Teuchos::null)
      dserror("expect coupling block matrix");

    s->Apply(*veln,*rhs);

    rhs->Scale(-1.*Dt());

    veln = StructureField()->Interface()->ExtractOtherVector(rhs); // only inner DOFs
    Extractor()->AddVector(*veln,0,*rhs_); // add inner structure contributions to 'f'

    veln = StructureField()->Interface()->ExtractFSICondVector(rhs); // only DOFs on interface
    veln = FluidField()->Interface()->InsertFSICondVector(StructureToFluidAtInterface(veln)); // convert to fluid map

    double scale     = FluidField()->ResidualScaling();

    veln->Scale((1.0-ftiparam)/((1.0-stiparam)*scale));

    Extractor()->AddVector(*veln,1,*rhs_); // add fluid contribution to 'f'

    LINALG::SparseMatrix& kig = k_fs->Matrix(0,1);
    LINALG::SparseMatrix& kgg = k_fs->Matrix(1,1);

    rhs = Teuchos::rcp(new Epetra_Vector(kig.RowMap()));
    kig.Apply(*fveln,*rhs);
    //veln = FluidField()->Interface()->InsertOtherVector(rhs);
    rhs->Scale(-1.*Dt());

    rhs = FluidField()->Interface()->InsertOtherVector(rhs);

    Extractor()->AddVector(*rhs,1,*rhs_);

    rhs = Teuchos::rcp(new Epetra_Vector(kgg.RowMap()));
    kgg.Apply(*fveln,*rhs);
    //FluidField()->Interface()->InsertFSICondVector(rhs,veln);

    rhs->Scale(-1.*Dt());

    rhs = FluidField()->Interface()->InsertFSICondVector(rhs);

    // veln->Scale(-1.*Dt());

    Extractor()->AddVector(*rhs,1,*rhs_);

    // Reset quantities for previous iteration step since they still store values from the last time step
    //ddiinc_ = LINALG::CreateVector(*StructureField()->Interface()->OtherMap(),true);
    //solipre_ = Teuchos::null;
    //ddginc_ = LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(),true);
    //solgpre_ = Teuchos::null;
    //duiinc_ = LINALG::CreateVector(*FluidField()->Interface()->OtherMap(),true);
    //solivelpre_ = Teuchos::null;
    //fgcur_ = LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(),true);
    //sgicur_ = Teuchos::null;
    //sggcur_ = Teuchos::null;
    //cgicur_ = Teuchos::null;
    //cggcur_ = Teuchos::null;
  }
#endif

  // store interface force onto the structure to know it in the next time step as previous force
  // in order to recover the Lagrange multiplier
  // fgpre_ = fgcur_;
  fgcur_ = StructureField()->Interface()->ExtractFSICondVector(StructureField()->RHS());
}


/*----------------------------------------------------------------------*
 |                                                         vuong 11/12  |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicStructureSplit::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicStructureSplit::SetupSystemMatrix");

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> s = StructureField()->BlockSystemMatrix();
  // Teuchos::RCP<LINALG::SparseMatrix> s = StructureField()->SystemMatrix();
  if (s == Teuchos::null) dserror("expect structure block matrix");
  Teuchos::RCP<LINALG::SparseMatrix> f = FluidField()->SystemMatrix();
  if (f == Teuchos::null) dserror("expect fluid matrix");

  // just to play it save ...
  mat.Zero();

  /*----------------------------------------------------------------------*/

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  /*----------------------------------------------------------------------*/
  // structural part k_sf (3nxn)
  // build mechanical-fluid block

  // create empty matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> k_sf = StructFluidCouplingBlockMatrix();
  if (k_sf == Teuchos::null) dserror("expect coupling block matrix");

  // call the element and calculate the matrix block
  ApplyStrCouplMatrix(k_sf);

  /*----------------------------------------------------------------------*/
  // fluid part k_fs ( (3n+1)x3n )
  // build fluid-mechanical block

  // create empty matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> k_fs = FluidStructCouplingBlockMatrix();
  if (k_fs == Teuchos::null) dserror("expect coupling block matrix");

  // call the element and calculate the matrix block
  ApplyFluidCouplMatrix(k_fs);

  /*----------------------------------------------------------------------*/

  k_fs->Complete();
  k_sf->Complete();

  f->UnComplete();

  /*----------------------------------------------------------------------*/
  if (evaluateinterface_)
  {
    double scale = FluidField()->ResidualScaling();
    double timescale = FluidField()->TimeScaling();

    // get time integration parameters of structure an fluid time integrators
    // to enable consistent time integration among the fields
    double stiparam = StructureField()->TimIntParam();
    double ftiparam = FluidField()->TimIntParam();

    (*sigtransform_)(s->FullRowMap(), s->FullColMap(), s->Matrix(0, 1), 1. / timescale,
        ADAPTER::CouplingMasterConverter(*icoupfs_),
        k_sf->Matrix(0, 1),  // k_sf->Matrix(0,1),mat.Matrix(0,1)
        true, true);

    (*sggtransform_)(s->Matrix(1, 1), (1.0 - ftiparam) / ((1.0 - stiparam) * scale * timescale),
        ADAPTER::CouplingMasterConverter(*icoupfs_), ADAPTER::CouplingMasterConverter(*icoupfs_),
        *f, true, true);

    (*sgitransform_)(s->Matrix(1, 0), (1.0 - ftiparam) / ((1.0 - stiparam) * scale),
        ADAPTER::CouplingMasterConverter(*icoupfs_),
        k_fs->Matrix(1, 0),  // k_fs->Matrix(1,0), mat.Matrix(1,0)
        true);

    (*cfggtransform_)(s->FullRowMap(),  // k_fs->FullRowMap(),
        s->FullColMap(),                // k_fs->FullColMap(),
        k_fs->Matrix(1, 1), 1. / timescale, ADAPTER::CouplingMasterConverter(*icoupfs_), *f, true,
        true);

    (*cfigtransform_)(s->FullRowMap(),  // k_fs->FullRowMap(),
        s->FullColMap(),                // k_fs->FullColMap(),
        k_fs->Matrix(0, 1), 1. / timescale, ADAPTER::CouplingMasterConverter(*icoupfs_), *f, true,
        true);

    (*csggtransform_)(k_sf->Matrix(1, 1), (1.0 - ftiparam) / ((1.0 - stiparam) * scale),
        ADAPTER::CouplingMasterConverter(*icoupfs_), *f, true);

    (*csgitransform_)(k_sf->Matrix(1, 0), (1.0 - ftiparam) / ((1.0 - stiparam) * scale),
        ADAPTER::CouplingMasterConverter(*icoupfs_), *f, true);
  }

  /*----------------------------------------------------------------------*/
  // pure structural part
  mat.Matrix(0, 0).Add(s->Matrix(0, 0), false, 1., 0.0);

  // structure coupling part
  mat.Matrix(0, 1).Add(k_sf->Matrix(0, 0), false, 1.0, 0.0);
  mat.Matrix(0, 1).Add(k_sf->Matrix(0, 1), false, 1.0, 1.0);

  // pure fluid part
  // uncomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  // this just once...
  // f->UnComplete();
  mat.Assign(1, 1, LINALG::View, *f);

  // fluid coupling part
  mat.Matrix(1, 0).Add(k_fs->Matrix(0, 0), false, 1.0, 0.0);
  mat.Matrix(1, 0).Add(k_fs->Matrix(1, 0), false, 1.0, 1.0);
  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();

  // store parts of structural matrix to know them in the next iteration as previous iteration
  // matrices
  // sgipre_ = sgicur_;
  // sggpre_ = sggcur_;
  // cgipre_ = cgicur_;
  // cggpre_ = cggcur_;
  sgicur_ = Teuchos::rcp(new LINALG::SparseMatrix(s->Matrix(1, 0)));
  sggcur_ = Teuchos::rcp(new LINALG::SparseMatrix(s->Matrix(1, 1)));
  cgicur_ = Teuchos::rcp(new LINALG::SparseMatrix(k_sf->Matrix(1, 0)));
  cggcur_ = Teuchos::rcp(new LINALG::SparseMatrix(k_sf->Matrix(1, 1)));
}

/*----------------------------------------------------------------------*
 |                                                         vuong 11/12  |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicStructureSplit::SetupVector(Epetra_Vector& f,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv, double fluidscale)
{
  // extract the inner and boundary dofs of all three fields

  Teuchos::RCP<Epetra_Vector> sov = StructureField()->Interface()->ExtractOtherVector(sv);

  if (fluidscale != 0.0)
  {
    // get time integration parameters of structure an fluid time integrators
    // to enable consistent time integration among the fields
    double stiparam = StructureField()->TimIntParam();
    double ftiparam = FluidField()->TimIntParam();

    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> scv = StructureField()->Interface()->ExtractFSICondVector(sv);
    Teuchos::RCP<Epetra_Vector> modfv =
        FluidField()->Interface()->InsertFSICondVector(StructureToFluidAtInterface(scv));
    modfv->Update(1.0, *fv, (1.0 - ftiparam) / ((1.0 - stiparam) * fluidscale));

    // add contribution of Lagrange multiplier from previous time step
    if (lambda_ != Teuchos::null)
      modfv->Update(-ftiparam + stiparam * (1.0 - ftiparam) / (1.0 - stiparam),
          *StructureToFluidAtInterface(lambda_), 1.0);

    Extractor()->InsertVector(*modfv, 1, f);
  }
  else
  {
    Extractor()->InsertVector(*fv, 1, f);
  }

  Extractor()->InsertVector(*sov, 0, f);
}

/*----------------------------------------------------------------------*
 | extract field vectors for calling Evaluate() of the       vuong 01/12|
 | single fields                                                        |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicStructureSplit::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicStructureSplit::ExtractFieldVectors");

  // process fluid unknowns of the second field
  fx = Extractor()->ExtractVector(x, 1);

  // process structure unknowns
  if (evaluateinterface_)
  {
    Teuchos::RCP<Epetra_Vector> fcx = FluidField()->Interface()->ExtractFSICondVector(fx);

    // if(firstcall)
    // FluidField()->VelocityToDisplacement(fcx);
    // else
    {
      double timescale = 1. / FluidField()->TimeScaling();
      fcx->Scale(timescale);
    }

    Teuchos::RCP<Epetra_Vector> scx = FluidToStructureAtInterface(fcx);
    Teuchos::RCP<const Epetra_Vector> sox = Extractor()->ExtractVector(x, 0);

    Teuchos::RCP<Epetra_Vector> s = StructureField()->Interface()->InsertOtherVector(sox);
    StructureField()->Interface()->InsertFSICondVector(scx, s);

    Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(new const Epetra_Vector(s->Map(), true));
    LINALG::ApplyDirichlettoSystem(s, zeros, *(StructureField()->GetDBCMapExtractor()->CondMap()));

    sx = s;

    Teuchos::RCP<Epetra_Vector> fox = FluidField()->Interface()->ExtractOtherVector(fx);

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

/*----------------------------------------------------------------------*/
/* Recover the Lagrange multiplier at the interface   vuong 01/12 */
/*----------------------------------------------------------------------*/
void POROELAST::MonolithicStructureSplit::RecoverLagrangeMultiplierAfterTimeStep()
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "POROELAST::MonolithicStructureSplit::RecoverLagrangeMultiplierAfterTimeStep");
  if (evaluateinterface_)
  {
    // get time integration parameter of structural time integrator
    // to enable consistent time integration among the fields
    double stiparam = StructureField()->TimIntParam();
    double timescale = FluidField()->TimeScaling();

    // store the product S_{\GammaI} \Delta d_I^{n+1} in here
    Teuchos::RCP<Epetra_Vector> sgiddi =
        LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    // sgipre_->Multiply(false, *ddiinc_, *sgiddi);
    sgicur_->Multiply(false, *ddiinc_, *sgiddi);
    //(sgipre_->EpetraMatrix())->Multiply(false, *ddiinc_, *sgiddi);

    // store the product C_{\GammaI} \Delta u_I^{n+1} in here
    Teuchos::RCP<Epetra_Vector> fgiddi =
        LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    // cgipre_->Multiply(false, *duiinc_, *fgiddi);
    cgicur_->Multiply(false, *duiinc_, *fgiddi);

    // store the product S_{\Gamma\Gamma} \Delta d_\Gamma^{n+1} in here
    Teuchos::RCP<Epetra_Vector> sggddg =
        LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    // sggpre_->Multiply(false, *ddginc_, *sggddg);
    sggcur_->Multiply(false, *ddginc_, *sggddg);

    // store the prodcut C_{\Gamma\Gamma} \Delta u_\Gamma^{n+1} in here
    Teuchos::RCP<Epetra_Vector> cggddg =
        LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    // cggpre_->Multiply(false, *ddginc_, *cggddg);
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
  return;
}
