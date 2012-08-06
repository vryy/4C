#include <Teuchos_TimeMonitor.hpp>

#include "../drt_adapter/adapter_coupling_mortar.H"
#include "../drt_adapter/adapter_coupling.H"

#include "fsi_mortarmonolithic_fluidsplit.H"
#include "fsi_debugwriter.H"
#include "fsi_statustest.H"
#include "fsi_overlapprec_fsiamg.H"
#include "fsi_monolithic_linearsystem.H"
#include "fsi_matrixtransform.H"
#include "fsi_utils.H"

#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_structure/stru_aux.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../drt_ale/ale_utils_mapextractor.H"

#include "../drt_constraint/constraint_manager.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#define FLUIDSPLITAMG

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MortarMonolithicFluidSplit::MortarMonolithicFluidSplit(const Epetra_Comm& comm,
                                                            const Teuchos::ParameterList& timeparams)
  : BlockMonolithic(comm,timeparams),
    comm_(comm)
{
  notsetup_ = true;

  coupsfm_ = Teuchos::rcp(new ADAPTER::CouplingMortar());
  icoupfa_ = Teuchos::rcp(new ADAPTER::Coupling());
  fscoupfa_ = Teuchos::rcp(new ADAPTER::Coupling());

  aigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmiitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplit::SetupSystem()
{
  if (notsetup_)
  {
    const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
    linearsolverstrategy_ = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");

    aleproj_ = DRT::INPUT::IntegralValue<INPAR::FSI::SlideALEProj>(fsidyn,"SLIDEALEPROJ");

    SetDefaultParameters(fsidyn,NOXParameterList());

    // we use non-matching meshes at the interface
    // mortar with: structure = master, fluid = slave

    const int ndim = DRT::Problem::Instance()->NDim();

    // structure to fluid

    coupsfm_->Setup(*StructureField()->Discretization(),
                    *FluidField().Discretization(),
                    *AleField().Discretization(),
                    comm_,false);

    // fluid to ale at the interface

    icoupfa_->SetupConditionCoupling(*FluidField().Discretization(),
                                     FluidField().Interface()->FSICondMap(),
                                     *AleField().Discretization(),
                                     AleField().Interface()->FSICondMap(),
                                     "FSICoupling",
                                     ndim);

    // we might have a free surface
    if (FluidField().Interface()->FSCondRelevant())
    {
      fscoupfa_->SetupConditionCoupling(*FluidField().Discretization(),
                                        FluidField().Interface()->FSCondMap(),
                                        *AleField().Discretization(),
                                        AleField().Interface()->FSCondMap(),
                                        "FREESURFCoupling",
                                        ndim);
    }

    ADAPTER::Coupling& coupfa = FluidAleCoupling();

    // the fluid-ale coupling always matches
    const Epetra_Map* fluidnodemap = FluidField().Discretization()->NodeRowMap();
    const Epetra_Map* alenodemap   = AleField().Discretization()->NodeRowMap();

    coupfa.SetupCoupling(*FluidField().Discretization(),
                         *AleField().Discretization(),
                         *fluidnodemap,
                         *alenodemap,
                          ndim);

    FluidField().SetMeshMap(coupfa.MasterDofMap());

    // create combined map

    std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;
    vecSpaces.push_back(StructureField()->DofRowMap());
  #ifdef FLUIDSPLITAMG
    vecSpaces.push_back(FluidField()    .DofRowMap());
  #else
    vecSpaces.push_back(FluidField()    .Interface().OtherMap());
  #endif
    vecSpaces.push_back(AleField()      .Interface()->OtherMap());

    if (vecSpaces[1]->NumGlobalElements()==0)
      dserror("No inner fluid equations. Splitting not possible.");

    SetDofRowMaps(vecSpaces);

    /*----------------------------------------------------------------------*/
    // Switch fluid to interface split block matrix
    FluidField().UseBlockMatrix(true);

    // build ale system matrix in splitted system
    AleField().BuildSystemMatrix(false);

    aleresidual_ = Teuchos::rcp(new Epetra_Vector(*AleField().Interface()->OtherMap()));

    vector<int> pciter;
    vector<double> pcomega;
    vector<int> spciter;
    vector<double> spcomega;
    vector<int> fpciter;
    vector<double> fpcomega;
    vector<int> apciter;
    vector<double> apcomega;
    vector<string> blocksmoother;
    vector<double> schuromega;
    {
      int    word1;
      double word2;
      {
        std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsidyn,"PCITER"));
        std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsidyn,"PCOMEGA"));
        while (pciterstream >> word1)
          pciter.push_back(word1);
        while (pcomegastream >> word2)
          pcomega.push_back(word2);
      }
      {
        std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsidyn,"STRUCTPCITER"));
        std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsidyn,"STRUCTPCOMEGA"));
        while (pciterstream >> word1)
          spciter.push_back(word1);
        while (pcomegastream >> word2)
          spcomega.push_back(word2);
      }
      {
        std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsidyn,"FLUIDPCITER"));
        std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsidyn,"FLUIDPCOMEGA"));
        while (pciterstream >> word1)
          fpciter.push_back(word1);
        while (pcomegastream >> word2)
          fpcomega.push_back(word2);
      }
      {
        std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsidyn,"ALEPCITER"));
        std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsidyn,"ALEPCOMEGA"));
        while (pciterstream >> word1)
          apciter.push_back(word1);
        while (pcomegastream >> word2)
          apcomega.push_back(word2);
      }
      {
        string word;
        std::istringstream blocksmootherstream(Teuchos::getNumericStringParameter(fsidyn,"BLOCKSMOOTHER"));
        while (blocksmootherstream >> word)
          blocksmoother.push_back(word);
      }
      {
        std::istringstream blocksmootherstream(Teuchos::getNumericStringParameter(fsidyn,"SCHUROMEGA"));
        while (blocksmootherstream >> word2)
          schuromega.push_back(word2);
      }
    }

    // enable debugging
    if (DRT::INPUT::IntegralValue<int>(fsidyn,"DEBUGOUTPUT") & 2)
    {
      pcdbg_ = Teuchos::rcp(new UTILS::MonolithicDebugWriter(*this));
    }

    // create block system matrix
    switch(linearsolverstrategy_)
    {
    case INPAR::FSI::PreconditionedKrylov:
    case INPAR::FSI::FSIAMG:
      systemmatrix_ = Teuchos::rcp(new OverlappingBlockMatrixFSIAMG(
                                     Extractor(),
                                     *StructureField(),
                                     FluidField(),
                                     AleField(),
                                     false,
                                     DRT::INPUT::IntegralValue<int>(fsidyn,"SYMMETRICPRECOND"),
                                     blocksmoother,
                                     schuromega,
                                     pcomega,
                                     pciter,
                                     spcomega,
                                     spciter,
                                     fpcomega,
                                     fpciter,
                                     apcomega,
                                     apciter,
                                     DRT::INPUT::IntegralValue<int>(fsidyn,"FSIAMGANALYZE"),
                                     linearsolverstrategy_,
                                     DRT::Problem::Instance()->ErrorFile()->Handle()));
      break;
    default:
      dserror("Unsupported type of monolithic solver");
      break;
    }

    if(aleproj_ != INPAR::FSI::ALEprojection_none)
    {
      // set up sliding ale utils
      slideale_ = rcp(new FSI::UTILS::SlideAleUtils(StructureField()->Discretization(),
                                                    FluidField().Discretization(),
                                                    *coupsfm_,
                                                    true,
                                                    aleproj_));

      iprojdispinc_ = Teuchos::rcp(new Epetra_Vector(*coupsfm_->SlaveDofRowMap(),true));
      iprojdisp_ = Teuchos::rcp(new Epetra_Vector(*coupsfm_->SlaveDofRowMap(),true));
    }
    notsetup_=false;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplit::SetupRHS(Epetra_Vector& f, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::SetupRHS");

  SetupVector(f,
              StructureField()->RHS(),
              FluidField().RHS(),
              AleField().RHS(),
              FluidField().ResidualScaling());

  // add additional ale residual
  Extractor().AddVector(*aleresidual_,2,f);

  if (firstcall)
  {
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockf = FluidField().BlockSystemMatrix();


    Teuchos::RCP<Epetra_Vector> rhs;

    if (aleproj_!= INPAR::FSI::ALEprojection_none)
    {
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField().BlockSystemMatrix();

      rhs = Teuchos::rcp(new Epetra_Vector(a->Matrix(0,1).RowMap()));
      a->Matrix(0,1).Apply(*icoupfa_->MasterToSlave(iprojdispinc_),*rhs);

      Extractor().AddVector(*rhs,2,f);

      Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField().ShapeDerivatives();
      if (mmm!=Teuchos::null)
      {
        LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);
        LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);

        rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RowMap()));

        fmgg.Apply(*iprojdispinc_,*rhs);

        Teuchos::RCP<LINALG::SparseMatrix> mortar = coupsfm_->GetMortarTrafo();

        Teuchos::RCP<Epetra_Vector> tmprhs = Teuchos::rcp(new Epetra_Vector(mortar->DomainMap()));
        mortar->Multiply(true,*rhs,*tmprhs);

        rhs = StructureField()->Interface()->InsertFSICondVector(tmprhs);

        Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(new const Epetra_Vector(rhs->Map(),true));
        LINALG::ApplyDirichlettoSystem(rhs,zeros,*(StructureField()->GetDBCMapExtractor()->CondMap()));

        if (StructureField()->GetSTCAlgo() == INPAR::STR::stc_currsym)
        {
          Teuchos::RCP<LINALG::SparseMatrix> stcmat = StructureField()->GetSTCMat();
          stcmat->Multiply(true,*rhs,*rhs);
        }

        Extractor().AddVector(*rhs,0,f);

        rhs = Teuchos::rcp(new Epetra_Vector(fmig.RowMap()));

        fmig.Apply(*iprojdispinc_,*rhs);

        #ifdef FLUIDSPLITAMG
          rhs = FluidField().Interface()->InsertOtherVector(rhs);
        #endif

        zeros = Teuchos::rcp(new const Epetra_Vector(rhs->Map(),true));
        LINALG::ApplyDirichlettoSystem(rhs,zeros,*(StructureField()->GetDBCMapExtractor()->CondMap()));

        Extractor().AddVector(*rhs,1,f);

      }
    }

    LINALG::SparseMatrix& fig = blockf->Matrix(0,1);
    LINALG::SparseMatrix& fgg = blockf->Matrix(1,1);

    Teuchos::RCP<Epetra_Vector> fveln = FluidField().ExtractInterfaceVeln();
    double timescale = FluidField().TimeScaling();
    double scale     = FluidField().ResidualScaling();

    rhs = Teuchos::rcp(new Epetra_Vector(fig.RowMap()));

    fig.Apply(*fveln,*rhs);
    rhs->Scale(timescale*Dt());

#ifdef FLUIDSPLITAMG
    rhs = FluidField().Interface()->InsertOtherVector(rhs);
#endif

    Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(new const Epetra_Vector(rhs->Map(),true));
    LINALG::ApplyDirichlettoSystem(rhs,zeros,*(StructureField()->GetDBCMapExtractor()->CondMap()));

    Extractor().AddVector(*rhs,1,f);

    rhs = Teuchos::rcp(new Epetra_Vector(fgg.RowMap()));

    fgg.Apply(*fveln,*rhs);
    rhs->Scale(scale*timescale*Dt());

    Teuchos::RCP<LINALG::SparseMatrix> mortar = coupsfm_->GetMortarTrafo();

    Teuchos::RCP<Epetra_Vector> tmprhs = Teuchos::rcp(new Epetra_Vector(mortar->DomainMap()));
    mortar->Multiply(true,*rhs,*tmprhs);

    rhs = StructureField()->Interface()->InsertFSICondVector(tmprhs);

    zeros = Teuchos::rcp(new const Epetra_Vector(rhs->Map(),true));
    LINALG::ApplyDirichlettoSystem(rhs,zeros,*(StructureField()->GetDBCMapExtractor()->CondMap()));

    if (StructureField()->GetSTCAlgo() == INPAR::STR::stc_currsym)
    {
      Teuchos::RCP<LINALG::SparseMatrix> stcmat = StructureField()->GetSTCMat();
      stcmat->Multiply(true,*rhs,*rhs);
    }

    Extractor().AddVector(*rhs,0,f);
  }

  // NOX expects a different sign here.
  f.Scale(-1.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplit::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::SetupSystemMatrix");

  // extract Jacobian matrices and put them into composite system
  // matrix W

  Teuchos::RCP<LINALG::SparseMatrix> mortar = coupsfm_->GetMortarTrafo();

  // get info about STC feature
  INPAR::STR::STC_Scale stcalgo = StructureField()->GetSTCAlgo();
  Teuchos::RCP<LINALG::SparseMatrix> stcmat = Teuchos::null;
  if (stcalgo != INPAR::STR::stc_none)
    stcmat = StructureField()->GetSTCMat();

  const ADAPTER::Coupling& coupfa = FluidAleCoupling();

  Teuchos::RCP<LINALG::SparseMatrix> s = StructureField()->SystemMatrix();

  // split fluid matrix

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> f = FluidField().BlockSystemMatrix();

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField().BlockSystemMatrix();

  if (a==Teuchos::null)
    dserror("expect ale block matrix");

  LINALG::SparseMatrix& aii = a->Matrix(0,0);
  LINALG::SparseMatrix& aig = a->Matrix(0,1);

  double scale     = FluidField().ResidualScaling();
  double timescale = FluidField().TimeScaling();

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  // uncomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  // this just once...
  s->UnComplete();

  RCP<LINALG::SparseMatrix> fgg = MLMultiply(f->Matrix(1,1),false,*mortar,false,false,false,true);
  fgg = MLMultiply(*mortar,true,*fgg,false,false,false,true);

  s->Add(*fgg,false,scale*timescale,1.0);

  RCP<LINALG::SparseMatrix> fgi = MLMultiply(*mortar,true,f->Matrix(1,0),false,false,false,true);
  RCP<LINALG::SparseMatrix> lfgi = rcp(new LINALG::SparseMatrix(s->RowMap(),81,false));

  lfgi->Add(*fgi,false,scale,0.0);
  lfgi->Complete(fgi->DomainMap(),s->RangeMap());

  lfgi->ApplyDirichlet( *(StructureField()->GetDBCMapExtractor()->CondMap()),false);

  if (stcalgo == INPAR::STR::stc_currsym)
    lfgi = LINALG::MLMultiply(*stcmat, true, *lfgi, false, true, true, true);

#ifdef FLUIDSPLITAMG
  mat.Matrix(0,1).UnComplete();
  mat.Matrix(0,1).Add(*lfgi,false,1.,0.0);
#else
  mat.Assign(0,1,View,*lfgi);
#endif

  RCP<LINALG::SparseMatrix> fig = MLMultiply(f->Matrix(0,1),false,*mortar,false,false,false,true);
  RCP<LINALG::SparseMatrix> lfig = rcp(new LINALG::SparseMatrix(fig->RowMap(),81,false));

  lfig->Add(*fig,false,timescale,0.0);
  lfig->Complete(s->DomainMap(),fig->RangeMap());

  lfig->ApplyDirichlet( *(FluidField().GetDBCMapExtractor()->CondMap()),false);

  if (stcalgo != INPAR::STR::stc_none)
  {
    lfig = LINALG::MLMultiply(*lfig,false,*stcmat, false, false, false,true);
  }

#ifdef FLUIDSPLITAMG
  mat.Matrix(1,0).UnComplete();
  mat.Matrix(1,0).Add(*lfig,false,1.,0.0);
#else
  mat.Assign(1,0,View,*lfig);
#endif

  LINALG::SparseMatrix& fii = f->Matrix(0,0);

#ifdef FLUIDSPLITAMG
  mat.Matrix(1,1).UnComplete();
  mat.Matrix(1,1).Add(fii,false,1.,0.0);
  Teuchos::RCP<LINALG::SparseMatrix> eye = LINALG::Eye(*FluidField().Interface()->FSICondMap());
  mat.Matrix(1,1).Add(*eye,false,1.,1.0);
#else
  mat.Assign(1,1,View,fii);
#endif

  RCP<LINALG::SparseMatrix> laig = rcp(new LINALG::SparseMatrix(aii.RowMap(),81,false));
  (*aigtransform_)(a->FullRowMap(),
                   a->FullColMap(),
                   aig,
                   1.,
                   ADAPTER::CouplingSlaveConverter(*icoupfa_),
                   *laig);

  laig->Complete(f->Matrix(1,1).DomainMap(),aii.RangeMap());
  RCP<LINALG::SparseMatrix> llaig = MLMultiply(*laig,false,*mortar,false,false,false,true);
  laig = rcp(new LINALG::SparseMatrix(llaig->RowMap(),81,false));

  laig->Add(*llaig,false,1.0,0.0);
  laig->Complete(s->DomainMap(),llaig->RangeMap());

  laig->ApplyDirichlet( *(AleField().GetDBCMapExtractor()->CondMap()),false);

  if (stcalgo != INPAR::STR::stc_none)
  {
    laig = LINALG::MLMultiply(*laig,false,*stcmat, false, false, false,true);
  }

  mat.Assign(2,0,View,*laig);

  mat.Assign(2,2,View,aii);

  /*----------------------------------------------------------------------*/
  // add optional fluid linearization with respect to mesh motion block

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField().ShapeDerivatives();
  if (mmm!=Teuchos::null)
  {
    LINALG::SparseMatrix& fmii = mmm->Matrix(0,0);
    LINALG::SparseMatrix& fmgi = mmm->Matrix(1,0);

    // reuse transform objects to add shape derivative matrices to structural blocks

    RCP<LINALG::SparseMatrix> fmgg = MLMultiply(mmm->Matrix(1,1),false,*mortar,false,false,false,true);
    fmgg = MLMultiply(*mortar,true,*fmgg,false,false,false,true);

    RCP<LINALG::SparseMatrix> lfmgg = rcp(new LINALG::SparseMatrix(fmgg->RowMap(),81,false));
    lfmgg->Add(*fmgg,false,1.0,0.0);
    lfmgg->Complete(s->DomainMap(),fmgg->RangeMap());

    s->Add(*lfmgg,false,scale,1.0);

    RCP<LINALG::SparseMatrix> fmig = MLMultiply(mmm->Matrix(0,1),false,*mortar,false,false,false,true);
    RCP<LINALG::SparseMatrix> lfmig = rcp(new LINALG::SparseMatrix(fmig->RowMap(),81,false));

    lfmig->Add(*fmig,false,1.0,0.0);
    lfmig->Complete(s->DomainMap(),fmig->RangeMap());

    lfmig->ApplyDirichlet( *(FluidField().GetDBCMapExtractor()->CondMap()),false);

    if (stcalgo != INPAR::STR::stc_none)
    {
      lfmig = LINALG::MLMultiply(*lfmig,false,*stcmat, false, false, false,true);
    }

    mat.Matrix(1,0).Add(*lfmig,false,1.0,1.0);

    // We cannot copy the pressure value. It is not used anyway. So no exact
    // match here.
    (*fmiitransform_)(mmm->FullRowMap(),
                      mmm->FullColMap(),
                      fmii,
                      1.,
                      ADAPTER::CouplingMasterConverter(coupfa),
                      mat.Matrix(1,2),
                      false);


    RCP<LINALG::SparseMatrix> lfmgi = rcp(new LINALG::SparseMatrix(fmgi.RowMap(),81,false));
    (*fmiitransform_)(mmm->FullRowMap(),
                      mmm->FullColMap(),
                      fmgi,
                      1.0,
                      ADAPTER::CouplingMasterConverter(coupfa),
                      *lfmgi,
                      false);

    lfmgi->Complete(aii.DomainMap(),mortar->RangeMap());
    RCP<LINALG::SparseMatrix> llfmgi = MLMultiply(*mortar,true,*lfmgi,false,false,false,true);
    lfmgi = rcp(new LINALG::SparseMatrix(s->RowMap(),81,false));

    lfmgi->Add(*llfmgi,false,scale,0.0);
    lfmgi->Complete(aii.DomainMap(),s->RangeMap());

    lfmgi->ApplyDirichlet( *(StructureField()->GetDBCMapExtractor()->CondMap()),false);
    if (stcalgo == INPAR::STR::stc_currsym)
      lfmgi = LINALG::MLMultiply(*stcmat, true, *lfmgi, false, true, true, false);
    mat.Assign(0,2,View,*lfmgi);

  }

  s->Complete();
  s->ApplyDirichlet( *(StructureField()->GetDBCMapExtractor()->CondMap()),true);

  if (stcalgo != INPAR::STR::stc_none)
  {
    s = LINALG::MLMultiply(*s, false, *stcmat, false, true, true, true);

    if (stcalgo == INPAR::STR::stc_currsym)
      s = LINALG::MLMultiply(*stcmat, true, *s, false, true, true, false);
  }
  else
  {
    s->UnComplete();
  }

  // finally assign structure matrix to block (0,0)
  mat.Assign(0,0,View,*s);

  // done. make sure all blocks are filled.
  mat.Complete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplit::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::InitialGuess");

  SetupVector(*ig,
              StructureField()->InitialGuess(),
              FluidField().InitialGuess(),
              AleField().InitialGuess(),
              0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplit::ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  const bool scaling_infnorm = (bool)DRT::INPUT::IntegralValue<int>(fsidyn,"INFNORMSCALING");

  if (scaling_infnorm)
  {
    // The matrices are modified here. Do we have to change them back later on?

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0,0).EpetraMatrix();
    srowsum_ = rcp(new Epetra_Vector(A->RowMap(),false));
    scolsum_ = rcp(new Epetra_Vector(A->RowMap(),false));
    A->InvRowSums(*srowsum_);
    A->InvColSums(*scolsum_);
    if (A->LeftScale(*srowsum_) or
        A->RightScale(*scolsum_) or
        mat.Matrix(0,1).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0,2).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(1,0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(2,0).EpetraMatrix()->RightScale(*scolsum_))
      dserror("structure scaling failed");

    A = mat.Matrix(2,2).EpetraMatrix();
    arowsum_ = rcp(new Epetra_Vector(A->RowMap(),false));
    acolsum_ = rcp(new Epetra_Vector(A->RowMap(),false));
    A->InvRowSums(*arowsum_);
    A->InvColSums(*acolsum_);
    if (A->LeftScale(*arowsum_) or
        A->RightScale(*acolsum_) or
        mat.Matrix(2,0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2,1).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0,2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(1,2).EpetraMatrix()->RightScale(*acolsum_))
      dserror("ale scaling failed");

    Teuchos::RCP<Epetra_Vector> sx = Extractor().ExtractVector(b,0);
    Teuchos::RCP<Epetra_Vector> ax = Extractor().ExtractVector(b,2);

    if (sx->Multiply(1.0, *srowsum_, *sx, 0.0))
      dserror("structure scaling failed");
    if (ax->Multiply(1.0, *arowsum_, *ax, 0.0))
      dserror("ale scaling failed");

    Extractor().InsertVector(*sx,0,b);
    Extractor().InsertVector(*ax,2,b);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplit::UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  const bool scaling_infnorm = (bool)DRT::INPUT::IntegralValue<int>(fsidyn,"INFNORMSCALING");

  if (scaling_infnorm)
  {
    Teuchos::RCP<Epetra_Vector> sy = Extractor().ExtractVector(x,0);
    Teuchos::RCP<Epetra_Vector> ay = Extractor().ExtractVector(x,2);

    if (sy->Multiply(1.0, *scolsum_, *sy, 0.0))
      dserror("structure scaling failed");
    if (ay->Multiply(1.0, *acolsum_, *ay, 0.0))
      dserror("ale scaling failed");

    // get info about STC feature and unscale solution if necessary
    INPAR::STR::STC_Scale stcalgo = StructureField()->GetSTCAlgo();
    if (stcalgo != INPAR::STR::stc_none)
    {
      StructureField()->GetSTCMat()->Multiply(false,*sy,*sy);
    }

    Extractor().InsertVector(*sy,0,x);
    Extractor().InsertVector(*ay,2,x);

    Teuchos::RCP<Epetra_Vector> sx = Extractor().ExtractVector(b,0);
    Teuchos::RCP<Epetra_Vector> ax = Extractor().ExtractVector(b,2);

    if (sx->ReciprocalMultiply(1.0, *srowsum_, *sx, 0.0))
      dserror("structure scaling failed");
    if (ax->ReciprocalMultiply(1.0, *arowsum_, *ax, 0.0))
      dserror("ale scaling failed");

    // get info about STC feature
    if (stcalgo != INPAR::STR::stc_none)
    {
      StructureField()->GetSTCMat()->Multiply(false,*sx,*sx);
    }

    Extractor().InsertVector(*sx,0,b);
    Extractor().InsertVector(*ax,2,b);

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0,0).EpetraMatrix();
    srowsum_->Reciprocal(*srowsum_);
    scolsum_->Reciprocal(*scolsum_);
    if (A->LeftScale(*srowsum_) or
        A->RightScale(*scolsum_) or
        mat.Matrix(0,1).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0,2).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(1,0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(2,0).EpetraMatrix()->RightScale(*scolsum_))
      dserror("structure scaling failed");

    A = mat.Matrix(2,2).EpetraMatrix();
    arowsum_->Reciprocal(*arowsum_);
    acolsum_->Reciprocal(*acolsum_);
    if (A->LeftScale(*arowsum_) or
        A->RightScale(*acolsum_) or
        mat.Matrix(2,0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2,1).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0,2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(1,2).EpetraMatrix()->RightScale(*acolsum_))
      dserror("ale scaling failed");
  }

  // very simple hack just to see the linear solution

  Epetra_Vector r(b.Map());
  mat.Apply(x,r);
  r.Update(1.,b,1.);

  Teuchos::RCP<Epetra_Vector> sr = Extractor().ExtractVector(r,0);
  Teuchos::RCP<Epetra_Vector> fr = Extractor().ExtractVector(r,1);
  Teuchos::RCP<Epetra_Vector> ar = Extractor().ExtractVector(r,2);

  // increment additional ale residual
  aleresidual_->Update(-1.,*ar,0.);

  ios_base::fmtflags flags = Utils()->out().flags();

  double n,ns,nf,na;
  r.Norm2(&n);
  sr->Norm2(&ns);
  fr->Norm2(&nf);
  ar->Norm2(&na);
  Utils()->out() << std::scientific
                 << "\nlinear solver quality:\n"
                 << "L_2-norms:\n"
                 << END_COLOR "   |r|=" YELLOW << n
                 << END_COLOR "   |rs|=" YELLOW << ns
                 << END_COLOR "   |rf|=" YELLOW << nf
                 << END_COLOR "   |ra|=" YELLOW << na
                 << END_COLOR "\n";
  r.NormInf(&n);
  sr->NormInf(&ns);
  fr->NormInf(&nf);
  ar->NormInf(&na);
  Utils()->out() << "L_inf-norms:\n"
                 << END_COLOR "   |r|=" YELLOW << n
                 << END_COLOR "   |rs|=" YELLOW << ns
                 << END_COLOR "   |rf|=" YELLOW << nf
                 << END_COLOR "   |ra|=" YELLOW << na
                 << END_COLOR "\n";

  Utils()->out().flags(flags);

  if (StructureField()->GetSTCAlgo() != INPAR::STR::stc_none)
    StructureField()->SystemMatrix()->Reset();

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplit::SetupVector(Epetra_Vector &f,
                                         Teuchos::RCP<const Epetra_Vector> sv,
                                         Teuchos::RCP<const Epetra_Vector> fv,
                                         Teuchos::RCP<const Epetra_Vector> av,
                                         double fluidscale)
{

  // extract the inner and boundary dofs of all three fields

  Teuchos::RCP<Epetra_Vector> fov = FluidField().Interface()->ExtractOtherVector(fv);
#ifdef FLUIDSPLITAMG
  fov = FluidField().Interface()->InsertOtherVector(fov);
#endif
  Teuchos::RCP<Epetra_Vector> aov = AleField().Interface()->ExtractOtherVector(av);

  if (fluidscale!=0)
  {
    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> fcv = FluidField().Interface()->ExtractFSICondVector(fv);
    Teuchos::RCP<Epetra_Vector> scv = LINALG::CreateVector(*StructureField()->Interface()->FSICondMap());

    Teuchos::RCP<LINALG::SparseMatrix> mortar = coupsfm_->GetMortarTrafo();
    mortar->Multiply(true,*fcv,*scv);

    Teuchos::RCP<Epetra_Vector> modsv = StructureField()->Interface()->InsertFSICondVector(scv);
    modsv->Update(1.0, *sv, fluidscale);

    Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(new const Epetra_Vector(modsv->Map(),true));
    LINALG::ApplyDirichlettoSystem(modsv,zeros,*(StructureField()->GetDBCMapExtractor()->CondMap()));

    if (StructureField()->GetSTCAlgo() == INPAR::STR::stc_currsym)
    {
      Teuchos::RCP<LINALG::SparseMatrix> stcmat = StructureField()->GetSTCMat();
      stcmat->Multiply(true,*modsv,*modsv);
    }

    Extractor().InsertVector(*modsv,0,f);

  }
  else
  {
    Teuchos::RCP<Epetra_Vector> modsv =  rcp(new Epetra_Vector(*sv));
    if (StructureField()->GetSTCAlgo() == INPAR::STR::stc_currsym)
    {
      Teuchos::RCP<LINALG::SparseMatrix> stcmat = StructureField()->GetSTCMat();
      stcmat->Multiply(true,*sv,*modsv);
    }
    Extractor().InsertVector(*modsv,0,f);
  }

  Extractor().InsertVector(*fov,1,f);
  Extractor().InsertVector(*aov,2,f);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem>
FSI::MortarMonolithicFluidSplit::CreateLinearSystem(ParameterList& nlParams,
                                           NOX::Epetra::Vector& noxSoln,
                                           Teuchos::RCP<NOX::Utils> utils)
{
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys;

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  NOX::Epetra::Interface::Jacobian* iJac = this;
  NOX::Epetra::Interface::Preconditioner* iPrec = this;
  const Teuchos::RCP< Epetra_Operator > J = systemmatrix_;
  const Teuchos::RCP< Epetra_Operator > M = systemmatrix_;

  switch (linearsolverstrategy_)
  {
  case INPAR::FSI::PreconditionedKrylov:
  case INPAR::FSI::FSIAMG:
    linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,
                                                               lsParams,
                                                               Teuchos::rcp(iJac,false),
                                                               J,
                                                               Teuchos::rcp(iPrec,false),
                                                               M,
                                                               noxSoln));
    break;
  default:
    dserror("unsupported linear block solver strategy: %d", linearsolverstrategy_);
  }

  return linSys;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo>
FSI::MortarMonolithicFluidSplit::CreateStatusTest(Teuchos::ParameterList& nlParams,
                                         Teuchos::RCP<NOX::Epetra::Group> grp)
{
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo       =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  Teuchos::RCP<NOX::StatusTest::Combo> converged   =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 100)));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv    =
      Teuchos::rcp(new NOX::StatusTest::FiniteValue);

  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // require one solve
  converged->addStatusTest(Teuchos::rcp(new NOX::FSI::MinIters(1)));

  // setup tests for structural displacements

  Teuchos::RCP<NOX::StatusTest::Combo> structcombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> structureDisp =
    Teuchos::rcp(new NOX::FSI::PartialNormF("displacement",
                                            Extractor(),0,
                                            nlParams.get("Norm abs disp", 1.0e-6),
                                            NOX::Abstract::Vector::TwoNorm,
                                            NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate =
    Teuchos::rcp(new NOX::FSI::PartialNormUpdate("displacement update",
                                                 Extractor(),0,
                                                 nlParams.get("Norm abs disp", 1.0e-6),
                                                 NOX::FSI::PartialNormUpdate::Scaled));

  AddStatusTest(structureDisp);
  structcombo->addStatusTest(structureDisp);
  //structcombo->addStatusTest(structureDispUpdate);

  converged->addStatusTest(structcombo);

  // setup tests for fluid velocities

  std::vector<Teuchos::RCP<const Epetra_Map> > fluidvel;
  fluidvel.push_back(FluidField().InnerVelocityRowMap());
  fluidvel.push_back(Teuchos::null);
  LINALG::MultiMapExtractor fluidvelextract(*DofRowMap(),fluidvel);

  Teuchos::RCP<NOX::StatusTest::Combo> fluidvelcombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> innerFluidVel =
    Teuchos::rcp(new NOX::FSI::PartialNormF("velocity",
                                            fluidvelextract,0,
                                            nlParams.get("Norm abs vel", 1.0e-6),
                                            NOX::Abstract::Vector::TwoNorm,
                                            NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> innerFluidVelUpdate =
    Teuchos::rcp(new NOX::FSI::PartialNormUpdate("velocity update",
                                                 fluidvelextract,0,
                                                 nlParams.get("Norm abs vel", 1.0e-6),
                                                 NOX::FSI::PartialNormUpdate::Scaled));

  AddStatusTest(innerFluidVel);
  fluidvelcombo->addStatusTest(innerFluidVel);
  //fluidvelcombo->addStatusTest(innerFluidVelUpdate);

  converged->addStatusTest(fluidvelcombo);

  // setup tests for fluid pressure

  std::vector<Teuchos::RCP<const Epetra_Map> > fluidpress;
  fluidpress.push_back(FluidField().PressureRowMap());
  fluidpress.push_back(Teuchos::null);
  LINALG::MultiMapExtractor fluidpressextract(*DofRowMap(),fluidpress);

  Teuchos::RCP<NOX::StatusTest::Combo> fluidpresscombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> fluidPress =
    Teuchos::rcp(new NOX::FSI::PartialNormF("pressure",
                                            fluidpressextract,0,
                                            nlParams.get("Norm abs pres", 1.0e-6),
                                            NOX::Abstract::Vector::TwoNorm,
                                            NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> fluidPressUpdate =
    Teuchos::rcp(new NOX::FSI::PartialNormUpdate("pressure update",
                                                 fluidpressextract,0,
                                                 nlParams.get("Norm abs pres", 1.0e-6),
                                                 NOX::FSI::PartialNormUpdate::Scaled));

  AddStatusTest(fluidPress);
  fluidpresscombo->addStatusTest(fluidPress);
  //fluidpresscombo->addStatusTest(fluidPressUpdate);

  converged->addStatusTest(fluidpresscombo);

  return combo;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplit::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                                 Teuchos::RCP<const Epetra_Vector>& sx,
                                                 Teuchos::RCP<const Epetra_Vector>& fx,
                                                 Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::ExtractFieldVectors");

  // We have overlap at the interface. Thus we need the interface part of the
  // structure vector and append it to the fluid and ale vector. (With the
  // right translation.)

  sx = Extractor().ExtractVector(x,0);
  Teuchos::RCP<const Epetra_Vector> scx = StructureField()->Interface()->ExtractFSICondVector(sx);

  // get fluid displacements

  Teuchos::RCP<const Epetra_Vector> fox = Extractor().ExtractVector(x,1);
#ifdef FLUIDSPLITAMG
  fox = FluidField().Interface()->ExtractOtherVector(fox);
#endif
  Teuchos::RCP<LINALG::SparseMatrix> mortar = coupsfm_->GetMortarTrafo();
  Teuchos::RCP<Epetra_Vector> fcx = LINALG::CreateVector(*FluidField().Interface()->FSICondMap());
  mortar->Apply(*scx,*fcx);

  // process ale unknowns based on fluid displacements
  Teuchos::RCP<Epetra_Vector> acx =  icoupfa_->MasterToSlave(fcx);

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x,2);

  Teuchos::RCP<Epetra_Vector> a = AleField().Interface()->InsertOtherVector(aox);
  AleField().Interface()->InsertFSICondVector(acx, a);
  ax = a;

  // process fluid unknowns
  FluidField().DisplacementToVelocity(fcx);

  Teuchos::RCP<Epetra_Vector> f = FluidField().Interface()->InsertOtherVector(fox);
  FluidField().Interface()->InsertFSICondVector(fcx, f);
  fx = f;

}


void FSI::MortarMonolithicFluidSplit::Update()
{

  // update history variabels for sliding ale
  if (aleproj_!= INPAR::FSI::ALEprojection_none)
  {
    iprojdisp_ = Teuchos::rcp(new Epetra_Vector(*coupsfm_->SlaveDofRowMap(),true));
    Teuchos::RCP<Epetra_Vector> idispale =
        icoupfa_->SlaveToMaster(AleField().Interface()->ExtractFSICondVector(AleField().ExtractDisplacement()));

    slideale_->Remeshing(*StructureField(),
                        FluidField().Discretization(),
                        idispale,
                        iprojdisp_,
                        *coupsfm_,
                        Comm());

    iprojdispinc_->Update(-1.0,*iprojdisp_,1.0,*idispale,0.0);

    slideale_->EvaluateMortar(StructureField()->ExtractInterfaceDispnp(), iprojdisp_, *coupsfm_);
    slideale_->EvaluateFluidMortar(idispale,iprojdisp_);

    RCP<Epetra_Vector> temp = rcp(new Epetra_Vector(*iprojdisp_));
    temp->ReplaceMap(idispale->Map());
    Teuchos::RCP<Epetra_Vector> acx = icoupfa_->MasterToSlave(temp);
    AleField().ApplyInterfaceDisplacements(acx);
    FluidField().ApplyMeshDisplacement(AleToFluid(AleField().ExtractDisplacement()));

    Teuchos::RCP<Epetra_Vector> unew = slideale_->InterpolateFluid(FluidField().ExtractInterfaceVelnp());
    FluidField().ApplyInterfaceVelocities(unew);
  }

  StructureField()->Update();
  FluidField().Update();
  AleField().Update();

}

void FSI::MortarMonolithicFluidSplit::Output()
{
  StructureField()->Output();
  FluidField().    Output();

  if (aleproj_!= INPAR::FSI::ALEprojection_none)
  {
    const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
    int uprestart = fsidyn.get<int>("RESTARTEVRY");
    if (uprestart != 0 && FluidField().Step() % uprestart == 0)
    {
      FluidField().DiscWriter()->WriteVector("slideALE", iprojdisp_);
      FluidField().DiscWriter()->WriteVector("slideALEincr", iprojdispinc_);
      slideale_->OutputRestart(*FluidField().DiscWriter());
    }
  }

  AleField().      Output();
  FluidField().LiftDrag();

  if (StructureField()->GetConstraintManager()->HaveMonitor())
  {
    StructureField()->GetConstraintManager()->ComputeMonitorValues(StructureField()->Dispnp());
    if(comm_.MyPID() == 0)
      StructureField()->GetConstraintManager()->PrintMonitorValues();
  }

}

void FSI::MortarMonolithicFluidSplit::ReadRestart(int step)
{
  StructureField()->ReadRestart(step);
  FluidField().ReadRestart(step);

  SetupSystem();

  if (aleproj_!= INPAR::FSI::ALEprojection_none)
  {
    IO::DiscretizationReader reader =
        IO::DiscretizationReader(FluidField().Discretization(),step);
    reader.ReadVector(iprojdisp_, "slideALE");
    reader.ReadVector(iprojdispinc_, "slideALEincr");
    slideale_->ReadRestart(reader);
  }

  AleField().ReadRestart(step);

  SetTimeStep(FluidField().Time(),FluidField().Step());

  if (aleproj_!= INPAR::FSI::ALEprojection_none)
    slideale_->EvaluateMortar(StructureField()->ExtractInterfaceDispn(), iprojdisp_, *coupsfm_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplit::PrepareTimeStep()
{
  IncrementTimeAndStep();

  PrintHeader();

  if (StructureField()->GetSTCAlgo() != INPAR::STR::stc_none)
      StructureField()->SystemMatrix()->Reset();
  StructureField()->PrepareTimeStep();
  FluidField().    PrepareTimeStep();
  AleField().      PrepareTimeStep();
}
