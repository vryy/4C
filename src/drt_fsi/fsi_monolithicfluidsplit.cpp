/*----------------------------------------------------------------------*/
/*!
\file fsi_monolithicfluidsplit.cpp

\brief Solve FSI problem with matching grids using a monolithic scheme
with condensed fluid interface velocities

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-15262
</pre>
*/

/*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "fsi_monolithicfluidsplit.H"
#include "fsi_matrixtransform.H"
#include "fsi_debugwriter.H"
#include "fsi_statustest.H"
#include "fsi_overlapprec_fsiamg.H"
#include "fsi_monolithic_linearsystem.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../linalg/linalg_utils.H"
#include "../drt_ale/ale_utils_mapextractor.H"

#include "../drt_constraint/constraint_manager.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#define FLUIDSPLITAMG

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicFluidSplit::MonolithicFluidSplit(const Epetra_Comm& comm,
                                                const Teuchos::ParameterList& timeparams)
  : BlockMonolithic(comm,timeparams)
{
  // Remove all interface DOFs holding a DBC from the fluid DBC map and give a 'warning'
  std::vector<Teuchos::RCP<const Epetra_Map> > intersectionmaps;
  intersectionmaps.push_back(FluidField().GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(FluidField().Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmap = LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

  if (intersectionmap->NumGlobalElements() != 0)
  {
    // remove interface DOFs from fluid DBC map
    FluidField().RemoveDirichCond(intersectionmap);

    // give a warning to the user that Dirichlet boundary conditions might not be correct
    if (comm.MyPID() == 0)
    {
      cout << "    +------------------------------------------------------------------------------------+" << endl;
      cout << "    |                                    PLEASE NOTE:                                    |" << endl;
      cout << "    +------------------------------------------------------------------------------------+" << endl;
      cout << "    | You run a monolithic fluid split scheme. Hence, there are no fluid interface DOFs. |" << endl;
      cout << "    | Fluid Dirichlet boundary conditions on the interface will be neglected.            |" << endl;
      cout << "    | Check whether you have prescribed appropriate DBCs on structural interface DOFs.   |" << endl;
      cout << "    +------------------------------------------------------------------------------------+" << endl;
    }
  }

#ifdef DEBUG
  // check if removing Dirichlet conditions was successful
  intersectionmaps.resize(0);
  intersectionmaps.push_back(FluidField().GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(FluidField().Interface()->FSICondMap());
  intersectionmap = LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);
  if (intersectionmap->NumGlobalElements() != 0)
    dserror("Removal of fluid interface Dirichlet conditions from fluid DBC map was not successful.");
#endif

  fggtransform_   = Teuchos::rcp(new UTILS::MatrixRowColTransform);
  fgitransform_   = Teuchos::rcp(new UTILS::MatrixRowTransform);
  figtransform_   = Teuchos::rcp(new UTILS::MatrixColTransform);
  aigtransform_   = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmiitransform_  = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmgitransform_  = Teuchos::rcp(new UTILS::MatrixRowColTransform);

  // Recovery of Lagrange multiplier happens on fluid field
  lambda_   = Teuchos::rcp(new Epetra_Vector(*FluidField().Interface()->FSICondMap(),true));
  fmgiprev_ = Teuchos::null;
  fmgicur_  = Teuchos::null;
  fmggprev_ = Teuchos::null;
  fmggcur_  = Teuchos::null;
  fgiprev_  = Teuchos::null;
  fgicur_   = Teuchos::null;
  fggprev_  = Teuchos::null;
  fggcur_   = Teuchos::null;

#ifdef DEBUG
  // check whether allocation was successful
  if (fggtransform_   == Teuchos::null) { dserror("Allocation of 'fggtransform_' failed."); }
  if (fgitransform_   == Teuchos::null) { dserror("Allocation of 'fgitransform_' failed."); }
  if (figtransform_   == Teuchos::null) { dserror("Allocation of 'figtransform_' failed."); }
  if (aigtransform_   == Teuchos::null) { dserror("Allocation of 'aigtransform_' failed."); }
  if (fmiitransform_  == Teuchos::null) { dserror("Allocation of 'fmiitransform_' failed."); }
  if (fmgitransform_  == Teuchos::null) { dserror("Allocation of 'fmgitransform_' failed."); }
  if (lambda_         == Teuchos::null) { dserror("Allocation of 'lambda_' failed."); }
#endif

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::SetupSystem()
{

  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  linearsolverstrategy_ = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");

  SetDefaultParameters(fsidyn,NOXParameterList());

  // call SetupSystem in base class
  FSI::Monolithic::SetupSystem();

  // create combined map

  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;
  vecSpaces.push_back(StructureField()->DofRowMap());
#ifdef FLUIDSPLITAMG
  vecSpaces.push_back(FluidField()    .DofRowMap());
#else
  vecSpaces.push_back(FluidField()    .Interface()->OtherMap());
#endif
  vecSpaces.push_back(AleField()      .Interface()->OtherMap());

  if (vecSpaces[1]->NumGlobalElements()==0)
    dserror("No inner fluid equations. Splitting not possible. Panic.");

  SetDofRowMaps(vecSpaces);

  /*----------------------------------------------------------------------*/
  // Switch fluid to interface split block matrix
  FluidField().UseBlockMatrix(true);

  // build ale system matrix in splitted system
  AleField().BuildSystemMatrix(false);

  aleresidual_ = Teuchos::rcp(new Epetra_Vector(*AleField().Interface()->OtherMap()));

  std::vector<int> pciter;
  std::vector<double> pcomega;
  std::vector<int> spciter;
  std::vector<double> spcomega;
  std::vector<int> fpciter;
  std::vector<double> fpcomega;
  std::vector<int> apciter;
  std::vector<double> apcomega;
  std::vector<std::string> blocksmoother;
  std::vector<double> schuromega;
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
      std::string word;
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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::SetupRHS(Epetra_Vector& f, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::SetupRHS");

  SetupVector(f,
              StructureField()->RHS(),
              FluidField().RHS(),
              AleField().RHS(),
              FluidField().ResidualScaling());

  // add additional ale residual
  Extractor().AddVector(*aleresidual_,2,f);

  firstcall_ = firstcall;

  // The following terms of rhs are only considered in the first Newton iteration
  // since we formulate our linear system in iteration increments. They transport
  // information from the last time step or the predictor into the Newton loop.
  if (firstcall)
  {
    // get time integration parameters of structure an fluid time integrators
    // to enable consistent time integration among the fields
    const double stiparam = StructureField()->TimIntParam();
    const double ftiparam = FluidField().TimIntParam();

    // some scaling factors for fluid
    const double timescale = FluidField().TimeScaling();
    const double scale     = FluidField().ResidualScaling();

    // old interface velocity of fluid field
    const Teuchos::RCP<const Epetra_Vector> fveln = FluidField().ExtractInterfaceVeln();
    
    // store structural interface displacement increment due to predictor
    // or inhomogeneous Dirichlet boundary conditions
    ddgpred_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->ExtractInterfaceDispnp()));
    ddgpred_->Update(-1.0, *StructureField()->ExtractInterfaceDispn(), 1.0);

    // store fluid interface velocity increment due to predictor
    // or inhomogeneous Dirichlet boundary conditions
    dugpred_ = Teuchos::rcp(new Epetra_Vector(*FluidField().ExtractInterfaceVelnp()));
    dugpred_->Update(-1.0, *fveln, 1.0);

    // get fluid matrix
    const Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockf = FluidField().BlockSystemMatrix();

    // get fluid shape derivatives matrix
    const Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField().ShapeDerivatives();

    // get ale matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocka = AleField().BlockSystemMatrix();

#ifdef DEBUG
    if (blockf==Teuchos::null)  { dserror("Expected Teuchos::rcp to fluid block matrix."); }
    if (blocka==Teuchos::null)  { dserror("Expected Teuchos::rcp to ale block matrix."); }
#endif

    // extract fluid and ale submatrices
    const LINALG::SparseMatrix& fig = blockf->Matrix(0,1); // F_{I\Gamma}
    const LINALG::SparseMatrix& fgg = blockf->Matrix(1,1); // F_{\Gamma\Gamma}
    const LINALG::SparseMatrix& aig = blocka->Matrix(0,1); // A_{I\Gamma}

    // some often re-used vectors
    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null;

    // Different contributions/terms to the rhs are separated by the following comment line
    // ---------- structural interface DOFs
    /* The following terms are added to the structural interface DOFs of right hand side:
     *
     * rhs_firstnewtonstep =
     *
     * (1)  + (1-stiparam)/(1-ftiparam) * dt / tau * F_{\Gamma\Gamma} * u^{n}_{\Gamma}
     *
     * (2)  - (1-stiparam)/(1-ftiparam) / tau * F_{\Gamma\Gamma} * \Delta d_{\Gamma,p}
     *
     * (3)  - (1-stiparam)/(1-ftiparam) * F^{G}_{\Gamma\Gamma} * \Delta d_{\Gamma,p}
     *
     * (4)  + (1-stiparam)/(1-ftiparam) * F_{\Gamma\Gamma} * \Delta u_{\Gamma,p}
     *
     * Remarks on all terms:
     * +  tau: time scaling factor for interface time integration (tau = 1/FluidField().TimeScaling())
     *
     */
    // ----------addressing term 1
    rhs = Teuchos::rcp(new Epetra_Vector(fgg.RangeMap(),true));

    fgg.Apply(*fveln,*rhs);

    rhs->Scale(scale * (1.-stiparam)/(1.-ftiparam) * Dt() * timescale);

    rhs = FluidToStruct(rhs);
    rhs = StructureField()->Interface()->InsertFSICondVector(rhs);

    if (StructureField()->GetSTCAlgo() == INPAR::STR::stc_currsym)
    {
      Teuchos::RCP<LINALG::SparseMatrix> stcmat = StructureField()->GetSTCMat();
      stcmat->Multiply(true,*rhs,*rhs);
    }

    Extractor().AddVector(*rhs,0,f);
    // ----------end of term 1

    // ----------addressing term 2:
    rhs = Teuchos::rcp(new Epetra_Vector(fgg.RangeMap(),true));

    fgg.Apply(*StructToFluid(ddgpred_), *rhs);

    rhs->Scale(-scale * (1.-stiparam) / (1.-ftiparam) * timescale);
    rhs = StructureField()->Interface()->InsertFSICondVector(FluidToStruct(rhs));

    Extractor().AddVector(*rhs,0,f);
    // ----------end of term 2

    // ----------addressing term 3:
    if (mmm != Teuchos::null)
    {
      // extract F^{G}_{\Gamma\Gamma}
      const LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);

      rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RangeMap(),true));

      fmgg.Apply(*StructToFluid(ddgpred_), *rhs);

      rhs->Scale(-(1.-stiparam) / (1.-ftiparam));
      rhs = StructureField()->Interface()->InsertFSICondVector(FluidToStruct(rhs));

      Extractor().AddVector(*rhs,0,f);
    }
    // ----------end of term 3

    // ----------addressing term 4:
    rhs = Teuchos::rcp(new Epetra_Vector(fgg.RowMap(),true));

    fgg.Apply(*dugpred_,*rhs);

    rhs->Scale((1.0-stiparam)/(1.0-ftiparam));
    rhs = StructureField()->Interface()->InsertFSICondVector(FluidToStruct(rhs));

    if (StructureField()->GetSTCAlgo() == INPAR::STR::stc_currsym)
    {
      Teuchos::RCP<LINALG::SparseMatrix> stcmat = StructureField()->GetSTCMat();
      stcmat->Multiply(true,*rhs,*rhs);
    }

    Extractor().AddVector(*rhs,0,f);
    // ----------end of term 4
    // ----------end of structural interface DOFs

    // ---------- inner fluid DOFs
    /* The following terms are added to the inner fluid DOFs of right hand side:
     *
     * rhs_firstnewtonstep =
     *
     * (1)  + dt / tau F_{I \Gamma} u^{n}_{\Gamma}
     *
     * (2)  - 1 / tau F_{I \Gamma} * \Delta d_{\Gamma,p}
     *
     * (3)  - F^{G}_{I \Gamma} * \Delta d_{\Gamma,p}
     *
     * (4)  + F_{I \Gamma} * \Delta u_{\Gamma,p}
     *
     * Remarks on all terms:
     * +  tau: time scaling factor for interface time integration (tau = 1/FluidField().TimeScaling())
     *
     */
    // ----------addressing term 1
    rhs = Teuchos::rcp(new Epetra_Vector(fig.RangeMap(),true));

    fig.Apply(*fveln,*rhs);

    rhs->Scale(Dt() * timescale);

#ifdef FLUIDSPLITAMG
    rhs = FluidField().Interface()->InsertOtherVector(rhs);
#endif

    Extractor().AddVector(*rhs,1,f);
    // ----------end of term 1

    // ----------addressing term 2
    rhs = Teuchos::rcp(new Epetra_Vector(fig.RangeMap(),true));

    fig.Apply(*StructToFluid(ddgpred_),*rhs);

    rhs->Scale(-timescale);

#ifdef FLUIDSPLITAMG
    rhs = FluidField().Interface()->InsertOtherVector(rhs);
#endif

    Extractor().AddVector(*rhs,1,f);
    // ----------end of term 2

    // ----------addressing term 3
    if(mmm != Teuchos::null)
    {
      // extract F^{G}_{I \Gamma}
      const LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);

      rhs = Teuchos::rcp(new Epetra_Vector(fmig.RangeMap(),true));

      fmig.Apply(*StructToFluid(ddgpred_),*rhs);

      rhs->Scale(-1.);

#ifdef FLUIDSPLITAMG
    rhs = FluidField().Interface()->InsertOtherVector(rhs);
#endif

      Extractor().AddVector(*rhs,1,f);
    }
    // ----------end of term 3

    // ----------addressing term 4
    rhs = Teuchos::rcp(new Epetra_Vector(fig.RangeMap(),true));

    fig.Apply(*dugpred_,*rhs);

#ifdef FLUIDSPLITAMG
    rhs = FluidField().Interface()->InsertOtherVector(rhs);
#endif

    Extractor().AddVector(*rhs,1,f);
    // ----------end of term 4
    // ----------end of inner fluid DOFs

    // ---------- inner ale DOFs
    /* The following terms are added to the inner ale DOFs of right hand side:
     *
     * rhs_firstnewtonstep =
     *
     * (1)  - A_{I \Gamma} * \Delta d_{\Gamma,p}
     *
     */
    // ----------addressing term 1
    rhs = Teuchos::rcp(new Epetra_Vector(aig.RangeMap(),true));

    aig.Apply(*StructToAle(ddgpred_),*rhs);
    rhs->Scale(-1.0);

    Extractor().AddVector(*rhs,2,f);
    // ----------end of term 1
    // ---------- end of inner ale DOFs

    // -----------------------------------------------------
    // Now, all contributions/terms to rhs in the first Newton iteration are added.

    // Apply Dirichlet boundary conditions
    // structure
    rhs = Extractor().ExtractVector(f,0);
    Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(rhs->Map(),true));
    LINALG::ApplyDirichlettoSystem(rhs,zeros,*(StructureField()->GetDBCMapExtractor()->CondMap()));
    Extractor().InsertVector(*rhs,0,f);

    // fluid
    rhs = Extractor().ExtractVector(f,1);
    zeros = Teuchos::rcp(new Epetra_Vector(rhs->Map(),true));
    LINALG::ApplyDirichlettoSystem(rhs,zeros,*(FluidField().GetDBCMapExtractor()->CondMap()));
    Extractor().InsertVector(*rhs,1,f);

    // ale
    rhs = Extractor().ExtractVector(f,2);
    zeros = Teuchos::rcp(new Epetra_Vector(rhs->Map(),true));
    LINALG::ApplyDirichlettoSystem(rhs,zeros,*(AleField().GetDBCMapExtractor()->CondMap()));
    Extractor().InsertVector(*rhs,2,f);
    // -----------------------------------------------------

    // Reset quantities of previous iteration step since they still store values from the last time step
    ddginc_ = LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(),true);
    duiinc_ = LINALG::CreateVector(*FluidField().Interface()->OtherMap(),true);
    ddialeinc_ = LINALG::CreateVector(*AleField().Interface()->OtherMap(),true);
    soliprev_ = Teuchos::null;
    solgprev_ = Teuchos::null;
    fgcur_ = LINALG::CreateVector(*FluidField().Interface()->FSICondMap(),true);
    fgicur_ = Teuchos::null;
    fggcur_ = Teuchos::null;
  }

  // NOX expects the 'positive' residual. The negative sign for the
  // linearized Newton system J*dx=-r is done internally by NOX.
  f.Scale(-1.);

  // store interface force onto the fluid to know it in the next time step as previous force
  // in order to recover the Lagrange multiplier
  fgprev_ = fgcur_;
  fgcur_ = FluidField().Interface()->ExtractFSICondVector(FluidField().RHS());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::SetupSystemMatrix");

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  const ADAPTER::Coupling& coupsa = StructureAleCoupling();
  const ADAPTER::Coupling& coupfa = FluidAleCoupling();

  // get info about STC feature
  INPAR::STR::STC_Scale stcalgo = StructureField()->GetSTCAlgo();
  Teuchos::RCP<LINALG::SparseMatrix> stcmat = Teuchos::null;
  // if STC is to be used, get STC matrix from structure field
  if (stcalgo != INPAR::STR::stc_none)
    stcmat = StructureField()->GetSTCMat();

  // get single field block matrices
  Teuchos::RCP<LINALG::SparseMatrix> s = StructureField()->SystemMatrix(); // can't be 'const' --> is modified by STC
  const Teuchos::RCP<LINALG::BlockSparseMatrixBase> f = FluidField().BlockSystemMatrix();
  const Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField().BlockSystemMatrix();

#ifdef DEBUG
  if (s==Teuchos::null) { dserror("expect structure block matrix"); }
  if (f==Teuchos::null) { dserror("expect fluid block matrix"); }
  if (a==Teuchos::null) { dserror("expect ale block matrix"); }
#endif

  // extract submatrices
  LINALG::SparseMatrix& fii = f->Matrix(0,0);
  LINALG::SparseMatrix& fig = f->Matrix(0,1);
  LINALG::SparseMatrix& fgi = f->Matrix(1,0);
  LINALG::SparseMatrix& fgg = f->Matrix(1,1);
  LINALG::SparseMatrix& aii = a->Matrix(0,0);
  LINALG::SparseMatrix& aig = a->Matrix(0,1);
  
  // scaling factors for fluid
  const double scale     = FluidField().ResidualScaling();
  const double timescale = FluidField().TimeScaling();

  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();
  const double ftiparam = FluidField().TimIntParam();

  /*----------------------------------------------------------------------*/

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  // uncomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  // this just once...
  s->UnComplete();

  (*fggtransform_)( fgg,
                    (1.0-stiparam)/(1.0-ftiparam)*scale*timescale,
                    ADAPTER::CouplingSlaveConverter(coupsf),
                    ADAPTER::CouplingSlaveConverter(coupsf),
                    *s,
                    true,
                    true);

  Teuchos::RCP<LINALG::SparseMatrix> lfgi = Teuchos::rcp(new LINALG::SparseMatrix(s->RowMap(),81,false));
  (*fgitransform_)(fgi,
                   (1.0-stiparam)/(1.0-ftiparam)*scale,
                   ADAPTER::CouplingSlaveConverter(coupsf),
                   *lfgi);

  lfgi->Complete(fgi.DomainMap(),s->RangeMap());
  lfgi->ApplyDirichlet( *(StructureField()->GetDBCMapExtractor()->CondMap()),false);

  if (stcalgo == INPAR::STR::stc_currsym)
    lfgi = LINALG::MLMultiply(*stcmat, true, *lfgi, false, true, true, true);

#ifdef FLUIDSPLITAMG
  mat.Matrix(0,1).UnComplete();
  mat.Matrix(0,1).Add(*lfgi,false,1.,0.0);
#else
  mat.Assign(0,1,View,*lfgi);
#endif

  if (stcalgo == INPAR::STR::stc_none)
  {
    Teuchos::RCP<LINALG::SparseMatrix> lfig = Teuchos::rcp(new LINALG::SparseMatrix(fig.RowMap(),81,false));
    (*figtransform_)(f->FullRowMap(),
                     f->FullColMap(),
                     fig,
                     timescale,
                     ADAPTER::CouplingSlaveConverter(coupsf),
                     mat.Matrix(1,0));
  }
  else
  {
    Teuchos::RCP<LINALG::SparseMatrix> lfig = Teuchos::rcp(new LINALG::SparseMatrix(fig.RowMap(),81,false));
    (*figtransform_)(f->FullRowMap(),
                     f->FullColMap(),
                     fig,
                     timescale,
                     ADAPTER::CouplingSlaveConverter(coupsf),
                     *lfig);

    lfig->Complete(s->DomainMap(),fig.RangeMap());

    lfig = LINALG::MLMultiply(*lfig,false,*stcmat, false, false, false,true);

#ifdef FLUIDSPLITAMG
    mat.Matrix(1,0).UnComplete();
    mat.Matrix(1,0).Add(*lfig,false,1.,0.0);
#else
    mat.Assign(1,0,View,*lfig);
#endif
  }

#ifdef FLUIDSPLITAMG
  mat.Matrix(1,1).Add(fii,false,1.,0.0);
  Teuchos::RCP<LINALG::SparseMatrix> eye = LINALG::Eye(*FluidField().Interface()->FSICondMap());
  mat.Matrix(1,1).Add(*eye,false,1.,1.0);
#else
  mat.Assign(1,1,View,fii);
#endif

  if (stcalgo == INPAR::STR::stc_none)
  {
    (*aigtransform_)( a->FullRowMap(),
                      a->FullColMap(),
                      aig,
                      1.,
                      ADAPTER::CouplingSlaveConverter(coupsa),
                      mat.Matrix(2,0));
  }
  else
  {
    Teuchos::RCP<LINALG::SparseMatrix> laig = Teuchos::rcp(new LINALG::SparseMatrix(aii.RowMap(),81,false));
    (*aigtransform_)( a->FullRowMap(),
                      a->FullColMap(),
                      aig,
                      1.,
                      ADAPTER::CouplingSlaveConverter(coupsa),
                      *laig);

    laig->Complete(s->DomainMap(),laig->RangeMap());
    laig->ApplyDirichlet( *(AleField().GetDBCMapExtractor()->CondMap()),false);

    if (stcalgo != INPAR::STR::stc_none)
    {
      laig = LINALG::MLMultiply(*laig,false,*stcmat, false, false, false,true);
    }

    mat.Assign(2,0,View,*laig);
  }

  mat.Assign(2,2,View,aii);

  /*----------------------------------------------------------------------*/
  // add optional fluid linearization with respect to mesh motion block

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField().ShapeDerivatives();
  if (mmm!=Teuchos::null)
  {
    LINALG::SparseMatrix& fmii = mmm->Matrix(0,0);
    LINALG::SparseMatrix& fmgi = mmm->Matrix(1,0);

    LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);
    LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);

    // reuse transform objects to add shape derivative matrices to structural blocks

    if (stcalgo == INPAR::STR::stc_none)
    {
      (*figtransform_)(f->FullRowMap(),
                       f->FullColMap(),
                       fmig,
                       1.,
                       ADAPTER::CouplingSlaveConverter(coupsf),
                       mat.Matrix(1,0),
                       false,
                       true);
    }
    else
    {
      Teuchos::RCP<LINALG::SparseMatrix> lfmig = Teuchos::rcp(new LINALG::SparseMatrix(fmig.RowMap(),81,false));
      (*figtransform_)(f->FullRowMap(),
                       f->FullColMap(),
                       fmig,
                       1.,
                       ADAPTER::CouplingSlaveConverter(coupsf),
                       *lfmig,
                       false,
                       true);


      lfmig->Complete(s->DomainMap(),fmig.RangeMap());

      lfmig->ApplyDirichlet( *(FluidField().GetDBCMapExtractor()->CondMap()),false);

      if (stcalgo != INPAR::STR::stc_none)
      {
        lfmig = LINALG::MLMultiply(*lfmig,false,*stcmat, false, false, false,true);
      }

      mat.Matrix(1,0).Add(*lfmig,false,1.0,1.0);
    }

    (*fggtransform_)(fmgg,
                     (1.0-stiparam)/(1.0-ftiparam)*scale,
                     ADAPTER::CouplingSlaveConverter(coupsf),
                     ADAPTER::CouplingSlaveConverter(coupsf),
                     *s,
                     false,
                     true);

    // We cannot copy the pressure value. It is not used anyway. So no exact
    // match here.
    (*fmiitransform_)(mmm->FullRowMap(),
                      mmm->FullColMap(),
                      fmii,
                      1.,
                      ADAPTER::CouplingMasterConverter(coupfa),
                      mat.Matrix(1,2),
                      false);

    {
      Teuchos::RCP<LINALG::SparseMatrix> lfmgi = Teuchos::rcp(new LINALG::SparseMatrix(s->RowMap(),81,false));
      (*fmgitransform_)(fmgi,
                        (1.0-stiparam)/(1.0-ftiparam)*scale,
                        ADAPTER::CouplingSlaveConverter(coupsf),
                        ADAPTER::CouplingMasterConverter(coupfa),
                        *lfmgi,
                        false,
                        false);

      lfmgi->Complete(aii.DomainMap(),s->RangeMap());
      lfmgi->ApplyDirichlet( *(StructureField()->GetDBCMapExtractor()->CondMap()),false);

      if (stcalgo == INPAR::STR::stc_currsym)
        lfmgi = LINALG::MLMultiply(*stcmat, true, *lfmgi, false, true, true, true);

#ifdef FLUIDSPLITAMG
      mat.Matrix(0,2).UnComplete();
      mat.Matrix(0,2).Add(*lfmgi,false,1.,0.0);
#else
      mat.Assign(0,2,View,*lfmgi);
#endif

    }
  }

  s->Complete();
  s->ApplyDirichlet( *(StructureField()->GetDBCMapExtractor()->CondMap()),true);

  if (stcalgo == INPAR::STR::stc_none)
  {
    s->UnComplete();
  }
  else  // apply STC matrix on block (0,0) if STC is used
  {
    s = LINALG::MLMultiply(*s, false, *stcmat, false, true, true, true);
    if (stcalgo == INPAR::STR::stc_currsym)
      s = LINALG::MLMultiply(*stcmat, true, *s, false, true, true, false);
  }

  // finally assign structure block
  mat.Matrix(0,0).Assign(View,*s);

  // done. make sure all blocks are filled.
  mat.Complete();

  // store parts of fluid matrix to know them in the next iteration as previous iteration matrices
  fgiprev_ = fgicur_;
  fggprev_ = fggcur_;
  fgicur_ = Teuchos::rcp(new LINALG::SparseMatrix(f->Matrix(1,0)));
  fggcur_ = Teuchos::rcp(new LINALG::SparseMatrix(f->Matrix(1,1)));

  // store parts of fluid shape derivative matrix to know them in the next iteration as previous iteration matrices
  fmgiprev_ = fmgicur_;
  fmggprev_ = fmggcur_;
  if (mmm!=Teuchos::null)
  {
    fmgicur_ = Teuchos::rcp(new LINALG::SparseMatrix(mmm->Matrix(1,0)));
    fmggcur_ = Teuchos::rcp(new LINALG::SparseMatrix(mmm->Matrix(1,1)));
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::InitialGuess");

  SetupVector(*ig,
              StructureField()->InitialGuess(),
              FluidField().InitialGuess(),
              AleField().InitialGuess(),
              0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  const bool scaling_infnorm = (bool)DRT::INPUT::IntegralValue<int>(fsidyn,"INFNORMSCALING");

  if (scaling_infnorm)
  {
    // The matrices are modified here. Do we have to change them back later on?

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0,0).EpetraMatrix();
    srowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
    scolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
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
    arowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
    acolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
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
void FSI::MonolithicFluidSplit::UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
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
void FSI::MonolithicFluidSplit::SetupVector(Epetra_Vector &f,
                                         Teuchos::RCP<const Epetra_Vector> sv,
                                         Teuchos::RCP<const Epetra_Vector> fv,
                                         Teuchos::RCP<const Epetra_Vector> av,
                                         const double fluidscale)
{
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();
  const double ftiparam = FluidField().TimIntParam();

  // extract inner dofs
  Teuchos::RCP<Epetra_Vector> fov = FluidField().Interface()->ExtractOtherVector(fv);
#ifdef FLUIDSPLITAMG
  fov = FluidField().Interface()->InsertOtherVector(fov);
#endif
  Teuchos::RCP<Epetra_Vector> aov = AleField().Interface()->ExtractOtherVector(av);

  if (fluidscale!=0)
  {
    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> fcv = FluidField().Interface()->ExtractFSICondVector(fv);
    Teuchos::RCP<Epetra_Vector> modsv = StructureField()->Interface()->InsertFSICondVector(FluidToStruct(fcv));
    modsv->Update(1.0, *sv, (1.0-stiparam)/(1.0-ftiparam)*fluidscale);

    // add contribution of Lagrange multiplier from previous time step
    if (lambda_ != Teuchos::null)
    {
      Teuchos::RCP<Epetra_Vector> lambdaglobal = StructureField()->Interface()->InsertFSICondVector(FluidToStruct(lambda_));
      modsv->Update(stiparam-(ftiparam*(1.0-stiparam))/(1.0-ftiparam), *lambdaglobal, 1.0);
    }

    Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(new const Epetra_Vector(modsv->Map(),true));
    LINALG::ApplyDirichlettoSystem(modsv,zeros,*(StructureField()->GetDBCMapExtractor()->CondMap()));

    if (StructureField()->GetSTCAlgo() == INPAR::STR::stc_currsym)
    {
      Teuchos::RCP<LINALG::SparseMatrix> stcmat = StructureField()->GetSTCMat();
      stcmat->Multiply(true,*modsv,*modsv);
    }

    Extractor().InsertVector(*modsv,0,f); // add structural contributions to 'f'
  }
  else
  {
    Teuchos::RCP<Epetra_Vector> modsv =  Teuchos::rcp(new Epetra_Vector(*sv));
    if (StructureField()->GetSTCAlgo() == INPAR::STR::stc_currsym)
    {
      Teuchos::RCP<LINALG::SparseMatrix> stcmat = StructureField()->GetSTCMat();
      stcmat->Multiply(true,*sv,*modsv);
    }

    Extractor().InsertVector(*modsv,0,f); // add structural contributions to 'f'
  }

  Extractor().InsertVector(*fov,1,f); // add fluid contributions to 'f'
  Extractor().InsertVector(*aov,2,f); // add ALE contributions to 'f'
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem>
FSI::MonolithicFluidSplit::CreateLinearSystem(ParameterList& nlParams,
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
    break;
  }

  return linSys;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo>
FSI::MonolithicFluidSplit::CreateStatusTest(Teuchos::ParameterList& nlParams,
                                         Teuchos::RCP<NOX::Epetra::Group> grp)
{
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo       = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  Teuchos::RCP<NOX::StatusTest::Combo> converged   = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  // Create some other plausibility tests
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 100)));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv    = Teuchos::rcp(new NOX::StatusTest::FiniteValue);

  Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
    Teuchos::rcp(new NOX::StatusTest::NormUpdate(nlParams.get("Norm Update", 1.0e-5)));

  // Add single tests to the test combo
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(update);
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

  // setup tests for interface
  std::vector<Teuchos::RCP<const Epetra_Map> > interface;
  interface.push_back(StructureField()->Interface()->FSICondMap());
  interface.push_back(Teuchos::null);
  LINALG::MultiMapExtractor interfaceextract(*DofRowMap(),interface);

  Teuchos::RCP<NOX::StatusTest::Combo> interfacecombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> interfaceTest =
    Teuchos::rcp(new NOX::FSI::PartialNormF("interface",
                                            interfaceextract,0,
                                            nlParams.get("Norm abs vel", 1.0e-6),
                                            NOX::Abstract::Vector::TwoNorm,
                                            NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> interfaceTestUpdate =
    Teuchos::rcp(new NOX::FSI::PartialNormUpdate("interface update",
                                                 interfaceextract,0,
                                                 nlParams.get("Norm abs vel", 1.0e-6),
                                                 NOX::FSI::PartialNormUpdate::Scaled));

  AddStatusTest(interfaceTest);
  interfacecombo->addStatusTest(interfaceTest);
  //interfacecombo->addStatusTest(interfaceTestUpdate);

  converged->addStatusTest(interfacecombo);

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
void FSI::MonolithicFluidSplit::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                                 Teuchos::RCP<const Epetra_Vector>& sx,
                                                 Teuchos::RCP<const Epetra_Vector>& fx,
                                                 Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::ExtractFieldVectors");

#ifdef DEBUG
  if(ddgpred_ == Teuchos::null) { dserror("Vector 'ddgpred_' has not been initialized properly."); }
  if(dugpred_ == Teuchos::null) { dserror("Vector 'dugpred_' has not been initialized properly."); }
#endif

  // We have overlap at the interface. Thus we need the interface part of the
  // structure vector and append it to the fluid and ale vector. (With the
  // right translation.)

  sx = Extractor().ExtractVector(x,0);
  Teuchos::RCP<Epetra_Vector> scx = StructureField()->Interface()->ExtractFSICondVector(sx);

  // process fluid unknowns
  Teuchos::RCP<const Epetra_Vector> fox = Extractor().ExtractVector(x,1);
#ifdef FLUIDSPLITAMG
  fox = FluidField().Interface()->ExtractOtherVector(fox);
#endif
  Teuchos::RCP<Epetra_Vector> fcx = StructToFluid(scx);

  // consider fluid predictor increment only in first Newton iteration
  if (firstcall_)
  {
    FluidField().DisplacementToVelocity(fcx,StructToFluid(ddgpred_),dugpred_);
  }
  else
  {
    Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(fcx->Map(),true));
    FluidField().DisplacementToVelocity(fcx,StructToFluid(ddgpred_),zeros);
//    FluidField().DisplacementToVelocity(fcx);
  }

  Teuchos::RCP<Epetra_Vector> f = FluidField().Interface()->InsertOtherVector(fox);
  FluidField().Interface()->InsertFSICondVector(fcx, f);
  fx = f;

  // process ale unknowns
  scx->Update(1.0,*ddgpred_,1.0);

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x,2);
  Teuchos::RCP<Epetra_Vector> acx = StructToAle(scx);

  Teuchos::RCP<Epetra_Vector> a = AleField().Interface()->InsertOtherVector(aox);
  AleField().Interface()->InsertFSICondVector(acx, a);
  ax = a;

  // Store field vectors to know them later on as previous quantities
  // inner ale displacement increment
  // interface structure displacement increment
  if (disgprev_ != Teuchos::null)
    ddginc_->Update(1.0, *scx, -1.0, *disgprev_, 0.0);        // compute current iteration increment
  else
    ddginc_ = Teuchos::rcp(new Epetra_Vector(*scx));                   // first iteration increment

  disgprev_ = scx;                                            // store current step increment
  // ------------------------------------

  // inner ale displacement increment
  if (solialeprev_ != Teuchos::null)
    ddialeinc_->Update(1.0, *aox, -1.0, *solialeprev_, 0.0);  // compute current iteration increment
  else
    ddialeinc_ = Teuchos::rcp(new Epetra_Vector(*aox));                // first iteration increment

  solialeprev_ = aox;                                         // store current step increment
  // ------------------------------------

  // fluid solution increment
  if (soliprev_ != Teuchos::null)                             // compute current iteration increment
    duiinc_->Update(1.0, *fox, -1.0, *soliprev_, 0.0);
  else                                                        // first iteration increment
    duiinc_ =  Teuchos::rcp(new Epetra_Vector(*fox));
                                                              // store current step increment
  soliprev_ = fox;
  // ------------------------------------
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::Output()
{
  StructureField()->Output();
  FluidField().    Output();

  // output Lagrange multiplier
  {
    /* 'lambda_' is only defined on the interface. So, insert 'lambda_' into
     * 'lambdafull' that is defined on the entire fluid field. Then, write
     * output or restart data.
     */
    Teuchos::RCP<Epetra_Vector> lambdafull = FluidField().Interface()->InsertFSICondVector(lambda_);
    const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
    const int uprestart = fsidyn.get<int>("RESTARTEVRY");
    const int upres = fsidyn.get<int>("UPRES");
    if ((uprestart != 0 && FluidField().Step() % uprestart == 0) || FluidField().Step() % upres == 0)
      FluidField().DiscWriter()->WriteVector("fsilambda", lambdafull);
  }

  AleField().      Output();
  FluidField().LiftDrag();

  if (StructureField()->GetConstraintManager()->HaveMonitor())
  {
    StructureField()->GetConstraintManager()->ComputeMonitorValues(StructureField()->Dispnp());
    if(Comm().MyPID() == 0)
      StructureField()->GetConstraintManager()->PrintMonitorValues();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::ReadRestart(int step)
{
  StructureField()->ReadRestart(step);
  FluidField().ReadRestart(step);

  // read Lagrange multiplier
  {
    Teuchos::RCP<Epetra_Vector> lambdafull = Teuchos::rcp(new Epetra_Vector(*FluidField().DofRowMap(),true));
    IO::DiscretizationReader reader = IO::DiscretizationReader(FluidField().Discretization(),step);
    reader.ReadVector(lambdafull, "fsilambda");
    lambda_ = FluidField().Interface()->ExtractFSICondVector(lambdafull);
  }

  AleField().ReadRestart(step);

  SetTimeStep(FluidField().Time(),FluidField().Step());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::PrepareTimeStep()
{
  IncrementTimeAndStep();

  PrintHeader();

  if (StructureField()->GetSTCAlgo() != INPAR::STR::stc_none)
    StructureField()->SystemMatrix()->Reset();
  StructureField()->PrepareTimeStep();
  FluidField().    PrepareTimeStep();
  AleField().      PrepareTimeStep();
}

/*----------------------------------------------------------------------*/
/* Recover the Lagrange multiplier at the interface   mayr.mt (03/2012) */
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::RecoverLagrangeMultiplier()
{
  // get time integration parameter of fluid time integrator
  // to enable consistent time integration among the fields
  const double ftiparam = FluidField().TimIntParam();

  // some scaling factors for fluid
  const double timescale  = FluidField().TimeScaling();
  const double scale      = FluidField().ResidualScaling();

  // get fluid shape derivative matrix
  const Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField().ShapeDerivatives();

  // some often re-used vectors
  Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::null;     // stores intermediate result of terms (3)-(8)
  Teuchos::RCP<Epetra_Vector> auxvec = Teuchos::null;     // just for convenience
  Teuchos::RCP<Epetra_Vector> auxauxvec = Teuchos::null;  // just for convenience

  /* Recovery of Lagrange multiplier \lambda_^{n+1} is done by the following
   * condensation expression:
   *
   * lambda_^{n+1} =
   *
   * (1)  - ftiparam / (1.-ftiparam) * lambda^{n}
   *
   * (2)  - 1. / (1.-ftiparam) * tmpvec
   *
   * with tmpvec =
   *
   * (3)    r_{\Gamma}^{F,n+1}
   *
   * (4)  + 1 / tau * F_{\Gamma\Gamma} * \Delta d_{\Gamma}^{S,n+1}
   *
   * (5)  + F_{\Gamma\Gamma}^{G} * \Delta d_{\Gamma}^{S,n+1}
   *
   * (6)  + F_{\Gamma I} * \Delta u_{I}^{F,n+1}
   *
   * (7)  + F_{\Gamma I}^{G} * \Delta d_{I}^{G,n+1}
   *
   * (8)  + dt / tau * F_{\Gamma\Gamma} * u_{\Gamma}^n]
   *
   * Remark on term (8):
   * Term (8) has to be considered only in the first Newton iteration.
   * Hence, it will usually not be computed since in general we need more
   * than one nonlinear iteration until convergence.
   *
   * Remarks on all terms:
   * +  Division by -(1.0 - ftiparam) will be done in the end
   *    since this is common to all terms
   * +  tau: time scaling factor for interface time integration (tau = 1/FluidField().TimeScaling())
   * +  neglecting terms (4)-(8) should not alter the results significantly
   *    since at the end of the time step the solution increments tend to zero.
   *
   *                                                 Matthias Mayr (10/2012)
   */

  // ---------Addressing term (1)
  lambda_->Update(ftiparam,*lambda_,0.0);
  // ---------End of term (1)

  // ---------Addressing term (3)
  Teuchos::RCP<Epetra_Vector> fluidresidual = FluidField().Interface()->ExtractFSICondVector(FluidField().RHS());
  fluidresidual->Scale(-1.0); // invert sign to obtain residual, not rhs
  tmpvec = Teuchos::rcp(new Epetra_Vector(*fluidresidual));
  // ---------End of term (3)

  // ---------Addressing term (4)
  auxvec = Teuchos::rcp(new Epetra_Vector(fggprev_->RangeMap(),true));
  fggprev_->Apply(*StructToFluid(ddginc_), *auxvec);
  tmpvec->Update(timescale,*auxvec,1.0);
  // ---------End of term (4)

  // ---------Addressing term (5)
  if (fmggprev_!=Teuchos::null)
  {
    auxvec = Teuchos::rcp(new Epetra_Vector(fmggprev_->RangeMap(),true));
    fmggprev_->Apply(*StructToFluid(ddginc_),*auxvec);
    tmpvec->Update(1.0,*auxvec,1.0);
  }
  // ---------End of term (5)

  // ---------Addressing term (6)
  auxvec = Teuchos::rcp(new Epetra_Vector(fgiprev_->RangeMap(),true));
  fgiprev_->Apply(*duiinc_,*auxvec);
  tmpvec->Update(1.0,*auxvec,1.0);
  // ---------End of term (6)

  // ---------Addressing term (7)
  if (fmgiprev_!=Teuchos::null)
  {
    /* For matrix-vector-product, the DomainMap() of the matrix and the Map() of the vector
     * have to match. DomainMap() contains inner velocity DOFs and all pressure DOFs.
     * The inner ale displacement increment is converted to the fluid map using AleToFluid().
     * This results in a map that contains all velocity but no pressure DOFs.
     *
     * We have to circumvent some trouble with Epetra_BlockMaps since we cannot split
     * an Epetra_BlockMap into inner and interface DOFs.
     *
     * We create a map extractor 'velothermap' in order to extract the inner velocity
     * DOFs after calling AleToFluid(). Afterwards, a second map extractor 'velotherpressuremapext'
     * is used to append pressure DOFs filled with zeros.
     *
     * Finally, maps match and matrix-vector-multiplication can be done.
     */

    // extract inner velocity DOFs after calling AleToFluid()
    Teuchos::RCP<Epetra_Map> velothermap = LINALG::SplitMap(*FluidField().VelocityRowMap(),*InterfaceFluidAleCoupling().MasterDofMap());
    LINALG::MapExtractor velothermapext = LINALG::MapExtractor(*FluidField().VelocityRowMap(),velothermap,false);
    auxvec = Teuchos::rcp(new Epetra_Vector(*velothermap, true));
    velothermapext.ExtractOtherVector(AleToFluid(AleField().Interface()->InsertOtherVector(ddialeinc_)),auxvec);

    // add pressure DOFs
    LINALG::MapExtractor velotherpressuremapext = LINALG::MapExtractor(fmgiprev_->DomainMap(),velothermap);
    auxauxvec = Teuchos::rcp(new Epetra_Vector(fmgiprev_->DomainMap(), true));
    velotherpressuremapext.InsertCondVector(auxvec,auxauxvec);

    // prepare vector to store result of matrix-vector-product
    auxvec = Teuchos::rcp(new Epetra_Vector(fmgiprev_->RangeMap(),true));

    // Now, do the actual matrix-vector-product
    fmgiprev_->Apply(*auxauxvec,*auxvec);
    tmpvec->Update(1.0,*auxvec,1.0);
  }
  // ---------End of term (7)

  // ---------Addressing term (8)
  if (firstcall_)
  {
    auxvec = Teuchos::rcp(new Epetra_Vector(fggprev_->RangeMap(),true));
    fggprev_->Apply(*FluidField().ExtractInterfaceVeln(),*auxvec);
    tmpvec->Update(Dt()*timescale,*auxvec,1.0);
  }
  // ---------End of term (8)

  // ---------Addressing term (2)
  lambda_->Update(scale,*tmpvec,1.0); // scale with ResidualScaling() to get [N/m^2]
  // ---------End of term (2)

  // Finally, divide by (1.0-ftiparam) which is common to all terms
  lambda_->Scale(-1.0/(1.0-ftiparam));

  // Finally, the Lagrange multiplier 'lambda_' is recovered here.
  // It represents nodal forces acting onto the structure.

  return;
}
