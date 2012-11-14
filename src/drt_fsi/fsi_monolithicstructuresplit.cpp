/*----------------------------------------------------------------------*/
/*!
\file fsi_mortarmonolithic_structuresplit.cpp

\brief Solve FSI problem with matching grids using a monolithic scheme
with condensed structure interface displacements

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-15262
</pre>
*/

/*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "fsi_monolithicstructuresplit.H"
#include "fsi_matrixtransform.H"
#include "fsi_overlapprec_fsiamg.H"
#include "fsi_statustest.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_adapter/ad_fld_fluid.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../linalg/linalg_utils.H"
#include "../drt_ale/ale_utils_mapextractor.H"

#include "fsi_nox_group.H"

#include "fsi_debugwriter.H"

#include "../drt_constraint/constraint_manager.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicStructureSplit::MonolithicStructureSplit(const Epetra_Comm& comm,
                                                        const Teuchos::ParameterList& timeparams)
  : BlockMonolithic(comm,timeparams)
{
  // Throw an error if there are DBCs on structural interface DOFs.
  std::vector<Teuchos::RCP<const Epetra_Map> > intersectionmaps;
  intersectionmaps.push_back(StructureField()->GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(StructureField()->Interface()->FSICondMap());
  Teuchos::RefCountPtr<Epetra_Map> intersectionmap = LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

  if (intersectionmap->NumGlobalElements() != 0)
  {
    // remove interface DOFs from structural DBC map
    StructureField()->RemoveDirichCond(intersectionmap);

    // give a warning to the user that Dirichlet boundary conditions might not be correct
    if (comm.MyPID() == 0)
    {
      cout << "  +---------------------------------------------------------------------------------------------+" << endl;
      cout << "  |                                        PLEASE NOTE:                                         |" << endl;
      cout << "  +---------------------------------------------------------------------------------------------+" << endl;
      cout << "  | You run a monolithic structure split scheme. Hence, there are no structural interface DOFs. |" << endl;
      cout << "  | Fluid Dirichlet boundary conditions on the interface will be neglected.                     |" << endl;
      cout << "  | Check whether you have prescribed appropriate DBCs on structural interface DOFs.            |" << endl;
      cout << "  +---------------------------------------------------------------------------------------------+" << endl;
    }
  }

#ifdef DEBUG
  // check if removing Dirichlet conditions was successful
  intersectionmaps.resize(0);
  intersectionmaps.push_back(StructureField()->GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(StructureField()->Interface()->FSICondMap());
  intersectionmap = LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);
  if (intersectionmap->NumGlobalElements() != 0)
    dserror("Could not remove structural interface Dirichlet conditions from structure DBC map.");
#endif

  sggtransform_  = Teuchos::rcp(new UTILS::MatrixRowColTransform);
  sgitransform_ = Teuchos::rcp(new UTILS::MatrixRowTransform);
  sigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  aigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  fmiitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmgitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  fsaigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fsmgitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  fscoupfa_ = Teuchos::rcp(new ADAPTER::Coupling());

  // Recovery of Lagrange multiplier happens on structure field
  lambda_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->FSICondMap(),true));
  ddiinc_ = Teuchos::null;
  soliprev_ = Teuchos::null;
  ddginc_ = Teuchos::null;
  duginc_ = Teuchos::null;
  disgprev_ = Teuchos::null;
  sgiprev_ = Teuchos::null;
  sggprev_ = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicStructureSplit::SetupSystem()
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  linearsolverstrategy_ = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");

  SetDefaultParameters(fsidyn,NOXParameterList());

  // call SetupSystem in base class
  FSI::Monolithic::SetupSystem();
  
  // we might have a free surface
  if (FluidField().Interface()->FSCondRelevant())
  {
    const int ndim = DRT::Problem::Instance()->NDim();

    fscoupfa_->SetupConditionCoupling(*FluidField().Discretization(),
                                       FluidField().Interface()->FSCondMap(),
                                      *AleField().Discretization(),
                                       AleField().Interface()->FSCondMap(),
                                      "FREESURFCoupling",
                                       ndim);
  }

  // create combined map

  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;
  vecSpaces.push_back(StructureField()->Interface()->OtherMap());
  vecSpaces.push_back(FluidField()    .DofRowMap());
  vecSpaces.push_back(AleField()      .Interface()->OtherMap());

  if (vecSpaces[0]->NumGlobalElements()==0)
    dserror("No inner structural equations. Splitting not possible. Panic.");

  SetDofRowMaps(vecSpaces);

  // Use normal matrix for fluid equations but build (splitted) mesh movement
  // linearization (if requested in the input file)
  FluidField().UseBlockMatrix(false);

  // Use splitted structure matrix
  StructureField()->UseBlockMatrix();

  // build ale system matrix in splitted system
  AleField().BuildSystemMatrix(false);

  aleresidual_ = Teuchos::rcp(new Epetra_Vector(*AleField().Interface()->OtherMap()));

  // get the PCITER from inputfile
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
                                                          true,
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
void FSI::MonolithicStructureSplit::SetupRHS(Epetra_Vector& f, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::SetupRHS");

  SetupVector(f,
              StructureField()->RHS(),
              FluidField().RHS(),
              AleField().RHS(),
              FluidField().ResidualScaling());

  // add additional ale residual
  Extractor().AddVector(*aleresidual_,2,f);

  firstcall_ = firstcall;

  // The following terms of rhs are only considered in the first Newton iteration
  // since we formulate our linear system in iteration increments.
  if (firstcall)
  {
    // get structure matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocks = StructureField()->BlockSystemMatrix();
    if (blocks==Teuchos::null)
      dserror("expect structure block matrix");

    // get fluid shape derivatives matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField().ShapeDerivatives();

    //get ale matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocka = AleField().BlockSystemMatrix();
    if (blocka==Teuchos::null)
      dserror("expect ale block matrix");

    // extract structure and ale submatrices
    LINALG::SparseMatrix& sig = blocks->Matrix(0,1);  // S_{I\Gamma}
    LINALG::SparseMatrix& sgg = blocks->Matrix(1,1);  // S_{\Gamma\Gamma}
    LINALG::SparseMatrix& aig = blocka->Matrix(0,1);  // A_{I\Gamma}

    // store structural interface displacement increment due to predictor
    // or inhomogeneous Dirichlet boundary conditions
    ddgpre_ = rcp(new Epetra_Vector(*StructureField()->ExtractInterfaceDispnp()));
    ddgpre_->Update(-1.0, *StructureField()->ExtractInterfaceDispn(), 1.0);

    // store fluid interface velocity increment due to predictor
    // or inhomogeneous Dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector> fveln = FluidField().ExtractInterfaceVeln();
    dugpre_ = rcp(new Epetra_Vector(*FluidField().ExtractInterfaceVelnp()));
    dugpre_->Update(-1.0, *fveln, 1.0);

    // get time integration parameters of structure an fluid time integrators
    // to enable consistent time integration among the fields
    const double stiparam = StructureField()->TimIntParam();
    const double ftiparam = FluidField().TimIntParam();

    // some scaling factors for fluid
    const double timescale = FluidField().TimeScaling();
    const double scale     = FluidField().ResidualScaling();
    const double dt        = FluidField().Dt();

    // some often re-used vectors
    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null;

    // ---------- inner structural DOFs
    /* The following terms are added to the inner structural DOFs of right hand side:
     *
     * rhs_firstnewtonstep =
     *
     * (1)  S_{I\Gamma} * \Delta d_{\Gamma,p}
     *
     * (2)  - dt * S_{I\Gamma} * u^{n}_{\Gamma}
     *
     * (3)  - 1/timescale * S_{I\Gamma} *  Delta u_{\Gamma,p}
     *
     */
    // ----------adressing term 1
    rhs = rcp(new Epetra_Vector(sig.RowMap(),true));
    sig.Apply(*ddgpre_,*rhs);

    Extractor().AddVector(*rhs,0,f);

    // ----------adressing term 2
    rhs = rcp(new Epetra_Vector(sig.RowMap(),true));
    sig.Apply(*FluidToStruct(fveln),*rhs);
    rhs->Scale(-dt);

    Extractor().AddVector(*rhs,0,f);

    // ----------adressing term 3
    rhs = rcp(new Epetra_Vector(sig.RowMap(),true));
    sig.Apply(*FluidToStruct(dugpre_),*rhs);
    rhs->Scale(-1.0/timescale);

    Extractor().AddVector(*rhs,0,f);
    // ---------- end of inner structural DOFs

    // ---------- inner fluid DOFs
    /* The following terms are added to the inner fluid DOFs of right hand side:
     *
     * rhs_firstnewtonstep =
     *
     * (1)  - dt * F^{G}_{I\Gamma} * u^{n}_{\Gamma}
     *
     * (2)  - 1/timescale * F^{G}_{I\Gamma} * \Delta u_{\Gamma,p}
     *
     */
    // ----------addressing term 1
    if (mmm!=Teuchos::null)
    {
      LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);
      rhs = Teuchos::rcp(new Epetra_Vector(fmig.RowMap(),true));
      fmig.Apply(*fveln,*rhs);
      rhs->Scale(-dt);

      rhs = FluidField().Interface()->InsertOtherVector(rhs);

      Extractor().AddVector(*rhs,1,f);
    }

    // ----------addressing term 2
    if (mmm!=Teuchos::null)
    {
      LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);
      rhs = Teuchos::rcp(new Epetra_Vector(fmig.RowMap(),true));
      fmig.Apply(*dugpre_,*rhs);
      rhs->Scale(-1.0/timescale);

      rhs = FluidField().Interface()->InsertOtherVector(rhs);

      Extractor().AddVector(*rhs,1,f);
    }
    // ---------- end of inner fluid DOFs

    // ---------- interface fluid DOFs
    /* The following terms are added to the interface fluid DOFs of right hand side:
     *
     * rhs_firstnewtonstep =
     *
     * (1)  - dt * F^{G}_{\Gamma\Gamma} * u^{n}_{\Gamma}
     *
     * (2)  - 1/timescale * F^{G}_{\Gamma\Gamma} * \Delta u_{\Gamma,p}
     *
     * (3)  - dt * (1-ftiparam)/(1-stiparam) * S_{\Gamma\Gamma} * u^{n}_{\Gamma}
     *
     * (4)  - 1/timescale * (1-ftiparam)/(1-stiparam) * S_{\Gamma\Gamma} * \Delta u_{\Gamma,p}
     *
     * (5)  + (1-ftiparam)/(1-stiparam) * S_{\Gamma\Gamma} * \Delta d_{\Gamma,p}
     */
    // ----------addressing term 1
    if (mmm!=Teuchos::null)
    {
      LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);
      rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RowMap(),true));
      fmgg.Apply(*fveln,*rhs);
      rhs->Scale(-dt);

      rhs = FluidField().Interface()->InsertFSICondVector(rhs);

      Extractor().AddVector(*rhs,1,f);
    }

    // ----------addressing term 2
    if (mmm!=Teuchos::null)
    {
      LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);
      rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RowMap(),true));
      fmgg.Apply(*dugpre_,*rhs);
      rhs->Scale(-1.0/timescale);

      rhs = FluidField().Interface()->InsertFSICondVector(rhs);

      Extractor().AddVector(*rhs,1,f);
    }

    // ----------addressing term 3
    rhs = rcp(new Epetra_Vector(sgg.RowMap(),true));
    sgg.Apply(*FluidToStruct(fveln),*rhs);
    rhs->Scale(-dt * (1.-ftiparam) / ((1-stiparam) * scale));

    rhs = StructToFluid(rhs);
    rhs = FluidField().Interface()->InsertFSICondVector(rhs);

    Extractor().AddVector(*rhs,1,f);

    // ----------addressing term 4
    rhs = rcp(new Epetra_Vector(sgg.RowMap(),true));
    sgg.Apply(*FluidToStruct(dugpre_),*rhs);
    rhs->Scale(-1.0 *  (1.-ftiparam) / ((1-stiparam) * timescale * scale));

    rhs = StructToFluid(rhs);
    rhs = FluidField().Interface()->InsertFSICondVector(rhs);

    Extractor().AddVector(*rhs,1,f);

    // ----------addressing term 5
    rhs = rcp(new Epetra_Vector(sgg.RowMap(),true));
    sgg.Apply(*ddgpre_,*rhs);
    rhs->Scale((1.-ftiparam) / ((1-stiparam) * scale));

    rhs = StructToFluid(rhs);
    rhs = FluidField().Interface()->InsertFSICondVector(rhs);

    Extractor().AddVector(*rhs,1,f);
    // ---------- end of interface fluid DOFs

    // ---------- inner ale DOFs
    /* The following terms are added to the inner ale DOFs of right hand side:
     *
     * rhs_firstnewtonstep =
     *
     * (1)  - dt * A_{I\Gamma} * u^{n}_{\Gamma}
     *
     * (2)  - 1/timescale * A_{I\Gamma} * \Delta u_{\Gamma,p}
     *
     */
    // ----------addressing term 1
    rhs = rcp(new Epetra_Vector(aig.RowMap(),true));
    aig.Apply(*FluidToAleInterface(fveln),*rhs);
    rhs->Scale(-dt);

    Extractor().AddVector(*rhs,2,f);

    // ----------addressing term 2
    rhs = rcp(new Epetra_Vector(aig.RowMap(),true));
    aig.Apply(*FluidToAleInterface(dugpre_),*rhs);
    rhs->Scale(-1.0/timescale);

    Extractor().AddVector(*rhs,2,f);
    // ---------- end of inner ale DOFs

    // -----------------------------------------------------
    // Now, all contributions/terms to rhs in the first Newton iteration are added.

    // Apply Dirichlet boundary conditions
    // structure
    rhs = Extractor().ExtractVector(f,0);
    Teuchos::RCP<const Epetra_Vector> zeros = rcp(new const Epetra_Vector(rhs->Map(),true));
    LINALG::ApplyDirichlettoSystem(rhs,zeros,*(StructureField()->GetDBCMapExtractor()->CondMap()));
    Extractor().InsertVector(*rhs,0,f);

    // fluid
    rhs = Extractor().ExtractVector(f,1);
    zeros = Teuchos::rcp(new const Epetra_Vector(rhs->Map(),true));
    LINALG::ApplyDirichlettoSystem(rhs,zeros,*(FluidField().GetDBCMapExtractor()->CondMap()));
    Extractor().InsertVector(*rhs,1,f);

    // ale
    rhs = Extractor().ExtractVector(f,2);
    zeros = Teuchos::rcp(new const Epetra_Vector(rhs->Map(),true));
    LINALG::ApplyDirichlettoSystem(rhs,zeros,*(AleField().GetDBCMapExtractor()->CondMap()));
    Extractor().InsertVector(*rhs,2,f);
    // -----------------------------------------------------

    // if there is a free surface
    if (FluidField().Interface()->FSCondRelevant())
    {
      // here we extract the free surface submatrices from position 2
      LINALG::SparseMatrix& aig = blocka->Matrix(0,2);

      // extract fluid free surface velocities.
      Teuchos::RCP<Epetra_Vector> fveln = FluidField().ExtractFreeSurfaceVeln();
      Teuchos::RCP<Epetra_Vector> aveln = InterfaceFluidAleCoupling().MasterToSlave(fveln);

      Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(aig.RowMap()));
      aig.Apply(*aveln,*rhs);

      rhs->Scale(-1.*Dt());

      Extractor().AddVector(*rhs,1,f);

      // shape derivatives
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField().ShapeDerivatives();
      if (mmm!=Teuchos::null)
      {
        // here we extract the free surface submatrices from position 2
        LINALG::SparseMatrix& fmig = mmm->Matrix(0,2);
        LINALG::SparseMatrix& fmgg = mmm->Matrix(2,2);

        rhs = Teuchos::rcp(new Epetra_Vector(fmig.RowMap()));
        fmig.Apply(*fveln,*rhs);
        Teuchos::RCP<Epetra_Vector> veln = FluidField().Interface()->InsertOtherVector(rhs);

        rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RowMap()));
        fmgg.Apply(*fveln,*rhs);
        FluidField().Interface()->InsertFSCondVector(rhs,veln);

        veln->Scale(-1.*Dt());

        Extractor().AddVector(*veln,0,f);
      }
    }

    // Reset quantities for previous iteration step since they still store values from the last time step
    ddiinc_ = LINALG::CreateVector(*StructureField()->Interface()->OtherMap(),true);
    soliprev_ = Teuchos::null;
    ddginc_ = LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(),true);
    duginc_ = LINALG::CreateVector(*FluidField().Interface()->FSICondMap(),true);
    disgprev_ = Teuchos::null;
    sgicur_ = Teuchos::null;
    sggcur_ = Teuchos::null;
  }

  // NOX expects the 'positive' residual. The negative sign for the
  // linearized Newton system J*dx=-r is done internally by NOX.
  f.Scale(-1.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicStructureSplit::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::SetupSystemMatrix");

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  const ADAPTER::Coupling& coupfa = FluidAleCoupling();
  const ADAPTER::Coupling& icoupfa = InterfaceFluidAleCoupling();

  // get single field block matrices
  const Teuchos::RCP<LINALG::BlockSparseMatrixBase> s = StructureField()->BlockSystemMatrix();
  const Teuchos::RCP<LINALG::SparseMatrix> f = FluidField().SystemMatrix();
  const Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField().BlockSystemMatrix();

#ifdef DEBUG
  if (s==Teuchos::null) { dserror("expect structure block matrix"); }
  if (f==Teuchos::null) { dserror("expect fluid block matrix"); }
  if (a==Teuchos::null) { dserror("expect ale block matrix"); }
#endif

  // extract submatrices
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

  // Uncomplete fluid matrix to be able to deal with slightly defective
  // interface meshes.
  f->UnComplete();

  mat.Assign(0,0,View,s->Matrix(0,0));

  (*sigtransform_)(s->FullRowMap(),
                   s->FullColMap(),
                   s->Matrix(0,1),
                   1./timescale,
                   ADAPTER::CouplingMasterConverter(coupsf),
                   mat.Matrix(0,1));
  (*sggtransform_)(s->Matrix(1,1),
                   (1.0-ftiparam)/((1.0-stiparam)*scale*timescale),
                   ADAPTER::CouplingMasterConverter(coupsf),
                   ADAPTER::CouplingMasterConverter(coupsf),
                   *f,
                   true,
                   true);

  RCP<LINALG::SparseMatrix> lsgi = rcp(new LINALG::SparseMatrix(f->RowMap(),81,false));
  (*sgitransform_)(s->Matrix(1,0),
                   (1.0-ftiparam)/((1.0-stiparam)*scale),
                   ADAPTER::CouplingMasterConverter(coupsf),
                   *lsgi);

  lsgi->Complete(s->Matrix(1,0).DomainMap(),f->RangeMap());
  lsgi->ApplyDirichlet( *(FluidField().GetDBCMapExtractor()->CondMap()),false);

  mat.Assign(1,0,View,*lsgi);

  (*aigtransform_)(a->FullRowMap(),
                   a->FullColMap(),
                   aig,
                   1./timescale,
                   ADAPTER::CouplingSlaveConverter(icoupfa),
                   mat.Matrix(2,1));
  mat.Assign(2,2,View,aii);

  /*----------------------------------------------------------------------*/
  // add optional fluid linearization with respect to mesh motion block

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField().ShapeDerivatives();
  if (mmm!=Teuchos::null)
  {
    LINALG::SparseMatrix& fmii = mmm->Matrix(0,0);
    LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);
    LINALG::SparseMatrix& fmgi = mmm->Matrix(1,0);
    LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);

    f->Add(fmgg,false,1./timescale,1.0);
    f->Add(fmig,false,1./timescale,1.0);

    RCP<LINALG::SparseMatrix> lfmgi = rcp(new LINALG::SparseMatrix(f->RowMap(),81,false));
    (*fmgitransform_)(mmm->FullRowMap(),
                      mmm->FullColMap(),
                      fmgi,
                      1.,
                      ADAPTER::CouplingMasterConverter(coupfa),
                      *lfmgi,
                      false,
                      false);

    (*fmiitransform_)(mmm->FullRowMap(),
                      mmm->FullColMap(),
                      fmii,
                      1.,
                      ADAPTER::CouplingMasterConverter(coupfa),
                      *lfmgi,
                      false,
                      true);

    lfmgi->Complete(aii.DomainMap(),f->RangeMap());
    lfmgi->ApplyDirichlet( *(FluidField().GetDBCMapExtractor()->CondMap()),false);

    mat.Assign(1,2,View,*lfmgi);
  }

  // if there is a free surface
  if (FluidField().Interface()->FSCondRelevant())
  {
    // here we extract the free surface submatrices from position 2
    LINALG::SparseMatrix& aig = a->Matrix(0,2);

    (*fsaigtransform_)(a->FullRowMap(),
                       a->FullColMap(),
                       aig,
                       1./timescale,
                       ADAPTER::CouplingSlaveConverter(*fscoupfa_),
                       mat.Matrix(2,1));

    if (mmm!=Teuchos::null)
    {
      // We assume there is some space between fsi interface and free
      // surface. Thus the matrices mmm->Matrix(1,2) and mmm->Matrix(2,1) are
      // zero.

      // here we extract the free surface submatrices from position 2
      LINALG::SparseMatrix& fmig = mmm->Matrix(0,2);
      LINALG::SparseMatrix& fmgi = mmm->Matrix(2,0);
      LINALG::SparseMatrix& fmgg = mmm->Matrix(2,2);

      mat.Matrix(1,1).Add(fmgg,false,1./timescale,1.0);
      mat.Matrix(1,1).Add(fmig,false,1./timescale,1.0);

      (*fsmgitransform_)(mmm->FullRowMap(),
                         mmm->FullColMap(),
                         fmgi,
                         1.,
                         ADAPTER::CouplingMasterConverter(coupfa),
                         mat.Matrix(1,2),
                         false,
                         false);
    }
  }

  f->Complete();
  f->ApplyDirichlet( *(FluidField().GetDBCMapExtractor()->CondMap()),true);

  // finally assign fluid block
  mat.Assign(1,1,View,*f);

//  //  print stiffness matrix to *.mtl-file that can be read in matlab
//  const std::string fname = "stiffmatrix.mtl";
//  cout << __LINE__ << __FILE__ << endl;
//  cout << "Printing stiffmatrix to file" << endl;
//  LINALG::PrintBlockMatrixInMatlabFormat(fname,mat);
//  exit(0);

  // done. make sure all blocks are filled.
  mat.Complete();

  // store parts of structural matrix to know them in the next iteration as previous iteration matrices
  sgiprev_ = sgicur_;
  sggprev_ = sggcur_;
  sgicur_ = rcp(new LINALG::SparseMatrix(s->Matrix(1,0)));
  sggcur_ = rcp(new LINALG::SparseMatrix(s->Matrix(1,1)));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicStructureSplit::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::InitialGuess");

  SetupVector(*ig,
              StructureField()->InitialGuess(),
              FluidField().InitialGuess(),
              AleField().InitialGuess(),
              0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicStructureSplit::ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
  //should we scale the system?
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
void FSI::MonolithicStructureSplit::UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
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

    Extractor().InsertVector(*sy,0,x);
    Extractor().InsertVector(*ay,2,x);

    Teuchos::RCP<Epetra_Vector> sx = Extractor().ExtractVector(b,0);
    Teuchos::RCP<Epetra_Vector> ax = Extractor().ExtractVector(b,2);

    if (sx->ReciprocalMultiply(1.0, *srowsum_, *sx, 0.0))
      dserror("structure scaling failed");
    if (ax->ReciprocalMultiply(1.0, *arowsum_, *ax, 0.0))
      dserror("ale scaling failed");

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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicStructureSplit::SetupVector(Epetra_Vector &f,
                                                Teuchos::RCP<const Epetra_Vector> sv,
                                                Teuchos::RCP<const Epetra_Vector> fv,
                                                Teuchos::RCP<const Epetra_Vector> av,
                                                double fluidscale)
{
  // get time integration parameters of structure an fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();
  const double ftiparam = FluidField().TimIntParam();

  // extract inner dofs
  Teuchos::RCP<Epetra_Vector> sov = StructureField()->Interface()->ExtractOtherVector(sv);
  Teuchos::RCP<Epetra_Vector> aov = AleField()      .Interface()->ExtractOtherVector(av);

  if (fluidscale!=0)
  {
    // add fluid interface values to structure vector
    // scv: structure interface DOFs
    Teuchos::RCP<Epetra_Vector> scv = StructureField()->Interface()->ExtractFSICondVector(sv);

    // modfv: whole fluid map but entries only at fsi DOFs
    Teuchos::RCP<Epetra_Vector> modfv = FluidField().Interface()->InsertFSICondVector(StructToFluid(scv));

    // add structure rhs on condensed structural interface DOFs to fluid DOFs
    int err = modfv->Update(1.0, *fv, (1.0-ftiparam)/((1.0-stiparam)*fluidscale));
    if (err != 0) { dserror("modfv->Update() returned err = %i.",err); }

    // add contribution of Lagrange multiplier from previous time step
    if (lambda_ != Teuchos::null)
    {
      Teuchos::RCP<Epetra_Vector> lambdaglobal = FluidField().Interface()->InsertFSICondVector(StructToFluid(lambda_));
      err = modfv->Update((-ftiparam+(stiparam*(1.0-ftiparam))/(1.0-stiparam))/fluidscale, *lambdaglobal, 1.0);
      if (err != 0) { dserror("modfv->Update() returned err = %i.",err); }
    }

    Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(new const Epetra_Vector(modfv->Map(),true));
    LINALG::ApplyDirichlettoSystem(modfv,zeros,*(FluidField().GetDBCMapExtractor()->CondMap()));

    Extractor().InsertVector(*modfv,1,f);
  }
  else
  {
    Extractor().InsertVector(*fv,1,f);
  }

  Extractor().InsertVector(*sov,0,f);
  Extractor().InsertVector(*aov,2,f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem>
FSI::MonolithicStructureSplit::CreateLinearSystem(ParameterList& nlParams,
                                                  NOX::Epetra::Vector& noxSoln,
                                                  Teuchos::RCP<NOX::Utils> utils)
{
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys;

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  Teuchos::ParameterList* lsParams = NULL;

  // in case of nonlinCG the linear solver list is somewhere else
  if (dirParams.get("Method","User Defined")=="User Defined")
    lsParams = &(newtonParams.sublist("Linear Solver"));
  else if (dirParams.get("Method","User Defined")=="NonlinearCG")
    lsParams = &(dirParams.sublist("Nonlinear CG").sublist("Linear Solver"));
  else dserror("Unknown nonlinear method");

  NOX::Epetra::Interface::Jacobian* iJac = this;
  NOX::Epetra::Interface::Preconditioner* iPrec = this;
  const Teuchos::RCP< Epetra_Operator > J = systemmatrix_;
  const Teuchos::RCP< Epetra_Operator > M = systemmatrix_;

  switch (linearsolverstrategy_)
  {
  case INPAR::FSI::PreconditionedKrylov:
  case INPAR::FSI::FSIAMG:
    linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,
                                                               *lsParams,
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
FSI::MonolithicStructureSplit::CreateStatusTest(Teuchos::ParameterList& nlParams,
                                                Teuchos::RCP<NOX::Epetra::Group> grp)
{
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo       = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  Teuchos::RCP<NOX::StatusTest::Combo> converged   = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 100)));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv    = Teuchos::rcp(new NOX::StatusTest::FiniteValue);

  Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
  Teuchos::rcp(new NOX::StatusTest::NormUpdate(nlParams.get("Norm Update", 1.0e-5)));

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
  interface.push_back(FluidField().Interface()->FSICondMap());
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
void FSI::MonolithicStructureSplit::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                                        Teuchos::RCP<const Epetra_Vector>& sx,
                                                        Teuchos::RCP<const Epetra_Vector>& fx,
                                                        Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::ExtractFieldVectors");

#ifdef DEBUG
  if(ddgpre_ == Teuchos::null) { dserror("Vector 'ddgpre_' has not been initialized properly."); }
  if(dugpre_ == Teuchos::null) { dserror("Vector 'dugpre_' has not been initialized properly."); }
#endif

  // process fluid unknowns
  fx = Extractor().ExtractVector(x,1);

  Teuchos::RCP<Epetra_Vector> fcx = FluidField().Interface()->ExtractFSICondVector(fx);

  // process structure unknowns
  FluidField().VelocityToDisplacement(fcx,StructToFluid(ddgpre_),dugpre_);

  Teuchos::RCP<const Epetra_Vector> sox = Extractor().ExtractVector(x,0);
  Teuchos::RCP<Epetra_Vector> scx = FluidToStruct(fcx);

  Teuchos::RCP<Epetra_Vector> s = StructureField()->Interface()->InsertOtherVector(sox);
  StructureField()->Interface()->InsertFSICondVector(scx, s);
  sx = s;

  // process ale unknowns
  Teuchos::RCP<Epetra_Vector> fcxforale = FluidField().Interface()->ExtractFSICondVector(fx);
  Teuchos::RCP<Epetra_Vector> zeros = rcp(new Epetra_Vector(fcxforale->Map(),true));
  FluidField().VelocityToDisplacement(fcxforale,zeros,dugpre_);
  Teuchos::RCP<Epetra_Vector> acx = FluidToStruct(fcxforale);
  acx = StructToAle(acx);

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x,2);
//  Teuchos::RCP<Epetra_Vector> acx = StructToAle(scx);

  Teuchos::RCP<Epetra_Vector> a = AleField().Interface()->InsertOtherVector(aox);
  AleField().Interface()->InsertFSICondVector(acx, a);

  // if there is a free surface
  if (FluidField().Interface()->FSCondRelevant())
  {
    Teuchos::RCP<Epetra_Vector> fcx = FluidField().Interface()->ExtractFSCondVector(fx);
    FluidField().FreeSurfVelocityToDisplacement(fcx);

    Teuchos::RCP<Epetra_Vector> acx = fscoupfa_->MasterToSlave(fcx);
    AleField().Interface()->InsertFSCondVector(acx, a);
  }
  ax = a;

  // Store field vectors to know them later on as previous quantities
  if (soliprev_ != Teuchos::null)
    ddiinc_->Update(1.0, *sox, -1.0, *soliprev_, 0.0);  // compute current iteration increment
  else
    ddiinc_ = rcp(new Epetra_Vector(*sox));           // first iteration increment

  soliprev_ = sox;                                      // store current step increment

  if (disgprev_ != Teuchos::null)
    ddginc_->Update(1.0, *scx, -1.0, *disgprev_, 0.0);  // compute current iteration increment
  else
    ddginc_ = rcp(new Epetra_Vector(*scx));           // first iteration increment

  disgprev_ = scx;                                      // store current step increment

  if (velgprev_ != Teuchos::null)
    duginc_->Update(1.0, *fcx, -1.0, *velgprev_, 0.0);  // compute current iteration increment
  else
    duginc_ = rcp(new Epetra_Vector(*fcx));             // first iteration increment

  velgprev_ = fcx;                                      // store current step increment
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicStructureSplit::Output()
{
  StructureField()->Output();

  // output Lagrange multiplier
  {
    Teuchos::RCP<Epetra_Vector> lambdafull = StructureField()->Interface()->InsertFSICondVector(lambda_);
    const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
    const int uprestart = fsidyn.get<int>("RESTARTEVRY");
    const int upres = fsidyn.get<int>("UPRES");
    if ((uprestart != 0 && FluidField().Step() % uprestart == 0) || FluidField().Step() % upres == 0)
      StructureField()->DiscWriter()->WriteVector("fsilambda", lambdafull);
  }

  FluidField().    Output();
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
void FSI::MonolithicStructureSplit::ReadRestart(int step)
{
  // read Lagrange multiplier
  {
    Teuchos::RCP<Epetra_Vector> lambdafull = rcp(new Epetra_Vector(*StructureField()->DofRowMap(),true));
    IO::DiscretizationReader reader = IO::DiscretizationReader(StructureField()->Discretization(),step);
    reader.ReadVector(lambdafull, "fsilambda");
    lambda_ = StructureField()->Interface()->ExtractFSICondVector(lambdafull);
  }

  StructureField()->ReadRestart(step);
  FluidField().ReadRestart(step);
  AleField().ReadRestart(step);

  SetTimeStep(FluidField().Time(),FluidField().Step());
}

/*----------------------------------------------------------------------*/
/* Recover the Lagrange multiplier at the interface   mayr.mt (03/2012) */
/*----------------------------------------------------------------------*/
void FSI::MonolithicStructureSplit::RecoverLagrangeMultiplier()
{
  // get time integration parameter of structural time integrator
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();

  // some scaling factors for fluid
  const double timescale = FluidField().TimeScaling();

  // some often re-used vectors
  Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::null;     // stores intermediate result of terms (3)-(8)
  Teuchos::RCP<Epetra_Vector> auxvec = Teuchos::null;     // just for convenience
  Teuchos::RCP<Epetra_Vector> auxauxvec = Teuchos::null;  // just for convenience

  /* Recovery of Lagrange multiplier lambda^{n+1} is done by the following
   * condensation expression:
   *
   * lambda^{n+1} =
   *
   * (1)  - stiparam / (1.-stiparam) * lambda^{n}
   *
   * (2)  - 1. / (1.-stiparam) * D^{-T} * tmpvec
   *
   * with tmpvec =
   *
   * (3)    r_{\Gamma}^{S,n+1}
   *
   * (4)  + S_{\Gamma I} * \Delta d_{I}^{S,n+1}
   *
   * (5)  + tau * S_{\Gamma\Gamma} * \Delta u_{\Gamma}^{F,n+1}
   *
   * (6)  + dt * S_{\Gamma\Gamma} * u_{\Gamma}^n]
   *
   * Remark on term (6):
   * Term (6) has to be considered only in the first Newton iteration.
   * Hence, it will usually not be computed since in general we need more
   * than one nonlinear iteration until convergence.
   *
   * Remarks on all terms:
   * +  Division by -(1.0 - stiparam) will be done in the end
   *    since this is common to all terms
   * +  tau: time scaling factor for interface time integration (tau = 1/FluidField().TimeScaling())
   * +  neglecting terms (4)-(6) should not alter the results significantly
   *    since at the end of the time step the solution increments tend to zero.
   *
   *                                                 Matthias Mayr (10/2012)
   */
  int err = 0;
  // ---------Addressing term (1)
  err = lambda_->Update(stiparam,*lambda_,0.0);
  if (err!=0) { dserror("Failed!"); }
  // ---------End of term (1)

  // ---------Addressing term (3)
  Teuchos::RCP<Epetra_Vector> structureresidual = StructureField()->Interface()->ExtractFSICondVector(StructureField()->RHS());
  tmpvec = rcp(new Epetra_Vector(*structureresidual));
  // ---------End of term (3)

  /* Commented out terms (4) to (6) sinve they tend to introdcue oscillations
   * in the Lagrange multiplier field for certain material properties of the
   * structure
   *                                                    Matthias Mayr 11/2012
  // ---------Addressing term (4)
  auxvec = rcp(new Epetra_Vector(sgiprev_->RangeMap(),true));
  err = sgiprev_->Apply(*ddiinc_,*auxvec);
  if (err!=0) { dserror("Failed!"); }
  err = tmpvec->Update(-1.0,*auxvec,1.0);
  if (err!=0) { dserror("Failed!"); }
  // ---------End of term (4)

  // ---------Addressing term (5)
  auxvec = rcp(new Epetra_Vector(sggprev_->RangeMap(),true));
  err = sggprev_->Apply(*FluidToStruct(duginc_),*auxvec);
  if (err!=0) { dserror("Failed!"); }
  err = tmpvec->Update(-1.0/timescale,*auxvec,1.0);
  if (err!=0) { dserror("Failed!"); }
  // ---------End of term (5)

  //---------Addressing term (6)
  if (firstcall_)
  {
    auxvec = rcp(new Epetra_Vector(sggprev_->RangeMap(),true));
    err = sggprev_->Apply(*FluidToStruct(FluidField().ExtractInterfaceVeln()),*auxvec);
    if (err!=0) { dserror("Failed!"); }
    err = tmpvec->Update(-Dt(),*auxvec,1.0);
    if (err!=0) { dserror("Failed!"); }
  }
  // ---------End of term (6)
   *
   */

  // ---------Addressing term (2)
  err = lambda_->Update(1.0,*tmpvec,1.0);
  if (err!=0) { dserror("Failed!"); }
  // ---------End of term (2)

  // finally, divide by -(1.-stiparam) which is common to all terms
  lambda_->Scale(-1./(1.0-stiparam));

  return;
}
