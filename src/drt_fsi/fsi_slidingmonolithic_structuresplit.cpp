/*--------------------------------------------------------------------------*/
/*!
\file fsi_slidingmonolithic_structuresplit.cpp

\brief Solve FSI problem with sliding grids using a monolithic scheme
with condensed structure interface displacements

\level 2

\maintainer Andy Wirtz
*/
/*--------------------------------------------------------------------------*/



#include <Teuchos_TimeMonitor.hpp>

#include "../drt_adapter/adapter_coupling_mortar.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"

#include "fsi_slidingmonolithic_structuresplit.H"
#include "fsi_debugwriter.H"
#include "fsi_overlapprec.H"
#include "fsi_overlapprec_fsiamg.H"
#include "fsi_statustest.H"
#include "fsi_monolithic_linearsystem.H"
#include "fsi_matrixtransform.H"
#include "fsi_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_ale/ale_utils_mapextractor.H"
#include "../drt_adapter/ad_ale_fsi.H"

#include "../drt_constraint/constraint_manager.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::SlidingMonolithicStructureSplit::SlidingMonolithicStructureSplit(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams)
  : BlockMonolithic(comm,timeparams),
    comm_(comm),
    lambda_(Teuchos::null),
    lambdaold_(Teuchos::null),
    energysum_(0.0)
{
  // ---------------------------------------------------------------------------
  // FSI specific check of Dirichlet boundary conditions
  // ---------------------------------------------------------------------------
  // Create intersection of slave DOFs that hold a Dirichlet boundary condition
  // and are located at the FSI interface
  std::vector<Teuchos::RCP<const Epetra_Map> > intersectionmaps;
  intersectionmaps.push_back(StructureField()->GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(StructureField()->Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmap = LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

  // Check whether the intersection is empty
  if (intersectionmap->NumGlobalElements() != 0)
  {
//    std::cout << "Slave interface nodes with Dirichlet boundary condition "
//                 "(input file numbering):" << std::endl;
//    for (int i=0; i < (int)FluidField()->Discretization()->NumMyRowNodes(); i++)
//    {
//      // get all nodes and add them
//      int gid = StructureField()->Discretization()->NodeRowMap()->GID(i);
//
//      // do only nodes that I have in my discretization
//      if (!StructureField()->Discretization()->NodeColMap()->MyGID(gid)) continue;
//      DRT::Node* node = StructureField()->Discretization()->gNode(gid);
//      if (!node) dserror("Cannot find node with gid %",gid);
//
//      std::vector<int> nodedofs = StructureField()->Discretization()->Dof(node);
//
//      for (int j=0; j < (int)nodedofs.size(); j++)
//      {
//        for (int k=0; k < (int)intersectionmap->NumGlobalElements(); k++)
//        {
//          if (nodedofs[j] == intersectionmap->GID(k))
//          {
//            std::cout << gid+1 << std::endl;
//            k = (int)intersectionmap->GID(k);
//            j = (int)nodedofs.size();
//          }
//        }
//      }
//    }

    // It is not allowed, that slave DOFs at the interface hold a Dirichlet
    // boundary condition. Thus --> Error message

    // We do not have to care whether ALE interface DOFs carry DBCs in the
    // input file since they do not occur in the monolithic system and, hence,
    // do not cause a conflict.

    std::stringstream errormsg;
    errormsg  << "  +---------------------------------------------------------------------------------------------+" << std::endl
              << "  |                DIRICHLET BOUNDARY CONDITIONS ON SLAVE SIDE OF FSI INTERFACE                 |" << std::endl
              << "  +---------------------------------------------------------------------------------------------+" << std::endl
              << "  | NOTE: The slave side of the interface is not allowed to carry Dirichlet boundary conditions.|" << std::endl
              << "  |                                                                                             |" << std::endl
              << "  | This is a structure split scheme. Hence, master and slave field are chosen as follows:      |" << std::endl
              << "  |     MASTER  = FLUID                                                                         |" << std::endl
              << "  |     SLAVE   = STRUCTURE                                                                     |" << std::endl
              << "  |                                                                                             |" << std::endl
              << "  | Dirichlet boundary conditions were detected on slave interface degrees of freedom. Please   |" << std::endl
              << "  | remove Dirichlet boundary conditions from the slave side of the FSI interface.              |" << std::endl
              << "  | Only the master side of the FSI interface is allowed to carry Dirichlet boundary conditions.|" << std::endl
              << "  +---------------------------------------------------------------------------------------------+" << std::endl;

    dserror(errormsg.str());
  }
  // ---------------------------------------------------------------------------

  notsetup_ = true;

  coupsfm_ = Teuchos::rcp(new ADAPTER::CouplingMortar());
  fscoupfa_ = Teuchos::rcp(new ADAPTER::Coupling());

  aigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmiitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmgitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fsaigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fsmgitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  SetLambda();
  ddiinc_ = Teuchos::null;
  disiprev_ = Teuchos::null;
  disgprev_ = Teuchos::null;
  sgiprev_ = Teuchos::null;
  sggprev_ = Teuchos::null;

#ifdef DEBUG
  if (coupsfm_ == Teuchos::null) { dserror("Allocation of 'coupsfm_' failed."); }
  if (fscoupfa_ == Teuchos::null) { dserror("Allocation of 'fscoupfa_' failed."); }
  if (aigtransform_ == Teuchos::null) { dserror("Allocation of 'aigtransform_' failed."); }
  if (fmiitransform_ == Teuchos::null) { dserror("Allocation of 'fmiitransform_' failed."); }
  if (fmgitransform_ == Teuchos::null) { dserror("Allocation of 'fmgitransform_' failed."); }
  if (fsaigtransform_ == Teuchos::null) { dserror("Allocation of 'fsaigtransform_' failed."); }
  if (fsmgitransform_ == Teuchos::null) { dserror("Allocation of 'fsmgitransform_' failed."); }
  if (lambda_ == Teuchos::null) { dserror("Allocation of 'lambda_' failed."); }
  if (lambdaold_ == Teuchos::null) { dserror("Allocation of 'lambdaold_' failed."); }
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::SetLambda()
{
  lambda_ = Teuchos::rcp(
      new Epetra_Vector(*StructureField()->Interface()->FSICondMap(), true));
  lambdaold_ = Teuchos::rcp(
      new Epetra_Vector(*StructureField()->Interface()->FSICondMap(), true));

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::SetupSystem()
{
  if (notsetup_)
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
    linearsolverstrategy_ = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

    aleproj_ = DRT::INPUT::IntegralValue<INPAR::FSI::SlideALEProj>(fsidyn, "SLIDEALEPROJ");

    SetDefaultParameters(fsidyn, NOXParameterList());

    // we use non-matching meshes at the interface
    // mortar with: structure = slave, fluid = master

    const int ndim = DRT::Problem::Instance()->NDim();

    // get coupling objects
    ADAPTER::Coupling& icoupfa = InterfaceFluidAleCoupling();

    /* structure to fluid
     * coupling condition at the fsi interface:
     * displacements (=number spatial dimensions) are coupled
     * e.g.: 3D: coupleddof = [1, 1, 1]
     */
    std::vector<int> coupleddof(ndim, 1);

    coupsfm_->Setup(FluidField()->Discretization(),
                    StructureField()->Discretization(),
                    AleField()->WriteAccessDiscretization(),
                    coupleddof,
                    "FSICoupling",
                    comm_,
                    false);

    // fluid to ale at the interface
    icoupfa.SetupConditionCoupling(*FluidField()->Discretization(),
                                   FluidField()->Interface()->FSICondMap(),
                                   *AleField()->Discretization(),
                                   AleField()->Interface()->FSICondMap(),
                                   "FSICoupling",
                                   ndim);

    // we might have a free surface
    if (FluidField()->Interface()->FSCondRelevant())
    {
      fscoupfa_->SetupConditionCoupling(*FluidField()->Discretization(),
                                        FluidField()->Interface()->FSCondMap(),
                                        *AleField()->Discretization(),
                                        AleField()->Interface()->FSCondMap(),
                                        "FREESURFCoupling",
                                        ndim);
    }

    // the fluid-ale coupling always matches
    const Epetra_Map* fluidnodemap = FluidField()->Discretization()->NodeRowMap();
    const Epetra_Map* alenodemap = AleField()->Discretization()->NodeRowMap();

    ADAPTER::Coupling& coupfa = FluidAleCoupling();

    coupfa.SetupCoupling(*FluidField()->Discretization(),
                         *AleField()->Discretization(),
                         *fluidnodemap,
                         *alenodemap,
                          ndim);

    FluidField()->SetMeshMap(coupfa.MasterDofMap());

    // create combined map
    CreateCombinedDofRowMap();

    // Use normal matrix for fluid equations but build (splitted) mesh movement
    // linearization (if requested in the input file)
    FluidField()->UseBlockMatrix(false);

    // Use splitted structure matrix
    StructureField()->UseBlockMatrix();

    // build ale system matrix in splitted system
    AleField()->CreateSystemMatrix(AleField()->Interface());

    aleresidual_ = Teuchos::rcp(new Epetra_Vector(*FsiAleField()->FsiInterface()->OtherMap()));

    // -------------------------------------------------------------------------
    // Build the global Dirichlet map extractor
    SetupDBCMapExtractor();
    // -------------------------------------------------------------------------

    // enable debugging
    if (DRT::INPUT::IntegralValue<int>(fsidyn, "DEBUGOUTPUT") & 2)
    {
      pcdbg_ = Teuchos::rcp(new UTILS::MonolithicDebugWriter(*this));
    }

    CreateSystemMatrix();

    // set up sliding ale if necessary
    if(aleproj_ != INPAR::FSI::ALEprojection_none)
    {
      // MeshInit possibly modifies reference configuration of slave side --> recompute element volume in InitializeElements()
      StructureField()->Discretization()->FillComplete(false,true,true);
      // set up sliding ale utils
      slideale_ = Teuchos::rcp(new FSI::UTILS::SlideAleUtils(StructureField()->Discretization(),
                                                    FluidField()->Discretization(),
                                                    *coupsfm_,
                                                    false,
                                                    aleproj_));

      iprojdisp_ = Teuchos::rcp(new Epetra_Vector(*coupsfm_->MasterDofMap(), true));
      iprojdispinc_ = Teuchos::rcp(new Epetra_Vector(*coupsfm_->MasterDofMap(), true));
    }
    notsetup_=false;
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::CreateCombinedDofRowMap()
{
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;
  vecSpaces.push_back(StructureField()->Interface()->OtherMap());
  vecSpaces.push_back(FluidField()->DofRowMap());
  vecSpaces.push_back(FsiAleField()->FsiInterface()->OtherMap());

  if (vecSpaces[0]->NumGlobalElements() == 0)
    dserror("No inner structural equations. Splitting not possible. Panic.");

  SetDofRowMaps(vecSpaces);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::SetupDBCMapExtractor()
{
  /* Dirichlet maps for structure and fluid do not intersect with interface map.
   * ALE Dirichlet map might intersect with interface map, but ALE interface
   * DOFs are not part of the final system of equations. Hence, we just need the
   * intersection of inner ALE DOFs with Dirichlet ALE DOFs.
   */
  std::vector<Teuchos::RCP<const Epetra_Map> > aleintersectionmaps;
  aleintersectionmaps.push_back(AleField()->GetDBCMapExtractor()->CondMap());
  aleintersectionmaps.push_back(FsiAleField()->FsiInterface()->OtherMap());
  Teuchos::RCP<Epetra_Map> aleintersectionmap =
      LINALG::MultiMapExtractor::IntersectMaps(aleintersectionmaps);

  // Merge Dirichlet maps of structure, fluid and ALE to global FSI Dirichlet map
  std::vector<Teuchos::RCP<const Epetra_Map> > dbcmaps;
  dbcmaps.push_back(StructureField()->GetDBCMapExtractor()->CondMap());
  dbcmaps.push_back(FluidField()->GetDBCMapExtractor()->CondMap());
  dbcmaps.push_back(aleintersectionmap);
  Teuchos::RCP<const Epetra_Map> dbcmap =
      LINALG::MultiMapExtractor::MergeMaps(dbcmaps);

  // Finally, create the global FSI Dirichlet map extractor
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor(*DofRowMap(),dbcmap,true));
  if (dbcmaps_ == Teuchos::null)
    dserror("Creation of FSI Dirichlet map extractor failed.");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase>
FSI::SlidingMonolithicStructureSplit::SystemMatrix() const
{
  return systemmatrix_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::SetupRHSResidual(Epetra_Vector& f)
{
  /* get time integration parameters of structure and fluid time integrators
   * to enable consistent time integration among the fields
   */
  const double stiparam = StructureField()->TimIntParam();
  const double ftiparam = FluidField()->TimIntParam();

  // some scaling factors for fluid
  const double fluidscale = FluidField()->ResidualScaling();

  // get the Mortar matrix M
  const Teuchos::RCP<LINALG::SparseMatrix> mortarp = coupsfm_->GetMortarTrafo();

  // get single field residuals
  Teuchos::RCP<const Epetra_Vector> sv = Teuchos::rcp(new Epetra_Vector(*StructureField()->RHS()));
  Teuchos::RCP<const Epetra_Vector> fv = Teuchos::rcp(new Epetra_Vector(*FluidField()->RHS()));
  Teuchos::RCP<const Epetra_Vector> av = Teuchos::rcp(new Epetra_Vector(*AleField()->RHS()));

  // extract only inner DOFs from structure (=slave) and ALE field
  Teuchos::RCP<const Epetra_Vector> sov = StructureField()->Interface()->ExtractOtherVector(sv);
  Teuchos::RCP<const Epetra_Vector> aov = FsiAleField()->FsiInterface()->ExtractOtherVector(av);

  // add structure interface residual to fluid interface residual considering temporal scaling
  Teuchos::RCP<const Epetra_Vector> scv = StructureField()->Interface()->ExtractFSICondVector(sv);
  Teuchos::RCP<Epetra_Vector> fcv = LINALG::CreateVector(*FluidField()->Interface()->FSICondMap(), true);
  mortarp->Multiply(true, *scv, *fcv);
  Teuchos::RCP<Epetra_Vector> modfv = FluidField()->Interface()->InsertFSICondVector(fcv);
  modfv->Update(1.0, *fv, (1.0 - ftiparam) / ((1.0 - stiparam) * fluidscale));

  // put the single field residuals together
  FSI::Monolithic::CombineFieldVectors(f, sov, modfv, aov);

  // add additional ale residual
  Extractor().AddVector(*aleresidual_, 2, f);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::SetupRHSLambda(Epetra_Vector& f)
{
  if (lambda_ != Teuchos::null)
  {
    // get time integration parameters of structure and fluid time integrators
    // to enable consistent time integration among the fields
    const double stiparam = StructureField()->TimIntParam();
    const double ftiparam = FluidField()->TimIntParam();

    // some scaling factors for fluid
    const double fluidscale = FluidField()->ResidualScaling();

    // get the Mortar matrix M
    Teuchos::RCP<const LINALG::SparseMatrix> mortarm = coupsfm_->GetMMatrix();

    // project Lagrange multiplier field onto the master interface DOFs and consider temporal scaling
    Teuchos::RCP<Epetra_Vector> lambda = Teuchos::rcp(new Epetra_Vector(mortarm->DomainMap(), true));
    mortarm->Multiply(true, *lambda_, *lambda);
    Teuchos::RCP<Epetra_Vector> lambdafull = FluidField()->Interface()->InsertFSICondVector(lambda);
    lambdafull->Scale((-ftiparam + (stiparam * (1.0 - ftiparam)) / (1.0 - stiparam)) / fluidscale);

    // add Lagrange multiplier
    Extractor().AddVector(*lambdafull, 1, f);
  }

  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::SetupRHSFirstiter(Epetra_Vector& f)
{
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();
  const double ftiparam = FluidField()->TimIntParam();

  // some scaling factors for fluid
  const double scale = FluidField()->ResidualScaling();

  // old interface velocity of fluid field
  const Teuchos::RCP<const Epetra_Vector> fveln = FluidField()->ExtractInterfaceVeln();

  // get the Mortar projection matrix P = D^{-1} * M
  const Teuchos::RCP<const LINALG::SparseMatrix> mortarp = coupsfm_->GetMortarTrafo();

  // get fluid shape derivatives matrix
  const Teuchos::RCP<const LINALG::BlockSparseMatrixBase> mmm = FluidField()->ShapeDerivatives();

  // get structure matrix
  const Teuchos::RCP<const LINALG::BlockSparseMatrixBase> blocks = StructureField()->BlockSystemMatrix();

  // get ale matrix
  const Teuchos::RCP<const LINALG::BlockSparseMatrixBase> blocka = AleField()->BlockSystemMatrix();

#ifdef DEBUG
  if (mortarp==Teuchos::null)
    dserror("Expected Teuchos::rcp to mortar matrix P.");
  if (blocks==Teuchos::null)
    dserror("Expected Teuchos::rcp to structure block matrix.");
  if (blocka==Teuchos::null)
    dserror("Expected Teuchos::rcp to ALE block matrix.");
#endif

  // extract submatrices
  const LINALG::SparseMatrix& sig = blocks->Matrix(0,1); // S_{I\Gamma}
  const LINALG::SparseMatrix& sgg = blocks->Matrix(1,1); // S_{\Gamma\Gamma}
  const LINALG::SparseMatrix& aig = blocka->Matrix(0,1); // A_{I\Gamma}

  // some often re-used vectors
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null; // right hand side of single set of DOFs
  Teuchos::RCP<Epetra_Vector> auxvec = Teuchos::null; // just for convenience
  Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::null; // just for convenience

  /* Different contributions/terms to the rhs are separated by the following
   * comment line */
  // ---------- inner structure DOFs
  /* The following terms are added to the inner structure DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  - dt * S_{I \Gamma} * P * u^{n}_{\Gamma}
   *
   * (2)  + S_{I \Gamma} * \Delta d_{\Gamma,p}
   *
   * Remarks on all terms:
   * +  tau: time scaling factor for interface time integration
   *         (tau = 1/FluidField()->TimeScaling())
   *
   */
  // ----------addressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(sig.RangeMap(), true));
  auxvec = Teuchos::rcp(new Epetra_Vector(mortarp->RangeMap(), true));

  mortarp->Apply(*fveln, *auxvec);
  sig.Apply(*auxvec, *rhs);

  rhs->Scale(-Dt());

  Extractor().AddVector(*rhs, 0, f);
  // ----------end of term 1

  // ----------addressing term 2
  rhs = Teuchos::rcp(new Epetra_Vector(sig.RangeMap(), true));

  sig.Apply(*ddgpred_, *rhs);

  Extractor().AddVector(*rhs, 0, f);
  // ----------end of term 2
  // ----------end of inner structure DOFs

  // ---------- inner fluid DOFs
  /* The following terms are added to the inner fluid DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  - dt * F^{G}_{I \Gamma} * u^{n}_{\Gamma}
   *
   */
  // ----------addressing term 1
  if (mmm != Teuchos::null)
  {
    // extract F^{G}_{I \Gamma}
    const LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);

    rhs = Teuchos::rcp(new Epetra_Vector(fmig.RangeMap(), true));

    fmig.Apply(*fveln, *rhs);

    rhs->Scale(-Dt());
    rhs = FsiFluidField()->FsiInterface()->InsertOtherVector(rhs);

    Extractor().AddVector(*rhs, 1, f);
  }
  // ----------end of term 1
  // ----------end of inner fluid DOFs

  // ---------- interface fluid DOFs
  /* The following terms are added to the interface fluid DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  - dt * F^{G}_{\Gamma\Gamma} * u^{n}_{\Gamma}
   *
   * (2)  - (1-ftiparam) / (1-stiparam) * dt * P^{T} * S_{\Gamma\Gamma} * P * u^{n}_{\Gamma}
   *
   * (3)  + (1-ftiparam) / (1-stiparam) * P^{T} * S_{\Gamma\Gamma} * \Delta d_{\Gamma,p}
   *
   * Remarks on all terms:
   * +  tau: time scaling factor for interface time integration (tau = 1/FluidField()->TimeScaling())
   *
   */
  // ----------addressing term 1
  if (mmm != Teuchos::null)
  {
    // extract F^{G}_{\Gamma\Gamma}
    const LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);

    rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RangeMap(), true));

    fmgg.Apply(*fveln, *rhs);

    rhs->Scale(-Dt());
    rhs = FluidField()->Interface()->InsertFSICondVector(rhs);

    Extractor().AddVector(*rhs, 1, f);
  }
  // ----------end of term 1

  // ----------addressing term 2
  rhs = Teuchos::rcp(new Epetra_Vector(mortarp->DomainMap(), true));
  auxvec = Teuchos::rcp(new Epetra_Vector(sgg.RangeMap(), true));
  tmpvec = Teuchos::rcp(new Epetra_Vector(mortarp->RangeMap(), true));

  mortarp->Apply(*fveln, *tmpvec);
  sgg.Apply(*tmpvec, *auxvec);
  mortarp->Multiply(true, *auxvec, *rhs);

  rhs->Scale(-(1. - ftiparam) / (1. - stiparam) * Dt() / scale);
  rhs = FluidField()->Interface()->InsertFSICondVector(rhs);

  Extractor().AddVector(*rhs, 1, f);
  // ----------end of term 2

  // ----------addressing term 3
  rhs = Teuchos::rcp(new Epetra_Vector(mortarp->DomainMap(), true));
  auxvec = Teuchos::rcp(new Epetra_Vector(sgg.RangeMap(), true));

  sgg.Apply(*ddgpred_, *auxvec);
  mortarp->Multiply(true, *auxvec, *rhs);

  rhs->Scale((1. - ftiparam) / (1. - stiparam) / scale);
  rhs = FluidField()->Interface()->InsertFSICondVector(rhs);

  Extractor().AddVector(*rhs, 1, f);
  // ----------end of term 3
  // ----------end of interface fluid DOFs

  // ---------- inner ALE DOFs
  /* The following terms are added to the inner ALE DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  - dt * A_{I \Gamma} * u^{n}_{\Gamma}
   *
   */
  // ----------addressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(aig.RangeMap(), true));

  aig.Apply(*FluidToAleInterface(fveln), *rhs);

  rhs->Scale(-Dt());

  Extractor().AddVector(*rhs, 2, f);
  // ----------end of term 1
  // ----------end of inner ALE DOFs

  // only if relative movement between ale and structure is possible
  if (aleproj_ != INPAR::FSI::ALEprojection_none)
  {
    rhs = Teuchos::rcp(new Epetra_Vector(aig.RowMap(), true));

    aig.Apply(*FluidToAleInterface(iprojdispinc_), *rhs);

    Extractor().AddVector(*rhs, 2, f);
  }

  // if there is a free surface
  if (FluidField()->Interface()->FSCondRelevant())
  {
    // here we extract the free surface submatrices from position 2
    const LINALG::SparseMatrix& afsig = blocka->Matrix(0, 2);

    // extract fluid free surface velocities.
    Teuchos::RCP<Epetra_Vector> fveln = FluidField()->ExtractFreeSurfaceVeln();
    Teuchos::RCP<Epetra_Vector> aveln = FluidToAleInterface(fveln);

    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(afsig.RowMap(), true));
    afsig.Apply(*aveln, *rhs);

    rhs->Scale(-1. * Dt());

    Extractor().AddVector(*rhs, 1, f);

    // shape derivatives
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField()->ShapeDerivatives();
    if (mmm!=Teuchos::null)
    {
      // here we extract the free surface submatrices from position 2
      LINALG::SparseMatrix& fmfsig = mmm->Matrix(0,2);
      LINALG::SparseMatrix& fmfsgg = mmm->Matrix(2,2);

      rhs = Teuchos::rcp(new Epetra_Vector(fmfsig.RowMap(), true));
      fmfsig.Apply(*fveln, *rhs);
      Teuchos::RCP<Epetra_Vector> veln = FsiFluidField()->FsiInterface()->InsertOtherVector(rhs);

      rhs = Teuchos::rcp(new Epetra_Vector(fmfsgg.RowMap(), true));
      fmfsgg.Apply(*fveln, *rhs);
      FluidField()->Interface()->InsertFSCondVector(rhs, veln);

      veln->Scale(-1. * Dt());

      Extractor().AddVector(*veln, 0, f);
    }
  }

  /* Reset quantities for previous iteration step since they still store values
   * from the last time step */
  ddiinc_   = LINALG::CreateVector(*StructureField()->Interface()->OtherMap(), true);
  disiprev_ = Teuchos::null;
  disgprev_ = Teuchos::null;
  sgicur_   = Teuchos::null;
  sggcur_   = Teuchos::null;

  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::SetupSystemMatrix(
    LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::SlidingMonolithicStructureSplit::SetupSystemMatrix");

  // get the Mortar projection matrix P = D^{-1} * M
  Teuchos::RCP<LINALG::SparseMatrix> mortarp = coupsfm_->GetMortarTrafo();

  // get single field block matrices
  const Teuchos::RCP<LINALG::BlockSparseMatrixBase> s = StructureField()->BlockSystemMatrix();
  const Teuchos::RCP<LINALG::SparseMatrix> f = FluidField()->SystemMatrix();
  const Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField()->BlockSystemMatrix();


#ifdef DEBUG
  // check whether allocation was successful
  if (mortarp == Teuchos::null)
    dserror("Expected Teuchos::rcp to mortar matrix P.");
  if (s == Teuchos::null) { dserror("expect structure block matrix"); }
  if (f == Teuchos::null) { dserror("expect fluid matrix"); }
  if (a == Teuchos::null) { dserror("expect ale block matrix"); }

  // some checks whether maps for matrix-matrix-multiplication do really match
  if (!s->Matrix(0,1).DomainMap().PointSameAs(mortarp->RangeMap()))
    dserror("Maps do not match.");
  if (!s->Matrix(1,0).RangeMap().PointSameAs(mortarp->RangeMap()))
    dserror("Maps do not match.");
  if (!s->Matrix(1,1).DomainMap().PointSameAs(mortarp->RangeMap()))
    dserror("Maps do not match.");
#endif

  // extract submatrices
  LINALG::SparseMatrix& aii = a->Matrix(0,0);
  const LINALG::SparseMatrix& aig = a->Matrix(0,1);

  // scaling factors for fluid
  const double scale = FluidField()->ResidualScaling();
  const double timescale = FluidField()->TimeScaling();

  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();
  const double ftiparam = FluidField()->TimIntParam();

  // Uncomplete fluid matrix to be able to deal with slightly defective
  // interface meshes.
  f->UnComplete();

  // ---------------------------------------------------------------------------
  // BEGIN building the global 4x4 system matrix
  // ---------------------------------------------------------------------------
  // Contributions to blocks in system matrix are listed separately.
  // Block numbering in comments ranges from (1,1) to (4,4).

  mat.Assign(0,0,LINALG::View,s->Matrix(0,0));

  // ----------Addressing contribution to block (1,3)
  Teuchos::RCP<LINALG::SparseMatrix> sig =
      MLMultiply(s->Matrix(0,1), false, *mortarp, false, false,  false, true);
  Teuchos::RCP<LINALG::SparseMatrix> lsig =
      Teuchos::rcp(new LINALG::SparseMatrix(sig->RowMap(), 81, false));

  lsig->Add(*sig, false, 1. / timescale, 0.0);
  lsig->Complete(f->DomainMap(), sig->RangeMap());

  mat.Assign(0, 1, LINALG::View, *lsig);

  // ----------Addressing contribution to block (3,1)
  Teuchos::RCP<LINALG::SparseMatrix> sgi =
      MLMultiply(*mortarp, true, s->Matrix(1,0), false, false, false, true);
  Teuchos::RCP<LINALG::SparseMatrix> lsgi =
      Teuchos::rcp(new LINALG::SparseMatrix(f->RowMap(), 81, false));

  lsgi->Add(*sgi, false, (1. - ftiparam) / ((1. - stiparam) * scale), 0.0);
  lsgi->Complete(sgi->DomainMap(), f->RangeMap());

  mat.Assign(1, 0, LINALG::View, *lsgi);

  // ----------Addressing contribution to block (3,3)
  Teuchos::RCP<LINALG::SparseMatrix> sgg =
      MLMultiply(s->Matrix(1,1), false, *mortarp, false, false, false, true);
  sgg = MLMultiply(*mortarp, true, *sgg, false, false, false, true);

  f->Add(*sgg, false, (1. - ftiparam) / ((1. - stiparam) * scale * timescale), 1.0);
  mat.Assign(1, 1, LINALG::View, *f);

  (*aigtransform_)(a->FullRowMap(),
                   a->FullColMap(),
                   aig,
                   1./timescale,
                   ADAPTER::CouplingSlaveConverter(InterfaceFluidAleCoupling()),
                   mat.Matrix(2,1));
  mat.Assign(2, 2, LINALG::View, aii);

  /*--------------------------------------------------------------------------*/
  // add optional fluid linearization with respect to mesh motion block

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField()->ShapeDerivatives();
  if (mmm != Teuchos::null)
  {
    // extract submatrices
    const LINALG::SparseMatrix& fmii = mmm->Matrix(0,0);
    const LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);
    const LINALG::SparseMatrix& fmgi = mmm->Matrix(1,0);
    const LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);

    // ----------Addressing contribution to block (3,3)
    mat.Matrix(1, 1).Add(fmgg, false, 1. / timescale, 1.0);

    // ----------Addressing contribution to block (2,3)
    mat.Matrix(1, 1).Add(fmig, false, 1. / timescale, 1.0);

    const ADAPTER::Coupling& coupfa = FluidAleCoupling();

    (*fmgitransform_)(mmm->FullRowMap(),
                      mmm->FullColMap(),
                      fmgi,
                      1.,
                      ADAPTER::CouplingMasterConverter(coupfa),
                      mat.Matrix(1,2),
                      false,
                      false);

    (*fmiitransform_)(mmm->FullRowMap(),
                      mmm->FullColMap(),
                      fmii,
                      1.,
                      ADAPTER::CouplingMasterConverter(coupfa),
                      mat.Matrix(1,2),
                      false,
                      true);
  }

  // if there is a free surface
  if (FluidField()->Interface()->FSCondRelevant())
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

      mat.Matrix(1, 1).Add(fmgg, false, 1. / timescale, 1.0);
      mat.Matrix(1, 1).Add(fmig, false, 1. / timescale, 1.0);

      const ADAPTER::Coupling& coupfa = FluidAleCoupling();

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

  // done. make sure all blocks are filled.
  mat.Complete();

  // Finally, take care of Dirichlet boundary conditions
  mat.ApplyDirichlet(*(dbcmaps_->CondMap()), true);
  //
  // ---------------------------------------------------------------------------
  // END building the global system matrix
  // ---------------------------------------------------------------------------

  // store parts of structural matrix to know them in the next iteration as previous iteration matrices
  sgiprev_ = sgicur_;
  sggprev_ = sggcur_;
  sgicur_ = Teuchos::rcp(new LINALG::SparseMatrix(s->Matrix(1,0)));
  sggcur_ = Teuchos::rcp(new LINALG::SparseMatrix(s->Matrix(1,1)));
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::Update()
{

  // update history variables for sliding ale
  if (aleproj_!= INPAR::FSI::ALEprojection_none)
  {
    iprojdisp_ = Teuchos::rcp(new Epetra_Vector(*coupsfm_->MasterDofMap(),true));
    Teuchos::RCP<Epetra_Vector> idispale =
        AleToFluidInterface(AleField()->Interface()->ExtractFSICondVector(AleField()->Dispnp()));

    slideale_->Remeshing(*StructureField(),
                         FluidField()->Discretization(),
                         idispale,
                         iprojdisp_,
                         *coupsfm_,
                         Comm());

    iprojdispinc_->Update(-1.0,*iprojdisp_,1.0,*idispale,0.0);

    slideale_->EvaluateMortar(StructureField()->ExtractInterfaceDispnp(), iprojdisp_, *coupsfm_);
    slideale_->EvaluateFluidMortar(idispale,iprojdisp_);

    Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*iprojdisp_));
    temp->ReplaceMap(idispale->Map());
    Teuchos::RCP<Epetra_Vector> acx = FluidToAleInterface(temp);
    AleField()->ApplyInterfaceDisplacements(acx);
    FluidField()->ApplyMeshDisplacement(AleToFluid(AleField()->Dispnp()));

    Teuchos::RCP<Epetra_Vector> unew =
        slideale_->InterpolateFluid(FluidField()->ExtractInterfaceVelnp());
    FluidField()->ApplyInterfaceVelocities(unew);
  }

  // call Update()-routine in base class to handle the single fields
  FSI::MonolithicBase::Update();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::ScaleSystem(
    LINALG::BlockSparseMatrixBase& mat,
    Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool)DRT::INPUT::IntegralValue<int>(fsimono,"INFNORMSCALING");

  if (scaling_infnorm)
  {
    // The matrices are modified here. Do we have to change them back later on?

    // do scaling of structure rows
    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0,0).EpetraMatrix();
    srowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    scolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*srowsum_);
    A->InvColSums(*scolsum_);
    if (A->LeftScale(*srowsum_) or
        A->RightScale(*scolsum_) or
        mat.Matrix(0,1).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0,2).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(1,0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(2,0).EpetraMatrix()->RightScale(*scolsum_))
      dserror("structure scaling failed");

    // do scaling of ale rows
    A = mat.Matrix(2,2).EpetraMatrix();
    arowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    acolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*arowsum_);
    A->InvColSums(*acolsum_);
    if (A->LeftScale(*arowsum_) or
        A->RightScale(*acolsum_) or
        mat.Matrix(2,0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2,1).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0,2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(1,2).EpetraMatrix()->RightScale(*acolsum_))
      dserror("ale scaling failed");

    // do scaling of structure and ale rhs vectors
    Teuchos::RCP<Epetra_Vector> sx = Extractor().ExtractVector(b,0);
    Teuchos::RCP<Epetra_Vector> ax = Extractor().ExtractVector(b,2);

    if (sx->Multiply(1.0, *srowsum_, *sx, 0.0))
      dserror("structure scaling failed");
    if (ax->Multiply(1.0, *arowsum_, *ax, 0.0))
      dserror("ale scaling failed");

    Extractor().InsertVector(*sx, 0, b);
    Extractor().InsertVector(*ax, 2, b);
  }
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::UnscaleSolution(
    LINALG::BlockSparseMatrixBase& mat,
    Epetra_Vector& x,
    Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool)DRT::INPUT::IntegralValue<int>(fsimono,"INFNORMSCALING");

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
  mat.Apply(x, r);
  r.Update(1., b, 1.);

  Teuchos::RCP<Epetra_Vector> sr = Extractor().ExtractVector(r,0);
  Teuchos::RCP<Epetra_Vector> fr = Extractor().ExtractVector(r,1);
  Teuchos::RCP<Epetra_Vector> ar = Extractor().ExtractVector(r,2);

  // increment additional ale residual
  aleresidual_->Update(-1., *ar, 0.);

  std::ios_base::fmtflags flags = Utils()->out().flags();

  double n,ns,nf,na;
  r.Norm2(&n);
  sr->Norm2(&ns);
  fr->Norm2(&nf);
  ar->Norm2(&na);
  Utils()->out() << std::scientific
                 << "\nlinear solver quality:\n"
                 << "L_2-norms:\n"
                 << "   |r|=" << n
                 << "   |rs|=" << ns
                 << "   |rf|=" << nf
                 << "   |ra|=" << na
                 << "\n";
  r.NormInf(&n);
  sr->NormInf(&ns);
  fr->NormInf(&nf);
  ar->NormInf(&na);
  Utils()->out() << "L_inf-norms:\n"
                 << "   |r|=" << n
                 << "   |rs|=" << ns
                 << "   |rf|=" << nf
                 << "   |ra|=" << na
                 << "\n";

  Utils()->out().flags(flags);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem>
FSI::SlidingMonolithicStructureSplit::CreateLinearSystem(
    Teuchos::ParameterList& nlParams,
    NOX::Epetra::Vector& noxSoln,
    Teuchos::RCP<NOX::Utils> utils)
{
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys;

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
//  Teuchos::ParameterList* lsParams = NULL;
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

//  // in case of nonlinCG the linear solver list is somewhere else
//  if (dirParams.get("Method","User Defined")=="User Defined")
//    lsParams = &(newtonParams.sublist("Linear Solver"));
//  else if (dirParams.get("Method","User Defined")=="NonlinearCG")
//    lsParams = &(dirParams.sublist("Nonlinear CG").sublist("Linear Solver"));
//  else dserror("Unknown nonlinear method");

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


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo>
FSI::SlidingMonolithicStructureSplit::CreateStatusTest(
    Teuchos::ParameterList& nlParams,
    Teuchos::RCP<NOX::Epetra::Group> grp)
{
  // ---------------------------------------------------------------------------
  // Setup the test framework
  // ---------------------------------------------------------------------------
  // Create the top-level test combo
  Teuchos::RCP<NOX::StatusTest::Combo> combo
    = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  // Create test combo for convergence of residuals and iterative increments
  Teuchos::RCP<NOX::StatusTest::Combo> converged
    = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  // Create some other plausibility tests
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters
    = Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get<int>("Max Iterations")));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv
    = Teuchos::rcp(new NOX::StatusTest::FiniteValue);

  // Add single tests to the top-level test combo
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Start filling the 'converged' combo here
  // require one solve
  converged->addStatusTest(Teuchos::rcp(new NOX::FSI::MinIters(1)));


  // ---------------------------------------------------------------------------
  // setup tests for structural displacement field
  // ---------------------------------------------------------------------------
  // create NOX::StatusTest::Combo for structural displacement field
  Teuchos::RCP<NOX::StatusTest::Combo> structcombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> structureDisp_L2 =
    Teuchos::rcp(new NOX::FSI::PartialNormF("DISPL residual",
                                            Extractor(),0,
                                            nlParams.get<double>("Tol dis res L2"),
                                            NOX::Abstract::Vector::TwoNorm,
                                            NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> structureDisp_inf =
    Teuchos::rcp(new NOX::FSI::PartialNormF("DISPL residual",
                                            Extractor(),0,
                                            nlParams.get<double>("Tol dis res Inf"),
                                            NOX::Abstract::Vector::MaxNorm,
                                            NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate_L2 =
    Teuchos::rcp(new NOX::FSI::PartialNormUpdate("DISPL update",
                                                 Extractor(),0,
                                                 nlParams.get<double>("Tol dis inc L2"),
                                                 NOX::Abstract::Vector::TwoNorm,
                                                 NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate_inf =
    Teuchos::rcp(new NOX::FSI::PartialNormUpdate("DISPL update",
                                                 Extractor(),0,
                                                 nlParams.get<double>("Tol dis inc Inf"),
                                                 NOX::Abstract::Vector::MaxNorm,
                                                 NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  AddStatusTest(structureDisp_L2);

  // add norm-tests to structural displacement NOX::StatusTest::Combo
  structcombo->addStatusTest(structureDisp_L2);
  structcombo->addStatusTest(structureDisp_inf);
  structcombo->addStatusTest(structureDispUpdate_L2);
  structcombo->addStatusTest(structureDispUpdate_inf);

  // add structural displacement test combo to top-level test combo
  converged->addStatusTest(structcombo);
  // ---------- end of structural displacement field tests

  // ---------------------------------------------------------------------------
  // setup tests for interface
  // ---------------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map> > interface;
  interface.push_back(FluidField()->Interface()->FSICondMap());
  interface.push_back(Teuchos::null);
  LINALG::MultiMapExtractor interfaceextract(*DofRowMap(),interface);

  // create NOX::StatusTest::Combo for interface
  Teuchos::RCP<NOX::StatusTest::Combo> interfacecombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> interfaceTest_L2 =
    Teuchos::rcp(new NOX::FSI::PartialNormF("GAMMA residual",
                                            interfaceextract,0,
                                            nlParams.get<double>("Tol fsi res L2"),
                                            NOX::Abstract::Vector::TwoNorm,
                                            NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> interfaceTest_inf =
    Teuchos::rcp(new NOX::FSI::PartialNormF("GAMMA residual",
                                            interfaceextract,0,
                                            nlParams.get<double>("Tol fsi res Inf"),
                                            NOX::Abstract::Vector::MaxNorm,
                                            NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> interfaceTestUpdate_L2 =
    Teuchos::rcp(new NOX::FSI::PartialNormUpdate("GAMMA update",
                                                 interfaceextract,0,
                                                 nlParams.get<double>("Tol fsi inc L2"),
                                                 NOX::Abstract::Vector::TwoNorm,
                                                 NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> interfaceTestUpdate_inf =
    Teuchos::rcp(new NOX::FSI::PartialNormUpdate("GAMMA update",
                                                 interfaceextract,0,
                                                 nlParams.get<double>("Tol fsi inc Inf"),
                                                 NOX::Abstract::Vector::MaxNorm,
                                                 NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  AddStatusTest(interfaceTest_L2);

  // add norm-tests to interface NOX::StatusTest::Combo
  interfacecombo->addStatusTest(interfaceTest_L2);
  interfacecombo->addStatusTest(interfaceTest_inf);
  interfacecombo->addStatusTest(interfaceTestUpdate_L2);
  interfacecombo->addStatusTest(interfaceTestUpdate_inf);

  // add interface test combo to top-level test combo
  converged->addStatusTest(interfacecombo);
  // ---------- end of interface tests

  // ---------------------------------------------------------------------------
  // setup tests for fluid velocity field
  // ---------------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map> > fluidvel;
  fluidvel.push_back(FluidField()->InnerVelocityRowMap());
  fluidvel.push_back(Teuchos::null);
  LINALG::MultiMapExtractor fluidvelextract(*DofRowMap(),fluidvel);

  // create NOX::StatusTest::Combo for fluid velocity field
  Teuchos::RCP<NOX::StatusTest::Combo> fluidvelcombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> innerFluidVel_L2 =
    Teuchos::rcp(new NOX::FSI::PartialNormF("VELOC residual",
                                            fluidvelextract,0,
                                            nlParams.get<double>("Tol vel res L2"),
                                            NOX::Abstract::Vector::TwoNorm,
                                            NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> innerFluidVel_inf =
    Teuchos::rcp(new NOX::FSI::PartialNormF("VELOC residual",
                                            fluidvelextract,0,
                                            nlParams.get<double>("Tol vel res Inf"),
                                            NOX::Abstract::Vector::MaxNorm,
                                            NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> innerFluidVelUpdate_L2 =
    Teuchos::rcp(new NOX::FSI::PartialNormUpdate("VELOC update",
                                                 fluidvelextract,0,
                                                 nlParams.get<double>("Tol vel inc L2"),
                                                 NOX::Abstract::Vector::TwoNorm,
                                                 NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> innerFluidVelUpdate_inf =
    Teuchos::rcp(new NOX::FSI::PartialNormUpdate("VELOC update",
                                                 fluidvelextract,0,
                                                 nlParams.get<double>("Tol vel inc Inf"),
                                                 NOX::Abstract::Vector::MaxNorm,
                                                 NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  AddStatusTest(innerFluidVel_L2);

  // add norm-tests to fluid velocity NOX::StatusTest::Combo
  fluidvelcombo->addStatusTest(innerFluidVel_L2);
  fluidvelcombo->addStatusTest(innerFluidVel_inf);
  fluidvelcombo->addStatusTest(innerFluidVelUpdate_L2);
  fluidvelcombo->addStatusTest(innerFluidVelUpdate_inf);

  // add fluid velocity test combo to top-level test combo
  converged->addStatusTest(fluidvelcombo);
  // ---------- end of fluid velocity field tests

  // ---------------------------------------------------------------------------
  // setup tests for fluid pressure field
  // ---------------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map> > fluidpress;
  fluidpress.push_back(FluidField()->PressureRowMap());
  fluidpress.push_back(Teuchos::null);
  LINALG::MultiMapExtractor fluidpressextract(*DofRowMap(),fluidpress);

  // create NOX::StatusTest::Combo for fluid pressure field
  Teuchos::RCP<NOX::StatusTest::Combo> fluidpresscombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> fluidPress_L2 =
    Teuchos::rcp(new NOX::FSI::PartialNormF("PRESS residual",
                                            fluidpressextract,0,
                                            nlParams.get<double>("Tol pre res L2"),
                                            NOX::Abstract::Vector::TwoNorm,
                                            NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> fluidPress_inf =
    Teuchos::rcp(new NOX::FSI::PartialNormF("PRESS residual",
                                            fluidpressextract,0,
                                            nlParams.get<double>("Tol pre res Inf"),
                                            NOX::Abstract::Vector::MaxNorm,
                                            NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> fluidPressUpdate_L2 =
    Teuchos::rcp(new NOX::FSI::PartialNormUpdate("PRESS update",
                                                 fluidpressextract,0,
                                                 nlParams.get<double>("Tol pre inc L2"),
                                                 NOX::Abstract::Vector::TwoNorm,
                                                 NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> fluidPressUpdate_inf =
    Teuchos::rcp(new NOX::FSI::PartialNormUpdate("PRESS update",
                                                 fluidpressextract,0,
                                                 nlParams.get<double>("Tol pre inc Inf"),
                                                 NOX::Abstract::Vector::MaxNorm,
                                                 NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  AddStatusTest(fluidPress_L2);

  // add norm-tests to fluid pressure NOX::StatusTest::Combo
  fluidpresscombo->addStatusTest(fluidPress_L2);
  fluidpresscombo->addStatusTest(fluidPress_inf);
  fluidpresscombo->addStatusTest(fluidPressUpdate_L2);
  fluidpresscombo->addStatusTest(fluidPressUpdate_inf);

  // add fluid pressure test combo to top-level test combo
  converged->addStatusTest(fluidpresscombo);
  // ---------- end of fluid pressure field tests

  // Finally, return the test combo
  return combo;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::ExtractFieldVectors(
    Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx,
    Teuchos::RCP<const Epetra_Vector>& fx,
    Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::SlidingMonolithicStructureSplit::ExtractFieldVectors");

#ifdef DEBUG
  if(ddgpred_ == Teuchos::null)
    dserror("Vector 'ddgpred_' has not been initialized properly.");
#endif

  // get the Mortar projection matrix P = D^{-1} * M
  const Teuchos::RCP<const LINALG::SparseMatrix> mortarp = coupsfm_->GetMortarTrafo();

  // ---------------------------------------------------------------------------
  // process fluid unknowns
  // ---------------------------------------------------------------------------
  // extract fluid solution increment from NOX increment
  Teuchos::RCP<Epetra_Vector> f = Extractor().ExtractVector(x,1);
  FluidField()->UpdateSlaveDOF(f);
  fx = f;

  // ---------------------------------------------------------------------------
  // process ale unknowns
  // ---------------------------------------------------------------------------
  // extract inner ALE solution increment from NOX increment
  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x,2);

  // convert fluid interface velocities into ALE interface displacements
  Teuchos::RCP<Epetra_Vector> fcx = FluidField()->Interface()->ExtractFSICondVector(fx);
  FluidField()->VelocityToDisplacement(fcx);
  Teuchos::RCP<Epetra_Vector> acx = FluidToAleInterface(fcx);

  // put inner and interface ALE solution increments together
  Teuchos::RCP<Epetra_Vector> a = FsiAleField()->FsiInterface()->InsertOtherVector(aox);
  AleField()->Interface()->InsertFSICondVector(acx, a);

  // if there is a free surface
  if (FluidField()->Interface()->FSCondRelevant())
  {
    Teuchos::RCP<Epetra_Vector> fcx = FluidField()->Interface()->ExtractFSCondVector(fx);
    FluidField()->FreeSurfVelocityToDisplacement(fcx);

    Teuchos::RCP<Epetra_Vector> acx = fscoupfa_->MasterToSlave(fcx);
    AleField()->Interface()->InsertFSCondVector(acx, a);
  }

  AleField()->UpdateSlaveDOF(a);
  ax = a;

  // ---------------------------------------------------------------------------
  // process structure unknowns
  // ---------------------------------------------------------------------------
  // extract inner structure solution increment from NOX increment
  Teuchos::RCP<const Epetra_Vector> sox = Extractor().ExtractVector(x,0);

  // convert ALE interface displacements to structure interface displacements
  Teuchos::RCP<Epetra_Vector> scx =
      LINALG::CreateVector(*StructureField()->Interface()->FSICondMap());
  acx = AleToFluidInterface(acx);
  mortarp->Apply(*acx,*scx);
  scx->Update(-1.0, *ddgpred_, 1.0);

  // put inner and interface structure solution increments together
  Teuchos::RCP<Epetra_Vector> s = StructureField()->Interface()->InsertOtherVector(sox);
  StructureField()->Interface()->InsertFSICondVector(scx, s);
  sx = s;

  // ---------------------------------------------------------------------------

  // Store field vectors to know them later on as previous quantities
  if (disiprev_ != Teuchos::null)
    ddiinc_->Update(1.0, *sox, -1.0, *disiprev_, 0.0);  // compute current iteration increment
  else
    ddiinc_ = Teuchos::rcp(new Epetra_Vector(*sox));    // first iteration increment

  disiprev_ = sox;                                      // store current step increment

  if (velgprev_ != Teuchos::null)
    duginc_->Update(1.0, *fcx, -1.0, *velgprev_, 0.0);  // compute current iteration increment
  else
    duginc_ = Teuchos::rcp(new Epetra_Vector(*fcx));    // first iteration increment

  velgprev_ = fcx;                                      // store current step increment
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::Output()
{
  StructureField()->Output();

  // output Lagrange multiplier
  OutputLambda();

  FluidField()->    Output();

  if (aleproj_!= INPAR::FSI::ALEprojection_none)
  {
    int uprestart = timeparams_.get<int>("RESTARTEVRY");
    if (uprestart != 0 && FluidField()->Step() % uprestart == 0)
    {
      FluidField()->DiscWriter()->WriteVector("slideALE", iprojdisp_);
      FluidField()->DiscWriter()->WriteVector("slideALEincr", iprojdispinc_);
      slideale_->OutputRestart(*FluidField()->DiscWriter());
    }
  }
  AleField()->Output();
  FluidField()->LiftDrag();

  if (StructureField()->GetConstraintManager()->HaveMonitor())
  {
    StructureField()->GetConstraintManager()->ComputeMonitorValues(StructureField()->Dispnp());
    if(comm_.MyPID() == 0)
      StructureField()->GetConstraintManager()->PrintMonitorValues();
  }

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::OutputLambda()
{
  /* 'lambda_' is only defined on the interface. So, insert 'lambda_' into

   * output or restart data.
   */
  Teuchos::RCP<Epetra_Vector> lambdafull = StructureField()->Interface()->InsertFSICondVector(lambda_);
  const int uprestart = timeparams_.get<int>("RESTARTEVRY");
  const int upres = timeparams_.get<int>("RESULTSEVRY");
  if ((uprestart != 0 && FluidField()->Step() % uprestart == 0) || FluidField()->Step() % upres == 0)
    StructureField()->DiscWriter()->WriteVector("fsilambda", lambdafull);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::ReadRestart(int step)
{
  // read Lagrange multiplier
  {
    Teuchos::RCP<Epetra_Vector> lambdafull = Teuchos::rcp(new Epetra_Vector(*StructureField()->DofRowMap(),true));
    IO::DiscretizationReader reader = IO::DiscretizationReader(StructureField()->Discretization(),step);
    reader.ReadVector(lambdafull, "fsilambda");
    lambda_ = StructureField()->Interface()->ExtractFSICondVector(lambdafull);
  }

  StructureField()->ReadRestart(step);
  FluidField()->ReadRestart(step);

  SetupSystem();

  if (aleproj_!= INPAR::FSI::ALEprojection_none)
  {
    IO::DiscretizationReader reader =
        IO::DiscretizationReader(FluidField()->Discretization(),step);
    reader.ReadVector(iprojdisp_, "slideALE");
    reader.ReadVector(iprojdispinc_, "slideALEincr");
    slideale_->ReadRestart(reader);
  }

  AleField()->ReadRestart(step);

  SetTimeStep(FluidField()->Time(),FluidField()->Step());

  if (aleproj_!= INPAR::FSI::ALEprojection_none)
    slideale_->EvaluateMortar(StructureField()->ExtractInterfaceDispn(), iprojdisp_, *coupsfm_);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::RecoverLagrangeMultiplier()
{
  // get time integration parameter of structural time integrator
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();

  // some scaling factors for fluid
//  const double timescale = FluidField()->TimeScaling();

  // get the Mortar projection matrix P = D^{-1} * M
  const Teuchos::RCP<LINALG::SparseMatrix> mortarp = coupsfm_->GetMortarTrafo();

  // get the inverted Mortar matrix D^{-1}
  const Teuchos::RCP<LINALG::SparseMatrix> mortardinv = coupsfm_->GetDinvMatrix();

#ifdef DEBUG
  if (mortarp == Teuchos::null)
    dserror("Expected Teuchos::rcp to mortar matrix P.");
  if (mortardinv == Teuchos::null)
    dserror("Expected Teuchos::rcp to mortar matrix D^{-1}.");
#endif

  // some often re-used vectors
  Teuchos::RCP<Epetra_Vector> tmpvec    = Teuchos::null;  // stores intermediate result of terms (3)-(8)
  Teuchos::RCP<Epetra_Vector> auxvec    = Teuchos::null;  // just for convenience
  Teuchos::RCP<Epetra_Vector> auxauxvec = Teuchos::null;  // just for convenience

  /* Recovery of Lagrange multiplier lambda^{n+1} is done by the following
   * condensation expression:
   *
   * lambda^{n+1} =
   *
   * (1)  - stiparam / (1.-stiparam) * lambda^{n}
   *
   * (2)  + 1. / (1.-stiparam) * D^{-T} * tmpvec
   *
   * with tmpvec =
   *
   * (3)    r_{\Gamma}^{S,n+1}
   *
   * (4)  + S_{\Gamma I} * \Delta d_{I}^{S,n+1}
   *
   * (5)  + tau * S_{\Gamma\Gamma} * P * \Delta u_{\Gamma}^{F,n+1}
   *
   * (6)  + dt * S_{\Gamma\Gamma} * P * u_{\Gamma}^n]
   *
   * Remark on term (6):
   * Term (6) has to be considered only in the first Newton iteration.
   * Hence, it will usually not be computed since in general we need more
   * than one nonlinear iteration until convergence.
   *
   * Remarks on all terms:
   * +  Division by -(1.0 - stiparam) will be done in the end
   *    since this is common to all terms
   * +  tau: time scaling factor for interface time integration (tau = 1/FluidField()->TimeScaling())
   * +  neglecting terms (4)-(6) should not alter the results significantly
   *    since at the end of the time step the solution increments tend to zero.
   *
   *                                                 Matthias Mayr (10/2012)
   */

  // ---------Addressing term (1)
  lambda_->Scale(-stiparam);
  // ---------End of term (1)

  // ---------Addressing term (3)
  Teuchos::RCP<Epetra_Vector> structureresidual
    = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->ExtractFSICondVector(StructureField()->RHS())));
  structureresidual->Scale(-1.0); // invert sign to obtain residual, not rhs
  tmpvec = Teuchos::rcp(new Epetra_Vector(*structureresidual));
  // ---------End of term (3)

  /* You might want to comment out terms (4) to (6) since they tend to
   * introduce oscillations in the Lagrange multiplier field for certain
   * material properties of the structure.
   *                                                    Matthias Mayr 11/2012
  // ---------Addressing term (4)
  auxvec = Teuchos::rcp(new Epetra_Vector(sgiprev_->RangeMap(),true));
  sgiprev_->Apply(*ddiinc_,*auxvec);
  tmpvec->Update(1.0,*auxvec,1.0);
  // ---------End of term (4)

  // ---------Addressing term (5)
  auxvec = Teuchos::rcp(new Epetra_Vector(mortarp->RangeMap(),true));
  mortarp->Apply(*duginc_,*auxvec);
  auxauxvec = Teuchos::rcp(new Epetra_Vector(sggprev_->RangeMap(),true));
  sggprev_->Apply(*auxvec,*auxauxvec);
  tmpvec->Update(1.0/timescale,*auxauxvec,1.0);
  // ---------End of term (5)

  // ---------Addressing term (6)
  if (firstcall_)
  {
    auxvec = Teuchos::rcp(new Epetra_Vector(mortarp->RangeMap(),true));
    mortarp->Apply(*FluidField()->ExtractInterfaceVeln(),*auxvec);
    auxauxvec = Teuchos::rcp(new Epetra_Vector(sggprev_->RangeMap(),true));
    sggprev_->Apply(*auxvec,*auxauxvec);
    tmpvec->Update(Dt(),*auxauxvec,1.0);
  }
  // ---------End of term (6)
   *
   */

  // ---------Addressing term (2)
  auxvec = Teuchos::rcp(new Epetra_Vector(mortardinv->DomainMap(),true));
  mortardinv->Multiply(true,*tmpvec,*auxvec);
  lambda_->Update(1.0,*auxvec,1.0);
  // ---------End of term (2)

  // finally, divide by -(1.-stiparam) which is common to all terms
  lambda_->Scale(1./(1.0-stiparam));

  /* Finally, the Lagrange multiplier lambda_ is recovered here. It has the
   * unit [N/m^2]. Actual nodal forces are obtained by multiplication with
   * mortar matrices M or D later on.
   */

//  CheckKinematicConstraint();
//  CheckDynamicEquilibrium();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::CalculateInterfaceEnergyIncrement()
{
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();
  const double ftiparam = FluidField()->TimIntParam();

  // interface traction weighted by time integration factors
  Teuchos::RCP<Epetra_Vector> tractionstructure =
      Teuchos::rcp(new Epetra_Vector(lambda_->Map(), true));
  tractionstructure->Update(stiparam-ftiparam, *lambdaold_, ftiparam-stiparam,
      *lambda_, 0.0);

  // displacement increment of this time step
  Teuchos::RCP<Epetra_Vector> deltad =
      Teuchos::rcp(new Epetra_Vector(*StructureField()->DofRowMap(), true));
  deltad->Update(1.0, *StructureField()->Dispnp(), -1.0,
      *StructureField()->Dispn(), 0.0);

  // calculate the energy increment
  double energy = 0.0;
  tractionstructure->Dot(*StructureField()->Interface()->ExtractFSICondVector(deltad), &energy);

  energysum_ += energy;

  WriteInterfaceEnergyFile(energy, energysum_);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::CheckKinematicConstraint()
{
  // some scaling factors for fluid
  const double timescale  = FluidField()->TimeScaling();

  // get the Mortar matrices D and M
  const Teuchos::RCP<LINALG::SparseMatrix> mortard = coupsfm_->GetDMatrix();
  const Teuchos::RCP<LINALG::SparseMatrix> mortarm = coupsfm_->GetMMatrix();

  // get interface displacements and velocities
  Teuchos::RCP<Epetra_Vector> disnp = StructureField()->ExtractInterfaceDispnp();
  Teuchos::RCP<Epetra_Vector> disn = StructureField()->ExtractInterfaceDispn();
  Teuchos::RCP<Epetra_Vector> velnp = FluidField()->ExtractInterfaceVelnp();
  Teuchos::RCP<Epetra_Vector> veln = FluidField()->ExtractInterfaceVeln();

  // prepare vectors for projected interface quantities
  Teuchos::RCP<Epetra_Vector> disnpproj = Teuchos::rcp(new Epetra_Vector(mortard->RangeMap(), true));
  Teuchos::RCP<Epetra_Vector> disnproj = Teuchos::rcp(new Epetra_Vector(mortard->RangeMap(), true));
  Teuchos::RCP<Epetra_Vector> velnpproj = Teuchos::rcp(new Epetra_Vector(mortarm->RangeMap(), true));
  Teuchos::RCP<Epetra_Vector> velnproj = Teuchos::rcp(new Epetra_Vector(mortarm->RangeMap(), true));

  // projection of interface displacements
  mortard->Apply(*disnp,*disnpproj);
  mortard->Apply(*disn,*disnproj);

  // projection of interface velocities
  mortarm->Apply(*velnp,*velnpproj);
  mortarm->Apply(*veln,*velnproj);

  // calculate violation of kinematic interface constraint
  Teuchos::RCP<Epetra_Vector> violation = Teuchos::rcp(new Epetra_Vector(*disnpproj));
  violation->Update(-1.0, *disnproj, 1.0);
  violation->Update(-1.0/timescale, *velnpproj, 1.0/timescale, *velnproj, 1.0);
  violation->Update(-Dt(), *velnproj, 1.0);

  // calculate some norms
  double violationl2 = 0.0;
  double violationinf = 0.0;
  violation->Norm2(&violationl2);
  violation->NormInf(&violationinf);

  // scale L2-Norm with length of vector
  violationl2 /= sqrt(violation->MyLength());

  // output to screen
  std::ios_base::fmtflags flags = Utils()->out().flags();

  Utils()->out() << std::scientific
                 << "\nViolation of kinematic interface constraint:\n"
                 << "L_2-norm: "
                 << violationl2
                 << "        L_inf-norm: "
                 << violationinf
                 << "\n";
  Utils()->out().flags(flags);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::CheckDynamicEquilibrium()
{
  // get the Mortar matrices D and M
  const Teuchos::RCP<LINALG::SparseMatrix> mortard = coupsfm_->GetDMatrix();
  const Teuchos::RCP<LINALG::SparseMatrix> mortarm = coupsfm_->GetMMatrix();

#ifdef DEBUG
  if (mortarm == Teuchos::null)
    dserror("Expected Teuchos::rcp to mortar matrix M.");
  if (mortard == Teuchos::null)
    dserror("Expected Teuchos::rcp to mortar matrix D.");
#endif

  // auxiliary vectors
  Teuchos::RCP<Epetra_Vector> tractionmaster =
      Teuchos::rcp(new Epetra_Vector(mortarm->DomainMap(),true));
  Teuchos::RCP<Epetra_Vector> tractionslave =
      Teuchos::rcp(new Epetra_Vector(mortard->DomainMap(),true));

  // calculate forces on master and slave side
  mortarm->Multiply(true,*lambda_,*tractionmaster);
  mortard->Multiply(true,*lambda_,*tractionslave);

  // calculate violation of dynamic equilibrium
  Teuchos::RCP<Epetra_Vector> violation =
      Teuchos::rcp(new Epetra_Vector(*tractionmaster));
  violation->Update(-1.0,*tractionslave,1.0);

  // calculate some norms
  double violationl2 = 0.0;
  double violationinf = 0.0;
  violation->Norm2(&violationl2);
  violation->NormInf(&violationinf);

  // scale L2-Norm with sqrt of length of interface vector
  violationl2 /= sqrt(StructureField()->Interface()->FSICondMap()->NumGlobalElements());

  // output to screen
  std::ios_base::fmtflags flags = Utils()->out().flags();

  Utils()->out() << std::scientific
                 << "\nViolation of dynamic interface equilibrium:\n"
                 << "L_2-norm: "
                 << violationl2
                 << "        L_inf-norm: "
                 << violationinf
                 << "\n";
  Utils()->out().flags(flags);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::CombineFieldVectors(Epetra_Vector& v,
    Teuchos::RCP<const Epetra_Vector> sv,
    Teuchos::RCP<const Epetra_Vector> fv,
    Teuchos::RCP<const Epetra_Vector> av,
    bool fullvectors)
{
  if (fullvectors)
  {
    // extract inner DOFs from slave vectors
    Teuchos::RCP<Epetra_Vector> sov =
        StructureField()->Interface()->ExtractOtherVector(sv);
    Teuchos::RCP<Epetra_Vector> aov =
        FsiAleField()->FsiInterface()->ExtractOtherVector(av);

    // put them together
    Extractor().AddVector(*sov,0,v);
    Extractor().AddVector(*fv,1,v);
    Extractor().AddVector(*aov,2,v);
  }
  else
    FSI::Monolithic::CombineFieldVectors(v,sv,fv,av);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double FSI::SlidingMonolithicStructureSplit::SelectDtErrorBased() const
{
  // get time step size suggestions based on some error norms
  const double dtfl = GetAdaFlDt(); // based on all fluid DOFs
  const double dtflfsi = GetAdaFlFSIDt(); // based on fluid FSI DOFs
  const double dtstrinner = GetAdaStrInnerDt(); // based on inner structural DOFs

  double dt = Dt();

  // select time step size based on error estimation
  if (IsAdaStructure() and IsAdaFluid())
    dt = std::min(std::min(dtfl, dtflfsi), dtstrinner);
  else if (IsAdaStructure() and (not IsAdaFluid()))
    dt = dtstrinner;
  else if((not IsAdaStructure()) and IsAdaFluid())
    dt = std::min(dtfl, dtflfsi);
  else
  {
    // no change in time step size based on structure or fluid field error estimation
  }

  return dt;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool FSI::SlidingMonolithicStructureSplit::SetAccepted() const
{
  // get error norms
  const double flnorm = GetAdaFlnorm(); // based on all fluid DOFs
  const double flfsinorm = GetAdaFlFSInorm(); // based on fluid FSI DOFs
  const double strinnernorm = GetAdaStrInnernorm(); // based on inner structural DOFs

  bool accepted = std::max(flnorm,flfsinorm) < errtolfl_ && strinnernorm < errtolstr_;

  // in case error estimation in the fluid field is turned off:
  if (not IsAdaFluid())
    accepted = strinnernorm < errtolstr_;

  // in case error estimation in the structure field is turned off:
  if (not IsAdaStructure())
    accepted = std::max(flnorm,flfsinorm) < errtolfl_;

  // no error based time adaptivity
  if ((not IsAdaStructure()) and (not IsAdaFluid()))
    accepted = true;

  return accepted;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::CreateSystemMatrix()
{
  FSI::BlockMonolithic::CreateSystemMatrix(systemmatrix_,true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::CreateNodeOwnerRelationship(
    std::map<int, int>* nodeOwner,
    std::map<int, std::list<int> >* inverseNodeOwner,
    std::map<int, DRT::Node*>* structurenodesPtr,
    std::map<int, DRT::Node*>* fluidgnodesPtr,
    Teuchos::RCP<DRT::Discretization> structuredis,
    Teuchos::RCP<DRT::Discretization> fluiddis,
    const INPAR::FSI::Redistribute domain)
{
  /*******************************************/
  /* distribute masternodes to future owners */
  /*******************************************/
  /* Idea:
   * - the P matrix maps dofs from fluid to structure (fluidsplit) and from structure to fluid (structuresplit)
   * - find "neighboring dofs" via entries in row of P matrix
   * - get owner of the node of the corresponding slave dof and global node
   *   id of the node of the corresponding master dof
   * - save this pair in map nodeOwner: global id of master node -> desired owner of this node
   */

  int numproc = comm_.NumProc();
  int myrank = comm_.MyPID();

  // get P matrix
  Teuchos::RCP<LINALG::SparseMatrix> P = coupsfm_->GetMortarTrafo(); // P matrix that couples structure dofs with fluid dofs
  Teuchos::RCP<Epetra_CrsMatrix> P_;
  P_ = P->EpetraMatrix();
  const Epetra_Map& P_Map = P_->RowMap();         // maps fluid dofs to procs

  int NumMyElements = P_Map.NumMyElements();
  int numFluidDofs = P_->NumGlobalCols();               // number of related fluid dofs on interface

  int NumMaxElements;
  comm_.MaxAll(&NumMyElements, &NumMaxElements, 1);

  int skip = DRT::Problem::Instance()->NDim();  // Only evaluate one dof per node, skip the other dofs related to the node. It is nsd_=skip.

  int re[2];            // return: re[0] = global node id, re[1] = owner of node
  int flnode;           // fluid
  int flowner;
  int stnode;           // structure
  int stowner = -1;


  std::map<int,int>::iterator nodeOwnerIt;

  // loop over fluid dofs, i.e. rows in the P matrix
  /* We loop over NumMaxElements because all procs have to stay in the loop until
   * the last proc checked each of its rows.
   */


  for (int stdofLID = 0; stdofLID < NumMaxElements; stdofLID = stdofLID + skip)
  {

    flnode = -1;
    flowner = -1;
    stnode = -1;
    stowner = -1;

    int NumEntries = 0;
    double Values[numFluidDofs];
    int fldofGID[numFluidDofs];

    if (NumMyElements > 0)
    {
      int stdofGID = P_Map.GID(stdofLID); // gid of structure dof
      // find related node and owner to stdofGID
      FindNodeRelatedToDof(structurenodesPtr, stdofGID, structuredis, re);   //fluid
      stnode = re[0];
      stowner = re[1];

      P_->ExtractGlobalRowCopy(stdofGID, numFluidDofs, NumEntries, Values, fldofGID);
    }

    // Loop over related structure dofs and get related nodes.

    int maxNumEntries;
    comm_.MaxAll(&NumEntries,&maxNumEntries,1);

    for (int j = 0; j < maxNumEntries; ++j){
      if (j < NumEntries){
        FindNodeRelatedToDof(fluidgnodesPtr, fldofGID[j], fluiddis, re);
        flnode = re[0];
        flowner = re[1];
      }

      // A processor can only find a neighbored node if it belongs to its domain (which in general is not the case).
      // Therefore, we have to search in the other processors' domains.

      int copy;
      int foundNode = -2;
      int sendNode = -2;
      int saveNode = -2;
      int foundOwner;
      int sendOwner;
      int saveOwner = -2;
      int dofid;

      for (int proc = 0; proc < numproc; ++proc){
        copy = flnode;
        comm_.Broadcast(&copy, 1, proc);    // send nodeid to check if the processor needs help to find its neighbor, code: copy=-2
        if (copy == -2){
          dofid = fldofGID[j];
          comm_.Broadcast(&dofid, 1, proc);
          FindNodeRelatedToDof(fluidgnodesPtr, dofid, fluiddis,re); // let each processor look for the node related to gstid
          foundNode = re[0];
          foundOwner = re[1];
          for (int j = 0; j < numproc; ++j){
            sendNode = foundNode;
            sendOwner = foundOwner;
            comm_.Broadcast(&sendNode, 1, j);
            comm_.Broadcast(&sendOwner, 1, j);
            if (sendNode != -2){ // check which processor found the node
              saveNode = sendNode;
              saveOwner = sendOwner;
              break;
            }
          }
          if (myrank == proc && saveNode != -2){ // save the nodegid on the respective processor
            stnode = saveNode;
            stowner = saveOwner;
          }
        }
      }

      if (domain == INPAR::FSI::Redistribute_structure && stnode != -1){   // map structure nodes to fluid owners
          (*nodeOwner)[stnode]=flowner;
      }
      else if (domain == INPAR::FSI::Redistribute_fluid)  // map fluid nodes to structure owners
        (*nodeOwner)[flnode]=stowner;

    } //for (int j = 0; j < maxNumEntries; ++j){

    if (NumMyElements != 0){
      NumMyElements -= skip;
    }

  } // end loop over structure dofs


//    std::map<int,int>::iterator nodeOwnerPrint;
//    for (nodeOwnerPrint = nodeOwner.begin(); nodeOwnerPrint != nodeOwner.end(); ++nodeOwnerPrint){
//      std::cout<<"\nNode: "<<nodeOwnerPrint->first<<" Owner: "<<nodeOwnerPrint->second<<" I am proc "<<myrank;
//    }


  // If the structure is redistributed, it might occur that one node is contained several times in the
  // list because several procs might have introduced it.
  for (int proc=0; proc<numproc; ++proc){
    //std::map<int,int>::iterator nodeOwnerIt;   already declared
    nodeOwnerIt = nodeOwner->begin();
    int nodeSize = (int)nodeOwner->size();
    comm_.Broadcast(&nodeSize,1,proc);
    int node;
    for (int i=0; i<nodeSize; ++i){
      node = nodeOwnerIt->first;
      comm_.Broadcast(&node,1,proc);
      if (myrank > proc){

        try{
           nodeOwner->at(node);
           nodeOwner->erase(node);
          }
         catch (std::exception& exc)
        {
        }
      }
      if (myrank==proc){
        nodeOwnerIt++;
      }
    }

  }

//  std::map<int,int>::iterator nodeOwnerPrint;
//  for (nodeOwnerPrint = nodeOwner.begin(); nodeOwnerPrint != nodeOwner.end(); ++nodeOwnerPrint){
//    std::cout<<"\nNode: "<<nodeOwnerPrint->first<<" Owner: "<<nodeOwnerPrint->second<<" I am proc "<<comm_.MyPID();
//  }

  std::list<int>::iterator listIt;
    int sendPair[2];
    nodeOwnerIt = nodeOwner->begin();    // already declared

    for (int proc=0; proc<numproc; ++proc){
      int numNodes = nodeOwner->size();
      comm_.Broadcast(&numNodes,1,proc);

      for (int i=0; i<numNodes; ++i){
        if (myrank==proc){
          sendPair[0]=nodeOwnerIt->first;
          sendPair[1]=nodeOwnerIt->second;
        }
        comm_.Broadcast(sendPair,2,proc);
        listIt = (*inverseNodeOwner)[sendPair[1]].begin();
        (*inverseNodeOwner)[sendPair[1]].insert(listIt,sendPair[0]);

        if (myrank==proc)
          nodeOwnerIt++;
      }
    }

}
