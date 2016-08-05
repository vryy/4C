/*--------------------------------------------------------------------------*/
/*!
\file fsi_slidingmonolithic_fluidsplit.cpp

\brief Solve FSI problem with sliding grids using a monolithic scheme
with condensed fluid interface velocities

\level 2

\maintainer Andy Wirtz
*/
/*--------------------------------------------------------------------------*/


#include <Teuchos_TimeMonitor.hpp>

#include "../drt_adapter/adapter_coupling_mortar.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

#include "fsi_slidingmonolithic_fluidsplit.H"
#include "fsi_debugwriter.H"
#include "fsi_statustest.H"
#include "fsi_overlapprec.H"
#include "fsi_overlapprec_fsiamg.H"
#include "fsi_monolithic_linearsystem.H"
#include "fsi_matrixtransform.H"
#include "fsi_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_structure/stru_aux.H"
#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../drt_ale/ale_utils_mapextractor.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"
#include "../drt_adapter/ad_ale_fsi.H"

#include "../drt_constraint/constraint_manager.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../linalg/linalg_mapextractor.H"

#include <math.h>

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
FSI::SlidingMonolithicFluidSplit::SlidingMonolithicFluidSplit(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams) :
    BlockMonolithic(comm, timeparams), comm_(comm), lambda_(Teuchos::null), lambdaold_(
        Teuchos::null), energysum_(0.0)
{
  // ---------------------------------------------------------------------------
  // FSI specific check of Dirichlet boundary conditions
  // ---------------------------------------------------------------------------
  // Create intersection of slave DOFs that hold a Dirichlet boundary condition
  // and are located at the FSI interface
  std::vector<Teuchos::RCP<const Epetra_Map> > intersectionmaps;
  intersectionmaps.push_back(FluidField()->GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(FluidField()->Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmap =
      LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

  // Check whether the intersection is empty
  if (intersectionmap->NumGlobalElements() != 0)
  {
//    std::cout << "Slave interface nodes with Dirichlet boundary condition "
//              "(input file numbering):" << std::endl;
//    for (int i=0; i < (int)FluidField()->Discretization()->NumMyRowNodes(); i++)
//    {
//      // get all nodes and add them
//      int gid = FluidField()->Discretization()->NodeRowMap()->GID(i);
//
//      // do only nodes that I have in my discretization
//      if (!FluidField()->Discretization()->NodeColMap()->MyGID(gid)) continue;
//      DRT::Node* node = FluidField()->Discretization()->gNode(gid);
//      if (!node) dserror("Cannot find node with gid %",gid);
//
//      std::vector<int> nodedofs = FluidField()->Discretization()->Dof(node);
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
    errormsg
        << "  +---------------------------------------------------------------------------------------------+"
        << std::endl
        << "  |                DIRICHLET BOUNDARY CONDITIONS ON SLAVE SIDE OF FSI INTERFACE                 |"
        << std::endl
        << "  +---------------------------------------------------------------------------------------------+"
        << std::endl
        << "  | NOTE: The slave side of the interface is not allowed to carry Dirichlet boundary conditions.|"
        << std::endl
        << "  |                                                                                             |"
        << std::endl
        << "  | This is a fluid split scheme. Hence, master and slave field are chosen as follows:          |"
        << std::endl
        << "  |     MASTER  = STRUCTURE                                                                     |"
        << std::endl
        << "  |     SLAVE   = FLUID                                                                         |"
        << std::endl
        << "  |                                                                                             |"
        << std::endl
        << "  | Dirichlet boundary conditions were detected on slave interface degrees of freedom. Please   |"
        << std::endl
        << "  | remove Dirichlet boundary conditions from the slave side of the FSI interface.              |"
        << std::endl
        << "  | Only the master side of the FSI interface is allowed to carry Dirichlet boundary conditions.|"
        << std::endl
        << "  +---------------------------------------------------------------------------------------------+"
        << std::endl;

    dserror(errormsg.str());
  }
  // ---------------------------------------------------------------------------

  notsetup_ = true;

  coupsfm_ = Teuchos::rcp(new ADAPTER::CouplingMortar());
  fscoupfa_ = Teuchos::rcp(new ADAPTER::Coupling());

  aigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmiitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  coupsfm_ = Teuchos::rcp(new ADAPTER::CouplingMortar());
  fscoupfa_ = Teuchos::rcp(new ADAPTER::Coupling());

  aigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmiitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  // Recovery of Lagrange multiplier happens on fluid field
  SetLambda();

  fmgiprev_ = Teuchos::null;
  fmgicur_ = Teuchos::null;
  fmggprev_ = Teuchos::null;
  fmggcur_ = Teuchos::null;
  fgiprev_ = Teuchos::null;
  fgicur_ = Teuchos::null;
  fggprev_ = Teuchos::null;
  fggcur_ = Teuchos::null;

#ifdef DEBUG
  if (coupsfm_ == Teuchos::null)
  { dserror("Allocation of 'coupsfm_' failed.");}
  if (fscoupfa_ == Teuchos::null)
  { dserror("Allocation of 'fscoupfa_' failed.");}
  if (aigtransform_ == Teuchos::null)
  { dserror("Allocation of 'aigtransform_' failed.");}
  if (fmiitransform_ == Teuchos::null)
  { dserror("Allocation of 'fmiitransform_' failed.");}
  if (lambda_ == Teuchos::null)
  { dserror("Allocation of 'lambda_' failed.");}
  if (lambdaold_ == Teuchos::null)
  { dserror("Allocation of 'lambdaold_' failed.");}
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::SetLambda()
{
  lambda_ = Teuchos::rcp(
      new Epetra_Vector(*FluidField()->Interface()->FSICondMap(), true));
  lambdaold_ = Teuchos::rcp(
      new Epetra_Vector(*FluidField()->Interface()->FSICondMap(), true));

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::SetupSystem()
{
  if (notsetup_)
  {
    const Teuchos::ParameterList& fsidyn =
        DRT::Problem::Instance()->FSIDynamicParams();
    const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
    linearsolverstrategy_ = DRT::INPUT::IntegralValue<
        INPAR::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

    aleproj_ = DRT::INPUT::IntegralValue<INPAR::FSI::SlideALEProj>(fsidyn,
        "SLIDEALEPROJ");

    SetDefaultParameters(fsidyn, NOXParameterList());

    // we use non-matching meshes at the interface
    // mortar with: structure = master, fluid = slave

    const int ndim = DRT::Problem::Instance()->NDim();

    // get coupling objects
    ADAPTER::Coupling& icoupfa = InterfaceFluidAleCoupling();

    /* structure to fluid
     * coupling condition at the fsi interface:
     * displacements (=number spacial dimensions) are coupled
     * e.g.: 3D: coupleddof = [1, 1, 1]
     */
    std::vector<int> coupleddof(ndim, 1);

    coupsfm_->Setup(StructureField()->Discretization(),
        FluidField()->Discretization(), AleField()->WriteAccessDiscretization(),
        coupleddof, "FSICoupling", comm_, true);

    // fluid to ale at the interface

    icoupfa.SetupConditionCoupling(*FluidField()->Discretization(),
        FluidField()->Interface()->FSICondMap(), *AleField()->Discretization(),
        AleField()->Interface()->FSICondMap(), "FSICoupling", ndim);

    // we might have a free surface
    if (FluidField()->Interface()->FSCondRelevant())
    {
      fscoupfa_->SetupConditionCoupling(*FluidField()->Discretization(),
          FluidField()->Interface()->FSCondMap(), *AleField()->Discretization(),
          AleField()->Interface()->FSCondMap(), "FREESURFCoupling", ndim);
    }

    ADAPTER::Coupling& coupfa = FluidAleCoupling();

    // the fluid-ale coupling always matches
    const Epetra_Map* fluidnodemap =
        FluidField()->Discretization()->NodeRowMap();
    const Epetra_Map* alenodemap = AleField()->Discretization()->NodeRowMap();

    coupfa.SetupCoupling(*FluidField()->Discretization(),
        *AleField()->Discretization(), *fluidnodemap, *alenodemap, ndim);

    FluidField()->SetMeshMap(coupfa.MasterDofMap());

    // create combined map
    CreateCombinedDofRowMap();

    /*------------------------------------------------------------------------*/
    // Switch fluid to interface split block matrix
    FluidField()->UseBlockMatrix(true);

    // build ale system matrix in splitted system
    AleField()->CreateSystemMatrix(AleField()->Interface());

    aleresidual_ = Teuchos::rcp(
        new Epetra_Vector(*FsiAleField()->FsiInterface()->OtherMap()));

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

    if (aleproj_ != INPAR::FSI::ALEprojection_none)
    {
      // set up sliding ale utils
      slideale_ = Teuchos::rcp(
          new FSI::UTILS::SlideAleUtils(StructureField()->Discretization(),
              FluidField()->Discretization(), *coupsfm_, true, aleproj_));

      iprojdispinc_ = Teuchos::rcp(
          new Epetra_Vector(*coupsfm_->SlaveDofMap(), true));
      iprojdisp_ = Teuchos::rcp(
          new Epetra_Vector(*coupsfm_->SlaveDofMap(), true));
    }
    notsetup_ = false;
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::CreateCombinedDofRowMap()
{
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;
  vecSpaces.push_back(StructureField()->DofRowMap());
  vecSpaces.push_back(FluidField()->DofRowMap());
  vecSpaces.push_back(FsiAleField()->FsiInterface()->OtherMap());

  if (vecSpaces[1]->NumGlobalElements() == 0)
    dserror("No inner fluid equations. Splitting not possible.");

  SetDofRowMaps(vecSpaces);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::SetupDBCMapExtractor()
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
  Teuchos::RCP<const Epetra_Map> dbcmap = LINALG::MultiMapExtractor::MergeMaps(
      dbcmaps);

  // Finally, create the global FSI Dirichlet map extractor
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor(*DofRowMap(), dbcmap, true));
  if (dbcmaps_ == Teuchos::null)
    dserror("Creation of FSI Dirichlet map extractor failed.");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase>
FSI::SlidingMonolithicFluidSplit::SystemMatrix() const
{
  return systemmatrix_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::SetupRHSResidual(Epetra_Vector& f)
{
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();
  const double ftiparam = FluidField()->TimIntParam();

  // some scaling factors for fluid
  const double fluidscale = FluidField()->ResidualScaling();

  // get the Mortar matrix M
  const Teuchos::RCP<LINALG::SparseMatrix> mortarp = coupsfm_->GetMortarTrafo();

  // get single field residuals
  Teuchos::RCP<const Epetra_Vector> sv = Teuchos::rcp(
      new Epetra_Vector(*StructureField()->RHS()));
  Teuchos::RCP<const Epetra_Vector> fv = Teuchos::rcp(
      new Epetra_Vector(*FluidField()->RHS()));
  Teuchos::RCP<const Epetra_Vector> av = Teuchos::rcp(
      new Epetra_Vector(*AleField()->RHS()));

  // extract only inner DOFs from fluid (=slave) and ALE field
  Teuchos::RCP<Epetra_Vector> fov =
      FsiFluidField()->FsiInterface()->ExtractOtherVector(fv);
  fov = FsiFluidField()->FsiInterface()->InsertOtherVector(fov);
  Teuchos::RCP<const Epetra_Vector> aov =
      FsiAleField()->FsiInterface()->ExtractOtherVector(av);

  /* add fluid interface residual to structure interface residual considering
   * temporal scaling
   */
  Teuchos::RCP<Epetra_Vector> fcv =
      FluidField()->Interface()->ExtractFSICondVector(fv);
  Teuchos::RCP<Epetra_Vector> scv = LINALG::CreateVector(
      *StructureField()->Interface()->FSICondMap(), true);
  mortarp->Multiply(true, *fcv, *scv);
  Teuchos::RCP<Epetra_Vector> modsv =
      StructureField()->Interface()->InsertFSICondVector(scv);
  modsv->Update(1.0, *sv, (1.0 - stiparam) / (1.0 - ftiparam) * fluidscale);

  if (StructureField()->GetSTCAlgo() == INPAR::STR::stc_currsym)
  {
    Teuchos::RCP<LINALG::SparseMatrix> stcmat = StructureField()->GetSTCMat();
    stcmat->Multiply(true, *modsv, *modsv);
  }

  // put the single field residuals together
  FSI::Monolithic::CombineFieldVectors(f, modsv, fov, aov);

  // add additional ale residual
  Extractor().AddVector(*aleresidual_, 2, f);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::SetupRHSLambda(Epetra_Vector& f)
{
  if (lambda_ != Teuchos::null)
  {
    // get time integration parameters of structure and fluid time integrators
    // to enable consistent time integration among the fields
    const double stiparam = StructureField()->TimIntParam();
    const double ftiparam = FluidField()->TimIntParam();

    // get the Mortar matrix M
    const Teuchos::RCP<LINALG::SparseMatrix> mortarm = coupsfm_->GetMMatrix();

    /* project Lagrange multiplier field onto the master interface DOFs and
     * consider temporal scaling */
    Teuchos::RCP<Epetra_Vector> lambda = Teuchos::rcp(
        new Epetra_Vector(mortarm->DomainMap(), true));
    mortarm->Multiply(true, *lambda_, *lambda);
    Teuchos::RCP<Epetra_Vector> lambdafull =
        StructureField()->Interface()->InsertFSICondVector(lambda);
    lambdafull->Scale(
        stiparam - (ftiparam * (1.0 - stiparam)) / (1.0 - ftiparam));

    // add Lagrange multiplier
    Extractor().AddVector(*lambdafull, 0, f);
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::SetupRHSFirstiter(Epetra_Vector& f)
{
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();
  const double ftiparam = FluidField()->TimIntParam();

  // some scaling factors for fluid
  const double timescale = FluidField()->TimeScaling();
  const double scale = FluidField()->ResidualScaling();

  // old interface velocity of fluid field
  const Teuchos::RCP<const Epetra_Vector> fveln =
      FluidField()->ExtractInterfaceVeln();

  // get the Mortar projection matrix P = D^{-1} * M
  const Teuchos::RCP<const LINALG::SparseMatrix> mortarp =
      coupsfm_->GetMortarTrafo();

  // get fluid matrix
  const Teuchos::RCP<const LINALG::BlockSparseMatrixBase> blockf =
      FluidField()->BlockSystemMatrix();

  // get fluid shape derivatives matrix
  const Teuchos::RCP<const LINALG::BlockSparseMatrixBase> mmm =
      FluidField()->ShapeDerivatives();

  // get ale matrix
  const Teuchos::RCP<const LINALG::BlockSparseMatrixBase> blocka =
      AleField()->BlockSystemMatrix();

#ifdef DEBUG
  if (mortarp ==Teuchos::null)
  dserror("Expected Teuchos::rcp to mortar matrix P.");
  if (blockf ==Teuchos::null)
  dserror("Expected Teuchos::rcp to fluid block matrix.");
  if (blocka ==Teuchos::null)
  dserror("Expected Teuchos::rcp to ale block matrix.");
#endif

  // extract fluid and ale submatrices
  const LINALG::SparseMatrix& fig = blockf->Matrix(0, 1); // F_{I\Gamma}
  const LINALG::SparseMatrix& fgg = blockf->Matrix(1, 1); // F_{\Gamma\Gamma}
  const LINALG::SparseMatrix& aig = blocka->Matrix(0, 1); // A_{I\Gamma}

  // some often re-used vectors
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null; // right hand side of single set of DOFs
  Teuchos::RCP<Epetra_Vector> auxvec = Teuchos::null; // just for convenience
  Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::null; // just for convenience

  /* Different contributions/terms to the rhs are separated by the following
   * comment line */
  // ---------- interface structure DOFs
  /* The following terms are added to the interface structure DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  + (1-stiparam)/(1-ftiparam) * dt / tau * P^{T} * F_{\Gamma\Gamma} * u^{n}_{\Gamma}
   *
   * (2)  - (1-stiparam)/(1-ftiparam) / tau * P^{T} * F_{\Gamma\Gamma} * P * \Delta d_{\Gamma,p}
   *
   * (3)  - (1-stiparam)/(1-ftiparam) * P^{T} * F^{G}_{\Gamma\Gamma} * P * \Delta d_{\Gamma,p}
   *
   * Remarks on all terms:
   * +  tau: time scaling factor for interface time integration (tau = 1/FluidField()->TimeScaling())
   *
   */
  // ----------addressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(mortarp->DomainMap(), true));
  auxvec = Teuchos::rcp(new Epetra_Vector(fgg.RowMap(), true));

  fgg.Apply(*fveln, *auxvec);
  mortarp->Multiply(true, *auxvec, *rhs);

  if (StructureField()->GetSTCAlgo() == INPAR::STR::stc_currsym)
  {
    Teuchos::RCP<LINALG::SparseMatrix> stcmat = StructureField()->GetSTCMat();
    stcmat->Multiply(true, *rhs, *rhs);
  }

  rhs->Scale(scale * (1. - stiparam) / (1. - ftiparam) * Dt() * timescale);
  rhs = StructureField()->Interface()->InsertFSICondVector(rhs);

  Extractor().AddVector(*rhs, 0, f);
  // ----------end of term 1

  // ----------addressing term 2
  rhs = Teuchos::rcp(new Epetra_Vector(mortarp->DomainMap(), true));
  auxvec = Teuchos::rcp(new Epetra_Vector(fgg.RangeMap(), true));
  tmpvec = Teuchos::rcp(new Epetra_Vector(mortarp->RangeMap(), true));

  mortarp->Apply(*ddgpred_, *tmpvec);
  fgg.Apply(*tmpvec, *auxvec);
  mortarp->Multiply(true, *auxvec, *rhs);

  rhs->Scale(-scale * (1. - stiparam) / (1. - ftiparam) * timescale);
  rhs = StructureField()->Interface()->InsertFSICondVector(rhs);

  Extractor().AddVector(*rhs, 0, f);
  // ----------end of term 2

  // ----------addressing term 3
  if (mmm != Teuchos::null)
  {
    // extract F^{G}_{\Gamma\Gamma}
    const LINALG::SparseMatrix& fmgg = mmm->Matrix(1, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(mortarp->DomainMap(), true));
    auxvec = Teuchos::rcp(new Epetra_Vector(fmgg.RangeMap(), true));
    tmpvec = Teuchos::rcp(new Epetra_Vector(mortarp->RangeMap(), true));

    mortarp->Apply(*ddgpred_, *tmpvec);
    fmgg.Apply(*tmpvec, *auxvec);
    mortarp->Multiply(true, *auxvec, *rhs);

    rhs->Scale(-(1. - stiparam) / (1. - ftiparam));
    rhs = StructureField()->Interface()->InsertFSICondVector(rhs);

    Extractor().AddVector(*rhs, 0, f);
  }
  // ----------end of term 3
  // ----------end of interface structure DOFs

  // ---------- inner fluid DOFs
  /* The following terms are added to the inner fluid DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  + dt / tau * F_{I \Gamma} * u^{n}_{\Gamma}
   *
   * (2)  - 1 / tau * F_{I \Gamma} * P * \Delta d_{\Gamma,p}
   *
   * (3)  - F^{G}_{I \Gamma} * P * \Delta d_{\Gamma,p}
   *
   * Remarks on all terms:
   * +  tau: time scaling factor for interface time integration
   *         (tau = 1/FluidField()->TimeScaling())
   *
   */
  // ----------addressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(fig.RowMap(), true));

  fig.Apply(*fveln, *rhs);

  rhs->Scale(Dt() * timescale);

  rhs = FsiFluidField()->FsiInterface()->InsertOtherVector(rhs);

  Extractor().AddVector(*rhs, 1, f);
  // ----------end of term 1

  // ----------addressing term 2
  rhs = Teuchos::rcp(new Epetra_Vector(fig.RangeMap(), true));
  auxvec = Teuchos::rcp(new Epetra_Vector(mortarp->RangeMap(), true));

  mortarp->Apply(*ddgpred_, *auxvec);
  fig.Apply(*auxvec, *rhs);

  rhs->Scale(-timescale);

  rhs = FsiFluidField()->FsiInterface()->InsertOtherVector(rhs);

  Extractor().AddVector(*rhs, 1, f);
  // ----------end of term 2

  // ----------addressing term 3
  if (mmm != Teuchos::null)
  {
    // extract F^{G}_{I \Gamma}
    const LINALG::SparseMatrix& fmig = mmm->Matrix(0, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(fmig.RangeMap(), true));
    auxvec = Teuchos::rcp(new Epetra_Vector(mortarp->RangeMap(), true));

    mortarp->Apply(*ddgpred_, *auxvec);
    fmig.Apply(*auxvec, *rhs);

    rhs->Scale(-1.);

    rhs = FsiFluidField()->FsiInterface()->InsertOtherVector(rhs);

    Extractor().AddVector(*rhs, 1, f);
  }
  // ----------end of term 3
  // ----------end of inner fluid DOFs

  // ---------- inner ALE DOFs
  /* The following terms are added to the inner ALE DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  - A_{I \Gamma} * P * \Delta d_{\Gamma,p}
   *
   */
  // ----------addressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(aig.RangeMap(), true));
  auxvec = Teuchos::rcp(new Epetra_Vector(mortarp->RangeMap(), true));

  mortarp->Apply(*ddgpred_, *auxvec);
  aig.Apply(*FluidToAleInterface(auxvec), *rhs);

  rhs->Scale(-1.0);

  Extractor().AddVector(*rhs, 2, f);
  // ----------end of term 1
  // ----------end of inner ALE DOFs

  // only if relative movement between ale and structure is possible
  if (aleproj_ != INPAR::FSI::ALEprojection_none)
  {
    // get block ale matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> a =
        AleField()->BlockSystemMatrix();
    if (a == Teuchos::null)
    {
      dserror("expect ale block matrix");
    }

    rhs = Teuchos::rcp(new Epetra_Vector(a->Matrix(0, 1).RowMap()));

    a->Matrix(0, 1).Apply(*FluidToAleInterface(iprojdispinc_), *rhs);

    Extractor().AddVector(*rhs, 2, f);

    // get fluid shape derivative matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm =
        FluidField()->ShapeDerivatives();
    if (mmm != Teuchos::null)
    {
      // extract submatrices
      LINALG::SparseMatrix& fmig = mmm->Matrix(0, 1);
      LINALG::SparseMatrix& fmgg = mmm->Matrix(1, 1);

      rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RowMap()));

      fmgg.Apply(*iprojdispinc_, *rhs);

      Teuchos::RCP<Epetra_Vector> tmprhs = Teuchos::rcp(
          new Epetra_Vector(mortarp->DomainMap()));
      mortarp->Multiply(true, *rhs, *tmprhs);

      rhs = StructureField()->Interface()->InsertFSICondVector(tmprhs);

      Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(
          new const Epetra_Vector(rhs->Map(), true));
      LINALG::ApplyDirichlettoSystem(rhs, zeros,
          *(StructureField()->GetDBCMapExtractor()->CondMap()));

      if (StructureField()->GetSTCAlgo() == INPAR::STR::stc_currsym)
      {
        Teuchos::RCP<LINALG::SparseMatrix> stcmat =
            StructureField()->GetSTCMat();
        stcmat->Multiply(true, *rhs, *rhs);
      }

      rhs->Scale(Dt() * timescale * (1. - stiparam) / (1. - ftiparam));
      Extractor().AddVector(*rhs, 0, f);

      rhs = Teuchos::rcp(new Epetra_Vector(fmig.RowMap()));

      fmig.Apply(*iprojdispinc_, *rhs);

      rhs = FsiFluidField()->FsiInterface()->InsertOtherVector(rhs);

      zeros = Teuchos::rcp(new const Epetra_Vector(rhs->Map(), true));
      LINALG::ApplyDirichlettoSystem(rhs, zeros,
          *(StructureField()->GetDBCMapExtractor()->CondMap()));

      rhs->Scale(-timescale * Dt());

      Extractor().AddVector(*rhs, 1, f);

    }
  }

  // Reset quantities of previous iteration step since they still store values from the last time step
  ddginc_ = LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(),
      true);
  duiinc_ = LINALG::CreateVector(*FsiFluidField()->FsiInterface()->OtherMap(), true);
  veliprev_ = Teuchos::null;
  velgprev_ = Teuchos::null;
  fgicur_ = Teuchos::null;
  fggcur_ = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::SetupSystemMatrix(
    LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::SetupSystemMatrix");

  // get the Mortar projection matrix P = D^{-1} * M
  const Teuchos::RCP<LINALG::SparseMatrix> mortarp = coupsfm_->GetMortarTrafo();

  // get info about STC feature
  INPAR::STR::STC_Scale stcalgo = StructureField()->GetSTCAlgo();
  Teuchos::RCP<LINALG::SparseMatrix> stcmat = Teuchos::null;
  if (stcalgo != INPAR::STR::stc_none)
    stcmat = StructureField()->GetSTCMat();

  const ADAPTER::Coupling& coupfa = FluidAleCoupling();

  // get single field block matrices
  Teuchos::RCP<LINALG::SparseMatrix> s = StructureField()->SystemMatrix(); // can't be 'const' --> is modified by STC
  const Teuchos::RCP<LINALG::BlockSparseMatrixBase> f =
      FluidField()->BlockSystemMatrix();
  const Teuchos::RCP<LINALG::BlockSparseMatrixBase> a =
      AleField()->BlockSystemMatrix();

#ifdef DEBUG
  // check whether allocation was successful
  if (mortarp == Teuchos::null)
  dserror("Expected Teuchos::rcp to mortar matrix P.");
  if (s == Teuchos::null)
  { dserror("expect structure block matrix");}
  if (f == Teuchos::null)
  { dserror("expect fluid block matrix");}
  if (a == Teuchos::null)
  { dserror("expect ale block matrix");}

  // some checks whether maps for matrix-matrix-multiplication do really match
  if (!f->Matrix(0,1).DomainMap().PointSameAs(mortarp->RangeMap()))
  dserror("Maps do not match.");
  if (!f->Matrix(1,0).RangeMap(). PointSameAs(mortarp->RangeMap()))
  dserror("Maps do not match.");
  if (!f->Matrix(1,1).DomainMap().PointSameAs(mortarp->RangeMap()))
  dserror("Maps do not match.");
#endif

  // extract submatrices
  LINALG::SparseMatrix& aii = a->Matrix(0, 0); // A_{II}
  LINALG::SparseMatrix& aig = a->Matrix(0, 1); // A_{I\Gamma}
  LINALG::SparseMatrix& fii = f->Matrix(0, 0); // F_{II}

  // scaling factors for fluid
  const double scale = FluidField()->ResidualScaling();
  const double timescale = FluidField()->TimeScaling();

  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();
  const double ftiparam = FluidField()->TimIntParam();

  // uncomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  // this just once...
  s->UnComplete();

  // ---------------------------------------------------------------------------
  // BEGIN building the global 4x4 system matrix
  // ---------------------------------------------------------------------------
  // Contributions to blocks in system matrix are listed separately.
  // Block numbering in comments ranges from (1,1) to (4,4).

  // ---------Addressing contribution to block (2,2)
  Teuchos::RCP<LINALG::SparseMatrix> fgg = MLMultiply(f->Matrix(1, 1), false,
      *mortarp, false, false, false, true);
  fgg = MLMultiply(*mortarp, true, *fgg, false, false, false, true);

  s->Add(*fgg, false, scale * timescale * (1. - stiparam) / (1. - ftiparam),
      1.0);

  // ---------Addressing contribution to block (2,3)
  Teuchos::RCP<LINALG::SparseMatrix> fgi = MLMultiply(*mortarp, true,
      f->Matrix(1, 0), false, false, false, true);
  Teuchos::RCP<LINALG::SparseMatrix> lfgi = Teuchos::rcp(
      new LINALG::SparseMatrix(s->RowMap(), 81, false));

  lfgi->Add(*fgi, false, scale, 0.0);
  lfgi->Complete(fgi->DomainMap(), s->RangeMap());

  if (stcalgo == INPAR::STR::stc_currsym)
    lfgi = LINALG::MLMultiply(*stcmat, true, *lfgi, false, true, true, true);

  mat.Matrix(0, 1).UnComplete();
  mat.Matrix(0, 1).Add(*lfgi, false, (1. - stiparam) / (1. - ftiparam), 0.0);

  // ---------Addressing contribution to block (3,2)
  Teuchos::RCP<LINALG::SparseMatrix> fig = MLMultiply(f->Matrix(0, 1), false,
      *mortarp, false, false, false, true);
  Teuchos::RCP<LINALG::SparseMatrix> lfig = Teuchos::rcp(
      new LINALG::SparseMatrix(fig->RowMap(), 81, false));

  lfig->Add(*fig, false, timescale, 0.0);
  lfig->Complete(s->DomainMap(), fig->RangeMap());

  if (stcalgo != INPAR::STR::stc_none)
  {
    lfig = LINALG::MLMultiply(*lfig, false, *stcmat, false, false, false, true);
  }

  mat.Matrix(1, 0).UnComplete();
  mat.Matrix(1, 0).Add(*lfig, false, 1., 0.0);

  // ---------Addressing contribution to block (3,3)
  mat.Matrix(1, 1).UnComplete();
  mat.Matrix(1, 1).Add(fii, false, 1., 0.0);
  Teuchos::RCP<LINALG::SparseMatrix> eye = LINALG::Eye(
      *FluidField()->Interface()->FSICondMap());
  mat.Matrix(1, 1).Add(*eye, false, 1., 1.0);

  // ---------Addressing contribution to block (4,2)
  Teuchos::RCP<LINALG::SparseMatrix> laig = Teuchos::rcp(
      new LINALG::SparseMatrix(aii.RowMap(), 81, false));
  (*aigtransform_)(a->FullRowMap(), a->FullColMap(), aig, 1.,
      ADAPTER::CouplingSlaveConverter(InterfaceFluidAleCoupling()), *laig);

  laig->Complete(f->Matrix(1, 1).DomainMap(), aii.RangeMap());
  Teuchos::RCP<LINALG::SparseMatrix> llaig = MLMultiply(*laig, false, *mortarp,
      false, false, false, true);
  laig = Teuchos::rcp(new LINALG::SparseMatrix(llaig->RowMap(), 81, false));

  laig->Add(*llaig, false, 1.0, 0.0);
  laig->Complete(s->DomainMap(), llaig->RangeMap());

  if (stcalgo != INPAR::STR::stc_none)
  {
    laig = LINALG::MLMultiply(*laig, false, *stcmat, false, false, false, true);
  }

  mat.Assign(2,0,LINALG::View,*laig);

  // ---------Addressing contribution to block (4,4)
  mat.Assign(2,2,LINALG::View,aii);

  /*--------------------------------------------------------------------------*/
  // add optional fluid linearization with respect to mesh motion block
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm =
      FluidField()->ShapeDerivatives();
  if (mmm != Teuchos::null)
  {
    // extract submatrices
    LINALG::SparseMatrix& fmii = mmm->Matrix(0, 0);
    LINALG::SparseMatrix& fmgi = mmm->Matrix(1, 0);

    // reuse transform objects to add shape derivative matrices to structural blocks

    // ---------Addressing contribution to block (2,2)
    Teuchos::RCP<LINALG::SparseMatrix> fmgg = MLMultiply(mmm->Matrix(1, 1),
        false, *mortarp, false, false, false, true);
    fmgg = MLMultiply(*mortarp, true, *fmgg, false, false, false, true);

    Teuchos::RCP<LINALG::SparseMatrix> lfmgg = Teuchos::rcp(
        new LINALG::SparseMatrix(fmgg->RowMap(), 81, false));
    lfmgg->Add(*fmgg, false, 1.0, 0.0);
    lfmgg->Complete(s->DomainMap(), fmgg->RangeMap());

    s->Add(*lfmgg, false, scale * (1. - stiparam) / (1. - ftiparam), 1.0);

    // ---------Addressing contribution to block (3,2)
    Teuchos::RCP<LINALG::SparseMatrix> fmig = MLMultiply(mmm->Matrix(0, 1),
        false, *mortarp, false, false, false, true);
    Teuchos::RCP<LINALG::SparseMatrix> lfmig = Teuchos::rcp(
        new LINALG::SparseMatrix(fmig->RowMap(), 81, false));

    lfmig->Add(*fmig, false, 1.0, 0.0);
    lfmig->Complete(s->DomainMap(), fmig->RangeMap());

    if (stcalgo != INPAR::STR::stc_none)
    {
      lfmig = LINALG::MLMultiply(*lfmig, false, *stcmat, false, false, false,
          true);
    }

    mat.Matrix(1, 0).Add(*lfmig, false, 1.0, 1.0);

    // We cannot copy the pressure value. It is not used anyway. So no exact
    // match here.
    (*fmiitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmii, 1.,
        ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false);

    Teuchos::RCP<LINALG::SparseMatrix> lfmgi = Teuchos::rcp(
        new LINALG::SparseMatrix(fmgi.RowMap(), 81, false));
    (*fmiitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmgi, 1.0,
        ADAPTER::CouplingMasterConverter(coupfa), *lfmgi, false);

    // ---------Addressing contribution to block (2,4)
    lfmgi->Complete(aii.DomainMap(), mortarp->RangeMap());
    Teuchos::RCP<LINALG::SparseMatrix> llfmgi = MLMultiply(*mortarp, true,
        *lfmgi, false, false, false, true);
    lfmgi = Teuchos::rcp(new LINALG::SparseMatrix(s->RowMap(), 81, false));

    lfmgi->Add(*llfmgi, false, scale, 0.0);
    lfmgi->Complete(aii.DomainMap(), s->RangeMap());

    if (stcalgo == INPAR::STR::stc_currsym)
      lfmgi = LINALG::MLMultiply(*stcmat, true, *lfmgi, false, true, true,
          false);
    lfmgi->Scale((1. - stiparam) / (1. - ftiparam));
    mat.Assign(0, 2, LINALG::View, *lfmgi);
  }

  s->Complete();

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
  mat.Assign(0,0,LINALG::View,*s);

  // done. make sure all blocks are filled.
  mat.Complete();

  // Finally, take care of Dirichlet boundary conditions
  mat.ApplyDirichlet(*(dbcmaps_->CondMap()), true);
  //
  // ---------------------------------------------------------------------------
  // END building the global system matrix
  // ---------------------------------------------------------------------------

//  Teuchos::RCP<Epetra_CrsMatrix> matrix = mat.Matrix(0,0).EpetraMatrix();
//  LINALG::PrintMatrixInMatlabFormat("mat.dat",*matrix,true);

//  LINALG::PrintBlockMatrixInMatlabFormat("mat.dat",mat);
//  std::cout<<"\nWROTE MATRIX!!";

  // ---------------------------------------------------------------------------
  // NOX related stuff needed for recovery of Lagrange multiplier
  // ---------------------------------------------------------------------------
  // store parts of fluid matrix to know them in the next iteration as previous
  // iteration matrices
  fgiprev_ = fgicur_;
  fggprev_ = fggcur_;
  fgicur_ = Teuchos::rcp(new LINALG::SparseMatrix(f->Matrix(1, 0)));
  fggcur_ = Teuchos::rcp(new LINALG::SparseMatrix(f->Matrix(1, 1)));

  // store parts of fluid shape derivative matrix to know them in the next
  // iteration as previous iteration matrices
  fmgiprev_ = fmgicur_;
  fmggprev_ = fmggcur_;
  if (mmm != Teuchos::null)
  {
    fmgicur_ = Teuchos::rcp(new LINALG::SparseMatrix(mmm->Matrix(1, 0)));
    fmggcur_ = Teuchos::rcp(new LINALG::SparseMatrix(mmm->Matrix(1, 1)));
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::ScaleSystem(
    LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn =
      DRT::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool) DRT::INPUT::IntegralValue<int>(fsimono,
      "INFNORMSCALING");

  if (scaling_infnorm)
  {
    // The matrices are modified here. Do we have to change them back later on?

    // do scaling of structure rows
    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0, 0).EpetraMatrix();
    srowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    scolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*srowsum_);
    A->InvColSums(*scolsum_);
    if (A->LeftScale(*srowsum_) or A->RightScale(*scolsum_)
        or mat.Matrix(0, 1).EpetraMatrix()->LeftScale(*srowsum_)
        or mat.Matrix(0, 2).EpetraMatrix()->LeftScale(*srowsum_)
        or mat.Matrix(1, 0).EpetraMatrix()->RightScale(*scolsum_)
        or mat.Matrix(2, 0).EpetraMatrix()->RightScale(*scolsum_))
      dserror("structure scaling failed");

    // do scaling of ale rows
    A = mat.Matrix(2, 2).EpetraMatrix();
    arowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    acolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*arowsum_);
    A->InvColSums(*acolsum_);
    if (A->LeftScale(*arowsum_) or A->RightScale(*acolsum_)
        or mat.Matrix(2, 0).EpetraMatrix()->LeftScale(*arowsum_)
        or mat.Matrix(2, 1).EpetraMatrix()->LeftScale(*arowsum_)
        or mat.Matrix(0, 2).EpetraMatrix()->RightScale(*acolsum_)
        or mat.Matrix(1, 2).EpetraMatrix()->RightScale(*acolsum_))
      dserror("ale scaling failed");

    // do scaling of structure and ale rhs vectors
    Teuchos::RCP<Epetra_Vector> sx = Extractor().ExtractVector(b, 0);
    Teuchos::RCP<Epetra_Vector> ax = Extractor().ExtractVector(b, 2);

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
void FSI::SlidingMonolithicFluidSplit::UnscaleSolution(
    LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn =
      DRT::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool) DRT::INPUT::IntegralValue<int>(fsimono,
      "INFNORMSCALING");

  if (scaling_infnorm)
  {
    Teuchos::RCP<Epetra_Vector> sy = Extractor().ExtractVector(x, 0);
    Teuchos::RCP<Epetra_Vector> ay = Extractor().ExtractVector(x, 2);

    if (sy->Multiply(1.0, *scolsum_, *sy, 0.0))
      dserror("structure scaling failed");
    if (ay->Multiply(1.0, *acolsum_, *ay, 0.0))
      dserror("ale scaling failed");

    // get info about STC feature and unscale solution if necessary
    INPAR::STR::STC_Scale stcalgo = StructureField()->GetSTCAlgo();
    if (stcalgo != INPAR::STR::stc_none)
    {
      StructureField()->GetSTCMat()->Multiply(false, *sy, *sy);
    }

    Extractor().InsertVector(*sy, 0, x);
    Extractor().InsertVector(*ay, 2, x);

    Teuchos::RCP<Epetra_Vector> sx = Extractor().ExtractVector(b, 0);
    Teuchos::RCP<Epetra_Vector> ax = Extractor().ExtractVector(b, 2);

    if (sx->ReciprocalMultiply(1.0, *srowsum_, *sx, 0.0))
      dserror("structure scaling failed");
    if (ax->ReciprocalMultiply(1.0, *arowsum_, *ax, 0.0))
      dserror("ale scaling failed");

    // get info about STC feature
    if (stcalgo != INPAR::STR::stc_none)
    {
      StructureField()->GetSTCMat()->Multiply(false, *sx, *sx);
    }

    Extractor().InsertVector(*sx, 0, b);
    Extractor().InsertVector(*ax, 2, b);

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0, 0).EpetraMatrix();
    srowsum_->Reciprocal(*srowsum_);
    scolsum_->Reciprocal(*scolsum_);
    if (A->LeftScale(*srowsum_) or A->RightScale(*scolsum_)
        or mat.Matrix(0, 1).EpetraMatrix()->LeftScale(*srowsum_)
        or mat.Matrix(0, 2).EpetraMatrix()->LeftScale(*srowsum_)
        or mat.Matrix(1, 0).EpetraMatrix()->RightScale(*scolsum_)
        or mat.Matrix(2, 0).EpetraMatrix()->RightScale(*scolsum_))
      dserror("structure scaling failed");

    A = mat.Matrix(2, 2).EpetraMatrix();
    arowsum_->Reciprocal(*arowsum_);
    acolsum_->Reciprocal(*acolsum_);
    if (A->LeftScale(*arowsum_) or A->RightScale(*acolsum_)
        or mat.Matrix(2, 0).EpetraMatrix()->LeftScale(*arowsum_)
        or mat.Matrix(2, 1).EpetraMatrix()->LeftScale(*arowsum_)
        or mat.Matrix(0, 2).EpetraMatrix()->RightScale(*acolsum_)
        or mat.Matrix(1, 2).EpetraMatrix()->RightScale(*acolsum_))
      dserror("ale scaling failed");
  }

  // very simple hack just to see the linear solution

  Epetra_Vector r(b.Map());
  mat.Apply(x, r);
  r.Update(1., b, 1.);

  Teuchos::RCP<Epetra_Vector> sr = Extractor().ExtractVector(r, 0);
  Teuchos::RCP<Epetra_Vector> fr = Extractor().ExtractVector(r, 1);
  Teuchos::RCP<Epetra_Vector> ar = Extractor().ExtractVector(r, 2);

  // increment additional ale residual
  aleresidual_->Update(-1., *ar, 0.);

  std::ios_base::fmtflags flags = Utils()->out().flags();

  double n, ns, nf, na;
  r.Norm2(&n);
  sr->Norm2(&ns);
  fr->Norm2(&nf);
  ar->Norm2(&na);
  Utils()->out() << std::scientific << "\nlinear solver quality:\n"
      << "L_2-norms:\n" << "   |r|=" << n << "   |rs|=" << ns << "   |rf|="
      << nf << "   |ra|=" << na << "\n";
  r.NormInf(&n);
  sr->NormInf(&ns);
  fr->NormInf(&nf);
  ar->NormInf(&na);
  Utils()->out() << "L_inf-norms:\n" << "   |r|=" << n << "   |rs|=" << ns
      << "   |rf|=" << nf << "   |ra|=" << na << "\n";

  Utils()->out().flags(flags);

  if (StructureField()->GetSTCAlgo() != INPAR::STR::stc_none)
    StructureField()->SystemMatrix()->Reset();

}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem> FSI::SlidingMonolithicFluidSplit::CreateLinearSystem(
    Teuchos::ParameterList& nlParams, NOX::Epetra::Vector& noxSoln,
    Teuchos::RCP<NOX::Utils> utils)
{
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys;

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  NOX::Epetra::Interface::Jacobian* iJac = this;
  NOX::Epetra::Interface::Preconditioner* iPrec = this;
  const Teuchos::RCP<Epetra_Operator> J = systemmatrix_;
  const Teuchos::RCP<Epetra_Operator> M = systemmatrix_;

  switch (linearsolverstrategy_)
  {
  case INPAR::FSI::PreconditionedKrylov:
  case INPAR::FSI::FSIAMG:
  case INPAR::FSI::HybridSchwarz: // ToDo (noll) Do we need a separate case for HybridSchwarz?
    linSys = Teuchos::rcp(
        new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
            Teuchos::rcp(iJac, false), J, Teuchos::rcp(iPrec, false), M,
            noxSoln));

    break;
  default:
    dserror("unsupported linear block solver strategy: %d",
        linearsolverstrategy_);
    break;
  }

  return linSys;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo> FSI::SlidingMonolithicFluidSplit::CreateStatusTest(
    Teuchos::ParameterList& nlParams, Teuchos::RCP<NOX::Epetra::Group> grp)
{
  // ---------------------------------------------------------------------------
  // Setup the test framework
  // ---------------------------------------------------------------------------
  // Create the top-level test combo
  Teuchos::RCP<NOX::StatusTest::Combo> combo = Teuchos::rcp(
      new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  // Create test combo for convergence of residuals and iterative increments
  Teuchos::RCP<NOX::StatusTest::Combo> converged = Teuchos::rcp(
      new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  // Create some other plausibility tests
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(
      new NOX::StatusTest::MaxIters(nlParams.get<int>("Max Iterations")));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv = Teuchos::rcp(
      new NOX::StatusTest::FiniteValue);

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
  Teuchos::RCP<NOX::StatusTest::Combo> structcombo = Teuchos::rcp(
      new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> structureDisp_L2 = Teuchos::rcp(
      new NOX::FSI::PartialNormF("DISPL residual", Extractor(), 0,
          nlParams.get<double>("Tol dis res L2"),
          NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> structureDisp_inf = Teuchos::rcp(
      new NOX::FSI::PartialNormF("DISPL residual", Extractor(), 0,
          nlParams.get<double>("Tol dis res Inf"),
          NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate_L2 =
      Teuchos::rcp(
          new NOX::FSI::PartialNormUpdate("DISPL update", Extractor(), 0,
              nlParams.get<double>("Tol dis inc L2"),
              NOX::Abstract::Vector::TwoNorm,
              NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate_inf =
      Teuchos::rcp(
          new NOX::FSI::PartialNormUpdate("DISPL update", Extractor(), 0,
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
  interface.push_back(StructureField()->Interface()->FSICondMap());
  interface.push_back(Teuchos::null);
  LINALG::MultiMapExtractor interfaceextract(*DofRowMap(), interface);

  // create NOX::StatusTest::Combo for interface
  Teuchos::RCP<NOX::StatusTest::Combo> interfacecombo = Teuchos::rcp(
      new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> interfaceTest_L2 = Teuchos::rcp(
      new NOX::FSI::PartialNormF("GAMMA residual", interfaceextract, 0,
          nlParams.get<double>("Tol fsi res L2"),
          NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> interfaceTest_inf = Teuchos::rcp(
      new NOX::FSI::PartialNormF("GAMMA residual", interfaceextract, 0,
          nlParams.get<double>("Tol fsi res Inf"),
          NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> interfaceTestUpdate_L2 =
      Teuchos::rcp(
          new NOX::FSI::PartialNormUpdate("GAMMA update", interfaceextract, 0,
              nlParams.get<double>("Tol fsi inc L2"),
              NOX::Abstract::Vector::TwoNorm,
              NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> interfaceTestUpdate_inf =
      Teuchos::rcp(
          new NOX::FSI::PartialNormUpdate("GAMMA update", interfaceextract, 0,
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
  LINALG::MultiMapExtractor fluidvelextract(*DofRowMap(), fluidvel);

  // create NOX::StatusTest::Combo for fluid velocity field
  Teuchos::RCP<NOX::StatusTest::Combo> fluidvelcombo = Teuchos::rcp(
      new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> innerFluidVel_L2 = Teuchos::rcp(
      new NOX::FSI::PartialNormF("VELOC residual", fluidvelextract, 0,
          nlParams.get<double>("Tol vel res L2"),
          NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> innerFluidVel_inf = Teuchos::rcp(
      new NOX::FSI::PartialNormF("VELOC residual", fluidvelextract, 0,
          nlParams.get<double>("Tol vel res Inf"),
          NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> innerFluidVelUpdate_L2 =
      Teuchos::rcp(
          new NOX::FSI::PartialNormUpdate("VELOC update", fluidvelextract, 0,
              nlParams.get<double>("Tol vel inc L2"),
              NOX::Abstract::Vector::TwoNorm,
              NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> innerFluidVelUpdate_inf =
      Teuchos::rcp(
          new NOX::FSI::PartialNormUpdate("VELOC update", fluidvelextract, 0,
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
  LINALG::MultiMapExtractor fluidpressextract(*DofRowMap(), fluidpress);

  // create NOX::StatusTest::Combo for fluid pressure field
  Teuchos::RCP<NOX::StatusTest::Combo> fluidpresscombo = Teuchos::rcp(
      new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> fluidPress_L2 = Teuchos::rcp(
      new NOX::FSI::PartialNormF("PRESS residual", fluidpressextract, 0,
          nlParams.get<double>("Tol pre res L2"),
          NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> fluidPress_inf = Teuchos::rcp(
      new NOX::FSI::PartialNormF("PRESS residual", fluidpressextract, 0,
          nlParams.get<double>("Tol pre res Inf"),
          NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> fluidPressUpdate_L2 = Teuchos::rcp(
      new NOX::FSI::PartialNormUpdate("PRESS update", fluidpressextract, 0,
          nlParams.get<double>("Tol pre inc L2"),
          NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> fluidPressUpdate_inf = Teuchos::rcp(
      new NOX::FSI::PartialNormUpdate("PRESS update", fluidpressextract, 0,
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
void FSI::SlidingMonolithicFluidSplit::ExtractFieldVectors(
    Teuchos::RCP<const Epetra_Vector> x, Teuchos::RCP<const Epetra_Vector>& sx,
    Teuchos::RCP<const Epetra_Vector>& fx,
    Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::ExtractFieldVectors");

#ifdef DEBUG
  if(ddgpred_ == Teuchos::null)
  dserror("Vector 'ddgpred_' has not been initialized properly.");
#endif

  // get the Mortar projection matrix P = D^{-1} * M
  const Teuchos::RCP<const LINALG::SparseMatrix> mortarp =
      coupsfm_->GetMortarTrafo();

  // ---------------------------------------------------------------------------
  // process structure unknowns
  // ---------------------------------------------------------------------------
  // extract structure solution increment from NOX increment
  sx = Extractor().ExtractVector(x, 0);

  // ---------------------------------------------------------------------------
  // process ale unknowns
  // ---------------------------------------------------------------------------
  // extract inner ALE solution increment from NOX increment
  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x, 2);

  // convert structure solution increment to ALE solution increment at the interface
  Teuchos::RCP<Epetra_Vector> scx =
      StructureField()->Interface()->ExtractFSICondVector(sx);
  scx->Update(1.0, *ddgpred_, 1.0);
  Teuchos::RCP<Epetra_Vector> acx = LINALG::CreateVector(
      *FluidField()->Interface()->FSICondMap());
  mortarp->Apply(*scx, *acx);
  acx = FluidToAleInterface(acx);

  // put inner and interface ALE solution increments together
  Teuchos::RCP<Epetra_Vector> a = FsiAleField()->FsiInterface()->InsertOtherVector(
      aox);
  AleField()->Interface()->InsertFSICondVector(acx, a);
  AleField()->UpdateSlaveDOF(a);
  ax = a;

  // ---------------------------------------------------------------------------
  // process fluid unknowns
  // ---------------------------------------------------------------------------
  // extract inner fluid solution increment from NOX increment
  Teuchos::RCP<const Epetra_Vector> fox = Extractor().ExtractVector(x, 1);
  fox = FsiFluidField()->FsiInterface()->ExtractOtherVector(fox);

  // convert ALE solution increment to fluid solution increment at the interface
  Teuchos::RCP<Epetra_Vector> fcx = AleToFluidInterface(acx);
  FluidField()->DisplacementToVelocity(fcx);

  // put inner and interface fluid solution increments together
  Teuchos::RCP<Epetra_Vector> f = FsiFluidField()->FsiInterface()->InsertOtherVector(
      fox);
  FluidField()->Interface()->InsertFSICondVector(fcx, f);
  FluidField()->UpdateSlaveDOF(f);
  fx = f;

  // ---------------------------------------------------------------------------

  // Store field vectors to know them later on as previous quantities
  // interface structure displacement increment
  if (disgprev_ != Teuchos::null)
    ddginc_->Update(1.0, *scx, -1.0, *disgprev_, 0.0); // compute current iteration increment
  else
    ddginc_ = Teuchos::rcp(new Epetra_Vector(*scx)); // first iteration increment

  disgprev_ = scx; // store current step increment
  // ------------------------------------

  // inner ale displacement increment
  if (aleiprev_ != Teuchos::null)
    ddialeinc_->Update(1.0, *aox, -1.0, *aleiprev_, 0.0); // compute current iteration increment
  else
    ddialeinc_ = Teuchos::rcp(new Epetra_Vector(*aox)); // first iteration increment

  aleiprev_ = aox; // store current step increment
  // ------------------------------------

  // inner fluid solution increment
  if (veliprev_ != Teuchos::null) // compute current iteration increment
    duiinc_->Update(1.0, *fox, -1.0, *veliprev_, 0.0);
  else
    // first iteration increment
    duiinc_ = Teuchos::rcp(new Epetra_Vector(*fox));
  // store current step increment
  veliprev_ = fox;
  // ------------------------------------
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::Update()
{

  // update history variabels for sliding ale
  if (aleproj_ != INPAR::FSI::ALEprojection_none)
  {
    iprojdisp_ = Teuchos::rcp(
        new Epetra_Vector(*coupsfm_->SlaveDofMap(), true));
    Teuchos::RCP<Epetra_Vector> idispale = AleToFluidInterface(
        AleField()->Interface()->ExtractFSICondVector(AleField()->Dispnp()));

    slideale_->Remeshing(*StructureField(), FluidField()->Discretization(),
        idispale, iprojdisp_, *coupsfm_, Comm());

    iprojdispinc_->Update(-1.0, *iprojdisp_, 1.0, *idispale, 0.0);

    slideale_->EvaluateMortar(StructureField()->ExtractInterfaceDispnp(),
        iprojdisp_, *coupsfm_);
    slideale_->EvaluateFluidMortar(idispale, iprojdisp_);

    Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(
        new Epetra_Vector(*iprojdisp_));
    temp->ReplaceMap(idispale->Map());
    Teuchos::RCP<Epetra_Vector> acx = FluidToAleInterface(temp);
    AleField()->ApplyInterfaceDisplacements(acx);
    FluidField()->ApplyMeshDisplacement(AleToFluid(AleField()->Dispnp()));

    Teuchos::RCP<Epetra_Vector> unew = slideale_->InterpolateFluid(
        FluidField()->ExtractInterfaceVelnp());
    FluidField()->ApplyInterfaceVelocities(unew);
  }

  // call Update()-routine in base class to handle the single fields
  FSI::MonolithicBase::Update();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::Output()
{
  StructureField()->Output();
  FluidField()->Output();

  if (aleproj_ != INPAR::FSI::ALEprojection_none)
  {
    int uprestart = timeparams_.get<int>("RESTARTEVRY");
    if (uprestart != 0 && FluidField()->Step() % uprestart == 0)
    {
      FluidField()->DiscWriter()->WriteVector("slideALE", iprojdisp_);
      FluidField()->DiscWriter()->WriteVector("slideALEincr", iprojdispinc_);
      slideale_->OutputRestart(*FluidField()->DiscWriter());
    }
  }

  // output Lagrange multiplier
  OutputLambda();

  AleField()->Output();

  if (StructureField()->GetConstraintManager()->HaveMonitor())
  {
    StructureField()->GetConstraintManager()->ComputeMonitorValues(
        StructureField()->Dispnp());
    if (comm_.MyPID() == 0)
      StructureField()->GetConstraintManager()->PrintMonitorValues();
  }

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::OutputLambda()
{
  /* 'lambda_' is only defined on the interface. So, insert 'lambda_' into
   * 'lambdafull' that is defined on the entire fluid field. Then, write
   * output or restart data.
   */
  Teuchos::RCP<Epetra_Vector> lambdafull =
      FluidField()->Interface()->InsertFSICondVector(lambda_);
  const int uprestart = timeparams_.get<int>("RESTARTEVRY");
  const int upres = timeparams_.get<int>("RESULTSEVRY");
  if ((uprestart != 0 && FluidField()->Step() % uprestart == 0)
      || FluidField()->Step() % upres == 0)
    FluidField()->DiscWriter()->WriteVector("fsilambda", lambdafull);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::ReadRestart(int step)
{
  StructureField()->ReadRestart(step);
  FluidField()->ReadRestart(step);

  // read Lagrange multiplier
  {
    Teuchos::RCP<Epetra_Vector> lambdafull = Teuchos::rcp(
        new Epetra_Vector(*FluidField()->DofRowMap(), true));
    IO::DiscretizationReader reader = IO::DiscretizationReader(
        FluidField()->Discretization(), step);
    reader.ReadVector(lambdafull, "fsilambda");
    lambda_ = FluidField()->Interface()->ExtractFSICondVector(lambdafull);
  }

  SetupSystem();

  if (aleproj_ != INPAR::FSI::ALEprojection_none)
  {
    IO::DiscretizationReader reader = IO::DiscretizationReader(
        FluidField()->Discretization(), step);
    reader.ReadVector(iprojdisp_, "slideALE");
    reader.ReadVector(iprojdispinc_, "slideALEincr");
    slideale_->ReadRestart(reader);
  }
  AleField()->ReadRestart(step);

  SetTimeStep(FluidField()->Time(), FluidField()->Step());

  if (aleproj_ != INPAR::FSI::ALEprojection_none)
    slideale_->EvaluateMortar(StructureField()->ExtractInterfaceDispn(),
        iprojdisp_, *coupsfm_);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::PrepareTimeStep()
{

  precondreusecount_ = 0;

  IncrementTimeAndStep();
  PrintHeader();

  PrepareTimeStepPreconditioner();

  if (StructureField()->GetSTCAlgo() != INPAR::STR::stc_none)
    StructureField()->SystemMatrix()->Reset();

  PrepareTimeStepFields();

  //Note: it's important to first prepare the single fields and than the fsi problem
  PrepareTimeStepFSI();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::RecoverLagrangeMultiplier()
{
  // store previous Lagrange multiplier for calculation of interface energy
  lambdaold_->Update(1.0, *lambda_, 0.0);

  // get time integration parameter of fluid time integrator
  // to enable consistent time integration among the fields
  const double ftiparam = FluidField()->TimIntParam();

  // some scaling factors for fluid
  const double timescale = FluidField()->TimeScaling();
  const double scale = FluidField()->ResidualScaling();

  // get the Mortar projection matrix P = D^{-1} * M
  const Teuchos::RCP<LINALG::SparseMatrix> mortarp = coupsfm_->GetMortarTrafo();

  // get the inverted Mortar matrix D^{-1}
  const Teuchos::RCP<LINALG::SparseMatrix> mortardinv =
      coupsfm_->GetDinvMatrix();

#ifdef DEBUG
  if (mortarp == Teuchos::null)
  dserror("Expected Teuchos::rcp to mortar matrix P.");
  if (mortardinv == Teuchos::null)
  dserror("Expected Teuchos::rcp to mortar matrix D^{-1}.");
#endif

  // get fluid shape derivative matrix
  const Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm =
      FluidField()->ShapeDerivatives();

  // some often re-used vectors
  Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::null; // stores intermediate result of terms (3)-(8)
  Teuchos::RCP<Epetra_Vector> auxvec = Teuchos::null; // just for convenience
  Teuchos::RCP<Epetra_Vector> auxauxvec = Teuchos::null; // just for convenience

  /* Recovery of Lagrange multiplier \lambda_^{n+1} is done by the following
   * condensation expression:
   *
   * lambda_^{n+1} =
   *
   * (1)  - ftiparam / (1.-ftiparam) * lambda^{n}
   *
   * (2)  - 1. / (1.-ftiparam) * D^{-T} * tmpvec
   *
   * with tmpvec =
   *
   * (3)    r_{\Gamma}^{F,n+1}
   *
   * (4)  + 1 / tau * F_{\Gamma\Gamma} * P * \Delta d_{\Gamma}^{S,n+1}
   *
   * (5)  + F_{\Gamma\Gamma}^{G} * P * \Delta d_{\Gamma}^{S,n+1}
   *
   * (6)  + F_{\Gamma I} * \Delta u_{I}^{F,n+1}
   *
   * (7)  + F_{\Gamma I}^{G} * \Delta d_{I}^{G,n+1}
   *
   * (8)  - dt / tau * F_{\Gamma\Gamma} * u_{\Gamma}^n]
   *
   * Remark on term (8):
   * Term (8) has to be considered only in the first Newton iteration.
   * Hence, it will usually not be computed since in general we need more
   * than one nonlinear iteration until convergence.
   *
   * Remarks on all terms:
   * +  Division by (1.0 - ftiparam) will be done in the end
   *    since this is common to all terms
   * +  tau: time scaling factor for interface time integration (tau = 1/FluidField()->TimeScaling())
   * +  neglecting terms (4)-(8) should not alter the results significantly
   *    since at the end of the time step the solution increments tend to zero.
   *
   * Terms arising from field specific predictors have to be considered only in
   * the first Newton iteration. Since we usually need more than on iteration,
   * these terms are not implemented, yet.
   *
   *                                                 Matthias Mayr (10/2012)
   */

  // ---------Addressing term (1)
  lambda_->Scale(ftiparam);
  // ---------End of term (1)

  // ---------Addressing term (3)
  Teuchos::RCP<Epetra_Vector> fluidresidual =
      FluidField()->Interface()->ExtractFSICondVector(FluidField()->RHS());
  fluidresidual->Scale(-1.0);
  tmpvec = Teuchos::rcp(new Epetra_Vector(*fluidresidual));
  // ---------End of term (3)

  // ---------Addressing term (4)
  auxvec = Teuchos::rcp(new Epetra_Vector(mortarp->RangeMap(), true));
  mortarp->Apply(*ddginc_, *auxvec);
  auxauxvec = Teuchos::rcp(new Epetra_Vector(fggprev_->RangeMap(), true));
  fggprev_->Apply(*auxvec, *auxauxvec);
  tmpvec->Update(timescale, *auxauxvec, 1.0);
  // ---------End of term (4)

  // ---------Addressing term (5)
  if (fmggprev_ != Teuchos::null)
  {
    auxvec = Teuchos::rcp(new Epetra_Vector(mortarp->RangeMap(), true));
    mortarp->Apply(*ddginc_, *auxvec);
    fmggprev_->Apply(*auxvec, *auxauxvec);
    tmpvec->Update(1.0, *auxauxvec, 1.0);
  }
  // ---------End of term (5)

  // ---------Addressing term (6)
  auxvec = Teuchos::rcp(new Epetra_Vector(fgiprev_->RangeMap(), true));
  fgiprev_->Apply(*duiinc_, *auxvec);
  tmpvec->Update(1.0, *auxvec, 1.0);
  // ---------End of term (6)

  // ---------Addressing term (7)
  if (fmgiprev_ != Teuchos::null)
  {
    /* For matrix-vector-product, the DomainMap() of the matrix and the Map()
     * of the vector have to match. DomaintMap() contains inner velocity DOFs
     * and all pressure DOFs. The inner ale displacement increment is converted
     * to the fluid map using AleToFluid(). This results in a map that contains
     * all velocity but no pressure DOFs.
     *
     * We have to circumvent some trouble with Epetra_BlockMaps since we cannot
     * split an Epetra_BlockMap into inner and interface DOFs.
     *
     * We create a map extractor 'velothermap' in order to extract the inner
     * velocity DOFs after calling AleToFluid(). Afterwards, a second map
     * extractor 'velotherpressuremapext' is used to append pressure DOFs filled
     * with zeros.
     *
     * Finally, maps match and matrix-vector-multiplication can be done.
     */

    // extract inner velocity DOFs after calling AleToFluid()
    Teuchos::RCP<Epetra_Map> velothermap = LINALG::SplitMap(
        *FluidField()->VelocityRowMap(),
        *InterfaceFluidAleCoupling().MasterDofMap());
    LINALG::MapExtractor velothermapext = LINALG::MapExtractor(
        *FluidField()->VelocityRowMap(), velothermap, false);
    auxvec = Teuchos::rcp(new Epetra_Vector(*velothermap, true));
    velothermapext.ExtractOtherVector(
        AleToFluid(FsiAleField()->FsiInterface()->InsertOtherVector(ddialeinc_)),
        auxvec);

    // add pressure DOFs
    LINALG::MapExtractor velotherpressuremapext = LINALG::MapExtractor(
        fmgiprev_->DomainMap(), velothermap);
    auxauxvec = Teuchos::rcp(new Epetra_Vector(fmgiprev_->DomainMap(), true));
    velotherpressuremapext.InsertCondVector(auxvec, auxauxvec);

    // prepare vector to store result of matrix-vector-product
    auxvec = Teuchos::rcp(new Epetra_Vector(fmgiprev_->RangeMap(), true));

    // Now, do the actual matrix-vector-product
    fmgiprev_->Apply(*auxauxvec, *auxvec);
    tmpvec->Update(1.0, *auxvec, 1.0);
  }
  // ---------End of term (7)

  // ---------Addressing term (8)
  if (firstcall_)
  {
    auxvec = Teuchos::rcp(new Epetra_Vector(fggprev_->RangeMap(), true));
    fggprev_->Apply(*FluidField()->ExtractInterfaceVeln(), *auxvec);
    tmpvec->Update(Dt() * timescale, *auxvec, 1.0);
  }
  // ---------End of term (8)

  // ---------Addressing term (2)
  auxvec = Teuchos::rcp(new Epetra_Vector(mortardinv->DomainMap(), true));
  mortardinv->Multiply(true, *tmpvec, *auxvec);
  lambda_->Update(scale, *auxvec, 1.0); // scale with ResidualScaling() to get [N/m^2]
  // ---------End of term (2)

  // Finally, divide by (1.0-ftiparam) which is common to all terms
  lambda_->Scale(-1.0 / (1.0 - ftiparam));

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
void FSI::SlidingMonolithicFluidSplit::CalculateInterfaceEnergyIncrement()
{
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();
  const double ftiparam = FluidField()->TimIntParam();

  // get the Mortar matrix M
  const Teuchos::RCP<LINALG::SparseMatrix> mortarm = coupsfm_->GetMMatrix();

  // interface traction weighted by time integration factors
  Teuchos::RCP<Epetra_Vector> tractionfluid = Teuchos::rcp(
      new Epetra_Vector(lambda_->Map(), true));
  Teuchos::RCP<Epetra_Vector> tractionstructure = Teuchos::rcp(
      new Epetra_Vector(*StructureField()->Interface()->FSICondMap(), true));
  tractionfluid->Update(stiparam - ftiparam, *lambdaold_, ftiparam - stiparam,
      *lambda_, 0.0);
  mortarm->Multiply(true, *tractionfluid, *tractionstructure);

  // displacement increment of this time step
  Teuchos::RCP<Epetra_Vector> deltad = Teuchos::rcp(
      new Epetra_Vector(*StructureField()->DofRowMap(), true));
  deltad->Update(1.0, *StructureField()->Dispnp(), -1.0,
      *StructureField()->Dispn(), 0.0);

  // calculate the energy increment
  double energy = 0.0;
  tractionstructure->Dot(
      *StructureField()->Interface()->ExtractFSICondVector(deltad), &energy);

  energysum_ += energy;

  WriteInterfaceEnergyFile(energy, energysum_);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::CheckKinematicConstraint()
{
  // some scaling factors for fluid
  const double timescale = FluidField()->TimeScaling();

  // get the Mortar matrices D and M
  const Teuchos::RCP<LINALG::SparseMatrix> mortard = coupsfm_->GetDMatrix();
  const Teuchos::RCP<LINALG::SparseMatrix> mortarm = coupsfm_->GetMMatrix();

#ifdef DEBUG
  if (mortarm == Teuchos::null)
  dserror("Expected Teuchos::rcp to mortar matrix M.");
  if (mortard == Teuchos::null)
  dserror("Expected Teuchos::rcp to mortar matrix D.");
#endif

  // get interface displacements and velocities
  const Teuchos::RCP<Epetra_Vector> disnp =
      StructureField()->ExtractInterfaceDispnp();
  const Teuchos::RCP<Epetra_Vector> disn =
      StructureField()->ExtractInterfaceDispn();
  const Teuchos::RCP<Epetra_Vector> velnp =
      FluidField()->ExtractInterfaceVelnp();
  const Teuchos::RCP<Epetra_Vector> veln = FluidField()->ExtractInterfaceVeln();

  // prepare vectors for projected interface quantities
  Teuchos::RCP<Epetra_Vector> disnpproj = Teuchos::rcp(
      new Epetra_Vector(mortarm->RangeMap(), true));
  Teuchos::RCP<Epetra_Vector> disnproj = Teuchos::rcp(
      new Epetra_Vector(mortarm->RangeMap(), true));
  Teuchos::RCP<Epetra_Vector> velnpproj = Teuchos::rcp(
      new Epetra_Vector(mortard->RangeMap(), true));
  Teuchos::RCP<Epetra_Vector> velnproj = Teuchos::rcp(
      new Epetra_Vector(mortard->RangeMap(), true));

  // projection of interface displacements
  mortarm->Apply(*disnp, *disnpproj);
  mortarm->Apply(*disn, *disnproj);

  // projection of interface velocities
  mortard->Apply(*velnp, *velnpproj);
  mortard->Apply(*veln, *velnproj);

  // calculate violation of kinematic interface constraint
  Teuchos::RCP<Epetra_Vector> violation = Teuchos::rcp(
      new Epetra_Vector(*disnpproj));
  violation->Update(-1.0, *disnproj, 1.0);
  violation->Update(-1.0 / timescale, *velnpproj, 1.0 / timescale, *velnproj,
      1.0);
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
      << "\nViolation of kinematic interface constraint:\n" << "L_2-norm: "
      << violationl2 << "        L_inf-norm: " << violationinf << "\n";
  Utils()->out().flags(flags);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::CheckDynamicEquilibrium()
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
  Teuchos::RCP<Epetra_Vector> tractionmaster = Teuchos::rcp(
      new Epetra_Vector(mortarm->DomainMap(), true));
  Teuchos::RCP<Epetra_Vector> tractionslave = Teuchos::rcp(
      new Epetra_Vector(mortard->DomainMap(), true));

  // calculate tractions on master and slave side
  mortarm->Multiply(true, *lambda_, *tractionmaster);
  mortard->Multiply(true, *lambda_, *tractionslave);

  // calculate violation of dynamic equilibrium
  Teuchos::RCP<Epetra_Vector> violation = Teuchos::rcp(
      new Epetra_Vector(*tractionmaster));
  violation->Update(-1.0, *tractionslave, 1.0);

  // calculate some norms
  double violationl2 = 0.0;
  double violationinf = 0.0;
  violation->Norm2(&violationl2);
  violation->NormInf(&violationinf);

  // scale L2-Norm with sqrt of length of interface vector
  violationl2 /= sqrt(
      FluidField()->Interface()->FSICondMap()->NumGlobalElements());

  // output to screen
  std::ios_base::fmtflags flags = Utils()->out().flags();

  Utils()->out() << std::scientific
      << "\nViolation of dynamic interface equilibrium:\n" << "L_2-norm: "
      << violationl2 << "        L_inf-norm: " << violationinf << "\n";
  Utils()->out().flags(flags);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::CombineFieldVectors(Epetra_Vector& v,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv,
    Teuchos::RCP<const Epetra_Vector> av, bool fullvectors)
{
  if (fullvectors)
  {
    // extract inner DOFs from slave vectors
    Teuchos::RCP<Epetra_Vector> fov =
        FsiFluidField()->FsiInterface()->ExtractOtherVector(fv);
    fov = FsiFluidField()->FsiInterface()->InsertOtherVector(fov);
    Teuchos::RCP<Epetra_Vector> aov =
        FsiAleField()->FsiInterface()->ExtractOtherVector(av);

    // put them together
    Extractor().AddVector(*sv, 0, v);
    Extractor().AddVector(*fov, 1, v);
    Extractor().AddVector(*aov, 2, v);
  }
  else
    FSI::Monolithic::CombineFieldVectors(v, sv, fv, av);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double FSI::SlidingMonolithicFluidSplit::SelectDtErrorBased() const
{
  // get time step size suggestions
  const double dtstr = GetAdaStrDt(); // based on all structure DOFs
  const double dtstrfsi = GetAdaStrFSIDt(); // based on structure FSI DOFs
  const double dtflinner = GetAdaFlInnerDt(); // based on inner fluid DOFs

  double dt = Dt();

  // select time step size based on error estimation
  if (IsAdaStructure() and IsAdaFluid())
    dt = std::min(std::min(dtstr, dtstrfsi), dtflinner);
  else if (IsAdaStructure() and (not IsAdaFluid()))
    dt = std::min(dtstr, dtstrfsi);
  else if ((not IsAdaStructure()) and IsAdaFluid())
    dt = dtflinner;
  else
  {
    // no change in time step size based on structure or fluid field error estimation
  }

  return dt;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool FSI::SlidingMonolithicFluidSplit::SetAccepted() const
{
  // get error norms
  const double strnorm = GetAdaStrnorm(); // based on all structure DOFs
  const double strfsinorm = GetAdaStrFSInorm(); // based on structure FSI DOFs
  const double flinnernorm = GetAdaFlInnerNorm(); // based on inner fluid DOFs

  bool accepted = std::max(strnorm, strfsinorm) < errtolstr_
      and flinnernorm < errtolfl_;

  // in case error estimation in the fluid field is turned off:
  if (not IsAdaFluid())
    accepted = std::max(strnorm, strfsinorm) < errtolstr_;

  // in case error estimation in the structure field is turned off:
  if (not IsAdaStructure())
    accepted = flinnernorm < errtolfl_;

  // no error based time adaptivity
  if ((not IsAdaStructure()) and (not IsAdaFluid()))
    accepted = true;

  return accepted;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::CreateSystemMatrix()
{
  FSI::BlockMonolithic::CreateSystemMatrix(systemmatrix_, false);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::CreateNodeOwnerRelationship(
    std::map<int, int>* nodeOwner,
    std::map<int, std::list<int> >* inverseNodeOwner,
    std::map<int, DRT::Node*>* fluidnodesPtr,
    std::map<int, DRT::Node*>* structuregnodesPtr,
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
  Teuchos::RCP<LINALG::SparseMatrix> P = coupsfm_->GetMortarTrafo(); // P matrix that couples fluid dofs with structure dofs
  Teuchos::RCP<Epetra_CrsMatrix> P_;
  P_ = P->EpetraMatrix();
  const Epetra_Map& P_Map = P_->RowMap();         // maps fluid dofs to procs

  int NumMyElements = P_Map.NumMyElements();
  int numStructureDofs = P_->NumGlobalCols();               // number of related structure dofs on interface

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


  for (int fldofLID = 0; fldofLID < NumMaxElements; fldofLID = fldofLID + skip)
  {

    flnode = -1;
    flowner = -1;
    stnode = -1;
    stowner = -1;

    int NumEntries = 0;
    double Values[numStructureDofs];
    int stdofGID[numStructureDofs];

    if (NumMyElements > 0)
    {
      int fldofGID = P_Map.GID(fldofLID); // gid of fluid dof
      // find related node and owner to fldofGID
      FindNodeRelatedToDof(fluidnodesPtr, fldofGID, fluiddis, re);   //fluid
      flnode = re[0];
      flowner = re[1];

      P_->ExtractGlobalRowCopy(fldofGID, numStructureDofs, NumEntries, Values, stdofGID);
    }

    // Loop over related structure dofs and get related nodes.

    int maxNumEntries;
    comm_.MaxAll(&NumEntries,&maxNumEntries,1);

    for (int j = 0; j < maxNumEntries; ++j){
      if (j < NumEntries){
        FindNodeRelatedToDof(structuregnodesPtr, stdofGID[j], structuredis, re);
        stnode = re[0];
        stowner = re[1];
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
        copy = stnode;
        comm_.Broadcast(&copy, 1, proc);    // send nodeid to check if the processor needs help to find its neighbor, code: copy=-2
        if (copy == -2){
          dofid = stdofGID[j];
          comm_.Broadcast(&dofid, 1, proc);
          FindNodeRelatedToDof(structuregnodesPtr, dofid, structuredis,re); // let each processor look for the node related to gstid
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

  } // end loop over fluid dofs


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

      if (myrank == proc)
        ++nodeOwnerIt;
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
