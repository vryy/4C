/*!----------------------------------------------------------------------
\file fsi_fluidfluidmonolithic_structuresplit.cpp
\brief Control routine for monolithic fluid-fluid-fsi using XFEM and NOX

<pre>
Maintainer:  Shadan Shahmiri
             shahmiri@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15265
</pre>

*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "fsi_fluidfluidmonolithic_structuresplit.H"
#include "fsi_matrixtransform.H"
#include "fsi_overlapprec.H"
#include "fsi_overlapprec_fsiamg.H"
#include "fsi_statustest.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_fluid_fluid_fsi.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../linalg/linalg_utils.H"
#include "../drt_ale/ale_utils_mapextractor.H"
#include "../drt_ale/ale.H"

#include "../linalg/linalg_solver.H"

#include "fsi_nox_group.H"

#include "fsi_debugwriter.H"

#include "../drt_inpar/inpar_xfem.H"

#include "../drt_constraint/constraint_manager.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_inpar/inpar_ale.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidFluidMonolithicStructureSplit::FluidFluidMonolithicStructureSplit(const Epetra_Comm& comm,
                                                                           const Teuchos::ParameterList& timeparams)
  : BlockMonolithic(comm,timeparams)
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
    // remove interface DOFs from structural DBC map
    StructureField()->RemoveDirichCond(intersectionmap);

    // give a warning to the user that Dirichlet boundary conditions might not be correct

   if (comm.MyPID() == 0)
    {
      cout << "  +---------------------------------------------------------------------------------------------+" << endl;
      cout << "  |                                        PLEASE NOTE:                                         |" << endl;
      cout << "  +---------------------------------------------------------------------------------------------+" << endl;
      cout << "  | You run a monolithic structure split scheme. Hence, there are no structural interface DOFs. |" << endl;
      cout << "  | Structure Dirichlet boundary conditions on the interface will be neglected.                 |" << endl;
      cout << "  | Check whether you have prescribed appropriate DBCs on structural interface DOFs.            |" << endl;
      cout << "  +---------------------------------------------------------------------------------------------+" << endl;
    }
  }

  // ---------------------------------------------------------------------------

  const Teuchos::ParameterList& adyn     = DRT::Problem::Instance()->AleDynamicParams();
  int aletype = DRT::INPUT::IntegralValue<int>(adyn,"ALE_TYPE");

  if ((aletype!=INPAR::ALE::incr_lin) and (monolithic_approach_!=INPAR::XFEM::XFFSI_Full_Newton))
    dserror("Relaxing Ale Aprooach is just posiible with Ale-incr-lin!");


  sggtransform_ = Teuchos::rcp(new UTILS::MatrixRowColTransform);
  sgitransform_ = Teuchos::rcp(new UTILS::MatrixRowTransform);
  sigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  aigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  fmiitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmgitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  fsaigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fsmgitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  fscoupfa_ = Teuchos::rcp(new ADAPTER::Coupling());

  // Recovery of Lagrange multiplier happens on structure field
  lambda_   = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->FSICondMap(),true));
  ddiinc_   = Teuchos::null;
  soliprev_ = Teuchos::null;
  ddginc_   = Teuchos::null;
  duginc_   = Teuchos::null;
  disgprev_ = Teuchos::null;
  sgiprev_  = Teuchos::null;
  sggprev_  = Teuchos::null;

#ifdef DEBUG
  // check whether allocation was successful
  if (sggtransform_   == Teuchos::null) { dserror("Allocation of 'sggtransform_' failed."); }
  if (sgitransform_   == Teuchos::null) { dserror("Allocation of 'sgitransform_' failed."); }
  if (sigtransform_   == Teuchos::null) { dserror("Allocation of 'sigtransform_' failed."); }
  if (aigtransform_   == Teuchos::null) { dserror("Allocation of 'aigtransform_' failed."); }
  if (fmiitransform_  == Teuchos::null) { dserror("Allocation of 'fmiitransform_' failed."); }
  if (fmgitransform_  == Teuchos::null) { dserror("Allocation of 'fmgitransform_' failed."); }
  if (fsaigtransform_ == Teuchos::null) { dserror("Allocation of 'fsaigtransform_' failed."); }
  if (fsmgitransform_ == Teuchos::null) { dserror("Allocation of 'fsmgitransform_' failed."); }
  if (fscoupfa_       == Teuchos::null) { dserror("Allocation of 'fscoupfa_' failed."); }
  if (lambda_         == Teuchos::null) { dserror("Allocation of 'lambda_' failed."); }
#endif

  const Teuchos::ParameterList& xfluiddyn  = DRT::Problem::Instance()->XFluidDynamicParams();
  monolithic_approach_ = DRT::INPUT::IntegralValue<INPAR::XFEM::Monolithic_xffsi_Approach>
                         (xfluiddyn.sublist("GENERAL"),"MONOLITHIC_XFFSI_APPROACH");

  if (monolithic_approach_ == INPAR::XFEM::XFFSI_Full_Newton)
	 dserror("NOX-based XFFSI Approach does not work with XFFSI_Full_Newton!");

  relaxing_ale_ = (bool)DRT::INPUT::IntegralValue<int>(xfluiddyn.sublist("GENERAL"),"RELAXING_ALE");
  relaxing_ale_every_ = xfluiddyn.sublist("GENERAL").get<int>("RELAXING_ALE_EVERY");

  return;
}
/*----------------------------------------------------------------------*
 |  map and extractor containing the dofs with Dirichlet BC
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> FSI::FluidFluidMonolithicStructureSplit::CombinedDBCMap()
{
  const Teuchos::RCP<const Epetra_Map > scondmap = StructureField()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map> ffcondmap = FluidField().FluidDirichMaps();
  const Teuchos::RCP<const Epetra_Map > acondmap = AleField().GetDBCMapExtractor()->CondMap();

  // this is a structure split so we leave the fluid map unchanged. It
  // means that the dirichlet dofs could also contain fsi dofs. So if
  // you want to prescribe any dirichlet values at the fsi-interface it
  // is the fluid field which decides.

  std::vector<Teuchos::RCP<const Epetra_Map> > vectoroverallfsimaps;
  vectoroverallfsimaps.push_back(scondmap);
  vectoroverallfsimaps.push_back(ffcondmap);
  vectoroverallfsimaps.push_back(acondmap);

  Teuchos::RCP<Epetra_Map> overallfsidbcmaps = LINALG::MultiMapExtractor::MergeMaps(vectoroverallfsimaps);

  //structure and ale maps should not have any fsiCondDofs, so we
  //throw them away from overallfsidbcmaps

  std::vector<int> otherdbcmapvector; //vector of dbc
  const int mylength = overallfsidbcmaps->NumMyElements(); //on each prossesor (lids)
  const int* mygids = overallfsidbcmaps->MyGlobalElements();
  for (int i=0; i<mylength; ++i)
  {
    int gid = mygids[i];
    int fullmaplid = fullmap_->LID(gid);
    // if it is not a fsi dof
    if (fullmaplid >= 0)
      otherdbcmapvector.push_back(gid);
  }

  Teuchos::RCP<Epetra_Map> otherdbcmap = Teuchos::rcp(new Epetra_Map(-1, otherdbcmapvector.size(), &otherdbcmapvector[0], 0, Comm()));

  // Finally, create the global FSI Dirichlet map extractor
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor(*DofRowMap(),otherdbcmap,true));

  return otherdbcmap;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::SetupSystem()
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  linearsolverstrategy_ = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");
  SetDefaultParameters(fsidyn,NOXParameterList());
  ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  ADAPTER::Coupling& coupsa = StructureAleCoupling();
  ADAPTER::Coupling& coupfa = FluidAleCoupling();
  ADAPTER::Coupling& icoupfa = InterfaceFluidAleCoupling();
  const int ndim = DRT::Problem::Instance()->NDim();

  // structure to fluid
  coupsf.SetupConditionCoupling(*StructureField()->Discretization(),
                                 StructureField()->Interface()->FSICondMap(),
                                *FluidField().Discretization(),
                                 FluidField().Interface()->FSICondMap(),
                                "FSICoupling",
                                 ndim);

  // structure to ale
  coupsa.SetupConditionCoupling(*StructureField()->Discretization(),
                                 StructureField()->Interface()->FSICondMap(),
                                *AleField().Discretization(),
                                 AleField().Interface()->FSICondMap(),
                                "FSICoupling",
                                 ndim);

  // fluid to ale at the interface
  icoupfa.SetupConditionCoupling(*FluidField().Discretization(),
                                 FluidField().Interface()->FSICondMap(),
                                 *AleField().Discretization(),
                                 AleField().Interface()->FSICondMap(),
                                 "FSICoupling",
                                 ndim);

  // In the following we assume that both couplings find the same dof
  // map at the structural side. This enables us to use just one
  // interface dof map for all fields and have just one transfer
  // operator from the interface map to the full field map.
  if (not coupsf.MasterDofMap()->SameAs(*coupsa.MasterDofMap()))
    dserror("structure interface dof maps do not match");

  if (coupsf.MasterDofMap()->NumGlobalElements()==0)
    dserror("No nodes in matching FSI interface. Empty FSI coupling condition?");

  // the fluid-ale coupling always matches
  const Epetra_Map* embfluidnodemap = FluidField().Discretization()->NodeRowMap();
  const Epetra_Map* alenodemap   = AleField().Discretization()->NodeRowMap();

  coupfa.SetupCoupling(*FluidField().Discretization(),
                       *AleField().Discretization(),
                       *embfluidnodemap,
                       *alenodemap,
                       ndim);

  FluidField().SetMeshMap(coupfa.MasterDofMap());

  // create combined map
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;
  vecSpaces.push_back(StructureField()->Interface()->OtherMap());
  vecSpaces.push_back(FluidField()    .DofRowMap());
  vecSpaces.push_back(AleField()      .Interface()->OtherMap());

  if (vecSpaces[0]->NumGlobalElements()==0)
    dserror("No inner structural equations. Splitting not possible. Panic.");

  SetDofRowMaps(vecSpaces);
  fullmap_ = LINALG::MultiMapExtractor::MergeMaps(vecSpaces);

  // Use normal matrix for fluid equations but build (splitted) mesh movement
  // linearization (if requested in the input file)
  FluidField().UseBlockMatrix(false);

  // Use splitted structure matrix
  StructureField()->UseBlockMatrix();

  // build ale system matrix in splitted system
  AleField().BuildSystemMatrix(false);

  aleresidual_ = Teuchos::rcp(new Epetra_Vector(*AleField().Interface()->OtherMap()));

  // ---------------------------------------------------------------------------

  // get the PCITER from inputfile
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
void FSI::FluidFluidMonolithicStructureSplit::SetupNewSystem()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::SetupNewSystem()");

  // create combined map
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;
  vecSpaces.push_back(StructureField()->Interface()->OtherMap());
  vecSpaces.push_back(FluidField()    .DofRowMap());
  vecSpaces.push_back(AleField()      .Interface()->OtherMap());

  if (vecSpaces[0]->NumGlobalElements()==0)
    dserror("No inner structural equations. Splitting not possible. Panic.");

  SetDofRowMaps(vecSpaces);
  fullmap_ = LINALG::MultiMapExtractor::MergeMaps(vecSpaces);
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();


  // get the PCITER from inputfile
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
Teuchos::RCP<LINALG::BlockSparseMatrixBase> FSI::FluidFluidMonolithicStructureSplit::SystemMatrix() const
{
  return systemmatrix_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::SetupRHS(Epetra_Vector& f, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::SetupRHS");

  firstcall_ = firstcall;

  // We want to add into a zero vector
  f.PutScalar(0.0);

  // contributions of single field residuals
  SetupRHSResidual(f);

  // contributions of Lagrange multiplier from last time step
  // is done in SetupRHSResidual()
  //SetupRHSLambda(f);

  // contributions of special "first nonlinear iteration" terms
  if (firstcall)
    SetupRHSFirstiter(f);

  // NOX expects the 'positive' residual. The negative sign for the
  // linearized Newton system J*dx=-r is done internally by NOX.
  // Since we assembled the right hand side, we have to invert the sign here.
  f.Scale(-1.);

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::SetupRHSResidual(Epetra_Vector& f)
{

  // build the combined dirichlet boundary maps.
  combined_dbcmaps_ = CombinedDBCMap();
  if (combined_dbcmaps_== Teuchos::null) { dserror("Creation of FSI Dirichlet map extractor failed."); }

  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();
  const double ftiparam = FluidField().TimIntParam();

  // some scaling factors for fluid
  const double fluidscale = FluidField().ResidualScaling();

  // get single field residuals
  Teuchos::RCP<const Epetra_Vector> sv = Teuchos::rcp(new Epetra_Vector(*StructureField()->RHS()));
  Teuchos::RCP<const Epetra_Vector> fv = Teuchos::rcp(new Epetra_Vector(*FluidField().RHS()));
  Teuchos::RCP<const Epetra_Vector> av = Teuchos::rcp(new Epetra_Vector(*AleField().RHS()));

//  // extract only inner DOFs from structure (=slave) and ALE field
  Teuchos::RCP<const Epetra_Vector> sov = StructureField()->Interface()->ExtractOtherVector(sv);
  Teuchos::RCP<const Epetra_Vector> aov = AleField().Interface()->ExtractOtherVector(av);

  // add structure interface residual to fluid interface residual considering temporal scaling
  Teuchos::RCP<Epetra_Vector> scv = StructureField()->Interface()->ExtractFSICondVector(sv);
  Teuchos::RCP<Epetra_Vector> modfv = FluidField().Interface()->InsertFSICondVector(StructToFluid(scv));
  // modfv = modfv * 1/fluidscale * d/b
  modfv->Scale( (1./fluidscale)*(1.0-ftiparam)/(1.0-stiparam));

  // add contribution of Lagrange multiplier from previous time step
  if (lambda_ != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> lambdaglobal = FluidField().Interface()->InsertFSICondVector(StructToFluid(lambda_));
    modfv->Update((-ftiparam+stiparam*(1.0-ftiparam)/(1.0-stiparam))/fluidscale, *lambdaglobal, 1.0);
  }

  // we need a temporary vector with the whole fluid dofs where we
  // can insert veln which has the embedded dofrowmap into it
  Teuchos::RCP<Epetra_Vector> fluidfluidtmp = LINALG::CreateVector(*FluidField().DofRowMap(),true);
  xfluidfluidsplitter_->InsertFluidVector(modfv,fluidfluidtmp);

  Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(new const Epetra_Vector(fluidfluidtmp->Map(),true));
  LINALG::ApplyDirichlettoSystem(fluidfluidtmp,zeros,*(FluidField().FluidDirichMaps()));

  // all fluid dofs
  Teuchos::RCP<Epetra_Vector> fvfluidfluid = Teuchos::rcp_const_cast< Epetra_Vector >(fv);

  // adding modfv to fvfluidfluid
  fvfluidfluid->Update(1.0,*fluidfluidtmp,1.0);

  Extractor().InsertVector(*fvfluidfluid,1,f);
  Extractor().InsertVector(*sov,0,f);
  Extractor().InsertVector(*aov,2,f);

  // put the single field residuals together
  //FSI::Monolithic::CombineFieldVectors(f,sov,fvfluidfluid,aov);

  // add additional ale residual
  Extractor().AddVector(*aleresidual_,2,f);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::SetupRHSLambda(Epetra_Vector& f)
{
  // dummy function, Adding lambda to right hand-side happens in SetupRHSResidual

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::SetupRHSFirstiter(Epetra_Vector& f)
{

	// get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();
  const double ftiparam = FluidField().TimIntParam();

  // some scaling factors for fluid
  const double scale = FluidField().ResidualScaling();

  // old interface velocity of fluid field
  const Teuchos::RCP<const Epetra_Vector> fveln = FluidField().ExtractInterfaceVeln();

  // get fluid shape derivatives matrix
  const Teuchos::RCP<const LINALG::BlockSparseMatrixBase> mmm = FluidField().ShapeDerivatives();

  // get structure matrix
  const Teuchos::RCP<const LINALG::BlockSparseMatrixBase> blocks = StructureField()->BlockSystemMatrix();

  // get ale matrix
  const Teuchos::RCP<const LINALG::BlockSparseMatrixBase> blocka = AleField().BlockSystemMatrix();

#ifdef DEBUG
  if (blocks==Teuchos::null) { dserror("Expected Teuchos::rcp to structure block matrix."); }
  if (blocka==Teuchos::null) { dserror("Expected Teuchos::rcp to ALE block matrix."); }
#endif

  // extract submatrices
  const LINALG::SparseMatrix& sig = blocks->Matrix(0,1);  // S_{I\Gamma}
  const LINALG::SparseMatrix& sgg = blocks->Matrix(1,1);  // S_{\Gamma\Gamma}
  const LINALG::SparseMatrix& aig = blocka->Matrix(0,1);  // A_{I\Gamma}

  // some often re-used vectors
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null;  // right hand side of single set of DOFs

  // Different contributions/terms to the rhs are separated by the following comment line
  // ---------- inner structure DOFs
  /* The following terms are added to the inner structure DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  + S_{I \Gamma} * \Delta d_{\Gamma,p}
   *
   * (2)  - dt * S_{I \Gamma} * u^{n}_{\Gamma}
   *
   * Remarks on all terms:
   * +  tau: time scaling factor for interface time integration (tau = 1/FluidField().TimeScaling())
   *
   */
  // ----------adressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(sig.RowMap(),true));
  sig.Apply(*ddgpred_,*rhs);

  Extractor().AddVector(*rhs,0,f);
  // ----------end of term 1

  // ----------adressing term 2
  rhs = Teuchos::rcp(new Epetra_Vector(sig.RowMap(),true));
  sig.Apply(*FluidToStruct(fveln),*rhs);
  rhs->Scale(-Dt());

  Extractor().AddVector(*rhs,0,f);
  // ----------end of term 2
  // ----------end of inner structural DOFs

  // we need a temporary vector with the whole fluid dofs where we
  // can insert the embedded dofrowmap into it
  Teuchos::RCP<Epetra_Vector> fluidfluidtmp = LINALG::CreateVector(*FluidField().DofRowMap(),true);


  // ---------- inner fluid DOFs
  /* The following terms are added to the inner fluid DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  - dt * F^{G}_{I\Gamma} * u^{n}_{\Gamma}
   *
   */
   // ----------addressing term 1
   if (mmm!=Teuchos::null)
  {
    const LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);
    rhs = Teuchos::rcp(new Epetra_Vector(fmig.RowMap(),true));
    fmig.Apply(*fveln,*rhs);
    rhs->Scale(-Dt());

    rhs = FluidField().Interface()->InsertOtherVector(rhs);
    xfluidfluidsplitter_->InsertFluidVector(rhs,fluidfluidtmp);

    Extractor().AddVector(*fluidfluidtmp,1,f);
  }
  // ----------end of term 1
  // ---------- end of inner fluid DOFs

  // ---------- interface fluid DOFs
  /* The following terms are added to the interface fluid DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  - dt * F^{G}_{\Gamma\Gamma} * u^{n}_{\Gamma}
   *
   * (2)  - dt * (1-ftiparam)/(1-stiparam) * S_{\Gamma\Gamma} * u^{n}_{\Gamma}
   *
   * (3)  + (1-ftiparam)/(1-stiparam) * S_{\Gamma\Gamma} * \Delta d_{\Gamma,p}
   *
   * Remarks on all terms:
   * +  tau: time scaling factor for interface time integration (tau = 1/FluidField().TimeScaling())
   *
   */
  // ----------addressing term 1
  if (mmm!=Teuchos::null)
  {
    const LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);
    rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RowMap(),true));
    fmgg.Apply(*fveln,*rhs);
    rhs->Scale(-Dt());

    rhs = FluidField().Interface()->InsertFSICondVector(rhs);

    xfluidfluidsplitter_->InsertFluidVector(rhs,fluidfluidtmp);
    fluidfluidtmp->PutScalar(0.0);

     Extractor().AddVector(*fluidfluidtmp,1,f);
 }
  // ----------end of term 1

  // ----------addressing term 2
  rhs = Teuchos::rcp(new Epetra_Vector(sgg.RowMap(),true));
  sgg.Apply(*FluidToStruct(fveln),*rhs);
  rhs->Scale(-Dt() * (1.-ftiparam) / ((1.-stiparam) * scale));

  rhs = StructToFluid(rhs);
  rhs = FluidField().Interface()->InsertFSICondVector(rhs);

  fluidfluidtmp->PutScalar(0.0);
  xfluidfluidsplitter_->InsertFluidVector(rhs,fluidfluidtmp);

  Extractor().AddVector(*fluidfluidtmp,1,f);
  // ----------end of term 2

  // ----------addressing term 3
  rhs = Teuchos::rcp(new Epetra_Vector(sgg.RowMap(),true));
  sgg.Apply(*ddgpred_,*rhs);
  rhs->Scale((1.-ftiparam) / ((1.-stiparam) * scale));

  rhs = StructToFluid(rhs);
  rhs = FluidField().Interface()->InsertFSICondVector(rhs);

  fluidfluidtmp->PutScalar(0.0);
  xfluidfluidsplitter_->InsertFluidVector(rhs,fluidfluidtmp);

  Extractor().AddVector(*fluidfluidtmp,1,f);
  // ----------end of term 3
  // ---------- end of interface fluid DOFs

  // ---------- inner ale DOFs
  /* The following terms are added to the inner ale DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  - dt * A_{I\Gamma} * u^{n}_{\Gamma}
   *
   */
   // ----------addressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(aig.RowMap(),true));
  aig.Apply(*FluidToAleInterface(fveln),*rhs);
  rhs->Scale(-Dt());
  Extractor().AddVector(*rhs,2,f);
  // ----------end of term 1
  // ---------- end of inner ale DOFs

  // -----------------------------------------------------
  // Now, all contributions/terms to rhs in the first Newton iteration are added.

  // Remark: in standard fsi apply Dirichlet is missing here. But for safety we apply
  // Dirichlet condition to RHS..

  // Apply Dirichlet boundary conditions
  // structure
  rhs = Extractor().ExtractVector(f,0);
  Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(new const Epetra_Vector(rhs->Map(),true));
  LINALG::ApplyDirichlettoSystem(rhs,zeros,*(StructureField()->GetDBCMapExtractor()->CondMap()));
  Extractor().InsertVector(*rhs,0,f);

  // fluid
  rhs = Extractor().ExtractVector(f,1);
  zeros = Teuchos::rcp(new const Epetra_Vector(rhs->Map(),true));
  LINALG::ApplyDirichlettoSystem(rhs,zeros,*(FluidField().FluidDirichMaps()));
  Extractor().InsertVector(*rhs,1,f);

  // ale
  rhs = Extractor().ExtractVector(f,2);
  zeros = Teuchos::rcp(new const Epetra_Vector(rhs->Map(),true));
  LINALG::ApplyDirichlettoSystem(rhs,zeros,*(AleField().GetDBCMapExtractor()->CondMap()));
  Extractor().InsertVector(*rhs,2,f);

  // Reset quantities for previous iteration step since they still store values from the last time step
  ddiinc_   = LINALG::CreateVector(*StructureField()->Interface()->OtherMap(),true);
  soliprev_ = Teuchos::null;
  ddginc_   = LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(),true);
  duginc_   = LINALG::CreateVector(*FluidField().Interface()->FSICondMap(),true);
  disgprev_ = Teuchos::null;
  sgicur_   = Teuchos::null;
  sggcur_   = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
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
  // check whether allocation was successful
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

  //----------------------------------------------------------------------

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

  Teuchos::RCP<LINALG::SparseMatrix> lsgi = Teuchos::rcp(new LINALG::SparseMatrix(f->RowMap(),81,false));
  (*sgitransform_)(s->Matrix(1,0),
		  	  	  ((1.0-ftiparam)/(1.0-stiparam))*(1./scale),
                   ADAPTER::CouplingMasterConverter(coupsf),
                   *lsgi);

  lsgi->Complete(s->Matrix(1,0).DomainMap(),f->RangeMap());
  lsgi->ApplyDirichlet( *(FluidField().FluidDirichMaps()),false);

  mat.Assign(1,0,View,*lsgi);

  (*aigtransform_)(a->FullRowMap(),
                   a->FullColMap(),
                   aig,
                   1./timescale,
                   ADAPTER::CouplingSlaveConverter(icoupfa),
                   mat.Matrix(2,1));
  mat.Assign(2,2,View,aii);

  //----------------------------------------------------------------------
  // add optional fluid linearization with respect to mesh motion block

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField().ShapeDerivatives();
  if (mmm!=Teuchos::null)
  {
    const LINALG::SparseMatrix& fmii = mmm->Matrix(0,0);
    const LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);
    const LINALG::SparseMatrix& fmgi = mmm->Matrix(1,0);
    const LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);

    mat.Matrix(1,1).Add(fmgg,false,1./timescale,1.0);
    mat.Matrix(1,1).Add(fmig,false,1./timescale,1.0);

    Teuchos::RCP<LINALG::SparseMatrix> lfmgi = Teuchos::rcp(new LINALG::SparseMatrix(f->RowMap(),81,false));
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
    lfmgi->ApplyDirichlet( *(FluidField().FluidDirichMaps()),false );

    mat.Assign(1,2,View,*lfmgi);
  }

  f->Complete();
  f->ApplyDirichlet( *(FluidField().FluidDirichMaps()),true);

  // finally assign fluid block
  mat.Assign(1,1,View,*f);

  // done. make sure all blocks are filled.
  mat.Complete();

  // Finally, take care of Dirichlet boundary conditions
  mat.ApplyDirichlet(*combined_dbcmaps_,true);

  // store parts of structural matrix to know them in the next iteration as previous iteration matrices
  sgiprev_ = sgicur_;
  sggprev_ = sggcur_;
  sgicur_ = Teuchos::rcp(new LINALG::SparseMatrix(s->Matrix(1,0)));
  sggcur_ = Teuchos::rcp(new LINALG::SparseMatrix(s->Matrix(1,1)));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  const bool scaling_infnorm = (bool)DRT::INPUT::IntegralValue<int>(fsidyn,"INFNORMSCALING");

  if (scaling_infnorm)
  {
    // The matrices are modified here. Do we have to change them back later on?

    // do scaling of structure rows
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

    // do scaling of ale rows
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

    // do scaling of structure and ale rhs vectors
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
void FSI::FluidFluidMonolithicStructureSplit::UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
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

  std::ios_base::fmtflags flags = Utils()->out().flags();

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
Teuchos::RCP<NOX::Epetra::LinearSystem>
FSI::FluidFluidMonolithicStructureSplit::CreateLinearSystem(ParameterList& nlParams,
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
FSI::FluidFluidMonolithicStructureSplit::CreateStatusTest(Teuchos::ParameterList& nlParams,
                                                Teuchos::RCP<NOX::Epetra::Group> grp)
{
  // --------------------------------------------------------------------
  // Setup the test framework
  // --------------------------------------------------------------------
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


  // --------------------------------------------------------------------
  // setup tests for structural displacement field
  // --------------------------------------------------------------------
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

  // --------------------------------------------------------------------
  // setup tests for interface
  // --------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map> > interface;
  interface.push_back(FluidField().Interface()->FSICondMap());
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

  // --------------------------------------------------------------------
  // setup tests for fluid velocity field
  // --------------------------------------------------------------------
  // build mapextractor
  // inner fluid velocity Dofs: inner velocity dofs without db-condition
  std::vector<Teuchos::RCP<const Epetra_Map> > fluidvel;
  fluidvel.push_back(FluidField().InnerVelocityRowMap());
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

  // --------------------------------------------------------------------
  // setup tests for fluid pressure field
  // --------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map> > fluidpress;
  fluidpress.push_back(FluidField().PressureRowMap());
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                                        Teuchos::RCP<const Epetra_Vector>& sx,
                                                        Teuchos::RCP<const Epetra_Vector>& fx,
                                                        Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::ExtractFieldVectors");

#ifdef DEBUG
  if(ddgpred_ == Teuchos::null) { dserror("Vector 'ddgpred_' has not been initialized properly."); }
#endif

  // ---------------------------------------------------------------------------
  // process fluid unknowns
  // ---------------------------------------------------------------------------
  // extract fluid solution increment from NOX increment
  fx = Extractor().ExtractVector(x,1);

  // extract embedded fluid vector
  Teuchos::RCP<Epetra_Vector> fx_emb = xfluidfluidsplitter_->ExtractFluidVector(fx);
  // extract fsi-fluid vector
  Teuchos::RCP<Epetra_Vector> fcx = FluidField().Interface()->ExtractFSICondVector(fx_emb);

  // ----------------------
  // process ale unknowns
  Teuchos::RCP<Epetra_Vector> fcxforale = FluidField().Interface()->ExtractFSICondVector(fx_emb);
  FluidField().VelocityToDisplacement(fcxforale);
  Teuchos::RCP<Epetra_Vector> acx = FluidToStruct(fcxforale);
  acx = StructToAle(acx);

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x,2);
  // Teuchos::RCP<Epetra_Vector> acx = StructToAle(scx);

  Teuchos::RCP<Epetra_Vector> a = AleField().Interface()->InsertOtherVector(aox);
  AleField().Interface()->InsertFSICondVector(acx, a);

  ax = a;

  // ---------------------------------------------------------------------------
  // process structure unknowns
  // ---------------------------------------------------------------------------
  // extract inner structure solution increment from NOX increment
  Teuchos::RCP<const Epetra_Vector> sox = Extractor().ExtractVector(x,0);

  // convert ALE interface displacements to structure interface displacements
  Teuchos::RCP<Epetra_Vector> scx = AleToStruct(acx);
  scx->Update(-1.0, *ddgpred_, 1.0);

  // put inner and interface structure solution increments together
  Teuchos::RCP<Epetra_Vector> s = StructureField()->Interface()->InsertOtherVector(sox);
  StructureField()->Interface()->InsertFSICondVector(scx, s);
  sx = s;

  // ---------------------------------------------------------------------------

  // Store field vectors to know them later on as previous quantities
  if (soliprev_ != Teuchos::null)
    ddiinc_->Update(1.0, *sox, -1.0, *soliprev_, 0.0);  // compute current iteration increment
  else
    ddiinc_ = Teuchos::rcp(new Epetra_Vector(*sox));    // first iteration increment

  soliprev_ = sox;                                      // store current step increment

  if (disgprev_ != Teuchos::null)
    ddginc_->Update(1.0, *scx, -1.0, *disgprev_, 0.0);  // compute current iteration increment
  else
    ddginc_ = Teuchos::rcp(new Epetra_Vector(*scx));    // first iteration increment

  disgprev_ = scx;                                      // store current step increment

  if (velgprev_ != Teuchos::null)
    duginc_->Update(1.0, *fcx, -1.0, *velgprev_, 0.0);  // compute current iteration increment
  else
    duginc_ = Teuchos::rcp(new Epetra_Vector(*fcx));    // first iteration increment

  velgprev_ = fcx;                                      // store current step increment
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::PrepareTimeStep()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::PrepareTimeStep");

  IncrementTimeAndStep();

  PrintHeader();

  StructureField()->PrepareTimeStep();
  FluidField().    PrepareTimeStep();
  AleField().      PrepareTimeStep();

  if (monolithic_approach_!=INPAR::XFEM::XFFSI_Full_Newton and (Step()!=1))
    SetupNewSystem();

  //xfluidfluid splitter
  xfluidfluidsplitter_ = FluidField().XFluidFluidMapExtractor();

  precondreusecount_ = 0;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::Output()
{

  StructureField()->Output();

  // output Lagrange multiplier
  {
    /* 'lambda_' is only defined on the interface. So, insert 'lambda_' into
    * 'lambdafull' that is defined on the entire structure field. Then, write
    * output or restart data.
    */
    Teuchos::RCP<Epetra_Vector> lambdafull = StructureField()->Interface()->InsertFSICondVector(lambda_);
    const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
    const int uprestart = fsidyn.get<int>("RESTARTEVRY");
    const int upres = fsidyn.get<int>("UPRES");
    if ((uprestart != 0 && FluidField().Step() % uprestart == 0) || FluidField().Step() % upres == 0)
      StructureField()->DiscWriter()->WriteVector("fsilambda", lambdafull);
  }

  FluidField().    Output();
  //FluidField().    OutputReducedD();
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
void FSI::FluidFluidMonolithicStructureSplit::ReadRestart(int step)
{
  // read Lagrange multiplier
  {
    Teuchos::RCP<Epetra_Vector> lambdafull = Teuchos::rcp(new Epetra_Vector(*StructureField()->DofRowMap(),true));
    IO::DiscretizationReader reader = IO::DiscretizationReader(StructureField()->Discretization(),step);
    reader.ReadVector(lambdafull, "fsilambda");
    lambda_ = StructureField()->Interface()->ExtractFSICondVector(lambdafull);
  }

  StructureField()->ReadRestart(step);
  FluidField().ReadRestart(step);
  FluidField().ReadRestartReducedD(step);
  AleField().ReadRestart(step);

  SetTimeStep(FluidField().Time(),FluidField().Step());
}

/*----------------------------------------------------------------------*/
/* Recover the Lagrange multiplier at the interface   mayr.mt (03/2012) */
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::RecoverLagrangeMultiplier()
{
  // get time integration parameter of structural time integrator
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();

  // some scaling factors for fluid
  //const double timescale = FluidField().TimeScaling();

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
   * (2)  + 1. / (1.-stiparam) * tmpvec
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
   * +  Division by (1.0 - stiparam) will be done in the end
   *    since this is common to all terms
   * +  tau: time scaling factor for interface time integration (tau = 1/FluidField().TimeScaling())
   * +  neglecting terms (4)-(6) should not alter the results significantly
   *    since at the end of the time step the solution increments tend to zero.
   *
   *                                                 Matthias Mayr (10/2012)
   */

  // ---------Addressing term (1)
  lambda_->Update(-stiparam,*lambda_,0.0);
  // ---------End of term (1)

  // ---------Addressing term (3)
  Teuchos::RCP<Epetra_Vector> structureresidual = StructureField()->Interface()->ExtractFSICondVector(StructureField()->RHS());
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
  auxvec = Teuchos::rcp(new Epetra_Vector(sggprev_->RangeMap(),true));
  sggprev_->Apply(*FluidToStruct(duginc_),*auxvec);
  tmpvec->Update(1.0/timescale,*auxvec,1.0);
  // ---------End of term (5)

  //---------Addressing term (6)
  if (firstcall_)
  {
    auxvec = Teuchos::rcp(new Epetra_Vector(sggprev_->RangeMap(),true));
    sggprev_->Apply(*FluidToStruct(FluidField().ExtractInterfaceVeln()),*auxvec);
    tmpvec->Update(Dt(),*auxvec,1.0);
  }
  // ---------End of term (6)
   *
   */

  // ---------Addressing term (2)
  lambda_->Update(1.0,*tmpvec,1.0);
  // ---------End of term (2)

  // finally, divide by -(1.-stiparam) which is common to all terms
  lambda_->Scale(1./(1.0-stiparam));

  // Finally, the Lagrange multiplier 'lambda_' is recovered here.
  // It represents nodal forces acting onto the structure.

  return;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::CombineFieldVectors(Epetra_Vector& v,
		                                      Teuchos::RCP<const Epetra_Vector> sv,
                                              Teuchos::RCP<const Epetra_Vector> fv,
                                              Teuchos::RCP<const Epetra_Vector> av,
                                              bool fullvectors)
{
  if (fullvectors)
  {
    // extract inner DOFs from slave vectors
    Teuchos::RCP<Epetra_Vector> sov = StructureField()->Interface()->ExtractOtherVector(sv);
    Teuchos::RCP<Epetra_Vector> aov = AleField().Interface()->ExtractOtherVector(av);

    // put them together
    Extractor().AddVector(*sov,0,v);
    Extractor().AddVector(*fv,1,v);
    Extractor().AddVector(*aov,2,v);
  }
  else
    FSI::Monolithic::CombineFieldVectors(v,sv,fv,av);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::Update()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::Update");

  //bool aleupdate = false;
  //if (relaxing_ale_ == true)
  //  if (currentstep_%relaxing_ale_every_==0) aleupdate = true;
  //--------------------------------------------------------
  // Ale Update in every time step
  bool aleupdate = true;


  if (monolithic_approach_!= INPAR::XFEM::XFFSI_Full_Newton)
  {
    FluidField().ApplyEmbFixedMeshDisplacement(AleToFluid(AleField().ExtractDispnp()));
  }

  // calling RecoverLagrangeMultiplier() in nonox inmplementaion happens here.
  // RecoverLagrangeMultiplier();


  if (monolithic_approach_!= INPAR::XFEM::XFFSI_Full_Newton and aleupdate)
  {
    if (Comm().MyPID() == 0)
      cout << "Relaxing Ale.." << endl;
    AleField().SolveAleXFluidFluidFSI();
    FluidField().ApplyMeshDisplacement(AleToFluid(AleField().ExtractDispnp()));
  }


  StructureField()->Update();
  FluidField().    Update();
  AleField().      Update();


  if (monolithic_approach_!= INPAR::XFEM::XFFSI_Full_Newton and aleupdate)
  {
    // build ale system matrix for the next time step. Here first we
    // update the vectors then we set the fluid-fluid dirichlet values
    // in buildsystemmatrix
    AleField().BuildSystemMatrix(false);
    aleresidual_ = Teuchos::rcp(new Epetra_Vector(*AleField().Interface()->OtherMap()));
  }
}
