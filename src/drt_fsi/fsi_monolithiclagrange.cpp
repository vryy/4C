#ifdef CCADISCRET

#include "fsi_monolithiclagrange.H"
#include "fsi_statustest.H"
#include "fsi_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicLagrange::MonolithicLagrange(Epetra_Comm& comm)
  : Monolithic(comm)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  SetDefaultParameters(fsidyn,NOXParameterList());

  // right now we use matching meshes at the interface

  ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  ADAPTER::Coupling& coupsa = StructureAleCoupling();
  ADAPTER::Coupling& coupfa = FluidAleCoupling();

  // structure to fluid

  coupsf.SetupConditionCoupling(*StructureField().Discretization(),
                                StructureField().Interface(),
                                *FluidField().Discretization(),
                                FluidField().Interface(),
                                "FSICoupling");

  // structure to ale

  coupsa.SetupConditionCoupling(*StructureField().Discretization(),
                                StructureField().Interface(),
                                *AleField().Discretization(),
                                AleField().Interface(),
                                "FSICoupling");

  // In the following we assume that both couplings find the same dof
  // map at the structural side. This enables us to use just one
  // interface dof map for all fields and have just one transfer
  // operator from the interface map to the full field map.
  if (not coupsf.MasterDofMap()->SameAs(*coupsa.MasterDofMap()))
    dserror("structure interface dof maps do not match");

  if (coupsf.MasterDofMap()->NumGlobalElements()==0)
    dserror("No nodes in matching FSI interface. Empty FSI coupling condition?");

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = FluidField().Discretization()->NodeRowMap();
  const Epetra_Map* alenodemap   = AleField().Discretization()->NodeRowMap();

  coupfa.SetupCoupling(*FluidField().Discretization(),
                       *AleField().Discretization(),
                       *fluidnodemap,
                       *alenodemap);

  FluidField().SetMeshMap(coupfa.MasterDofMap());

  // create combined map

  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;
  vecSpaces.push_back(StructureField().DofRowMap());
  vecSpaces.push_back(FluidField()    .Interface().OtherMap());
  vecSpaces.push_back(FluidField()    .Interface().CondMap());
  vecSpaces.push_back(AleField()      .Interface().OtherMap());
  vecSpaces.push_back(AleField()      .Interface().CondMap());
  vecSpaces.push_back(UTILS::ShiftMap(StructureField().Interface().CondMap(),vecSpaces));
  vecSpaces.push_back(UTILS::ShiftMap(StructureField().Interface().CondMap(),vecSpaces));

  SetDofRowMaps(vecSpaces);

  // create block system matrix

  systemmatrix_ = Teuchos::rcp(new LagrangianBlockMatrix(Extractor(),
                                                         StructureField().LinearSolver(),
                                                         FluidField().LinearSolver(),
                                                         AleField().LinearSolver()));

  /*----------------------------------------------------------------------*/
  // Switch fluid to interface split block matrix
  FluidField().UseBlockMatrix(FluidField().Interface(),
                              FluidField().Interface());

  // build ale system matrix in splitted system
  AleField().BuildSystemMatrix(false);

  // setup Lagrangian coupling terms

  coupsf.SetupCouplingMatrices(*vecSpaces[5],*StructureField().DofRowMap(),*FluidField().Interface().CondMap());
  coupsa.SetupCouplingMatrices(*vecSpaces[6],*StructureField().DofRowMap(),*AleField()  .Interface().CondMap());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicLagrange::SetupRHS(Epetra_Vector& f)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicLagrange::SetupRHS");

  if (Comm().MyPID()==0)
    std::cout << "FSI::MonolithicLagrange::SetupRHS\n";

  SetupVector(f,
              StructureField().RHS(),
              FluidField().RHS(),
              AleField().RHS(),
              FluidField().ResidualScaling());

  // NOX expects a different sign here.
  f.Scale(-1.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicLagrange::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicLagrange::SetupSystemMatrix");

  if (Comm().MyPID()==0)
    std::cout << "FSI::MonolithicLagrange::SetupSystemMatrix\n";

  // extract Jacobian matrices and put them into composite system
  // matrix W

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  const ADAPTER::Coupling& coupsa = StructureAleCoupling();
  //const ADAPTER::Coupling& coupfa = FluidAleCoupling();

  Teuchos::RCP<LINALG::SparseMatrix> s = StructureField().SystemMatrix();
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> f = FluidField().BlockSystemMatrix();

  LINALG::SparseMatrix& fii = f->Matrix(0,0);
  LINALG::SparseMatrix& fig = f->Matrix(0,1);
  LINALG::SparseMatrix& fgi = f->Matrix(1,0);
  LINALG::SparseMatrix& fgg = f->Matrix(1,1);

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField().BlockSystemMatrix();
  LINALG::SparseMatrix& aii = a->Matrix(0,0);
  LINALG::SparseMatrix& aig = a->Matrix(0,1);
  LINALG::SparseMatrix& agi = a->Matrix(1,0);
  LINALG::SparseMatrix& agg = a->Matrix(1,1);

  mat.Assign(0,0,View,*s);

  mat.Assign(1,1,View,fii);
  mat.Assign(1,2,View,fig);
  mat.Assign(2,1,View,fgi);
  mat.Assign(2,2,View,fgg);

  mat.Assign(3,3,View,aii);
  mat.Assign(3,4,View,aig);
  mat.Assign(4,3,View,agi);
  mat.Assign(4,4,View,agg);

  // note: these scaling factors are not guaranteed to be available until the
  // elements have been evaluated.

  double timescale = FluidField().TimeScaling();
  double scale     = FluidField().ResidualScaling();

  mat.Matrix(5,0).Add(*coupsf.MasterToMasterMat(),false,1.0,0.0);
  mat.Matrix(5,2).Add(*coupsf.SlaveToMasterMat() ,false,-timescale,0.0);

  mat.Matrix(0,5).Add(*coupsf.MasterToMasterMatTrans(),false,1.0,0.0);
  mat.Matrix(2,5).Add(*coupsf.SlaveToMasterMatTrans() ,false,-1./scale,0.0);

  mat.Matrix(6,0).Add(*coupsa.MasterToMasterMat(),false,1.0,0.0);
  mat.Matrix(6,4).Add(*coupsa.SlaveToMasterMat() ,false,-1.0,0.0);

  mat.Matrix(4,6).Add(*coupsa.SlaveToMasterMatTrans(),false,-1.0,0.0);

  // done. make sure all blocks are filled.
  mat.Complete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicLagrange::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicLagrange::InitialGuess");

  SetupVector(*ig,
              StructureField().InitialGuess(),
              FluidField().InitialGuess(),
              AleField().InitialGuess(),
              0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicLagrange::SetupVector(Epetra_Vector &f,
                                          Teuchos::RCP<const Epetra_Vector> sv,
                                          Teuchos::RCP<const Epetra_Vector> fv,
                                          Teuchos::RCP<const Epetra_Vector> av,
                                          double fluidscale)
{
  f.PutScalar(0.);

  Extractor().InsertVector(*sv,0,f);

  Teuchos::RCP<Epetra_Vector> fov = FluidField()    .Interface().ExtractOtherVector(fv);
  Teuchos::RCP<Epetra_Vector> fcv = FluidField()    .Interface().ExtractCondVector (fv);
  Teuchos::RCP<Epetra_Vector> aov = AleField()      .Interface().ExtractOtherVector(av);
  Teuchos::RCP<Epetra_Vector> acv = AleField()      .Interface().ExtractCondVector (av);

  Extractor().InsertVector(*fov,1,f);
  Extractor().InsertVector(*fcv,2,f);

  Extractor().InsertVector(*aov,3,f);
  Extractor().InsertVector(*acv,4,f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem>
FSI::MonolithicLagrange::CreateLinearSystem(ParameterList& nlParams,
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

  linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,
                                                             lsParams,
                                                             Teuchos::rcp(iJac,false),
                                                             J,
                                                             Teuchos::rcp(iPrec,false),
                                                             M,
                                                             noxSoln));

  return linSys;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo>
FSI::MonolithicLagrange::CreateStatusTest(Teuchos::ParameterList& nlParams,
                                          Teuchos::RCP<NOX::Epetra::Group> grp)
{
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo       = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  Teuchos::RCP<NOX::StatusTest::Combo> converged   = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 100)));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv    = Teuchos::rcp(new NOX::StatusTest::FiniteValue);

  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // require one solve
  converged->addStatusTest(Teuchos::rcp(new NOX::FSI::MinIters(1)));

  // setup tests for structural displacements

  std::vector<Teuchos::RCP<const Epetra_Map> > structdisp;
  structdisp.push_back(StructureField().Interface().OtherMap());
  structdisp.push_back(Teuchos::null);
  LINALG::MultiMapExtractor structdispextract(*DofRowMap(),structdisp);

  Teuchos::RCP<NOX::StatusTest::Combo> structcombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> structureDisp =
    Teuchos::rcp(new NOX::FSI::PartialNormF("displacement",
                                            structdispextract,0,
                                            nlParams.get("Norm abs disp", 1.0e-6),
                                            NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate =
    Teuchos::rcp(new NOX::FSI::PartialNormUpdate("displacement update",
                                                 structdispextract,0,
                                                 nlParams.get("Norm abs disp", 1.0e-6),
                                                 NOX::FSI::PartialNormUpdate::Scaled));

  AddStatusTest(structureDisp);
  structcombo->addStatusTest(structureDisp);
  //structcombo->addStatusTest(structureDispUpdate);

  converged->addStatusTest(structcombo);

  // test for interface forces

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  Teuchos::RCP<ADAPTER::Coupling::Converter> converter =
    Teuchos::rcp(new ADAPTER::Coupling::SlaveConverter(coupsf));

  Teuchos::RCP<NOX::StatusTest::Combo> interfacecombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialSumNormF> interfaceForce =
    Teuchos::rcp(new NOX::FSI::PartialSumNormF("interface",
                                               LINALG::MapExtractor(*DofRowMap(),
                                                                    StructureField().Interface().CondMap(),
                                                                    Teuchos::null),
                                               1.0,
                                               LINALG::MapExtractor(*DofRowMap(),
                                                                    FluidField().Interface().CondMap(),
                                                                    Teuchos::null),
                                               FluidField().ResidualScaling(),
                                               converter,
                                               nlParams.get("Norm abs disp", 1.0e-6),
                                               NOX::FSI::PartialNormF::Scaled));

  AddStatusTest(interfaceForce);
  interfacecombo->addStatusTest(interfaceForce);
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
void FSI::MonolithicLagrange::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                                  Teuchos::RCP<const Epetra_Vector>& sx,
                                                  Teuchos::RCP<const Epetra_Vector>& fx,
                                                  Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicLagrange::ExtractFieldVectors");

  sx = Extractor().ExtractVector(x,0);

  Teuchos::RCP<const Epetra_Vector> fox = Extractor().ExtractVector(x,1);
  Teuchos::RCP<const Epetra_Vector> fcx = Extractor().ExtractVector(x,2);

  Teuchos::RCP<Epetra_Vector> f = FluidField().Interface().InsertOtherVector(fox);
  FluidField().Interface().InsertCondVector(fcx, f);
  fx = f;

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x,3);
  Teuchos::RCP<const Epetra_Vector> acx = Extractor().ExtractVector(x,4);

  Teuchos::RCP<Epetra_Vector> a = AleField().Interface().InsertOtherVector(aox);
  AleField().Interface().InsertCondVector(acx, a);
  ax = a;
}


#endif
