#ifdef CCADISCRET

#include "fsi_partitionedmonolithic.H"
#include "fsi_statustest.H"

#include "fsi_nox_linearsystem_partitioned.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::PartitionedMonolithic::PartitionedMonolithic(Epetra_Comm& comm)
  : Monolithic(comm)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::PartitionedMonolithic::SetupSystem()
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  linearsolverstrategy_ = Teuchos::getIntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");

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

  // fluid to ale at the interface

  icoupfa_.SetupConditionCoupling(*FluidField().Discretization(),
                                  FluidField().Interface(),
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
  vecSpaces.push_back(AleField()      .Interface().OtherMap());

  if (vecSpaces[1]->NumGlobalElements()==0)
    dserror("No inner fluid equations. Splitting not possible. Panic.");

  SetDofRowMaps(vecSpaces);

  // build ale system matrix in splitted system
  AleField().BuildSystemMatrix(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::PartitionedMonolithic::SetupRHS(Epetra_Vector& f, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::PartitionedMonolithic::SetupRHS");

  SetupVector(f,
              StructureField().RHS(),
              FluidField().RHS(),
              AleField().RHS(),
              FluidField().ResidualScaling());

#if 0
  if (firstcall)
  {
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockf = FluidField().BlockSystemMatrix();

    LINALG::SparseMatrix& fig = blockf->Matrix(0,1);
    LINALG::SparseMatrix& fgg = blockf->Matrix(1,1);

    Teuchos::RCP<Epetra_Vector> fveln = FluidField().ExtractInterfaceVeln();
    double timescale = FluidField().TimeScaling();
    double scale     = FluidField().ResidualScaling();

    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(fig.RowMap()));

    fig.Apply(*fveln,*rhs);
    rhs->Scale(timescale*Dt());

    Extractor().AddVector(*rhs,1,f);

    rhs = Teuchos::rcp(new Epetra_Vector(fgg.RowMap()));

    fgg.Apply(*fveln,*rhs);
    rhs->Scale(scale*timescale*Dt());

    rhs = FluidToStruct(rhs);
    rhs = StructureField().Interface().InsertCondVector(rhs);
    Extractor().AddVector(*rhs,0,f);
  }
#endif

  // NOX expects a different sign here.
  f.Scale(-1.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::PartitionedMonolithic::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::PartitionedMonolithic::InitialGuess");

  SetupVector(*ig,
              StructureField().InitialGuess(),
              FluidField().InitialGuess(),
              AleField().InitialGuess(),
              0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem>
FSI::PartitionedMonolithic::CreateLinearSystem(Teuchos::ParameterList& nlParams,
                                               NOX::Epetra::Vector& noxSoln,
                                               Teuchos::RCP<NOX::Utils> utils)
{
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys;

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  linSys = Teuchos::rcp(new NOX::FSI::LinearPartitionedSolver(printParams,
                                                              lsParams,
                                                              Extractor(),
                                                              *this,
                                                              StructureField(),
                                                              FluidField(),
                                                              AleField(),
                                                              linearsolverstrategy_));

  return linSys;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo> FSI::PartitionedMonolithic::CreateStatusTest(Teuchos::ParameterList& nlParams,
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
void FSI::PartitionedMonolithic::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                                     Teuchos::RCP<const Epetra_Vector>& sx,
                                                     Teuchos::RCP<const Epetra_Vector>& fx,
                                                     Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::PartitionedMonolithic::ExtractFieldVectors");

  // We have overlap at the interface. Thus we need the interface part of the
  // structure vector and append it to the fluid and ale vector. (With the
  // right translation.)

  sx = Extractor().ExtractVector(x,0);
  Teuchos::RCP<const Epetra_Vector> scx = StructureField().Interface().ExtractCondVector(sx);

  // process fluid unknowns

  Teuchos::RCP<const Epetra_Vector> fox = Extractor().ExtractVector(x,1);
  Teuchos::RCP<Epetra_Vector> fcx = StructToFluid(scx);

  FluidField().DisplacementToVelocity(fcx);

  Teuchos::RCP<Epetra_Vector> f = FluidField().Interface().InsertOtherVector(fox);
  FluidField().Interface().InsertCondVector(fcx, f);
  fx = f;

  // process ale unknowns

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x,2);
  Teuchos::RCP<Epetra_Vector> acx = StructToAle(scx);

  Teuchos::RCP<Epetra_Vector> a = AleField().Interface().InsertOtherVector(aox);
  AleField().Interface().InsertCondVector(acx, a);
  ax = a;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::PartitionedMonolithic::SetupVector(Epetra_Vector &f,
                                             Teuchos::RCP<const Epetra_Vector> sv,
                                             Teuchos::RCP<const Epetra_Vector> fv,
                                             Teuchos::RCP<const Epetra_Vector> av,
                                             double fluidscale)
{

  // extract the inner and boundary dofs of all three fields

  Teuchos::RCP<Epetra_Vector> fov = FluidField()    .Interface().ExtractOtherVector(fv);
  Teuchos::RCP<Epetra_Vector> aov = AleField()      .Interface().ExtractOtherVector(av);

  if (fluidscale!=0)
  {
    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> fcv = FluidField().Interface().ExtractCondVector(fv);
    Teuchos::RCP<Epetra_Vector> modsv = StructureField().Interface().InsertCondVector(FluidToStruct(fcv));
    modsv->Update(1.0, *sv, fluidscale);

    Extractor().InsertVector(*modsv,0,f);
  }
  else
  {
    Extractor().InsertVector(*sv,0,f);
  }

  Extractor().InsertVector(*fov,1,f);
  Extractor().InsertVector(*aov,2,f);
}

#endif
