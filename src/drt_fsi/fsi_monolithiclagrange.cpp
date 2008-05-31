#ifdef CCADISCRET

#include "fsi_monolithiclagrange.H"
#include "fsi_statustest.H"

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
  vecSpaces.push_back(AleField()      .Interface().OtherMap());

  SetDofRowMaps(vecSpaces);

  // create block system matrix

  systemmatrix_ = Teuchos::rcp(new LagrangianBlockMatrix(Extractor(),
                                                         StructureField().LinearSolver(),
                                                         FluidField().LinearSolver(),
                                                         AleField().LinearSolver()));

  // build ale system matrix nonsplitted
  AleField().BuildSystemMatrix(true);
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
  const ADAPTER::Coupling& coupfa = FluidAleCoupling();
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
}


#endif
