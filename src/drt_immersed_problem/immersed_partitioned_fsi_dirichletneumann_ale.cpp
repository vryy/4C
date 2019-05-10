/*!----------------------------------------------------------------------

\brief partitioned immersed fsi algorithm for neumann-neumann like coupling (volume force coupling)

\level 2

\maintainer Jonas Eichinger

*----------------------------------------------------------------------*/
#include "immersed_partitioned_fsi_dirichletneumann_ale.H"

#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"
#include "../drt_adapter/ad_fld_fluid_ale_immersed.H"

#include "../drt_structure/stru_aux.H"

// aitken
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_fsi/fsi_nox_fixpoint.H"
#include "../drt_fsi/fsi_nox_aitken_immersed_ale.H"

#include "immersed_field_exchange_manager.H"


IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::ImmersedPartitionedFSIDirichletNeumannALE(
    const Epetra_Comm& comm)
    : ImmersedPartitionedFSIDirichletNeumann(comm),
      combined_newstate_(Teuchos::null),
      idispnp_(Teuchos::null),
      ivelnp_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::Setup()
{
  // call setup of base class
  IMMERSED::ImmersedPartitionedFSIDirichletNeumann::Setup();

  /// Immersed+FSI interface vector of new state to use in residual calculation
  combined_newstate_ =
      Teuchos::rcp(new Epetra_Vector(*(immersedstructure_->CombinedInterface()->FullMap()), true));

  DRT::ImmersedFieldExchangeManager::Instance()->SetAdapter(immersedstructure_);

  // check whether deformable background mesh is requested
  isALE_ =
      (globalproblem_->ImmersedMethodParams().get<std::string>("DEFORM_BACKGROUND_MESH") == "yes");
  if (isALE_ and myrank_ == 0) std::cout << " Deformable background mesh is used." << std::endl;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::SetupCoupling(
    const Teuchos::ParameterList& fsidyn, const Epetra_Comm& comm)
{
  if (myrank_ == 0)
    std::cout
        << "\n Setup ALE coupling at FSI Interface for ImmersedPartitionedFSIDirichletNeumannALE :"
        << std::endl;
  FSI::Partitioned::SetupCoupling(fsidyn, comm);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::SetStatesFluidOP()
{
  // for FluidOP
  structdis_->SetState(0, "displacement", immersedstructure_->Dispnp());
  structdis_->SetState(0, "velocity", immersedstructure_->Velnp());
  fluiddis_->SetState(0, "dispnp", MBFluidField()->FluidField()->Dispnp());


  // get displacements of ALE FSI interface
  idispnp_ = immersedstructure_->ExtractInterfaceDispnp();

  // store original input information
  std::string store_input = globalproblem_->FSIDynamicParams().get<std::string>("COUPALGO");

  // update current position of fluid nodes (only ALE is to compute)
  globalproblem_->getNonconstParameterList()
      .
      operator*()
      .sublist("FSI DYNAMIC")
      .set<std::string>("COUPALGO", "pseudo_structure");
  MBFluidField()->NonlinearSolve(StructToFluid(idispnp_), Teuchos::null);

  // substitute with original input information
  globalproblem_->getNonconstParameterList()
      .
      operator*()
      .sublist("FSI DYNAMIC")
      .set<std::string>("COUPALGO", store_input);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::SetStatesVelocityCorrection()
{
  fluiddis_->SetState(0, "velnp", MBFluidField()->FluidField()->Velnp());
  fluiddis_->SetState(0, "dispnp", MBFluidField()->FluidField()->Dispnp());

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::SetStatesStructOP()
{
  // for StructOP
  fluiddis_->SetState(0, "velnp", MBFluidField()->FluidField()->Velnp());
  fluiddis_->SetState(0, "dispnp", MBFluidField()->FluidField()->Dispnp());

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::InitialGuess()
{
  if (displacementcoupling_)
    return immersedstructure_->PredictFullInterfaceDispnp();
  else
  {
    Teuchos::RCP<Epetra_Vector> combined_traction = Teuchos::rcp(
        new Epetra_Vector(*(immersedstructure_->CombinedInterface()->FullMap()), true));

    // insert Immersed and FSI vector to combined vector (vector of whole interface)
    immersedstructure_->CombinedInterface()->InsertOtherVector(
        immersedstructure_->Interface()->ExtractIMMERSEDCondVector(struct_bdry_traction_),
        combined_traction);
    immersedstructure_->CombinedInterface()->InsertCondVector(
        FluidToStruct(MBFluidField()->ExtractInterfaceForces()), combined_traction);

    return combined_traction;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::FluidOp(
    Teuchos::RCP<Epetra_Vector> bforce, const FillType fillFlag)
{
  // the displacement -> velocity conversion at the ALE FSI interface
  idispnp_ = immersedstructure_->ExtractInterfaceDispnp();
  ivelnp_ = InterfaceVelocity(idispnp_);

  // immersed part
  IMMERSED::ImmersedPartitionedFSIDirichletNeumann::FluidOp(bforce, fillFlag);

  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::SolveFluid()
{
  MBFluidField()->NonlinearSolve(StructToFluid(idispnp_), StructToFluid(ivelnp_));

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::StructOp(
    Teuchos::RCP<Epetra_Vector> struct_bdry_traction, const FillType fillFlag)
{
  if (fillFlag == User)
  {
    dserror("fillFlag==User not implemented");
    return Teuchos::null;
  }
  else
  {
    // immersed part
    IMMERSED::ImmersedPartitionedFSIDirichletNeumann::StructOp(struct_bdry_traction, fillFlag);

    // add FSI part
    immersedstructure_->CombinedInterface()->AddCondVector(
        FluidToStruct(MBFluidField()->ExtractInterfaceForces()), struct_bdry_traction);

    // solve the structure
    IMMERSED::ImmersedPartitionedFSIDirichletNeumann::SolveStruct();

    return ExtractInterfaceDispnp();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::UpdateCurrentPositionsFluidNodes()
{
  // set and get displacement state (export happens inside)
  fluiddis_->SetState("dispnp", fluid_->FluidField()->Dispnp());
  Teuchos::RCP<const Epetra_Vector> displacements = fluiddis_->GetState("dispnp");

  // get number of column nodes on this proc
  int nummycolnodes = fluiddis_->NumMyColNodes();

  // update positions
  for (int lid = 0; lid < nummycolnodes; ++lid)
  {
    const DRT::Node* node = fluiddis_->lColNode(lid);
    LINALG::Matrix<3, 1> currpos;
    std::vector<int> dofstoextract(4);
    std::vector<double> mydisp(4);

    // get the current displacement
    fluiddis_->Dof(node, 0, dofstoextract);
    DRT::UTILS::ExtractMyValues(*displacements, mydisp, dofstoextract);

    currpos(0) = node->X()[0] + mydisp.at(0);
    currpos(1) = node->X()[1] + mydisp.at(1);
    currpos(2) = node->X()[2] + mydisp.at(2);

    currpositions_fluid_[node->Id()] = currpos;
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::ExtractInterfaceDispnp()
{
  return immersedstructure_->ExtractFullInterfaceDispnp();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::ApplyInterfaceForces(
    Teuchos::RCP<Epetra_Vector> full_traction_vec)
{
  immersedstructure_->ApplyImmersedInterfaceForces(
      FluidToStruct(MBFluidField()->ExtractInterfaceForces()),
      immersedstructure_->CombinedInterface()->ExtractOtherVector(*full_traction_vec));

  return;
}

void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::ExtractPreviousInterfaceSolution()
{
  FSI::Partitioned::ExtractPreviousInterfaceSolution();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::AddDirichCond()
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FluidAleImmersed>(MBFluidField())
      ->AddDirichCond(dbcmap_immersed_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::RemoveDirichCond()
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FluidAleImmersed>(MBFluidField())
      ->RemoveDirichCond(dbcmap_immersed_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::CalcResidual(Epetra_Vector& F,
    const Teuchos::RCP<Epetra_Vector> newstate, const Teuchos::RCP<Epetra_Vector> oldstate)
{
  // insert Immersed and FSI vector to combined vector (vector of whole interface)
  immersedstructure_->CombinedInterface()->InsertOtherVector(
      immersedstructure_->Interface()->ExtractIMMERSEDCondVector(newstate), combined_newstate_);
  immersedstructure_->CombinedInterface()->InsertCondVector(
      immersedstructure_->Interface()->ExtractFSICondVector(newstate), combined_newstate_);

  int err = F.Update(1.0, *combined_newstate_, -1.0, *oldstate, 0.0);

  double immersed_length =
      immersedstructure_->Interface()->ExtractIMMERSEDCondVector(newstate)->GlobalLength();
  double fsi_length =
      immersedstructure_->Interface()->ExtractFSICondVector(newstate)->GlobalLength();

  double immersed_norm = 0.0;
  immersedstructure_->CombinedInterface()->ExtractOtherVector(F)->Norm2(&immersed_norm);
  double fsi_norm = -1234.0;
  immersedstructure_->CombinedInterface()->ExtractCondVector(F)->Norm2(&fsi_norm);
  std::cout << "Immersed Residual = " << std::setprecision(14)
            << immersed_norm / sqrt(immersed_length) << " (length=" << immersed_length << ")"
            << std::endl;
  std::cout << "ALE-FSI Residual  = " << std::setprecision(14) << fsi_norm / sqrt(fsi_length)
            << " (length=" << fsi_length << ")" << std::endl;

  return err;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::SetDefaultParameters(
    const Teuchos::ParameterList& fsidyn, Teuchos::ParameterList& list)
{
  if (myrank_ == 0)
    std::cout << "\n SetDefaultParameters in ImmersedPartitionedFSIDirichletNeumannALE ..."
              << std::endl;

  // extract sublist with settings for partitioned solver
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = list;

  nlParams.set("Nonlinear Solver", "Line Search Based");
  nlParams.set("Preconditioner", "None");
  nlParams.set("Norm abs F", fsipart.get<double>("CONVTOL"));
  nlParams.set("Max Iterations", fsipart.get<int>("ITEMAX"));

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& lineSearchParams = nlParams.sublist("Line Search");

  //
  // Set parameters for NOX to chose the solver direction and line
  // search step.
  //
  switch (DRT::INPUT::IntegralValue<int>(fsidyn, "COUPALGO"))
  {
    case fsi_iter_stagg_AITKEN_rel_param:
    {
      std::cout << " Set AitkenFactoryImmersedAle ... " << std::endl;
      // fixed-point solver with Aitken relaxation parameter
      SetMethod("ITERATIVE STAGGERED SCHEME WITH RELAXATION PARAMETER VIA AITKEN ITERATION");

      nlParams.set("Jacobian", "None");

      dirParams.set("Method", "User Defined");
      Teuchos::RCP<NOX::Direction::UserDefinedFactory> fixpointfactory =
          Teuchos::rcp(new NOX::FSI::FixPointFactory());
      dirParams.set("User Defined Direction Factory", fixpointfactory);

      Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> aitkenfactory =
          Teuchos::rcp(new NOX::IMMERSED::AitkenFactoryImmersedAle());

      lineSearchParams.set("Method", "User Defined");
      lineSearchParams.set("User Defined Line Search Factory", aitkenfactory);

      lineSearchParams.sublist("Aitken").set("max step size", fsipart.get<double>("MAXOMEGA"));
      lineSearchParams.sublist("Aitken").set("min step size", -0.1);
      break;
    }
    default:
    {
      FSI::Partitioned::SetDefaultParameters(fsidyn, list);
    }
  }

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", Comm().MyPID());

  // set default output flag to no output
  // The field solver will output a lot, anyway.
  printParams.get("Output Information",
      ::NOX::Utils::Warning | ::NOX::Utils::OuterIteration | ::NOX::Utils::OuterIterationStatusTest
      // ::NOX::Utils::Parameters
  );

  Teuchos::ParameterList& solverOptions = nlParams.sublist("Solver Options");
  solverOptions.set<std::string>("Status Test Check Type", "Complete");

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::CalcFluidTractionsOnStructure()
{
  // immersed part
  IMMERSED::ImmersedPartitionedFSIDirichletNeumann::CalcFluidTractionsOnStructure();

  // ale-fsi part
  Teuchos::RCP<Epetra_Vector> ALEFSIForce = FluidToStruct(MBFluidField()->ExtractInterfaceForces());
  immersedstructure_->Interface()->AddFSICondVector(ALEFSIForce, struct_bdry_traction_);

  return;
}
