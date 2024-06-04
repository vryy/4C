/*----------------------------------------------------------------------*/
/*! \file
\file fsi_dirichletneumann_vel.cpp

\brief Solve FSI problems using a Dirichlet-Neumann partitioned approach


\level 3
*/
/*----------------------------------------------------------------------*/


#include "4C_fsi_dirichletneumann_vel.hpp"

#include "4C_adapter_fld_fbi_movingboundary.hpp"
#include "4C_adapter_str_fbiwrapper.hpp"
#include "4C_binstrategy.hpp"
#include "4C_fbi_adapter_constraintbridge.hpp"
#include "4C_fbi_beam_to_fluid_meshtying_output_writer.hpp"
#include "4C_fbi_beam_to_fluid_meshtying_params.hpp"
#include "4C_fbi_constraintenforcer.hpp"
#include "4C_fbi_constraintenforcer_factory.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fbi.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_io_control.hpp"
#include "4C_io_visualization_parameters.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <iostream>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::DirichletNeumannVel::DirichletNeumannVel(const Epetra_Comm& comm)
    : DirichletNeumann(comm),
      constraint_manager_(ADAPTER::ConstraintEnforcerFactory::CreateEnforcer(
          GLOBAL::Problem::Instance()->FSIDynamicParams(),
          GLOBAL::Problem::Instance()->FBIParams()))
{
  // empty constructor
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannVel::Setup()
{
  // call setup of base class
  FSI::DirichletNeumann::Setup();
  const Teuchos::ParameterList& fsidyn = GLOBAL::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
  if (CORE::UTILS::IntegralValue<int>(fsipart, "COUPVARIABLE") == INPAR::FSI::CoupVarPart::disp)
    FOUR_C_THROW("Please set the fsi coupling variable to Velocity or Force!\n");
  set_kinematic_coupling(
      CORE::UTILS::IntegralValue<int>(fsipart, "COUPVARIABLE") == INPAR::FSI::CoupVarPart::vel);
  if (Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(structure_field(), true) ==
      Teuchos::null)
  {
    FOUR_C_THROW("Something went very wrong here! You should have a FBIStructureWrapper!\n");
  }
  if (Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(MBFluidField(), true) == Teuchos::null)
  {
    FOUR_C_THROW("Something went very wrong here! You should have a FBIFluidMB adapter!\n");
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannVel::fluid_op(
    Teuchos::RCP<Epetra_Vector> ivel, const FillType fillFlag)
{
  const Teuchos::ParameterList& fbi = GLOBAL::Problem::Instance()->FBIParams();

  FSI::Partitioned::fluid_op(ivel, fillFlag);

  if (fillFlag == User)
  {
    FOUR_C_THROW("Not implemented!\n");
    return fluid_to_struct(MBFluidField()->RelaxationSolve(Teuchos::null, Dt()));
  }
  else
  {
    // A rather simple hack. We need something better!
    const int itemax = MBFluidField()->Itemax();

    if (fillFlag == MF_Res and mfresitemax_ > 0) MBFluidField()->SetItemax(mfresitemax_ + 1);

    MBFluidField()->nonlinear_solve(Teuchos::null, Teuchos::null);

    MBFluidField()->SetItemax(itemax);

    if (Teuchos::getIntegralValue<INPAR::FBI::BeamToFluidCoupling>(fbi, "COUPLING") !=
            INPAR::FBI::BeamToFluidCoupling::fluid &&
        fbi.get<int>("STARTSTEP") < Step())
    {
      constraint_manager_->recompute_coupling_without_pair_creation();
    }

    return fluid_to_struct(ivel);
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannVel::struct_op(
    Teuchos::RCP<Epetra_Vector> iforce, const FillType fillFlag)
{
  FSI::Partitioned::struct_op(iforce, fillFlag);

  const Teuchos::ParameterList& fbi = GLOBAL::Problem::Instance()->FBIParams();
  if (!(Teuchos::getIntegralValue<INPAR::FBI::BeamToFluidCoupling>(fbi, "COUPLING") ==
          INPAR::FBI::BeamToFluidCoupling::fluid) &&
      fbi.get<int>("STARTSTEP") < Step())
  {
    if (not use_old_structure_)
      structure_field()->apply_interface_forces(iforce);
    else
      FOUR_C_THROW(
          "Fluid beam interaction is not tested with the old structure time. You should not be "
          "here! Change the IntStrategy in your Input file to Standard.\n");
  }

  structure_field()->Solve();
  structure_field()->write_gmsh_struc_output_step();

  if (fbi.get<int>("STARTSTEP") < Step())
  {
    constraint_manager_->PrepareFluidSolve();
    constraint_manager_->Evaluate();
    if (!(Teuchos::getIntegralValue<INPAR::FBI::BeamToFluidCoupling>(fbi, "COUPLING") ==
            INPAR::FBI::BeamToFluidCoupling::solid))
      Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(MBFluidField(), true)->ResetExternalForces();
  }

  return struct_to_fluid(iforce);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannVel::initial_guess()
{
  if (get_kinematic_coupling())
  {
    // predict velocity
    // For now, we want to use no predictor. This function returns the current interface velocity.
    if (Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(structure_field(), true) !=
        Teuchos::null)
      return Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(structure_field(), true)
          ->predict_interface_velnp();
  }
  else
  {
    const Teuchos::ParameterList& fsidyn = GLOBAL::Problem::Instance()->FSIDynamicParams();
    const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
    if (CORE::UTILS::IntegralValue<int>(fsipart, "PREDICTOR") != 1)
    {
      FOUR_C_THROW(
          "unknown interface force predictor '%s'", fsipart.get<std::string>("PREDICTOR").c_str());
    }
    return interface_force();
  }
  FOUR_C_THROW("Something went very wrong here!\n");
  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannVel::fluid_to_struct(
    Teuchos::RCP<Epetra_Vector> iv)
{
  return constraint_manager_->FluidToStructure();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannVel::struct_to_fluid(
    Teuchos::RCP<Epetra_Vector> iv)
{
  return constraint_manager_->StructureToFluid(Step());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void FSI::DirichletNeumannVel::output()
{
  FSI::DirichletNeumann::output();
  constraint_manager_->Output(Time(), Step());
  visualization_output_writer_->write_output_runtime(constraint_manager_, Step(), Time());
  structure_field()->discretization()->Writer()->ClearMapCache();
  MBFluidField()->discretization()->Writer()->ClearMapCache();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void FSI::DirichletNeumannVel::Timeloop(
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& interface)
{
  constraint_manager_->Setup(structure_field(), MBFluidField());
  if (get_kinematic_coupling()) constraint_manager_->PrepareFluidSolve();
  visualization_output_writer_ =
      Teuchos::rcp(new BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter());
  visualization_output_writer_->Init();
  visualization_output_writer_->Setup(
      CORE::IO::VisualizationParametersFactory(
          GLOBAL::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"),
          *GLOBAL::Problem::Instance()->OutputControlFile(), Time()),
      Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(structure_field(), true)->GetIOData(),
      constraint_manager_->GetBridge()->GetParams()->get_visualization_ouput_params_ptr());
  constraint_manager_->Evaluate();
  if (get_kinematic_coupling()) struct_to_fluid(Teuchos::null);

  FSI::Partitioned::Timeloop(interface);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void FSI::DirichletNeumannVel::SetBinning(Teuchos::RCP<BINSTRATEGY::BinningStrategy> binning)
{
  constraint_manager_->SetBinning(binning);
};

FOUR_C_NAMESPACE_CLOSE
