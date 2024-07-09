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
      constraint_manager_(Adapter::ConstraintEnforcerFactory::create_enforcer(
          Global::Problem::instance()->fsi_dynamic_params(),
          Global::Problem::instance()->fbi_params()))
{
  // empty constructor
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannVel::setup()
{
  // call setup of base class
  FSI::DirichletNeumann::setup();
  const Teuchos::ParameterList& fsidyn = Global::Problem::instance()->fsi_dynamic_params();
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
  if (Core::UTILS::IntegralValue<int>(fsipart, "COUPVARIABLE") == Inpar::FSI::CoupVarPart::disp)
    FOUR_C_THROW("Please set the fsi coupling variable to Velocity or Force!\n");
  set_kinematic_coupling(
      Core::UTILS::IntegralValue<int>(fsipart, "COUPVARIABLE") == Inpar::FSI::CoupVarPart::vel);
  if (Teuchos::rcp_dynamic_cast<Adapter::FBIStructureWrapper>(structure_field(), true) ==
      Teuchos::null)
  {
    FOUR_C_THROW("Something went very wrong here! You should have a FBIStructureWrapper!\n");
  }
  if (Teuchos::rcp_dynamic_cast<Adapter::FBIFluidMB>(mb_fluid_field(), true) == Teuchos::null)
  {
    FOUR_C_THROW("Something went very wrong here! You should have a FBIFluidMB adapter!\n");
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannVel::fluid_op(
    Teuchos::RCP<Epetra_Vector> ivel, const FillType fillFlag)
{
  const Teuchos::ParameterList& fbi = Global::Problem::instance()->fbi_params();

  FSI::Partitioned::fluid_op(ivel, fillFlag);

  if (fillFlag == User)
  {
    FOUR_C_THROW("Not implemented!\n");
    return fluid_to_struct(mb_fluid_field()->relaxation_solve(Teuchos::null, dt()));
  }
  else
  {
    // A rather simple hack. We need something better!
    const int itemax = mb_fluid_field()->itemax();

    if (fillFlag == MF_Res and mfresitemax_ > 0) mb_fluid_field()->set_itemax(mfresitemax_ + 1);

    mb_fluid_field()->nonlinear_solve(Teuchos::null, Teuchos::null);

    mb_fluid_field()->set_itemax(itemax);

    if (Teuchos::getIntegralValue<Inpar::FBI::BeamToFluidCoupling>(fbi, "COUPLING") !=
            Inpar::FBI::BeamToFluidCoupling::fluid &&
        fbi.get<int>("STARTSTEP") < step())
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

  const Teuchos::ParameterList& fbi = Global::Problem::instance()->fbi_params();
  if (!(Teuchos::getIntegralValue<Inpar::FBI::BeamToFluidCoupling>(fbi, "COUPLING") ==
          Inpar::FBI::BeamToFluidCoupling::fluid) &&
      fbi.get<int>("STARTSTEP") < step())
  {
    if (not use_old_structure_)
      structure_field()->apply_interface_forces(iforce);
    else
      FOUR_C_THROW(
          "Fluid beam interaction is not tested with the old structure time. You should not be "
          "here! Change the IntStrategy in your Input file to Standard.\n");
  }

  structure_field()->solve();
  structure_field()->write_gmsh_struc_output_step();

  if (fbi.get<int>("STARTSTEP") < step())
  {
    constraint_manager_->prepare_fluid_solve();
    constraint_manager_->evaluate();
    if (!(Teuchos::getIntegralValue<Inpar::FBI::BeamToFluidCoupling>(fbi, "COUPLING") ==
            Inpar::FBI::BeamToFluidCoupling::solid))
      Teuchos::rcp_dynamic_cast<Adapter::FBIFluidMB>(mb_fluid_field(), true)
          ->reset_external_forces();
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
    if (Teuchos::rcp_dynamic_cast<Adapter::FBIStructureWrapper>(structure_field(), true) !=
        Teuchos::null)
      return Teuchos::rcp_dynamic_cast<Adapter::FBIStructureWrapper>(structure_field(), true)
          ->predict_interface_velnp();
  }
  else
  {
    const Teuchos::ParameterList& fsidyn = Global::Problem::instance()->fsi_dynamic_params();
    const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
    if (Core::UTILS::IntegralValue<int>(fsipart, "PREDICTOR") != 1)
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
  return constraint_manager_->fluid_to_structure();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannVel::struct_to_fluid(
    Teuchos::RCP<Epetra_Vector> iv)
{
  return constraint_manager_->structure_to_fluid(step());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void FSI::DirichletNeumannVel::output()
{
  FSI::DirichletNeumann::output();
  constraint_manager_->output(time(), step());
  visualization_output_writer_->write_output_runtime(constraint_manager_, step(), time());
  structure_field()->discretization()->writer()->clear_map_cache();
  mb_fluid_field()->discretization()->writer()->clear_map_cache();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void FSI::DirichletNeumannVel::timeloop(
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& interface)
{
  constraint_manager_->setup(structure_field(), mb_fluid_field());
  if (get_kinematic_coupling()) constraint_manager_->prepare_fluid_solve();
  visualization_output_writer_ =
      Teuchos::rcp(new BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter());
  visualization_output_writer_->init();
  visualization_output_writer_->setup(
      Core::IO::VisualizationParametersFactory(
          Global::Problem::instance()->io_params().sublist("RUNTIME VTK OUTPUT"),
          *Global::Problem::instance()->output_control_file(), time()),
      Teuchos::rcp_dynamic_cast<Adapter::FBIStructureWrapper>(structure_field(), true)
          ->get_io_data(),
      constraint_manager_->get_bridge()->get_params()->get_visualization_ouput_params_ptr());
  constraint_manager_->evaluate();
  if (get_kinematic_coupling()) struct_to_fluid(Teuchos::null);

  FSI::Partitioned::timeloop(interface);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void FSI::DirichletNeumannVel::set_binning(Teuchos::RCP<Core::Binstrategy::BinningStrategy> binning)
{
  constraint_manager_->set_binning(binning);
};

FOUR_C_NAMESPACE_CLOSE
