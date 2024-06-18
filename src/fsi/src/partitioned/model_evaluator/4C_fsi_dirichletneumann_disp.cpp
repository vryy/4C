/*----------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problems using a Dirichlet-Neumann partitioned approach
based on the interface displacements


\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_fsi_dirichletneumann_disp.hpp"

#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_fsi_debugwriter.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fsi.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::DirichletNeumannDisp::DirichletNeumannDisp(const Epetra_Comm& comm) : DirichletNeumann(comm)
{
  // empty constructor
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannDisp::setup()
{
  // call setup of base class
  FSI::DirichletNeumann::setup();
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
  set_kinematic_coupling(
      Core::UTILS::IntegralValue<int>(fsipart, "COUPVARIABLE") == Inpar::FSI::CoupVarPart::disp);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannDisp::fluid_op(
    Teuchos::RCP<Epetra_Vector> idisp, const FillType fillFlag)
{
  FSI::Partitioned::fluid_op(idisp, fillFlag);

  if (fillFlag == User)
  {
    // SD relaxation calculation
    return fluid_to_struct(MBFluidField()->RelaxationSolve(struct_to_fluid(idisp), Dt()));
  }
  else
  {
    // normal fluid solve

    // the displacement -> velocity conversion at the interface
    const Teuchos::RCP<Epetra_Vector> ivel = interface_velocity(idisp);

    // A rather simple hack. We need something better!
    const int itemax = MBFluidField()->Itemax();
    if (fillFlag == MF_Res and mfresitemax_ > 0) MBFluidField()->SetItemax(mfresitemax_ + 1);

    MBFluidField()->nonlinear_solve(struct_to_fluid(idisp), struct_to_fluid(ivel));

    MBFluidField()->SetItemax(itemax);

    return fluid_to_struct(MBFluidField()->extract_interface_forces());
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannDisp::struct_op(
    Teuchos::RCP<Epetra_Vector> iforce, const FillType fillFlag)
{
  FSI::Partitioned::struct_op(iforce, fillFlag);

  if (fillFlag == User)
  {
    // SD relaxation calculation
    return structure_field()->RelaxationSolve(iforce);
  }
  else
  {
    // normal structure solve
    if (not use_old_structure_)
      structure_field()->apply_interface_forces(iforce);
    else
      structure_field()->apply_interface_forces_temporary_deprecated(
          iforce);  // todo remove this line as soon as possible!
    structure_field()->Solve();
    structure_field()->write_gmsh_struc_output_step();
    return structure_field()->extract_interface_dispnp();
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannDisp::initial_guess()
{
  if (get_kinematic_coupling())
  {
    // predict displacement
    return structure_field()->predict_interface_dispnp();
  }
  else
  {
    const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
    const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
    if (Core::UTILS::IntegralValue<int>(fsipart, "PREDICTOR") != 1)
    {
      FOUR_C_THROW(
          "unknown interface force predictor '%s'", fsipart.get<std::string>("PREDICTOR").c_str());
    }
    return interface_force();
  }
}

FOUR_C_NAMESPACE_CLOSE
