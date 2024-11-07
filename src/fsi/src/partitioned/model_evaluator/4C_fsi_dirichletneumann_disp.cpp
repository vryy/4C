// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

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
  const Teuchos::ParameterList& fsidyn = Global::Problem::instance()->fsi_dynamic_params();
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
  set_kinematic_coupling(Teuchos::getIntegralValue<Inpar::FSI::CoupVarPart>(
                             fsipart, "COUPVARIABLE") == Inpar::FSI::CoupVarPart::disp);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> FSI::DirichletNeumannDisp::fluid_op(
    std::shared_ptr<Core::LinAlg::Vector<double>> idisp, const FillType fillFlag)
{
  FSI::Partitioned::fluid_op(idisp, fillFlag);

  if (fillFlag == User)
  {
    // SD relaxation calculation
    return fluid_to_struct(mb_fluid_field()->relaxation_solve(struct_to_fluid(idisp), dt()));
  }
  else
  {
    // normal fluid solve

    // the displacement -> velocity conversion at the interface
    const std::shared_ptr<Core::LinAlg::Vector<double>> ivel = interface_velocity(*idisp);

    // A rather simple hack. We need something better!
    const int itemax = mb_fluid_field()->itemax();
    if (fillFlag == MF_Res and mfresitemax_ > 0) mb_fluid_field()->set_itemax(mfresitemax_ + 1);

    mb_fluid_field()->nonlinear_solve(struct_to_fluid(idisp), struct_to_fluid(ivel));

    mb_fluid_field()->set_itemax(itemax);

    return fluid_to_struct(mb_fluid_field()->extract_interface_forces());
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> FSI::DirichletNeumannDisp::struct_op(
    std::shared_ptr<Core::LinAlg::Vector<double>> iforce, const FillType fillFlag)
{
  FSI::Partitioned::struct_op(iforce, fillFlag);

  if (fillFlag == User)
  {
    // SD relaxation calculation
    return structure_field()->relaxation_solve(iforce);
  }
  else
  {
    // normal structure solve
    if (not use_old_structure_)
      structure_field()->apply_interface_forces(iforce);
    else
      structure_field()->apply_interface_forces_temporary_deprecated(
          iforce);  // todo remove this line as soon as possible!
    structure_field()->solve();
    structure_field()->write_gmsh_struc_output_step();
    return structure_field()->extract_interface_dispnp();
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> FSI::DirichletNeumannDisp::initial_guess()
{
  if (get_kinematic_coupling())
  {
    // predict displacement
    return structure_field()->predict_interface_dispnp();
  }
  else
  {
    const Teuchos::ParameterList& fsidyn = Global::Problem::instance()->fsi_dynamic_params();
    const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
    if (fsipart.get<std::string>("PREDICTOR") != "d(n)")
    {
      FOUR_C_THROW(
          "unknown interface force predictor '%s'", fsipart.get<std::string>("PREDICTOR").c_str());
    }
    return interface_force();
  }
}

FOUR_C_NAMESPACE_CLOSE
