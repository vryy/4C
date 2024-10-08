/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for xfem-ale-fluids with moving boundaries

\level 2

*/
/*----------------------------------------------------------------------*/
#include "4C_adapter_fld_fluid_ale_xfem.hpp"

#include "4C_adapter_ale_fluid.hpp"
#include "4C_adapter_fld_fluid_xfsi.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_inpar_validparameters.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::FluidAleXFEM::FluidAleXFEM(const Teuchos::ParameterList& prbdyn, std::string condname)
    : FluidAle(prbdyn, condname)
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::FE::Discretization> Adapter::FluidAleXFEM::boundary_discretization()
{
  // returns the boundary discretization
  // REMARK:
  // the returned discretization has to match the structure discretization at the interface coupling
  // (see FSI::Partitioned::Partitioned(const Epetra_Comm& comm) ) therefore return the boundary dis
  // this is similar to the matching of fluid dis and ale dis in case of Adapter::FluidALE

  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(fluid_field(), true);

  return xfluid->boundary_discretization();
}


/*----------------------------------------------------------------------*/
/// communication object at the struct interface
/*----------------------------------------------------------------------*/
Teuchos::RCP<FLD::UTILS::MapExtractor> const& Adapter::FluidAleXFEM::struct_interface()
{
  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(fluid_field(), true);

  return xfluid->struct_interface();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidAleXFEM::nonlinear_solve(Teuchos::RCP<Core::LinAlg::Vector<double>> idisp,
    Teuchos::RCP<Core::LinAlg::Vector<double>> ivel)
{
  // if we have values at the interface we need to apply them

  // REMARK: for XFLUID idisp = Teuchos::null, ivel = Teuchos::null (called by fsi_fluid_xfem with
  // default Teuchos::null)
  //         for XFSI   idisp != Teuchos::null

  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(fluid_field(), true);

  // set idispnp in Xfluid
  if (idisp != Teuchos::null) xfluid->apply_struct_mesh_displacement(idisp);

  // set ivelnp in Xfluid
  if (ivel != Teuchos::null) xfluid->apply_struct_interface_velocities(ivel);

  // Update the ale update part
  if (fluid_field()->interface()->au_cond_relevant())
  {
    Teuchos::RCP<const Core::LinAlg::Vector<double>> dispnp = fluid_field()->dispnp();
    Teuchos::RCP<Core::LinAlg::Vector<double>> audispnp =
        fluid_field()->interface()->extract_au_cond_vector(*dispnp);
    ale_field()->apply_ale_update_displacements(aucoupfa_->master_to_slave(audispnp));
  }

  // Note: We do not look for moving ale boundaries (outside the coupling
  // interface) on the fluid side. Thus if you prescribe time variable ale
  // Dirichlet conditions the according fluid Dirichlet conditions will not
  // notice.


  ale_field()->solve();
  Teuchos::RCP<Core::LinAlg::Vector<double>> fluiddisp = ale_to_fluid_field(ale_field()->dispnp());
  fluid_field()->apply_mesh_displacement(fluiddisp);
  fluid_field()->solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Vector<double>> Adapter::FluidAleXFEM::relaxation_solve(
    Teuchos::RCP<Core::LinAlg::Vector<double>> idisp, double dt)
{
  FOUR_C_THROW("RelaxationSolve for XFEM useful?");
  std::cout << "WARNING: RelaxationSolve for XFEM useful?" << std::endl;

  // the displacement -> velocity conversion at the interface
  idisp->Scale(1. / dt);

  return fluid_field()->relaxation_solve(idisp);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Vector<double>> Adapter::FluidAleXFEM::extract_interface_forces()
{
  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(fluid_field(), true);
  return xfluid->extract_struct_interface_forces();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Vector<double>> Adapter::FluidAleXFEM::extract_interface_velnp()
{
  FOUR_C_THROW("Robin stuff");
  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(fluid_field(), true);
  return xfluid->extract_struct_interface_velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Vector<double>> Adapter::FluidAleXFEM::extract_interface_veln()
{
  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(fluid_field(), true);
  return xfluid->extract_struct_interface_veln();
}

FOUR_C_NAMESPACE_CLOSE
