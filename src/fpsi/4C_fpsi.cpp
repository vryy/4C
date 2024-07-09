/*----------------------------------------------------------------------*/
/*! \file

\level 2

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  rauch 12/12 |
 *----------------------------------------------------------------------*/
#include "4C_fpsi.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fpsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_poroelast_utils.hpp"

FOUR_C_NAMESPACE_OPEN

FPSI::FpsiBase::FpsiBase(const Epetra_Comm& comm, const Teuchos::ParameterList& fpsidynparams)
    : AlgorithmBase(comm, fpsidynparams)
{
  // nothing to do ... so far
}


/*----------------------------------------------------------------------*
 | redistribute the FPSI interface                           thon 11/14 |
 *----------------------------------------------------------------------*/
void FPSI::FpsiBase::redistribute_interface()
{
  Global::Problem* problem = Global::Problem::instance();
  const Epetra_Comm& comm = problem->get_dis("structure")->get_comm();
  Teuchos::RCP<FPSI::Utils> FPSI_UTILS = FPSI::Utils::instance();

  if (comm.NumProc() >
      1)  // if we have more than one processor, we need to redistribute at the FPSI interface
  {
    Teuchos::RCP<std::map<int, int>> Fluid_PoroFluid_InterfaceMap =
        FPSI_UTILS->get_fluid_poro_fluid_interface_map();
    Teuchos::RCP<std::map<int, int>> PoroFluid_Fluid_InterfaceMap =
        FPSI_UTILS->get_poro_fluid_fluid_interface_map();

    FPSI_UTILS->redistribute_interface(problem->get_dis("fluid"), problem->get_dis("porofluid"),
        "fpsi_coupling", *PoroFluid_Fluid_InterfaceMap);
    FPSI_UTILS->redistribute_interface(problem->get_dis("ale"), problem->get_dis("porofluid"),
        "fpsi_coupling", *PoroFluid_Fluid_InterfaceMap);
    FPSI_UTILS->redistribute_interface(problem->get_dis("porofluid"), problem->get_dis("fluid"),
        "fpsi_coupling", *Fluid_PoroFluid_InterfaceMap);
    FPSI_UTILS->redistribute_interface(problem->get_dis("structure"), problem->get_dis("fluid"),
        "fpsi_coupling", *Fluid_PoroFluid_InterfaceMap);

    // Material pointers need to be reset after redistribution.
    PoroElast::UTILS::SetMaterialPointersMatchingGrid(
        problem->get_dis("structure"), problem->get_dis("porofluid"));
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
