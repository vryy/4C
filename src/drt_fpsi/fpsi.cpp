/*!----------------------------------------------------------------------

\level 2

\maintainer  Christoph Ager
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  rauch 12/12 |
 *----------------------------------------------------------------------*/
#include "fpsi.H"

#include "fpsi_utils.H"

// drt_lib includes
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

// POROELAST includes
#include "../drt_poroelast/poroelast_utils.H"

FPSI::FPSI_Base::FPSI_Base(const Epetra_Comm& comm, const Teuchos::ParameterList& fpsidynparams)
    : AlgorithmBase(comm, fpsidynparams)
{
  // nothing to do ... so far
}


/*----------------------------------------------------------------------*
 | redistribute the FPSI interface                           thon 11/14 |
 *----------------------------------------------------------------------*/
void FPSI::FPSI_Base::RedistributeInterface()
{
  DRT::Problem* problem = DRT::Problem::Instance();
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();
  Teuchos::RCP<FPSI::Utils> FPSI_UTILS = FPSI::Utils::Instance();

  if (comm.NumProc() >
      1)  // if we have more than one processor, we need to redistribute at the FPSI interface
  {
    Teuchos::RCP<std::map<int, int>> Fluid_PoroFluid_InterfaceMap =
        FPSI_UTILS->Get_Fluid_PoroFluid_InterfaceMap();
    Teuchos::RCP<std::map<int, int>> PoroFluid_Fluid_InterfaceMap =
        FPSI_UTILS->Get_PoroFluid_Fluid_InterfaceMap();

    FPSI_UTILS->RedistributeInterface(problem->GetDis("fluid"), problem->GetDis("porofluid"),
        "FPSICoupling", *PoroFluid_Fluid_InterfaceMap);
    FPSI_UTILS->RedistributeInterface(problem->GetDis("ale"), problem->GetDis("porofluid"),
        "FPSICoupling", *PoroFluid_Fluid_InterfaceMap);
    FPSI_UTILS->RedistributeInterface(problem->GetDis("porofluid"), problem->GetDis("fluid"),
        "FPSICoupling", *Fluid_PoroFluid_InterfaceMap);
    FPSI_UTILS->RedistributeInterface(problem->GetDis("structure"), problem->GetDis("fluid"),
        "FPSICoupling", *Fluid_PoroFluid_InterfaceMap);

    // Material pointers need to be reset after redistribution.
    POROELAST::UTILS::SetMaterialPointersMatchingGrid(
        problem->GetDis("structure"), problem->GetDis("porofluid"));
  }

  return;
}
