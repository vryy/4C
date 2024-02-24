/*---------------------------------------------------------------------*/
/*! \file
\brief Control routine for structural dynamics (outsourced to adapter layer)

\level 1

*/
/*---------------------------------------------------------------------*/

#include "baci_structure_dyn_nln_drt.hpp"

#include "baci_adapter_str_factory.hpp"
#include "baci_adapter_str_structure.hpp"
#include "baci_adapter_str_structure_new.hpp"
#include "baci_comm_utils.hpp"
#include "baci_global_data.hpp"
#include "baci_inpar_structure.hpp"
#include "baci_io.hpp"
#include "baci_io_control.hpp"
#include "baci_lib_discret.hpp"
#include "baci_lib_periodicbc.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"
#include "baci_linear_solver_method_linalg.hpp"
#include "baci_structure_resulttest.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void caldyn_drt()
{
  // get input lists
  const Teuchos::ParameterList& sdyn = GLOBAL::Problem::Instance()->StructuralDynamicParams();
  // major switch to different time integrators
  switch (CORE::UTILS::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP"))
  {
    case INPAR::STR::dyna_statics:
    case INPAR::STR::dyna_genalpha:
    case INPAR::STR::dyna_genalpha_liegroup:
    case INPAR::STR::dyna_onesteptheta:
    case INPAR::STR::dyna_gemm:
    case INPAR::STR::dyna_expleuler:
    case INPAR::STR::dyna_centrdiff:
    case INPAR::STR::dyna_ab2:
    case INPAR::STR::dyna_ab4:
    case INPAR::STR::dyna_euma:
    case INPAR::STR::dyna_euimsto:
      dyn_nlnstructural_drt();
      break;
    default:
      dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
      break;
  }

  return;
}


/*----------------------------------------------------------------------*
 | structural nonlinear dynamics                                        |
 *----------------------------------------------------------------------*/
void dyn_nlnstructural_drt()
{
  // get input lists
  const Teuchos::ParameterList& sdyn = GLOBAL::Problem::Instance()->StructuralDynamicParams();
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = GLOBAL::Problem::Instance()->GetDis("structure");

  // connect degrees of freedom for periodic boundary conditions
  {
    PeriodicBoundaryConditions pbc_struct(structdis);

    if (pbc_struct.HasPBC())
    {
      pbc_struct.UpdateDofsForPeriodicBoundaryConditions();
    }
  }

  // create an adapterbase and adapter
  Teuchos::RCP<ADAPTER::Structure> structadapter = Teuchos::null;
  // FixMe The following switch is just a temporal hack, such we can jump between the new and the
  // old structure implementation. Has to be deleted after the clean-up has been finished!
  const enum INPAR::STR::IntegrationStrategy intstrat =
      CORE::UTILS::IntegralValue<INPAR::STR::IntegrationStrategy>(sdyn, "INT_STRATEGY");
  switch (intstrat)
  {
    // -------------------------------------------------------------------
    // old implementation
    // -------------------------------------------------------------------
    case INPAR::STR::int_old:
    {
      Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> adapterbase_old_ptr =
          Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(
              sdyn, const_cast<Teuchos::ParameterList&>(sdyn), structdis));
      structadapter = adapterbase_old_ptr->StructureField();
      structadapter->Setup();
      break;
    }
    // -------------------------------------------------------------------
    // new implementation
    // -------------------------------------------------------------------
    default:
    {
      Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> adapterbase_ptr =
          ADAPTER::BuildStructureAlgorithm(sdyn);
      adapterbase_ptr->Init(sdyn, const_cast<Teuchos::ParameterList&>(sdyn), structdis);
      adapterbase_ptr->Setup();
      structadapter = adapterbase_ptr->StructureField();
      break;
    }
  }

  const bool write_initial_state = CORE::UTILS::IntegralValue<int>(
      GLOBAL::Problem::Instance()->IOParams(), "WRITE_INITIAL_STATE");
  const bool write_final_state =
      CORE::UTILS::IntegralValue<int>(GLOBAL::Problem::Instance()->IOParams(), "WRITE_FINAL_STATE");

  // do restart
  const int restart = GLOBAL::Problem::Instance()->Restart();
  if (restart)
  {
    structadapter->ReadRestart(restart);
  }
  // write output at beginnning of calc
  else
  {
    if (write_initial_state)
    {
      constexpr bool force_prepare = true;
      structadapter->PrepareOutput(force_prepare);
      structadapter->Output();
      structadapter->PostOutput();
    }
  }

  // run time integration
  structadapter->Integrate();

  if (write_final_state && !structadapter->HasFinalStateBeenWritten())
  {
    constexpr bool forceWriteRestart = true;
    constexpr bool force_prepare = true;
    structadapter->PrepareOutput(force_prepare);
    structadapter->Output(forceWriteRestart);
    structadapter->PostOutput();
  }

  // test results
  GLOBAL::Problem::Instance()->AddFieldTest(structadapter->CreateFieldTest());
  GLOBAL::Problem::Instance()->TestAll(structadapter->DofRowMap()->Comm());

  // print monitoring of time consumption
  Teuchos::RCP<const Teuchos::Comm<int>> TeuchosComm =
      CORE::COMM::toTeuchosComm<int>(structdis->Comm());
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, true);

  // time to go home...
  return;

}  // end of dyn_nlnstructural_drt()

BACI_NAMESPACE_CLOSE
