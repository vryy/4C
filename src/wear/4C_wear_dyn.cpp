/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for structure with ale problems.

\level 2


*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  headers                                                  mgit 04/11 |
 *----------------------------------------------------------------------*/
#include "4C_wear_dyn.hpp"

#include "4C_ale_utils_clonestrategy.hpp"
#include "4C_discretization_condition_utils.hpp"
#include "4C_discretization_fem_general_utils_createdis.hpp"
#include "4C_fsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_wear.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_wear_partitioned.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | entry point for structure ale in DRT                      mgit 04/11 |
 *----------------------------------------------------------------------*/
void wear_dyn_drt(int restart)
{
  // create a communicator
  const Epetra_Comm& comm = GLOBAL::Problem::Instance()->GetDis("structure")->Comm();

  // ***********************************************************
  // Setup the problem
  // ***********************************************************
  const Teuchos::ParameterList& sdyn = GLOBAL::Problem::Instance()->structural_dynamic_params();
  const Teuchos::ParameterList& wearpara = GLOBAL::Problem::Instance()->WearParams();

  // check if quasistatic analysis
  if (sdyn.get<std::string>("DYNAMICTYP") != "Statics")
  {
    std::cout << "WARNING: wear without dynamic effects!!!" << std::endl;
    // FOUR_C_THROW ("ERROR: Structure with ale only for quasistatic analysis so in new sti so
    // far.");
  }

  // access the structure discretization, make sure it is filled
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
  structdis = GLOBAL::Problem::Instance()->GetDis("structure");
  // set degrees of freedom in the discretization
  if (!structdis->Filled() or !structdis->HaveDofs()) structdis->fill_complete();

  // access the ale discretization
  Teuchos::RCP<DRT::Discretization> aledis = Teuchos::null;
  aledis = GLOBAL::Problem::Instance()->GetDis("ale");
  if (!aledis->Filled()) aledis->fill_complete();

  // we use the structure discretization as layout for the ale discretization
  if (structdis->NumGlobalNodes() == 0) FOUR_C_THROW("Structure discretization is empty!");

  // clone ale mesh from structure discretization
  if (aledis->NumGlobalNodes() == 0)
  {
    CORE::FE::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(
        structdis, aledis, GLOBAL::Problem::Instance()->CloningMaterialMap());
    aledis->fill_complete();
    // setup material in every ALE element
    Teuchos::ParameterList params;
    params.set<std::string>("action", "setup_material");
    aledis->Evaluate(params);
  }
  else  // filled ale discretization
  {
    // if we have non-matching meshes:
    if (!CORE::UTILS::IntegralValue<bool>(
            GLOBAL::Problem::Instance()->WearParams(), "MATCHINGGRID"))
    {
      // create vector of discr.
      std::vector<Teuchos::RCP<DRT::Discretization>> dis;
      dis.push_back(structdis);
      dis.push_back(aledis);

      CORE::REBALANCE::RebalanceDiscretizationsByBinning(dis, false);
    }
  }
  // ***********************************************************

  Teuchos::RCP<WEAR::Algorithm> stru_ale = Teuchos::null;

  // structure ale object
  if (CORE::UTILS::IntegralValue<INPAR::WEAR::WearCoupAlgo>(wearpara, "WEAR_COUPALGO") ==
      INPAR::WEAR::wear_stagg)
  {
    stru_ale = Teuchos::rcp(new WEAR::Partitioned(comm));
  }
  else if (CORE::UTILS::IntegralValue<INPAR::WEAR::WearCoupAlgo>(wearpara, "WEAR_COUPALGO") ==
           INPAR::WEAR::wear_iterstagg)
  {
    stru_ale = Teuchos::rcp(new WEAR::Partitioned(comm));
  }
  else
  {
    FOUR_C_THROW("Chosen algorithm not supported");
  }

  // read restart before joining the time loop
  if (restart != 0) stru_ale->read_restart(restart);

  // solve the whole problem
  stru_ale->TimeLoop();

  // summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform the result test
  GLOBAL::Problem::Instance()->AddFieldTest(stru_ale->structure_field()->CreateFieldTest());
  GLOBAL::Problem::Instance()->TestAll(comm);

  return;
}  // wear_dyn_drt()

/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
