/*----------------------------------------------------------------------------*/
/** \file

\brief Adaption for the xstructure case

\maintainer Matthias Mayr


\level 3

*/
/*----------------------------------------------------------------------------*/


#include "ad_str_xstructure_algorithm.H"

#include "../drt_structure_new/str_timint_factory.H"
#include "../drt_structure_new/str_timint_base.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::XStructureAlgorithm::SetGlobalState(
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& dataglobalstate,
    const Teuchos::RCP<const STR::TIMINT::BaseDataSDyn>& datasdyn)
{
  dataglobalstate = STR::TIMINT::BuildDataGlobalState();
  dataglobalstate->Init(actdis_, *sdyn_, datasdyn);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::XStructureAlgorithm::SetTimeIntegrationStrategy(
    Teuchos::RCP<STR::TIMINT::Base>& ti_strategy,
    const Teuchos::RCP<STR::TIMINT::BaseDataIO>& dataio,
    const Teuchos::RCP<STR::TIMINT::BaseDataSDyn>& datasdyn,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& dataglobalstate, const int& restart)
{
  ti_strategy = STR::TIMINT::BuildStrategy(*sdyn_);
  ti_strategy->Init(dataio, datasdyn, dataglobalstate);
}
