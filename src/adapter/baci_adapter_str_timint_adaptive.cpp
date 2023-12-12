/*----------------------------------------------------------------------*/
/*! \file

\brief Structure field adapter for time step size adaptivity

\level 2


*/

/*----------------------------------------------------------------------*/
/* headers */
#include "baci_adapter_str_timint_adaptive.H"

#include "baci_structure_timada.H"
#include "baci_structure_timint.H"
#include "baci_structure_timint_create.H"

BACI_NAMESPACE_OPEN


/*======================================================================*/
/* constructor */
ADAPTER::StructureTimIntAda::StructureTimIntAda(
    Teuchos::RCP<STR::TimAda> sta, Teuchos::RCP<Structure> sti)
    : StructureWrapper(sti), sta_(sta)
{
  // make sure
  if (sta_ == Teuchos::null) dserror("Failed to create structural integrator");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::StructureTimIntAda::Integrate() { return sta_->Integrate(); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimIntAda::PrepareOutput(bool force_prepare)
{
  sta_->PrepareOutputPeriod(force_prepare);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimIntAda::Output() { sta_->OutputPeriod(); }


/*----------------------------------------------------------------------*/

BACI_NAMESPACE_CLOSE
