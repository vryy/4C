/*----------------------------------------------------------------------*/
/*!

\brief Structure field adapter for time step size adaptivity

\level 2

\maintainer Matthias Mayr

*/

/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_structure/strtimint_create.H"
#include "../drt_structure/strtimada.H"
#include "../drt_structure/strtimint.H"
#include "ad_str_timint_adaptive.H"


/*======================================================================*/
/* constructor */
ADAPTER::StructureTimIntAda::StructureTimIntAda(
    Teuchos::RCP<STR::TimAda> sta, Teuchos::RCP<Structure> sti)
    : StructureWrapper(sti), structure_(sta)
{
  // make sure
  if (structure_ == Teuchos::null) dserror("Failed to create structural integrator");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::StructureTimIntAda::Integrate() { return structure_->Integrate(); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimIntAda::PrepareOutput() { structure_->PrepareOutputPeriod(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimIntAda::Output() { structure_->OutputPeriod(); }


/*----------------------------------------------------------------------*/
