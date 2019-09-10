/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for SSI problems.

\maintainer Christoph Schmidt

\level 1
*/


#include "ad_str_ssiwrapper.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::SSIStructureWrapper::SSIStructureWrapper(Teuchos::RCP<Structure> structure)
    : StructureWrapper(structure)
{
}
