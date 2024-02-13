/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for SSI problems.


\level 1
*/


#include "baci_adapter_str_ssiwrapper.hpp"

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::SSIStructureWrapper::SSIStructureWrapper(Teuchos::RCP<Structure> structure)
    : StructureWrapper(structure)
{
}

BACI_NAMESPACE_CLOSE
