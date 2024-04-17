/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for monolithic porous multiphase flow through elastic medium problems

   \level 3

 *----------------------------------------------------------------------*/


#include "baci_poromultiphase_monolithic.hpp"

#include "baci_global_data.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
POROMULTIPHASE::PoroMultiPhaseMonolithic::PoroMultiPhaseMonolithic(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : PoroMultiPhaseBase(comm, globaltimeparams)
{
}

FOUR_C_NAMESPACE_CLOSE
