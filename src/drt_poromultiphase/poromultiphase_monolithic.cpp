/*----------------------------------------------------------------------*/
/*!
 \brief base class for monolithic porous multiphase flow through elastic medium problems

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/


#include "poromultiphase_monolithic.H"

#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
POROMULTIPHASE::PoroMultiPhaseMonolithic::PoroMultiPhaseMonolithic(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : PoroMultiPhaseBase(comm, globaltimeparams)
{
}
