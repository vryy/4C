/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for partitioned porous multiphase flow through elastic medium problems

   \level 3

 *----------------------------------------------------------------------*/


#include "baci_poromultiphase_partitioned.H"

#include "baci_lib_globalproblem.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
POROMULTIPHASE::PoroMultiPhasePartitioned::PoroMultiPhasePartitioned(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : PoroMultiPhaseBase(comm, globaltimeparams)
{
}

BACI_NAMESPACE_CLOSE
