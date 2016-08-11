/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_partitioned.cpp

 \brief base class for partitioned porous multiphase flow through elastic medium problems

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/


#include "poromultiphase_partitioned.H"

#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
POROMULTIPHASE::PoroMultiPhasePartitioned::PoroMultiPhasePartitioned(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams):
    PoroMultiPhaseBase(comm, globaltimeparams)
{

}


