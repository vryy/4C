/*!----------------------------------------------------------------------
\file mortar_strategy_base.cpp

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialComm.h"
#include "mortar_strategy_base.H"
#include "mortar_defines.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_mortar.H"


/*----------------------------------------------------------------------*
 | ctor (public)                                             popp 01/10 |
 *----------------------------------------------------------------------*/
MORTAR::StrategyBase::StrategyBase(
    const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap,
    Teuchos::ParameterList params,
    int dim,
    Teuchos::RCP<Epetra_Comm> comm,
    double alphaf,
    int maxdof) :
probdofs_(Teuchos::rcp(new Epetra_Map(*(DofRowMap)))),
probnodes_(Teuchos::rcp(new Epetra_Map(*(NodeRowMap)))),
comm_(comm),
scontact_(params),
dim_(dim),
alphaf_(alphaf),
parredist_(false),
maxdof_(maxdof)

{
  // empty constructor
  // (this is an abstract base class)
}
