/*!----------------------------------------------------------------------
\file mortar_abstract_strategy.cpp

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
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_inpar/inpar_mortar.H"


/*----------------------------------------------------------------------*
 | ctor (public)                                             popp 01/10 |
 *----------------------------------------------------------------------*/
MORTAR::StrategyBase::StrategyBase(DRT::Discretization& probdiscret,
                                   Teuchos::ParameterList params,
                                   int dim, Teuchos::RCP<Epetra_Comm> comm,
                                   double alphaf, int maxdof) :
probdiscret_(probdiscret),
comm_(comm),
scontact_(params),
dim_(dim),
alphaf_(alphaf),
maxdof_(maxdof)

{
  // create and store Epetra_Maps for problem dof, node, ele row maps
  probdofs_ =  Teuchos::rcp(new Epetra_Map(*(probdiscret.DofRowMap())));
  probnodes_ = Teuchos::rcp(new Epetra_Map(*(probdiscret.NodeRowMap())));
  probeles_ =  Teuchos::rcp(new Epetra_Map(*(probdiscret.ElementRowMap())));

  // empty constructor
  // (this is an abstract base class)
}

