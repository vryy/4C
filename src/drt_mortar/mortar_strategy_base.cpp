/*!----------------------------------------------------------------------
\file mortar_abstract_strategy.cpp

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

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
#include "../linalg/linalg_ana.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../linalg/linalg_utils.H"


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

