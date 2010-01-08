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
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_SerialComm.h"
#include "mortar_strategy_base.H"
#include "mortar_defines.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_sparsematrix.H"
#include "../drt_lib/linalg_ana.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_lib/linalg_utils.H"

using namespace std;
using namespace Teuchos;

/*----------------------------------------------------------------------*
 | ctor (public)                                             popp 01/10 |
 *----------------------------------------------------------------------*/
MORTAR::StrategyBase::StrategyBase(RCP<Epetra_Map> problemrowmap,
                                   Teuchos::ParameterList params,
                                   int dim, RCP<Epetra_Comm> comm,
                                   double alphaf) :
problemrowmap_(problemrowmap),
scontact_(params),
dim_(dim),
comm_(comm),
alphaf_(alphaf)

{
  // empty constructor
  // (this is an abstract base class)
}

#endif // CCADISCRET
