/*!----------------------------------------------------------------------
\file drt_dirichletextractor.cpp
\brief Implementation

<pre>
\brief Implementation
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "drt_dirichletextractor.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::DirichletExtractor::Setup(const DRT::Discretization& dis)
{
  DRT::UTILS::MultiConditionSelector mcs;
  mcs.AddSelector( Teuchos::rcp( new DRT::UTILS::DirichletSelector( dis ) ) );
  mcs.SetupExtractor(dis,*dis.DofRowMap(),*this);
}

void DRT::DirichletExtractor::ZeroDirichlets( Teuchos::RCP<Epetra_Vector> residual ) const
{
  DirichletPutScalar( *residual, 0.0 );
}
