/*!----------------------------------------------------------------------
\file cardiovascular0d_dofset.cpp

\brief A set of degrees of freedom

\level 2

<pre>
\maintainer Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>

*----------------------------------------------------------------------*/

#include "../drt_cardiovascular0d/cardiovascular0d_dofset.H"

#include <iostream>
#include <algorithm>
#include <numeric>


/*----------------------------------------------------------------------*
 |  ctor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
UTILS::Cardiovascular0DDofSet::Cardiovascular0DDofSet()
: DRT::DofSet()
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
UTILS::Cardiovascular0DDofSet::~Cardiovascular0DDofSet()
{
  return;
}


/*----------------------------------------------------------------------*
 |  reset everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
void UTILS::Cardiovascular0DDofSet::Reset()
{
  dofrowmap_ = Teuchos::null;
  dofcolmap_ = Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  setup everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
int UTILS::Cardiovascular0DDofSet::AssignDegreesOfFreedom
(
    const Teuchos::RCP<DRT::Discretization> dis,
    const int ndofs,
    const int start,
    const Teuchos::RCP<UTILS::MOR> mor
)
{
  // A definite offset is currently not supported.
  if (start!=0)
    dserror("right now user specified dof offsets are not supported");

  // Add DofSets in order of assignment to list. Once it is there it has its
  // place and will get its starting id from the previous DofSet.
  AddDofSettoList();

  // We assume that all dof sets before this one have been set up. Otherwise
  // we'd have to reorder the list.
  //
  // There is no test anymore to make sure that all prior dof sets have been
  // assigned. It seems people like to manipulate dof sets. People do create
  // dof sets that do not contain any dofs (on its first assignment), people
  // even shift dof set numbers to create overlapping dof sets. This is
  // perfectly fine.
  //
  // However if you rely on non-overlapping dof sets, you have to
  // FillComplete() your discretizations in the order of their creation. This
  // is guaranteed for all discretizations read from the input file since the
  // input reader calls FillComplete(). If you create your own discretizations
  // try to understand what you do.

  // Get highest GID used so far and add one
  // (In case of POD-MOR, the highest GID will be projmatrix->NumVectors()-1 because indexing starts from 0. Therefore, there is no need to add anything.)
  int count;
  if (mor == Teuchos::null or not mor->HaveMOR())
    count = MaxGIDinList(dis->Comm()) + 1;
  else
    count = mor->GetRedDim();

  // dofrowmap with index base = count, which is undesired
  Teuchos::RCP<Epetra_Map> dofrowmap = Teuchos::rcp(new Epetra_Map(ndofs,count,dis->Comm()));

  std::vector<int> gids;
  for(int i=0;i<dofrowmap->NumMyElements();i++)
    gids.push_back(dofrowmap->GID(i));

  // dofrowmap with index base = 0
  dofrowmap_ = Teuchos::rcp(new Epetra_Map(-1,gids.size(),&gids[0],0,dis->Comm()));

  return count;
}


