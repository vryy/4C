/*----------------------------------------------------------------------*/
/*! \file
\brief A set of degrees of freedom for constraint problems
\level 2

\maintainer Matthias Mayr

*----------------------------------------------------------------------*/

#include <iostream>
#include <algorithm>
#include <numeric>

#include "constraintdofset.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
UTILS::ConstraintDofSet::ConstraintDofSet() : DRT::DofSet() { return; }


/*----------------------------------------------------------------------*
 |  dtor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
UTILS::ConstraintDofSet::~ConstraintDofSet() { return; }


/*----------------------------------------------------------------------*
 |  reset everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
void UTILS::ConstraintDofSet::Reset()
{
  dofrowmap_ = Teuchos::null;
  dofcolmap_ = Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  setup everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
int UTILS::ConstraintDofSet::AssignDegreesOfFreedom(
    const Teuchos::RCP<DRT::Discretization> dis, const int ndofs, const int start)
{
  // A definite offset is currently not supported.
  if (start != 0) dserror("right now user specified dof offsets are not supported");

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
  const int count = MaxGIDinList(dis->Comm()) + 1;

  // dofrowmap with index base = count, which is undesired
  Teuchos::RCP<Epetra_Map> dofrowmap = Teuchos::rcp(new Epetra_Map(ndofs, count, dis->Comm()));

  std::vector<int> gids;
  for (int i = 0; i < dofrowmap->NumMyElements(); i++) gids.push_back(dofrowmap->GID(i));

  // dofrowmap with index base = 0
  dofrowmap_ = Teuchos::rcp(new Epetra_Map(-1, gids.size(), &gids[0], 0, dis->Comm()));

  return count;
}
