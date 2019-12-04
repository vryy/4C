/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of algebraic creation methods for namespace LINALG

\level 0
\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/

#include "../headers/compiler_definitions.h" /* access to fortran routines */
#include "linalg_utils_sparse_algebra_create.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 |  create a Epetra_CrsMatrix                                mwgee 12/06|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix> LINALG::CreateMatrix(const Epetra_Map& rowmap, const int npr)
{
  if (!rowmap.UniqueGIDs()) dserror("Row map is not unique");
  return Teuchos::rcp(new Epetra_CrsMatrix(::Copy, rowmap, npr, false));
}
/*----------------------------------------------------------------------*
 |  create a Epetra_Vector                                   mwgee 12/06|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> LINALG::CreateVector(const Epetra_BlockMap& rowmap, const bool init)
{
  return Teuchos::rcp(new Epetra_Vector(rowmap, init));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> LINALG::CreateMultiVector(
    const Epetra_BlockMap& rowmap, const int numrows, const bool init)
{
  return Teuchos::rcp(new Epetra_MultiVector(rowmap, numrows, init));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> LINALG::CreateMap(const std::set<int>& gids, const Epetra_Comm& comm)
{
  std::vector<int> mapvec;
  mapvec.reserve(gids.size());
  mapvec.assign(gids.begin(), gids.end());
  Teuchos::RCP<Epetra_Map> map =
      Teuchos::rcp(new Epetra_Map(-1, mapvec.size(), &mapvec[0], 0, comm));
  mapvec.clear();
  return map;
}

/*----------------------------------------------------------------------*
 | create epetra_map with out-of-bound check                 farah 06/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> LINALG::CreateMap(const std::vector<int>& gids, const Epetra_Comm& comm)
{
  Teuchos::RCP<Epetra_Map> map;

  if ((int)gids.size() > 0)
    map = Teuchos::rcp(new Epetra_Map(-1, gids.size(), &gids[0], 0, comm));
  else
    map = Teuchos::rcp(new Epetra_Map(-1, gids.size(), 0, 0, comm));

  return map;
}
