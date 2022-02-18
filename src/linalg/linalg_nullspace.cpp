/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of calculation routines related to Nullspaces

\level 0

*----------------------------------------------------------------------*/

#include "linalg_nullspace.H"
#include "linalg_utils_nullspace.H"
#include "linalg_utils_sparse_algebra_manipulation.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Nullspace::ComputeNullSpace(
    DRT::Discretization& dis, Teuchos::RCP<std::vector<double>> nullspace, int numdf, int dimns)
{
  // get number of dofs to check or allocate nullspace storage
  const int lrows = dis.DofRowMap(0)->NumMyElements();

  // if nullspace is null then do everything from scratch
  if (nullspace->size() == 0)
  {
    int nv, np = 0;  // dummy - not needed here

    // downwinding needs nodal block information, compute it
    if (dis.NumMyRowElements())
    {
      // We assume that all elements are of equal type
      DRT::Element* dwele = dis.lRowElement(0);
      dwele->ElementType().NodalBlockInformation(dwele, numdf, dimns, nv, np);
    }

    // communicate data to procs without row element
    int ldata[2] = {numdf, dimns};
    int gdata[2] = {0, 0};
    dis.Comm().MaxAll(&ldata[0], &gdata[0], 2);
    numdf = gdata[0];
    dimns = gdata[1];

    // allocate storage for vector holding the nullspace data
    nullspace->resize(dimns * lrows, 0.0);
  }
  // check if nullspace has correct size
  else if (((ssize_t)(nullspace->size()) != (dimns * lrows) or numdf == 0))
  {
    dserror("nullspace does not have correct size or numdf or dimns are zero");
  }

  // check whether numdf and dimns are consistent among all procs
  const int numproc = dis.Comm().NumProc();
  int sumnumdf;
  dis.Comm().SumAll(&numdf, &sumnumdf, 1);
  if (sumnumdf != numdf * numproc) dserror("numdf not consistent among procs");
  int sumdimns;
  dis.Comm().SumAll(&dimns, &sumdimns, 1);
  if (sumdimns != dimns * numproc) dserror("dimns not consistent among procs");

  // check if dimns is possible
  if (dimns > 10) dserror("Nullspace size only upto 10 supported");

  // compute nullspace for simple case: vector of ones
  if (dimns == 1 && numdf == 1)
  {
    for (int i = 0; i < lrows; ++i) (*nullspace)[i] = 1.0;
    return;
  }

  // for rigid body rotations
  // compute nodal center of the discretization
  double x0send[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < dis.NumMyRowNodes(); ++i)
    for (int j = 0; j < 3; ++j) x0send[j] += dis.lRowNode(i)->X()[j];
  double x0[3];
  dis.Comm().SumAll(x0send, x0, 3);
  for (int i = 0; i < 3; ++i) x0[i] /= dis.NumGlobalNodes();

  // let the elementtype compute the nullspace
  DRT::Element* dwele = dis.lRowElement(0);
  dwele->ElementType().ComputeNullSpace(dis, *nullspace, x0, numdf, dimns);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Nullspace::FixNullSpace(std::string field, const Epetra_Map& oldmap,
    const Epetra_Map& newmap, Teuchos::ParameterList& solveparams)
{
  // there is no ML list, do nothing
  if (!solveparams.isSublist("ML Parameters") && !solveparams.isSublist("MueLu Parameters")) return;

  // find the ML or MueLu list
  Teuchos::RCP<Teuchos::ParameterList> params_ptr = Teuchos::null;
  if (solveparams.isSublist("ML Parameters"))
    params_ptr = Teuchos::rcp(&(solveparams.sublist("ML Parameters")), false);
  else if (solveparams.isSublist("MueLu Parameters"))
    params_ptr = Teuchos::rcp(&(solveparams.sublist("MueLu Parameters")), false);
  else
    return;
  Teuchos::ParameterList& params = *params_ptr;

  const int ndim = params.get("null space: dimension", -1);
  if (ndim == -1) dserror("List does not contain nullspace dimension");

  Teuchos::RCP<Epetra_MultiVector> nullspace =
      params.get<Teuchos::RCP<Epetra_MultiVector>>("nullspace", Teuchos::null);
  if (nullspace == Teuchos::null) dserror("List does not contain nullspace");

  double* ons = nullspace->Values();

  Teuchos::RCP<std::vector<double>> ns =
      Teuchos::rcp(new std::vector<double>(nullspace->MyLength() * nullspace->NumVectors()));
  LINALG::EpetraMultiVectorToStdVector(nullspace, ns, ndim);

  const int olength = (int)ns->size() / ndim;

  const int nlength = newmap.NumMyElements();

  if (olength == nlength) return;  // everything should be ok, do nothing

  if (olength != oldmap.NumMyElements())
    dserror("Nullspace does not match old map length, %d %d", olength, oldmap.NumMyElements());

  if (nlength > olength)
    dserror("New problem size larger than old - full rebuild of nullspace neccessary");

  // Allocate a new nullspace and fill it
  Teuchos::RCP<std::vector<double>> nsnew =
      Teuchos::rcp(new std::vector<double>(nlength * ndim, 0.0));
  double* nns = &((*nsnew)[0]);

  for (int i = 0; i < nlength; ++i)
  {
    int gid = newmap.GID(i);
    int olid = oldmap.LID(gid);
    if (olid == -1) continue;

    // transfer entries for this dof to new nullspace vector
    for (int j = 0; j < ndim; ++j) nns[j * nlength + i] = ons[j * olength + olid];
  }

  // put new nullspace in parameter list
  // this print message can go away at some point
  if (!oldmap.Comm().MyPID()) printf("Fixing %s ML Nullspace\n", field.c_str());

  Teuchos::RCP<Epetra_MultiVector> nullspacenew =
      Teuchos::rcp(new Epetra_MultiVector(newmap, ndim, true));
  LINALG::StdVectorToEpetraMultiVector(nsnew, nullspacenew, ndim);

  params.set<Teuchos::RCP<Epetra_MultiVector>>("nullspace", nullspacenew);
  params.set("null space: vectors", nullspacenew->Values());

  return;
}
