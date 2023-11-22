/*---------------------------------------------------------------------*/
/*! \file

\brief Auxiliary functions of the main discretization class

\level 0


*/
/*---------------------------------------------------------------------*/

#include "baci_comm_exporter.H"
#include "baci_lib_discret.H"
#include "baci_lib_elementtype.H"
#include "baci_linalg_utils_sparse_algebra_manipulation.H"
#include "baci_linear_solver_method_parameters.H"
#include "baci_utils_exceptions.H"

/*----------------------------------------------------------------------*
 |  compute nullspace of system (public)                     mwgee 02/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ComputeNullSpaceIfNecessary(
    Teuchos::ParameterList& solveparams, bool recompute)
{
  // see whether we have a list for an iterative solver
  if (!solveparams.isSublist("Belos Parameters") || solveparams.isSublist("IFPACK Parameters"))
  {
    return;
  }

  int numdf = 1;  // default value for no. of degrees of freedom per node
  int dimns = 1;  // default value for no. of nullspace vectors
  int nv = 0;     // default value for no. of velocity dofs
  int np = 0;     // default value for no. of pressure dofs

  // downwinding needs nodal block information, compute it
  if (NumMyRowElements())
  {
    // We assume that all elements are of equal type
    DRT::Element* dwele = lRowElement(0);
    dwele->ElementType().NodalBlockInformation(dwele, numdf, dimns, nv, np);
  }

  // communicate data to procs without row element
  std::array<int, 4> ldata = {numdf, dimns, nv, np};
  std::array<int, 4> gdata = {0, 0, 0, 0};
  Comm().MaxAll(ldata.data(), gdata.data(), 4);
  numdf = gdata[0];
  dimns = gdata[1];
  nv = gdata[2];
  np = gdata[3];

  if (!(nv + np)) dserror("Cannot determine nodal block size");

  // store nv and np at unique location in solver parameter list
  solveparams.sublist("NodalBlockInformation").set("number of momentum dofs", nv);
  solveparams.sublist("NodalBlockInformation").set("number of constraint dofs", np);
  solveparams.sublist("NodalBlockInformation").set("number of dofs per node", numdf);
  solveparams.sublist("NodalBlockInformation").set("nullspace dimension", dimns);

  // adapt multigrid settings (if a multigrid preconditioner is used)
  // see whether we have a sublist indicating usage of Trilinos::ML or Trilinos::MueLu
  if (!solveparams.isSublist("ML Parameters") && !solveparams.isSublist("MueLu Parameters") &&
      !solveparams.isSublist("MueLu (Contact) Parameters") &&
      !solveparams.isSublist("MueLu (Fluid) Parameters") &&
      !solveparams.isSublist("MueLu (TSI) Parameters") &&
      !solveparams.isSublist("MueLu (BeamSolid) Parameters") &&
      !solveparams.isSublist("MueLu (FSI) Parameters"))
    return;
  Teuchos::ParameterList* mllist_ptr = nullptr;
  if (solveparams.isSublist("ML Parameters"))
    mllist_ptr = &(solveparams.sublist("ML Parameters"));
  else if (solveparams.isSublist("MueLu Parameters"))
    mllist_ptr = &(solveparams.sublist("MueLu Parameters"));
  else if (solveparams.isSublist("MueLu (Contact) Parameters"))
    mllist_ptr = &(solveparams.sublist("MueLu (Contact) Parameters"));
  else if (solveparams.isSublist("MueLu (Fluid) Parameters"))
    mllist_ptr = &(solveparams.sublist("MueLu (Fluid) Parameters"));
  else if (solveparams.isSublist("MueLu (TSI) Parameters"))
    mllist_ptr = &(solveparams);
  else if (solveparams.isSublist("MueLu (BeamSolid) Parameters"))
    mllist_ptr = &(solveparams);
  else if (solveparams.isSublist("MueLu (FSI) Parameters"))
    mllist_ptr = &(solveparams);
  else
    return;

  // see whether we have previously computed the nullspace
  // and recomputation is enforced
  Teuchos::ParameterList& mllist = *mllist_ptr;
  Teuchos::RCP<Epetra_MultiVector> ns =
      mllist.get<Teuchos::RCP<Epetra_MultiVector>>("nullspace", Teuchos::null);
  if (ns != Teuchos::null && !recompute) return;

  // no, we have not previously computed the nullspace
  // or want to recompute it anyway
  // -> compute nullspace
  // do the usual tests
  if (!Filled()) dserror("FillComplete was not called on discretization");
  if (!HaveDofs()) dserror("Discretization has no dofs assigned");

  // compute solver parameters and set them into list
  CORE::LINEAR_SOLVER::Parameters::ComputeSolverParameters(*this, mllist);
}

/*----------------------------------------------------------------------*
 |  SetState surrogate for node based vectors                  (public) |
 |                                                            gjb 06/09 |
 *----------------------------------------------------------------------*/
void DRT::Discretization::AddMultiVectorToParameterList(
    Teuchos::ParameterList& p, const std::string name, Teuchos::RCP<const Epetra_MultiVector> vec)
{
  // provide data in node-based multi-vector for usage on element level
  // -> export to column map is necessary for parallel evaluation
  // SetState cannot be used since this multi-vector is nodebased and not dofbased!
  if (vec != Teuchos::null)
  {
    const Epetra_Map* nodecolmap = NodeColMap();
    const int numcol = vec->NumVectors();

    // if it's already in column map just copy it
    // This is a rough test, but it might be ok at this place.
    if (vec->Map().PointSameAs(*nodecolmap))
    {
      // make a copy as in parallel such that no additional RCP points to the state vector
      Teuchos::RCP<Epetra_MultiVector> tmp =
          Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, numcol));
      tmp->Update(1.0, *vec, 0.0);
      p.set(name, tmp);
    }
    else  // if it's not in column map export and allocate
    {
      Teuchos::RCP<Epetra_MultiVector> tmp =
          Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, numcol));
      CORE::LINALG::Export(*vec, *tmp);
      p.set(name, tmp);
    }
  }
  else
    p.set(name, Teuchos::null);

  return;
}
