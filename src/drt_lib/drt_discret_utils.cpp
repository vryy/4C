/*!----------------------------------------------------------------------
\file drt_discret_utils.cpp

\brief Auxiliary functions of the main discretization class

\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*----------------------------------------------------------------------*/

#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "drt_elementtype.H"

/*----------------------------------------------------------------------*
 |  compute nullspace of system (public)                     mwgee 02/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ComputeNullSpaceIfNecessary(
                                              Teuchos::ParameterList& solveparams,
                                              bool recompute)
{
  // see whether we have an aztec list
  if (!solveparams.isSublist("Aztec Parameters") &&
      !solveparams.isSublist("Belos Parameters") &&
      !solveparams.isSublist("Stratimikos Parameters")) return;

  int numdf = 1; // default value for no. of degrees of freedom per node
  int dimns = 1; // default value for no. of nullspace vectors
  int nv=0; // default value for no. of velocity dofs
  int np=0; // default value for no. of pressure dofs

  // downwinding needs nodal block information, compute it
  if (NumMyRowElements())
  {
    // We assume that all elements are of equal type
    DRT::Element* dwele = lRowElement(0);
    dwele->ElementType().NodalBlockInformation( dwele, numdf, dimns, nv, np );
  }

  // communicate data to procs without row element
  int ldata[4] = {numdf, dimns, nv, np};
  int gdata[4] = {0, 0, 0 ,0};
  Comm().MaxAll(&ldata[0], &gdata[0], 4);
  numdf = gdata[0];
  dimns = gdata[1];
  nv    = gdata[2];
  np    = gdata[3];

  if (!(nv+np)) dserror("Cannot determine nodal block size");

  // store nv and np at unique location in solver parameter list
  solveparams.sublist("NodalBlockInformation").set("nv",nv); // TODO improve names
  solveparams.sublist("NodalBlockInformation").set("np",np);
  solveparams.sublist("NodalBlockInformation").set("numdf",numdf);
  solveparams.sublist("NodalBlockInformation").set("dimns",dimns);

  if(solveparams.isSublist("Aztec Parameters"))
  {

    // get the aztec list and see whether we use downwinding
    Teuchos::ParameterList& azlist = solveparams.sublist("Aztec Parameters");

    azlist.set<int>("downwinding nv",nv);
    azlist.set<int>("downwinding np",np);
  }
  else if(solveparams.isSublist("Belos Parameters"))
  {
    // get the belos list and see whether we use downwinding
    Teuchos::ParameterList& beloslist = solveparams.sublist("Belos Parameters");

    beloslist.set<int>("downwinding nv",nv);
    beloslist.set<int>("downwinding np",np);
  }
  else if(solveparams.isSublist("Stratimikos Parameters"))
  {
    // no up and downwinding supported within Stratimikos...
  }
  else
  {
    dserror("no Aztec and no Belos list");
  }

  // adapt ML settings (if ML preconditioner is used)
  // see whether we have a sublist indicating usage of Trilinos::ML
  if (!solveparams.isSublist("ML Parameters") &&
      !solveparams.isSublist("MueLu Parameters") &&
      !solveparams.isSublist("MueLu (Contact) Parameters") &&
      !solveparams.isSublist("MueLu (Contact2) Parameters") &&
      !solveparams.isSublist("MueLu (Contact3) Parameters") &&
      !solveparams.isSublist("MueLu (PenaltyContact) Parameters") &&
      !solveparams.isSublist("Stratimikos Parameters")) return;
  Teuchos::ParameterList* mllist_ptr = NULL;
  if (solveparams.isSublist("Stratimikos Parameters"))
  {
    // TODO: what about MueLu?
    if (solveparams.sublist("Stratimikos Parameters").get<std::string>("Preconditioner Type") != "ML")
        return;
    else
      mllist_ptr = &(solveparams.sublist("Stratimikos Parameters").sublist("Preconditioner Types").sublist("ML").sublist("ML Settings"));
  }
  else if (solveparams.isSublist("ML Parameters"))
    mllist_ptr = &(solveparams.sublist("ML Parameters"));
  else if (solveparams.isSublist("MueLu Parameters"))
    mllist_ptr = &(solveparams.sublist("MueLu Parameters"));
  else if (solveparams.isSublist("MueLu (Contact) Parameters"))
    mllist_ptr = &(solveparams.sublist("MueLu (Contact) Parameters"));
  else if (solveparams.isSublist("MueLu (Contact2) Parameters"))
    mllist_ptr = &(solveparams.sublist("MueLu (Contact2) Parameters"));
  else if (solveparams.isSublist("MueLu (Contact3) Parameters"))
    mllist_ptr = &(solveparams.sublist("MueLu (Contact3) Parameters"));
  else if (solveparams.isSublist("MueLu (PenaltyContact) Parameters"))
    mllist_ptr = &(solveparams.sublist("MueLu (PenaltyContact) Parameters"));
  else return;

  // check whether all procs have at least one row element. ML cannot handle this!
  if (!NumMyRowElements()) dserror("Proc does not have any elements. ML is not working in this case, use ILU,...");

  // see whether we have previously computed the nullspace
  // and recomputation is enforced
  Teuchos::ParameterList& mllist = *mllist_ptr;
  Teuchos::RCP<std::vector<double> > ns = mllist.get<Teuchos::RCP<std::vector<double> > >("nullspace",Teuchos::null);
  if (ns != Teuchos::null && !recompute) return;

  // no, we have not previously computed the nullspace
  // or want to recompute it anyway
  // -> compute nullspace
  // do the usual tests
  if (!Filled()) dserror("FillComplete was not called on discretization");
  if (!HaveDofs()) dserror("Discretization has no dofs assigned");

  // compute nullspace and fill it into the ML parameter list
  ComputeNullSpaceML(mllist,numdf,dimns);
}

/*--------------------------------------------------------------------------*
 |  directly compute nullspace (for Krylov projection)   (public) nis Feb13 |
 *--------------------------------------------------------------------------*/
void DRT::Discretization::ComputeNullSpaceML(
    Teuchos::ParameterList& mllist,
    const int numdf,
    const int dimns
    )
{
  mllist.set("PDE equations",numdf);
  mllist.set("null space: dimension",dimns);
  mllist.set("null space: type","pre-computed");
  mllist.set("null space: add default vectors",false);

  // allocate dimns times the local length of the rowmap
  const int lrows = DofRowMap(0)->NumMyElements();
  Teuchos::RCP<std::vector<double> > ns = Teuchos::rcp(new std::vector<double>(dimns*lrows));
  double* nullsp = &((*ns)[0]);
  mllist.set<Teuchos::RCP<std::vector<double> > >("nullspace",ns);
  mllist.set("null space: vectors",nullsp);
  mllist.set<bool>("ML validate parameter list",false); // otherwise, ML would not tolerate the Teuchos::RCP pointer to the null space in its list

  // compute null space directly. that will call eletypes.
  ComputeNullSpace(ns, numdf, dimns);
}

/*--------------------------------------------------------------------------*
 |  directly compute nullspace (for Krylov projection)   (public) nis Feb13 |
 *--------------------------------------------------------------------------*/
void DRT::Discretization::ComputeNullSpace(
  Teuchos::RCP<std::vector<double> > nullspace,
  int numdf,
  int dimns
  )
{
  // get number of dofs to check or allocate nullspace storage
  const int lrows = DofRowMap(0)->NumMyElements();

  // if nullspace is null then do everything from scratch
  if (nullspace->size()==0)
  {
    int nv, np = 0; // dummy - not needed here

    // downwinding needs nodal block information, compute it
    if (NumMyRowElements())
    {
      // We assume that all elements are of equal type
      DRT::Element* dwele = lRowElement(0);
      dwele->ElementType().NodalBlockInformation( dwele, numdf, dimns, nv, np );
    }

    // communicate data to procs without row element
    int ldata[2] = {numdf, dimns};
    int gdata[2] = {0, 0};
    Comm().MaxAll(&ldata[0], &gdata[0], 2);
    numdf = gdata[0];
    dimns = gdata[1];

    // allocate storage for vector holding the nullspace data
    nullspace->resize(dimns*lrows,0.0);
  }
  // check if nullspace has correct size
  else if( ((ssize_t)(nullspace->size())!=(dimns*lrows) or numdf==0) )
  {
    dserror("nullspace does not have correct size or numdf or dimns are zero");
  }

  // check whether numdf and dimns are consistent among all procs
  const int numproc = Comm().NumProc();
  int sumnumdf;
  Comm().SumAll(&numdf,&sumnumdf,1);
  if (sumnumdf != numdf*numproc) dserror("numdf not consistent among procs");
  int sumdimns;
  Comm().SumAll(&dimns,&sumdimns,1);
  if (sumdimns != dimns*numproc) dserror("dimns not consistent among procs");

  // check if dimns is possible
  if (dimns>6) dserror("Nullspace size only upto 6 supported");

  // compute nullspace for simple case: vector of ones
  if (dimns==1 && numdf==1)
  {
    for (int i=0; i<lrows; ++i)
      (*nullspace)[i] = 1.0;
    return;
  }

  // for rotational degrees of freedom:
  // compute nodal center of the discretization
  double x0send[3] = {0.0,0.0,0.0};
  for (int i=0; i<NumMyRowNodes(); ++i)
    for (int j=0; j<3; ++j)
      x0send[j] += lRowNode(i)->X()[j];
  double x0[3];
  Comm().SumAll(x0send,x0,3);
  for (int i=0; i<3; ++i)
    x0[i] /= NumGlobalNodes();

  // let the elementtype compute the nullspace
  DRT::Element* dwele = lRowElement(0);
  dwele->ElementType().ComputeNullSpace( *this, *nullspace, x0, numdf, dimns );
}

/*----------------------------------------------------------------------*
 |  SetState surrogate for node based vectors                  (public) |
 |                                                            gjb 06/09 |
 *----------------------------------------------------------------------*/
void DRT::Discretization::AddMultiVectorToParameterList(Teuchos::ParameterList& p,
                                                        const std::string name,
                                                        Teuchos::RCP<const Epetra_MultiVector> vec)
{
  //provide data in node-based multi-vector for usage on element level
  // -> export to column map is necessary for parallel evaluation
  //SetState cannot be used since this multi-vector is nodebased and not dofbased!
  if (vec != Teuchos::null)
  {
    const Epetra_Map* nodecolmap = NodeColMap();
    const int numcol = vec->NumVectors();

    // if it's already in column map just copy it
    // This is a rough test, but it might be ok at this place.
    if (vec->Map().PointSameAs(*nodecolmap))
    {
      // make a copy as in parallel such that no additional RCP points to the state vector
      Teuchos::RCP<Epetra_MultiVector> tmp = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,numcol));
      tmp->Update(1.0, *vec, 0.0);
      p.set(name,tmp);
    }
    else // if it's not in column map export and allocate
    {
      Teuchos::RCP<Epetra_MultiVector> tmp = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,numcol));
      LINALG::Export(*vec,*tmp);
      p.set(name,tmp);
    }
  }
  else
    p.set(name,Teuchos::null);

  return;
}
