/*!----------------------------------------------------------------------
\file simpler_operator.cpp

\class LINALG::SIMPLER_Operator

\brief An approximate block factorization preconditioner based on the 
       SIMPLE family of methods

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#include "simpler_operator.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 02/08|
 *----------------------------------------------------------------------*/
LINALG::SIMPLER_Operator::SIMPLER_Operator(RCP<Epetra_CrsMatrix> A,
                                           const ParameterList& velocitylist,
                                           const ParameterList& pressurelist) 
  : vlist_(velocitylist),
    plist_(pressurelist)
{
  // remove the SIMPLER sublist from the vlist_, 
  // otherwise it will try to recursively create a SIMPLE
  // preconditioner when we do the block solvers
  if (vlist_.isSublist("SIMPLER")) vlist_.remove("SIMPLER");
  
  Setup(A);
  
  return;
}


/*----------------------------------------------------------------------*
 |  apply const operator (public)                            mwgee 02/08|
 *----------------------------------------------------------------------*/
int LINALG::SIMPLER_Operator::ApplyInverse(const Epetra_MultiVector& X, 
                                           Epetra_MultiVector& Y) const
{
  // do not do anything for testing
  Y.Update(1.0,X,0.0);
  return 0;

  // wrap incoming and outgoing vectors as MLAPI::MultiVector
  // note: Aztec might pass X and Y as physically identical objects, 
  // so we deep copy here
  
  
  return 0;
} 

/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 02/08|
 *----------------------------------------------------------------------*/
void LINALG::SIMPLER_Operator::Setup(RCP<Epetra_CrsMatrix> A)
{
  // see whether velocity and pressure solver where configured as ML
  bool visml = vlist_.isSublist("ML Parameters");
  bool pisml = plist_.isSublist("ML Parameters");
  if (!visml) dserror("SIMPLER only works with ML-AMG for velocity");
  
  // get # dofs per node from vlist_ and split row map
  const int ndofpernode = vlist_.sublist("ML Parameters").get<int>("PDE equations",0);
  if (ndofpernode != 4 && ndofpernode !=3) dserror("You should have either 3 or 4 dofs per node at this point");
  const Epetra_Map& fullmap = A->RowMap();
  const int length = fullmap.NumMyElements();
  const int nv     = ndofpernode-1;
  const int nlnode = length / ndofpernode;
  vector<int> vgid(nlnode*nv);
  vector<int> pgid(nlnode);
  int vcount=0;
  for (int i=0; i<nlnode; ++i)
  {
    for (int j=0; j<ndofpernode-1; ++j)
      vgid[vcount++] = fullmap.GID(i*ndofpernode+j);
    pgid[i] = fullmap.GID(i*ndofpernode+ndofpernode-1);
  }
  vector<RCP<Epetra_Map> > maps(2);
  maps[0] = rcp(new Epetra_Map(-1,nlnode*nv,&vgid[0],0,fullmap.Comm()));
  maps[1] = rcp(new Epetra_Map(-1,nlnode,&pgid[0],0,fullmap.Comm()));
  vgid.clear(); pgid.clear();
  mmex_.Setup(fullmap,maps);
  
  // wrap matrix in SparseMatrix and split it into 2x2 BlockMatrix
  SparseMatrix fullmatrix(A);
  A_ = fullmatrix.Split<LINALG::DefaultBlockMatrixStrategy>(mmex_,mmex_);
  A_->Complete();
  
  // split nullspace into velocity and pressure subproblem
  if (visml)
  {
    vlist_.sublist("ML Parameters").set("PDE equations",nv);
    vlist_.sublist("ML Parameters").set("null space: dimension",nv);
    const int vlength = (*A_)(0,0).RowMap().NumMyElements();
    RCP<vector<double> > newns = rcp(new vector<double>(nv*vlength,0.0));
    for (int i=0; i<nlnode; ++i)
    {
      (*newns)[i*nv] = 1.0;
      (*newns)[vlength+i*nv+1] = 1.0;
      (*newns)[2*vlength+i*nv+2] = 1.0;
    }
    vlist_.sublist("ML Parameters").set("null space: vectors",&((*newns)[0]));
    vlist_.sublist("ML Parameters").set<RCP<vector<double> > >("nullspace",newns);
    cout << "vlist_" << vlist_;
  }
  if (pisml)
  {
    plist_.sublist("ML Parameters").set("PDE equations",1);
    plist_.sublist("ML Parameters").set("null space: dimension",1);
    const int plength = (*A_)(1,1).RowMap().NumMyElements();
    RCP<vector<double> > newns = rcp(new vector<double>(plength,1.0));
    plist_.sublist("ML Parameters").set("null space: vectors",&((*newns)[0]));
    plist_.sublist("ML Parameters").set<RCP<vector<double> > >("nullspace",newns);
    cout << "plist_" << plist_;
  }
  
  // Allocate solver for pressure and velocity
  
  
  return;
} 


#endif  // #ifdef CCADISCRET
