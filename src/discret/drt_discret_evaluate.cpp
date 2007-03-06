/*!----------------------------------------------------------------------
\file drt_discret_fillcomplete.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "drt_discret.H"
#include "drt_dserror.H"
#include "linalg_utils.H"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"




/*----------------------------------------------------------------------*
 |  evaluate (public)                                        mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(
                              ParameterList&                params, 
                              RefCountPtr<Epetra_CrsMatrix> systemmatrix1,
                              RefCountPtr<Epetra_CrsMatrix> systemmatrix2,
                              RefCountPtr<Epetra_Vector>    systemvector1,
                              RefCountPtr<Epetra_Vector>    systemvector2,
                              RefCountPtr<Epetra_Vector>    systemvector3)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  
  // see what we have for input
  bool havesysmatrix1 = false;
  bool havesysmatrix2 = false;
  bool havesysvector1 = false;
  bool havesysvector2 = false;
  bool havesysvector3 = false;
  if (systemmatrix1!=null) havesysmatrix1 = true;
  if (systemmatrix2!=null) havesysmatrix2 = true;
  if (systemvector1!=null) havesysvector1 = true;
  if (systemvector2!=null) havesysvector2 = true;
  if (systemvector3!=null) havesysvector3 = true;
  
  // see what we want to assemble (default is no assembly)
  const bool assemblemat1 = params.get("assemble matrix 1",false);
  const bool assemblemat2 = params.get("assemble matrix 2",false);
  const bool assemblevec1 = params.get("assemble vector 1",false);
  const bool assemblevec2 = params.get("assemble vector 2",false);
  const bool assemblevec3 = params.get("assemble vector 3",false);
  // check whether we have system matrices and vectors supplied to do this
  if (assemblemat1 && !havesysmatrix1) dserror("Do not have system matrix 1 for assembly");
  if (assemblemat2 && !havesysmatrix2) dserror("Do not have system matrix 2 for assembly");
  if (assemblevec1 && !havesysvector1) dserror("Do not have system vector 1 for assembly");
  if (assemblevec2 && !havesysvector2) dserror("Do not have system vector 2 for assembly");
  if (assemblevec3 && !havesysvector3) dserror("Do not have system vector 3 for assembly");

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;
  
  // loop over column elements
  const int numcolele = NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* actele = lColElement(i);
    
    // get element location vector, dirichlet flags and ownerships
    vector<int> lm;
    vector<int> lmowner;
    actele->LocationVector(lm,lmowner);
    
    // get dimension of element matrices and vectors
    // Reshape element matrices and vectors and init to zero
    const int eledim = (int)lm.size();
    elematrix1.Shape(eledim,eledim); 
    elematrix2.Shape(eledim,eledim); 
    elevector1.Size(eledim);
    elevector2.Size(eledim);
    elevector3.Size(eledim);
    
    // call the element evaluate method
    int err = actele->Evaluate(params,*this,lm,elematrix1,elematrix2,
                               elevector1,elevector2,elevector3);
    if (err) dserror("Proc %d: Element %d returned err=%d",Comm().MyPID(),actele->Id(),err);

    if (assemblemat1) LINALG::Assemble(*systemmatrix1,elematrix1,lm,lmowner);
    if (assemblemat2) LINALG::Assemble(*systemmatrix2,elematrix2,lm,lmowner);
    if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,lm,lmowner);
    if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,lm,lmowner);
    if (assemblevec3) LINALG::Assemble(*systemvector3,elevector3,lm,lmowner);

    
  } // for (int i=0; i<numcolele; ++i)
  return;
}

extern "C"
{
  void dyn_facfromcurve(int actcurve,double T,double *fac);  
}
/*----------------------------------------------------------------------*
 |  evaluate Neumann conditions (public)                     mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateNeumann(ParameterList& params, Epetra_Vector& systemvector)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  
  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  multimap<string,RefCountPtr<Condition> >::iterator fool;
  //--------------------------------------------------------
  // loop through Point Neumann conditions and evaluate them
  //--------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != (string)"PointNeumann") continue;
    DRT::Condition& cond = *(fool->second);
    const vector<int>* nodeids = cond.Get<vector<int> >("Node Ids");
    if (!nodeids) dserror("PointNeumann condition does not have nodal cloud");
    const int nnode = (*nodeids).size();
    vector<int>*    curve  = cond.Get<vector<int> >("curve");
    vector<int>*    onoff  = cond.Get<vector<int> >("onoff");
    vector<double>* val    = cond.Get<vector<double> >("val");
    // Neumann BCs for some historic reason only have one curve
    int curvenum = -1;
    if (curve) curvenum = (*curve)[0]; 
    double curvefac = 1.0;
      if (curvenum>=0 && usetime)
        dyn_facfromcurve(curvenum,time,&curvefac); 
    for (int i=0; i<nnode; ++i)
    {
      // do only nodes in my row map
      if (!NodeRowMap()->MyGID((*nodeids)[i])) continue;
      DRT::Node* actnode = gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);
      const int  numdf = actnode->Dof().NumDof();
      const int* dofs  = actnode->Dof().Dofs();
      for (int j=0; j<numdf; ++j)
      {
        if ((*onoff)[j]==0) continue;
        const int gid = dofs[j];
        double value  = (*val)[j];
        value *= curvefac;
        const int lid = systemvector.Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
        systemvector[lid] += value;
      }
    }
  }
  //--------------------------------------------------------
  // loop through line/surface/volume Neumann BCs and evaluate them
  //--------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
    if (fool->first == (string)"LineNeumann" ||
        fool->first == (string)"SurfaceNeumann" ||
        fool->first == (string)"VolumeNeumann"
       )
    {
      DRT::Condition& cond = *(fool->second);
      map<int,RefCountPtr<DRT::Element> >& geom = cond.Geometry();
      map<int,RefCountPtr<DRT::Element> >::iterator curr;
      Epetra_SerialDenseVector elevector;
      for (curr=geom.begin(); curr!=geom.end(); ++curr)
      {
        // get element location vector, dirichlet flags and ownerships
        vector<int> lm;
        vector<int> lmowner;
        curr->second->LocationVector(lm,lmowner);
        elevector.Size((int)lm.size());
        curr->second->EvaluateNeumann(params,*this,cond,lm,elevector);
        LINALG::Assemble(systemvector,elevector,lm,lmowner);
      }
    }
  return;
}

static void DoDirichletCondition(DRT::Condition&      cond,
                                 DRT::Discretization& dis,
                                 const bool           usetime,
                                 const double         time,
                                 Epetra_Vector&       systemvector,
                                 Epetra_Vector&       toggle);
/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateDirichlet(ParameterList& params, 
                                            Epetra_Vector& systemvector,
                                            Epetra_Vector& toggle)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  
  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  multimap<string,RefCountPtr<Condition> >::iterator fool;
  //--------------------------------------------------------
  // loop through Dirichlet conditions and evaluate them
  //--------------------------------------------------------
  // Note that this method does not sum up but 'sets' values in systemvector.
  // For this reason, Dirichlet BCs are evaluated hierarchical meaning
  // in this order:
  //                VolumeDirichlet
  //                SurfaceDirichlet
  //                LineDirichlet
  //                PointDirichlet
  // This way, lower entities override higher ones which is
  // equivalent to inheritance of dirichlet BCs as done in the old
  // ccarat discretization with design          (mgee 1/07)
  
  // Do VolumeDirichlet first
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != (string)"Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::VolumeDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,systemvector,toggle);
  }
  // Do SurfaceDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != (string)"Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::SurfaceDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,systemvector,toggle);
  }
  // Do LineDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != (string)"Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::LineDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,systemvector,toggle);
  }
  // Do LineDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != (string)"Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::PointDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,systemvector,toggle);
  }
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   mwgee 01/07|
 *----------------------------------------------------------------------*/
void DoDirichletCondition(DRT::Condition&      cond,
                          DRT::Discretization& dis,
                          const bool           usetime,
                          const double         time,
                          Epetra_Vector&       systemvector,
                          Epetra_Vector&       toggle)
{
  const vector<int>* nodeids = cond.Get<vector<int> >("Node Ids");
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
  const int nnode = (*nodeids).size();
  vector<int>*    curve  = cond.Get<vector<int> >("curve");
  vector<int>*    onoff  = cond.Get<vector<int> >("onoff");
  vector<double>* val    = cond.Get<vector<double> >("val");
  for (int i=0; i<nnode; ++i)
  {
    // do only nodes in my row map
    if (!dis.NodeRowMap()->MyGID((*nodeids)[i])) continue;
    DRT::Node* actnode = dis.gNode((*nodeids)[i]);
    if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);
    const int  numdf = actnode->Dof().NumDof();
    const int* dofs  = actnode->Dof().Dofs();
    for (int j=0; j<numdf; ++j)
    {
      if ((*onoff)[j]==0) continue;
      const int gid = dofs[j];
      double value  = (*val)[j];
      double curvefac = 1.0;
      int    curvenum = -1;
      if (curve) curvenum = (*curve)[j];
      if (curvenum>=0 && usetime)
        dyn_facfromcurve(curvenum,time,&curvefac);
      //cout << "Dirichlet value " << value << " curvefac " <<  curvefac << endl;
      value *= curvefac;
      const int lid = systemvector.Map().LID(gid);
      if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
      systemvector[lid] = value;
      toggle[lid] = 1.0;
    }
  }
  return;
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
