/*!----------------------------------------------------------------------
\file discret_fillcomplete.cpp
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
    vector<int> lmdirich;
    vector<int> lmowner;
    actele->LocationVector(lm,lmdirich,lmowner);
    
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
  systemvector.PutScalar(0.0);
  
  // loop through column nodes and set flags indicating whether Neumann
  // condition was already evaluated to false
  for (int i=0; i<NumMyColNodes(); ++i)
  {
    DRT::Node* actnode = lColNode(i);
    // get the Neumann conditions and set evaluated flag to false
    const int zero = 0;
    vector<DRT::Condition*> conds;
    actnode->GetCondition("PointNeumann",conds);
    for (int j=0; j<(int)conds.size(); ++j)
      conds[j]->Add("isevaluated",&zero,1);
    actnode->GetCondition("LineNeumann",conds);
    for (int j=0; j<(int)conds.size(); ++j)
      conds[j]->Add("isevaluated",&zero,1);
    actnode->GetCondition("SurfaceNeumann",conds);
    for (int j=0; j<(int)conds.size(); ++j)
      conds[j]->Add("isevaluated",&zero,1);
    actnode->GetCondition("VolumeNeumann",conds);
    for (int j=0; j<(int)conds.size(); ++j)
      conds[j]->Add("isevaluated",&zero,1);
  } // for (int i=0; i<NumMyColNodes(); ++i)

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  //--------------------------------------------------------
  // loop through Point Neumann conditions and evaluate them
  //--------------------------------------------------------
  for (int i=0; i<NumMyRowNodes(); ++i)
  {
    DRT::Node* actnode = lRowNode(i);
    vector<DRT::Condition*> conds;
    actnode->GetCondition("PointNeumann",conds);
    const int numcond = (int)conds.size();
    if (!numcond) continue;
    const int numdf = actnode->Dof().NumDof();
    const int* dofs = actnode->Dof().Dofs();
    for (int j=0; j<numcond; ++j)
    {
      vector<int>*    iseval = conds[j]->GetVector<int>("isevaluated");
      if ((*iseval)[0]) continue;
      vector<int>*    curve  = conds[j]->GetVector<int>("curve");
      vector<int>*    onoff  = conds[j]->GetVector<int>("onoff");
      vector<double>* val    = conds[j]->GetVector<double>("val");
      // Neumann BCs for some historic reason only have one curve
      int curvenum = -1;
      if (curve) curvenum = (*curve)[0]; 
      double curvefac = 1.0;
      if (curvenum>=0 && usetime)
        dyn_facfromcurve(curvenum,time,&curvefac); 
      for (int k=0; k<numdf; ++k)
      {
        if ((*onoff)[k]==0) continue;
        const int gid   = dofs[k];
        double    value = (*val)[k];
        value *= curvefac;
        int lid = systemvector.Map().LID(gid);
        systemvector[lid] += value;
      }
      // set flag indicating that this Point Neumann condition has been evaluated
      (*iseval)[0] = 1;
    }
  } // for (int i=0; i<NumMyRowNodes(); ++i)   

  //--------------------------------------------------------
  // evaluate line/surface/volume/Neumann conditions
  //--------------------------------------------------------
  RefCountPtr<Epetra_Vector> sysvec = rcp(&systemvector);
  sysvec.release();
  Evaluate(params,null,null,sysvec,null,null);

  cout << systemvector;
  exit(0);


  return;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
