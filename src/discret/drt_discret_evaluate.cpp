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
    
  } // for (int i=0; i<numcolele; ++i)
  
  
  
  
  
  
  return;
}




#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
