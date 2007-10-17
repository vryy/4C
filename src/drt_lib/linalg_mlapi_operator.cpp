/*!----------------------------------------------------------------------
\file linalg_mlapi_operator.cpp

\class LINALG::AMG_Operator

\brief A multipurpose experimental multigrid operator

This operator based on the ml advanced programming interface is a
multipurpose development object for amg ideas that shall be tested in
the baci framework

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#include "linalg_mlapi_operator.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
LINALG::AMG_Operator::AMG_Operator(RCP<Epetra_CrsMatrix> A, 
                                   ParameterList& params, 
                                   const bool compute) :
Epetra_Operator(),
label_("LINALG::AMG_Operator"),
params_(params),
A_(A)
{
  return;
}


/*----------------------------------------------------------------------*
 |  return label of operator                                 mwgee 10/07|
 *----------------------------------------------------------------------*/
const char* LINALG::AMG_Operator::Label() const
{
  return &label_[0];
}


/*----------------------------------------------------------------------*
 |  setup phase (private)                                    mwgee 10/07|
 *----------------------------------------------------------------------*/
void LINALG::AMG_Operator::Compute()
{
  return;
} 


/*----------------------------------------------------------------------*
 |  apply operator (public)                                  mwgee 10/07|
 *----------------------------------------------------------------------*/
int LINALG::AMG_Operator::ApplyInverse(const Epetra_MultiVector& X, 
                                             Epetra_MultiVector& Y) const
{
  // do not do anything for testing
  Y.Update(1.0,X,0.0);
  
  return 0;
} 




#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
