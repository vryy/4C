/*!----------------------------------------------------------------------
\file linalg_downwindmatrix.cpp

\class LINALG::DownwindMatrix

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
#include "linalg_downwindmatrix.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 03/08|
 *----------------------------------------------------------------------*/
LINALG::DownwindMatrix::DownwindMatrix(RCP<Epetra_CrsMatrix> A) :
A_(A)
{
  Setup(A);

  return;
}


/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 02/08|
 *----------------------------------------------------------------------*/
void LINALG::DownwindMatrix::Setup(RCP<Epetra_Operator> A)
{

  return;
}














#endif  // #ifdef CCADISCRET
