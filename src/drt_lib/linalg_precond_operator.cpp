/*!----------------------------------------------------------------------
\file  linalg_precond_operator.cpp

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "linalg_precond_operator.H"

/* --------------------------------------------------------------------
                          Constructor
   -------------------------------------------------------------------- */
LINALG::LinalgPrecondOperator::LinalgPrecondOperator(
  Teuchos::RCP<Epetra_Operator> precond,
  bool                          project) :
  project_(project),
  precond_(precond)
{
  return;
} // LINALG::LinalgPrecondOperator::LinalgPrecondOperator

/* --------------------------------------------------------------------
                          Destructor
   -------------------------------------------------------------------- */
LINALG::LinalgPrecondOperator::~LinalgPrecondOperator()
{
  return;
} // LINALG::LinalgPrecondOperator::~LinalgPrecondOperator

/* --------------------------------------------------------------------
                    (Modified) ApplyInverse call
   -------------------------------------------------------------------- */
int LINALG::LinalgPrecondOperator::ApplyInverse(
  const Epetra_MultiVector &X, 
  Epetra_MultiVector &Y
  ) const
{
  int ierr=precond_->ApplyInverse(X,Y);

  if(project_)
  {
    const int kerneldim = w_->NumVectors();
    const int numsolvecs= Y.NumVectors();

    // define a C vector of 1.0s
    Epetra_Vector c(*((*w_)(0)));
    c.PutScalar(1.0);

    // loop all solution vectors
    for(int sv=0;sv<numsolvecs;++sv)
    {
    
      // loop all basis vectors of kernel and orthogonalize against them
      for(int rr=0;rr<kerneldim;++rr)
      {
        double wTc=0.0;
        double cTY=0.0;

        /*
                   T
                  w * c
        */
        c.Dot(*((*w_)(rr)),&wTc);

        /*
                   T
                  c * sol
        */
        c.Dot(*(Y(sv)),&cTY);

        /*
                                  T
                       T         c * Y
                      P Y = Y - ------- * w
                                  T
                                 w * c
        */
        (Y(sv))->Update(-cTY/wTc,*((*w_)(rr)),1.0);
      }
    }
  }

  return(ierr);
} // LINALG::LinalgPrecondOperator::ApplyInverse


#endif // CCADISCRET
