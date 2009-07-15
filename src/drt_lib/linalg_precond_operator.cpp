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
  precond_(precond),
  cTw_    (0,0)
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
  Epetra_MultiVector       &Y
  ) const
{

  int ierr=0;

  // if necessary, project out matrix kernel to maintain well-posedness
  // of problem
  if(project_)
  {
    // Apply the inverse preconditioner to get new basis vector for the
    // Krylov space
    ierr=precond_->ApplyInverse(X,Y);

    // check for vectors for matrix kernel and weighted basis mean vector
    if(c_ == Teuchos::null || w_ == Teuchos::null)
    {
      dserror("no c_ and w_ supplied");
    }

    // there is only one solution vector --- so solution
    // vector index is zero
    if(Y.NumVectors()!=1)
    {
      dserror("expecting only one solution vector during AZTEC Apply call\n");
    }

    int v=0;

    // loop all weight vectors and orthogonalize against them
    for(int rr=0;rr<w_->NumVectors();++rr)
    {

      /*
                   T
                  w * Y
      */
      double wTY=0.0;

      (Y(v))->Dot(*((*w_)(rr)),&wTY);

      // loop all basis vectors of kernel
      for(int mm=0;mm<c_->NumVectors();++mm)
      {

        /*
                   T
                  w * c
        */
        double cTw=cTw_(mm,rr);

        /*

                                    T
                                   x * w
                        P x = x - ------- c
                                    T
                                   w * c

        */
        (Y(v))->Update(-wTY/cTw,*((*c_)(mm)),1.0);
      } // loop all weight vectors
    } // loop kernel basis vectors
  }
  else
  {
    ierr=precond_->ApplyInverse(X,Y);
  }

  return(ierr);
} // LINALG::LinalgPrecondOperator::ApplyInverse


#endif // CCADISCRET
