/*!----------------------------------------------------------------------
\file  linalg_projected_operator.cpp

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "linalg_projected_operator.H"

/* --------------------------------------------------------------------
                          Constructor
   -------------------------------------------------------------------- */
LINALG::LinalgProjectedOperator::LinalgProjectedOperator(
  Teuchos::RCP<Epetra_Operator> A      ,
  bool                          project) :
  project_(project),
  A_      (A)
{
  return;
} // LINALG::LinalgProjectedOperator::LinalgProjectedOperator

/* --------------------------------------------------------------------
                           Destructor
   -------------------------------------------------------------------- */
LINALG::LinalgProjectedOperator::~LinalgProjectedOperator()
{
  return;
} // LINALG::LinalgProjectedOperator::~LinalgProjectedOperator

/* --------------------------------------------------------------------
                      (Modified) Apply call
   -------------------------------------------------------------------- */
int LINALG::LinalgProjectedOperator::Apply(
  const Epetra_MultiVector &X, 
  Epetra_MultiVector       &Y
  ) const
{ 
  int ierr=0;

  // if necessary, project out matrix kernel
  if(project_)
  {
    // safety checks
    if(w_ == Teuchos::null)
    {
      dserror("weight vector has not been set\n");
    }
    if(c_ == Teuchos::null)
    {
      dserror("kernel vector has not been set\n");
    }

    // get a copy of the input vector to apply projector P
    Epetra_MultiVector PX(X);

    // loop all solution vectors
    for(int v=0;v<PX.NumVectors();++v)
    {
      // loop all basis vectors of kernel
      for(int mm=0;mm<c_->NumVectors();++mm)
      {
        // loop all weight vectors and orthogonalize against them
        for(int rr=0;rr<w_->NumVectors();++rr)
        {
          /*
                   T
                  w * c
          */
          double wTc=0.0;

          ((*c_)(mm))->Dot(*((*w_)(rr)),&wTc);

          if(fabs(wTc)<1e-14)
          {
            dserror("weight vector must not be orthogonal to c");
          }

          /*
                   T
                  w * X
          */
          double wTX=0.0;
 
          (PX(v))->Dot(*((*w_)(rr)),&wTX);

          /*

                                    T   
                                   x * w
                        P x = x - ------- c
                                    T
                                   w * c

          */
          (PX(v))->Update(-wTX/wTc,*((*c_)(mm)),1.0);
        } // rr
      } // mm
    } // v

    // Apply the operator to the projected input vector in order 
    // to get new basis vector for the Krylov space
    ierr=A_->Apply(PX,Y);
  }
  else
  {
    ierr=A_->Apply(X,Y);
  }

  return(ierr);
} // LINALG::LinalgProjectedOperator::Apply


#endif // CCADISCRET
