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
  A_      (A)      ,
  cTw_    (0,0)
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

    // Apply the operator
    ierr=A_->Apply(X,Y);

    // if necessary, orthogonalize to matrix kernel in 
    // order to get suitable new basis vectors for the 
    // restricted Krylov space

    // there is only one solution vector --- so solution 
    // vector index is zero
    if(Y.NumVectors()!=1)
    {
      dserror("expecting only one solution vector during AZTEC Apply call\n");
    }

    int sv=0;
    
    // loop all basis vectors of kernel and orthogonalize 
    // against them
    for(int mm=0;mm<c_->NumVectors();++mm)
    {
      /*
                   T
                  c * Y
      */
      double cTY=0.0;
          
      ((*c_)(mm))->Dot(*(Y(sv)),&cTY);

        // loop all weight vectors 
        for(int rr=0;rr<w_->NumVectors();++rr)
        {
          /*
                   T
                  w * c
          */
          double cTw=cTw_(mm,rr);

          /*
                                  T
                       T         c * Y
                      P Y = Y - ------- * w
                                  T
                                 w * c
          */
          (Y(sv))->Update(-cTY/cTw,*((*w_)(rr)),1.0);
          
        } // loop all weight vectors
      } // loop kernel basis vectors
  }
  else
  {
    ierr=A_->Apply(X,Y);
  }

  return(ierr);
} // LINALG::LinalgProjectedOperator::Apply

#endif // CCADISCRET
