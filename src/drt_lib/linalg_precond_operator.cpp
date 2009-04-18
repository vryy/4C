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
  Epetra_MultiVector       &Y
  ) const
{ 

  int ierr=0;

  // if necessary, project out matrix kernel
  if(project_)
  {
    const int kerneldim = w_->NumVectors();
    Epetra_MultiVector PX(X);

    const int num= PX.NumVectors();
  
    // loop all solution vectors
    for(int v=0;v<num;++v)
    {
      // loop all basis vectors of kernel and orthogonalize against them
      for(int rr=0;rr<kerneldim;++rr)
      {
        /*
                   T
                  w * c
        */
        double wTc=0.0;

        c_->Dot(*((*w_)(rr)),&wTc);

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

        (PX(v))->Update(-wTX/wTc,*c_,1.0);
      }
      

    }   

    // Apply the inverse preconditioner to get new basis vector for the
    // Krylov space
    ierr=precond_->ApplyInverse(PX,Y);

    const int numsolvecs= Y.NumVectors();

    // check for vectors for matrix kernel and weighted basis mean vector
    if(c_ == Teuchos::null || w_ == Teuchos::null)
    {
      dserror("no c_ and w_ supplied");
    }

    // loop all solution vectors
    for(int sv=0;sv<numsolvecs;++sv)
    {
   
      // loop all basis vectors of kernel and orthogonalize against them
      for(int rr=0;rr<kerneldim;++rr)
      {
        /*
                   T
                  w * c
        */
        double wTc=0.0;

        c_->Dot(*((*w_)(rr)),&wTc);

	if(fabs(wTc)<1e-14)
	{
	  dserror("weight vector must not be orthogonal to c");
	}

        /*
                   T
                  c * Y
        */
        double cTY=0.0;

        c_->Dot(*(Y(sv)),&cTY);

        /*
                                  T
                       T         c * Y
                      P Y = Y - ------- * w
                                  T
                                 w * c
        */
        (Y(sv))->Update(-cTY/wTc,*((*w_)(rr)),1.0);

        /*
                   T
                  w * Y 
        */
        double wTY=0.0;
 
        (Y(sv))->Dot(*((*w_)(rr)),&wTY);

        //        (Y(sv))->Update(-wTY/wTc,*c_,1.0);

      } // loop kernel basis vectors
    } // loop all solution vectors

  }
  else
  {
    ierr=precond_->ApplyInverse(X,Y);
  }

  // if necessary, project out matrix kernel
  if(project_)
  {
  } // if (project_)

  return(ierr);
} // LINALG::LinalgPrecondOperator::ApplyInverse


#endif // CCADISCRET
