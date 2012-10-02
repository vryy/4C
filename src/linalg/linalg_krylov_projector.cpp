/*!----------------------------------------------------------------------
\file  linalg_krylov_projector.cpp

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "linalg_krylov_projector.H"

/* --------------------------------------------------------------------
                          Constructor
   -------------------------------------------------------------------- */
LINALG::KrylovProjector::KrylovProjector(
  const bool                       project          ,
  Teuchos::RCP<Epetra_MultiVector> w                ,
  Teuchos::RCP<Epetra_MultiVector> c                ,
  Teuchos::RCP<Epetra_Operator>    A
  ) :
  project_(project),
  nsdim_(0),
  cTw_(0)
{
  w_=w;
  c_=c;

  if(project_)
  {
    if(w_==Teuchos::null || c_==Teuchos::null)
    {
      dserror("no kernel supplied for projection (but projection flag was set)");
    }

    nsdim_ = c_->NumVectors();
    if(w_->NumVectors() != nsdim_)
    {
      dserror("number of basis and weight vectors are not the same");
    }
    cTw_.Resize(nsdim_);

    // loop all kernel basis vectors
    for(int mm=0;mm<nsdim_;++mm)
    {

      // for each kernel vector provided check A*c=0
      if(A!=Teuchos::null)
      {
        Epetra_Vector result(((*c_)(mm))->Map(),false);

        A->Apply(*((*c_)(mm)),result);

        double norm=1e9;

        result.Norm2(&norm);

        if(norm>1e-12)
        {
          std::cout << "#####################################################" << std::endl;
          std::cout << "Krylov projection failed!                            " << std::endl;
          std::cout << "This might be caused by:                             " << std::endl;
          std::cout << " - you don't have pure Dirichlet boundary conditions " << std::endl;
          std::cout << "   or pbcs -> check your inputfile                   " << std::endl;
          std::cout << " - you don't integrate the pressure exactly and the  " << std::endl;
          std::cout << "   given kernel is not a kernel of your system -> to " << std::endl;
          std::cout << "   check this, use more gauss points (often problem  " << std::endl;
          std::cout << "   with nurbs)                                       " << std::endl;
          std::cout << " - there is indeed a problem with the Krylov projection " << std::endl;
          std::cout << "#####################################################" << std::endl;
          dserror("krylov projection failed, Ac returned %12.5e for kernel basis vector %d",norm,mm);
        }
      }

      // loop all weight vectors
      for(int rr=0;rr<nsdim_;++rr)
      {
        /*
                 T
                w * c
         */
        double cTw;
        ((*c_)(mm))->Dot(*((*w_)(rr)),&cTw);

        if (rr==mm)
        {
          if(fabs(cTw)<1e-14)
          {
            dserror("weight vector w_%i must not be orthogonal to c_%i",rr,mm);
          }
          else
            cTw_(mm) = cTw; // store value of scalar product
        }
        else
        {
          if(fabs(cTw)>1e-14)
          {
            dserror("weight vector w_%i must be orthogonal to c_%i",rr,mm);
          }
        }

      }
    }
  }
  return;
}
// LINALG::KrylovProjector::KrylovProjector

/* --------------------------------------------------------------------
                          Destructor
   -------------------------------------------------------------------- */
LINALG::KrylovProjector::~KrylovProjector()
{
  return;
} // LINALG::KrylovProjector::~KrylovProjector

/* --------------------------------------------------------------------
                    Apply projector P
   -------------------------------------------------------------------- */
int LINALG::KrylovProjector::ApplyP(Epetra_MultiVector& Y) const
{

  int ierr=0;

  // if necessary, project out matrix kernel to maintain well-posedness
  // of problem
  if(project_)
  {
    // there is only one solution vector --- so solution
    // vector index is zero
    if(Y.NumVectors()!=1)
    {
      dserror("expecting only one solution vector during AZTEC Apply call\n");
    }

    int v=0;

    // loop kernel basis/weight vector pairs
    for(int rr=0;rr<nsdim_;++rr)
    {

      /*
                   T
                  w * Y
      */
      double wTY=0.0;

      (Y(v))->Dot(*((*w_)(rr)),&wTY);

        /*
                   T
                  w * c
        */
        double cTw=cTw_(rr);

        /*

                                    T
                                   x * w
                        P x = x - ------- c
                                    T
                                   w * c

        */
        (Y(v))->Update(-wTY/cTw,*((*c_)(rr)),1.0);

    } // loop kernel basis/weight vector pairs
  }

  return(ierr);
} // LINALG::KrylovProjector::ApplyP


/* --------------------------------------------------------------------
                    Apply projector P^T
   -------------------------------------------------------------------- */
int LINALG::KrylovProjector::ApplyPT(Epetra_MultiVector& Y) const
{
  int ierr=0;

  // if necessary, project out matrix kernel
  if(project_)
  {
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

    // loop kernel basis/weight vector pairs and orthogonalize against them
    for(int rr=0;rr<nsdim_;++rr)
    {
      /*
                   T
                  c * Y
      */
      double cTY=0.0;

      ((*c_)(rr))->Dot(*(Y(sv)),&cTY);

          /*
                   T
                  w * c
          */
          double cTw=cTw_(rr);

          /*
                                  T
                       T         c * Y
                      P Y = Y - ------- * w
                                  T
                                 w * c
          */
          (Y(sv))->Update(-cTY/cTw,*((*w_)(rr)),1.0);

      } // loop kernel basis/weight vector pairs
  }

  return(ierr);
} // LINALG::KrylovProjector::ApplyPT


