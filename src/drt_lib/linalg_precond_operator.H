/*!----------------------------------------------------------------------
\file  linalg_precond_operator.H

\brief A common interface for ifpack, ml and simpler preconditioners.
       This interface allows a modification of the vector returned
       by the ApplyInverse call, which is necessary to do a solution on
       a Krylov space orthogonalized to certain (for example rigid body
       or zero pressure) modes.

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef LINALG_PRECOND_OPERATOR_H_
#define LINALG_PRECOND_OPERATOR_H_


#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "drt_dserror.H"


namespace LINALG
{
  /*!

  A common interface for ifpack, ml and simpler preconditioners.
  This interface allows a modification of the vector returned
  by the ApplyInverse call, which is necessary to do a solution on
  a Krylov space orthogonalized to certain (for example rigid body)
  modes.

  The linalg preconditioner interface class holds a pointer (rcp)
  to the actual preconditioner. All methods implemented to support
  the Epetra_Operator interface just call the corresponding functions
  of the actual preconditioner.

  Only the ApplyInverse method is modified.

  See linalg_projected_operator.H for related docu and code.

  \author gammi
  \date 03/09
  */
class LinalgPrecondOperator : public Epetra_Operator
{
public:

  /*!
  \brief Standard Constructor

  */
  LinalgPrecondOperator(Teuchos::RCP<Epetra_Operator> precond,bool project);


  /*!
  \brief Destructor

  */
  ~LinalgPrecondOperator();

  //! @name Special atribute set methods
  /*!

  \brief Set kernel space and weight vectors.

  The projector P is completely defined by w and c.

  \param w vector(s) of weights
  \param c kernel basis vector(s)

  */
  void SetKernelSpace(
    Teuchos::RCP<Epetra_MultiVector> w,
    Teuchos::RCP<Epetra_MultiVector> c
    )
  {
    w_=w;
    c_=c;

    if(project_)
    {
      if(w_==Teuchos::null || c_==Teuchos::null)
      {
        dserror("no kernel supplied for projection (but projection flag was set)");
      }

      cTw_.Shape(c_->NumVectors(),w_->NumVectors());

      for(int mm=0;mm<c_->NumVectors();++mm)
      {
        // loop all weight vectors
        for(int rr=0;rr<w_->NumVectors();++rr)
        {

          /*
                   T
                  w * c
          */
          ((*c_)(mm))->Dot(*((*w_)(rr)),&(cTw_(mm,rr)));

          if(fabs(cTw_(mm,rr))<1e-14)
          {
            dserror("weight vector must not be orthogonal to c");
          }
        }
      }
    }
    return;
  }
  //! @}

  //! @name Atribute set methods required to support the Epetra_Operator interface

  virtual int SetUseTranspose (bool UseTranspose)
  {
    return(precond_->SetUseTranspose(UseTranspose));
  }
  //! @}

  //! @name Mathematical functions required to support the Epetra_Operator interface (pass through)
  virtual int Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
  {
    return(precond_->Apply(X,Y));
  }

  virtual double NormInf () const
  {
    return(precond_->NormInf());
  }
  //! @}

  //! @name Atribute access functions required to support the Epetra_Operator interface
  virtual const char * Label () const
  {
    return(precond_->Label());
  }

  virtual bool UseTranspose () const
  {
    return(precond_->UseTranspose());
  }

  virtual bool HasNormInf () const
  {
    return(precond_->HasNormInf());
  }

  virtual const Epetra_Comm & Comm () const
  {
    return(precond_->Comm());
  }

  virtual const Epetra_Map & OperatorDomainMap () const
  {
    return(precond_->OperatorDomainMap());
  }

  virtual const Epetra_Map & OperatorRangeMap () const
  {
    return(precond_->OperatorRangeMap());
  }

  //! @}

  //! @name Mathematical functions required to support the Epetra_Operator interface (modified)
  /*
    (Modified) ApplyInverse call

    This method calls ApplyInverse on the actual preconditioner and, the
    solution is orthogonalized against a set of weight vectors provided in a
    multivector.

    This is done using a projector P defined by

                                    T
                                   x * w
                        P  x = x - ------ c
                                    T
                                   w * c

    w is the vector of weights, c a vector of ones (in the dofs under
    consideration) corresponding to the matrix kernel.

    The system we are solving with this procedure is not Au=b for u (since A
    might be singular), but we are solving

                        / T \         T
                       | P A | P u = P b ,
                        \   /

    for the projection of the solution Pu, i.e. in the preconditioned case


                                                          -+
           / T   \     /      -1 \          T              |
          | P * A | * |  P * M    | * xi = P  * b          |
           \     /     \         /                         |
                                                -1         |
                                       x = P * M  * xi     |
                                                          -+


    Hence, P is always associated with the apply inverse call of the
    preconditioner (the right bracket) and always called after the call
    to ApplyInverse.


    Properties of P are:

    1) c defines the kernel of P, i.e. P projects out the matrix kernel

                          T
                         c * w
              P c = c - ------- c = c - c = 0
                          T
                         w * c

    2) The space spanned by P x is orthogonal to the weight vector

                       /      T      \              T
     T   /   \     T  |      x * w    |    T       x * w     T       T       T
    w * | P x | = w * | x - ------- c | = w * x - ------- * w * c = w * x - w * x = 0
         \   /        |       T       |             T
                       \     w * c   /             w * c


    This modified Apply call is for singular matrices A when c is
    a vector defining A's nullspace. The preceding projection
    operator ensures
                            |           |
                           -+-         -+-T
                  A u = A u     where u    * c =0,

    even if A*c != 0 (for numerical inaccuracies during the computation
    of A)

    See the following article for further reading:

    @article{1055401,
     author = {Bochev,, Pavel and Lehoucq,, R. B.},
     title = {On the Finite Element Solution of the Pure Neumann Problem},
     journal = {SIAM Rev.},
     volume = {47},
     number = {1},
     year = {2005},
     issn = {0036-1445},
     pages = {50--66},
     doi = {http://dx.doi.org/10.1137/S0036144503426074},
     publisher = {Society for Industrial and Applied Mathematics},
     address = {Philadelphia, PA, USA},
     }

  */
  virtual int ApplyInverse (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;


private:
  //! flag whether to do a projection or just pass through
  bool                               project_;

  //! the actual preconditioner
  Teuchos::RCP<Epetra_Operator>      precond_;

  //! a set of vectors defining weighted (basis integral) vector for the projector
  Teuchos::RCP<Epetra_MultiVector>   w_;

  //! a set of vectors defining the vectors of ones (in the respective components)
  //! for the matrix kernel
  Teuchos::RCP<Epetra_MultiVector>   c_;

  //! values for scalar products c_i o w_j, computed once after setting c_ and w_
  Epetra_SerialDenseMatrix           cTw_;

};

} // namespace LINALG

#endif // LINALG_PRECOND_OPERATOR_H_

#endif // CCADISCRET
