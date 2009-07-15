/*!----------------------------------------------------------------------
\file  linalg_projected_operator.H

\brief

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef LINALG_PROJECTED_OPERATOR_H_
#define LINALG_PROJECTED_OPERATOR_H_

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "drt_dserror.H"

#include "linalg_sparsematrix.H"

namespace LINALG
{
  /*!
  \brief A slightly modified Epetra_Operator class allowing the
         application of a 'rank 1' projector P^T after applying the
         underlying operator

  Only the Apply method is modified. All other calls are 'pass-through'.
  The modification of the apply call is such that for singular matrices
  the result after the apply call is orthogonal to the matrix kernel c.

  System of equations to solve:

                         A * x = b

  (Right-) preconditioned system:
                                                -+
                              -1                 |
                         A * M  xi = b           |
                                                 |
                                      -1         |
                                 x = M  * xi     |
                                                -+

  A is singular, so use projectors for singular matrices; for the right
  preconditioned system we solve:

                                                -+
           T           -1         T              |
          P * A * P * M   * xi = P  * b          |
                                                 |
                                      -1         |
                             x = P * M  * xi     |
                                                -+

  P^T makes resulting Krylov vectors orthogonal to the kernel c of A.
  P ensures A * P * c = 0 even if for numerical errors A * c != 0 and
  defines a solution x with zero mean.

  The whole operator is split in matrix and preconditioner part:

           / T   \     /      -1 \
          | P * A | * |  P * M    |
           \     /     \         /

  This class deals with the first brackets via a modified Apply call, the
  right bracket is done in the corresponding preconditioner class.

  Abbreviations:

  A     : operator (singular matrix)
  M^{-1}: inverse preconditioner
  P, P^T: projectors

  \author gammi
  \date 04/09
  */
class LinalgProjectedOperator : public Epetra_Operator
{
public:

  /*!
  \brief Standard Constructor

  */
  LinalgProjectedOperator(Teuchos::RCP<Epetra_Operator> A,bool project);


  /*!
  \brief Destructor

  */
  ~LinalgProjectedOperator();

  //! @name Special atribute set methods
  /*!

  \brief Set kernel space and weight vectors.

  The projector P^T is completely defined by w and c.

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
        dserror("no kernel/weight supplied for projection (but projection flag was set)");
      }

      // resize
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
    return(A_->SetUseTranspose(UseTranspose));
  }
  //! @}

  //! @name Mathematical functions required to support the Epetra_Operator interface (modified)
  /*!
    \brief (Modified) Apply call

    Given the underlying Epetra_Operator A, we apply P^T*A using

                                        T
                         T         w * c
                        P  x = x - ------ x
                                    T
                                   w * c                                    T

    instead of A if projetion_ flag is set. If projection_ is false,
    this is just the standard Apply call.


    See the following article for further reading

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

     and the documentation to LinalgPrecondOperator as well as the source
     to LINALG::Solver::Solve_aztec.

  */
  //! @}
  virtual int Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

  virtual double NormInf () const
  {
    return(A_->NormInf());
  }
  //! @}

  //! @name Atribute access functions required to support the Epetra_Operator interface
  virtual const char * Label () const
  {
    return(A_->Label());
  }

  virtual bool UseTranspose () const
  {
    return(A_->UseTranspose());
  }

  virtual bool HasNormInf () const
  {
    return(A_->HasNormInf());
  }

  virtual const Epetra_Comm & Comm () const
  {
    return(A_->Comm());
  }

  virtual const Epetra_Map & OperatorDomainMap () const
  {
    return(A_->OperatorDomainMap());
  }

  virtual const Epetra_Map & OperatorRangeMap () const
  {
    return(A_->OperatorRangeMap());
  }

  virtual int ApplyInverse (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
  {
    return(A_->ApplyInverse(X,Y));
  }

  Teuchos::RCP<Epetra_Operator> UnprojectedOperator ()
  {
    return(A_);
  }

private:
  //! flag whether to do a projection or just pass through
  bool                               project_;

  //! the actual unprojected operator
  Teuchos::RCP<Epetra_Operator>      A_;

  //! a set of vectors defining weighted basis mean vector for the projector
  Teuchos::RCP<Epetra_MultiVector>   w_;

  //! a set of vectors defining the vectors of ones (in the respective components)
  //! for the matrix kernel
  Teuchos::RCP<Epetra_MultiVector>   c_;

  //! values for scalar products c_i o w_j, computed once after setting c_ and w_
  Epetra_SerialDenseMatrix           cTw_;

};

} // namespace LINALG

#endif // LINALG_PROJECTED_OPERATOR_H_

#endif // CCADISCRET
