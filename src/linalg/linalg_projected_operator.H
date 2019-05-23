/*!----------------------------------------------------------------------

\brief Declaration

\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
*/
/*----------------------------------------------------------------------*/

#ifndef LINALG_PROJECTED_OPERATOR_H_
#define LINALG_PROJECTED_OPERATOR_H_

#include "Teuchos_RCP.hpp"
#include "Epetra_Operator.h"

// forward declarations
class Epetra_MultiVector;

namespace LINALG
{
  class KrylovProjector;

  /*! @name

  \brief A slightly modified Epetra_Operator class allowing the
         application of a 'rank 1' projector P^T after applying the
         underlying operator

  Only the Apply method is modified. All other calls are 'pass-through'.
  The modification of the apply call is such that for singular matrices
  the result after the apply call is krylov to the matrix kernel c.

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

  P^T makes resulting Krylov vectors krylov to the kernel c of A.
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
    LinalgProjectedOperator(Teuchos::RCP<Epetra_Operator> A, bool project,
        Teuchos::RCP<LINALG::KrylovProjector> projector);

    /*!
    \brief Destructor

    */
    ~LinalgProjectedOperator();


    //! @name Atribute set methods required to support the Epetra_Operator interface

    virtual int SetUseTranspose(bool UseTranspose) { return (A_->SetUseTranspose(UseTranspose)); }
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
    virtual int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    virtual double NormInf() const { return (A_->NormInf()); }
    //! @}

    //! @name Atribute access functions required to support the Epetra_Operator interface
    virtual const char *Label() const { return (A_->Label()); }

    virtual bool UseTranspose() const { return (A_->UseTranspose()); }

    virtual bool HasNormInf() const { return (A_->HasNormInf()); }

    virtual const Epetra_Comm &Comm() const { return (A_->Comm()); }

    virtual const Epetra_Map &OperatorDomainMap() const { return (A_->OperatorDomainMap()); }

    virtual const Epetra_Map &OperatorRangeMap() const { return (A_->OperatorRangeMap()); }

    virtual int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
    {
      return (A_->ApplyInverse(X, Y));
    }

    Teuchos::RCP<Epetra_Operator> UnprojectedOperator() { return (A_); }

   private:
    //! flag whether to do a projection or just pass through
    bool project_;

    //! the actual unprojected operator
    Teuchos::RCP<Epetra_Operator> A_;

    //! Krylov space projector
    Teuchos::RCP<LINALG::KrylovProjector> projector_;
  };

}  // namespace LINALG

#endif  // LINALG_PROJECTED_OPERATOR_H_
