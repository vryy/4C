/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_LINALG_PROJECTED_OPERATOR_HPP
#define FOUR_C_LINALG_PROJECTED_OPERATOR_HPP

#include "4C_config.hpp"

#include <Epetra_Operator.h>
#include <Teuchos_RCP.hpp>

// forward declarations
class Epetra_MultiVector;

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
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
        Teuchos::RCP<Core::LinAlg::KrylovProjector> projector);



    //! @name Atribute set methods required to support the Epetra_Operator interface

    int SetUseTranspose(bool UseTranspose) override { return (a_->SetUseTranspose(UseTranspose)); }
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
       to Core::LinAlg::Solver.

    */
    //! @}
    int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const override;

    double NormInf() const override { return (a_->NormInf()); }
    //! @}

    //! @name Atribute access functions required to support the Epetra_Operator interface
    const char *Label() const override { return (a_->Label()); }

    bool UseTranspose() const override { return (a_->UseTranspose()); }

    bool HasNormInf() const override { return (a_->HasNormInf()); }

    const Epetra_Comm &Comm() const override { return (a_->Comm()); }

    const Epetra_Map &OperatorDomainMap() const override { return (a_->OperatorDomainMap()); }

    const Epetra_Map &OperatorRangeMap() const override { return (a_->OperatorRangeMap()); }

    int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const override
    {
      return (a_->ApplyInverse(X, Y));
    }

    Teuchos::RCP<Epetra_Operator> UnprojectedOperator() { return (a_); }

   private:
    //! flag whether to do a projection or just pass through
    bool project_;

    //! the actual unprojected operator
    Teuchos::RCP<Epetra_Operator> a_;

    //! Krylov space projector
    Teuchos::RCP<Core::LinAlg::KrylovProjector> projector_;
  };

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
