/*-----------------------------------------------------------------------*/
/*! \file
\brief Class performing coupling (condensation/recovery) for dual mortar
       methods in (volume) monolithic multi-physics applications, i.e. in
       block matrix systems. This also accounts for the correct condensation
       in the off-diagonal matrix blocks

\level 2

*/
/*-----------------------------------------------------------------------*/
#ifndef FOUR_C_MORTAR_MULTIFIELD_COUPLING_HPP
#define FOUR_C_MORTAR_MULTIFIELD_COUPLING_HPP

#include "4C_config.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Core::LinAlg
{
  class SparseMatrix;
  class BlockSparseMatrixBase;
}  // namespace Core::LinAlg

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Mortar
{
  class MultiFieldCoupling
  {
   public:
    /// c-tor
    MultiFieldCoupling(){};


    /// add a new discretization to perform coupling on
    void PushBackCoupling(const Teuchos::RCP<Core::FE::Discretization>& dis,  ///< discretization
        const int nodeset,                                                    ///< nodeset to couple
        const std::vector<int> dofs_to_couple                                 ///< dofs to couple
    );

    /// Perform condensation in all blocks of the matrix
    void CondenseMatrix(Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>& mat);

    /// Perform condensation in the right-hand side
    void CondenseRhs(Teuchos::RCP<Epetra_Vector>& rhs);

    /// recover condensed primal slave-sided dofs
    void RecoverIncr(Teuchos::RCP<Epetra_Vector>& incr);

   private:
    std::vector<Teuchos::RCP<Core::LinAlg::SparseMatrix>> p_;
  };
}  // namespace Mortar



FOUR_C_NAMESPACE_CLOSE

#endif
