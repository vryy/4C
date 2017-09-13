/*!----------------------------------------------------------------------
\file mortar_multifield_coupling.H

\brief Class performing coupling (condensation/recovery) for dual mortar
       methods in (volume) monolithic multi-physics applications, i.e. in
       block matrix systems. This also accounts for the correct condensation
       in the off-diagonal matrix blocks

\level 2

\maintainer Alexander Seitz

*-----------------------------------------------------------------------*/
#ifndef MORTAR_MULTIFIELD_COUPLING_H
#define MORTAR_MULTIFIELD_COUPLING_H

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

// forward declaration
namespace LINALG
{
class SparseMatrix;
class BlockSparseMatrixBase;
}
namespace DRT {
class Discretization;
}

namespace MORTAR
{
class MultiFieldCoupling
{
  public:
  /// c-tor
  MultiFieldCoupling(){};

  /// d-tor
  ~MultiFieldCoupling(){};

  /// add a new discretization to perform coupling on
  void PushBackCoupling(
      const Teuchos::RCP<DRT::Discretization>&   dis,      ///< discretization
      const int nodeset,                                   ///< nodeset to couple
      const std::vector<int> dofs_to_couple                ///< dofs to couple
      );

  /// Perform condensation in all blocks of the matrix
  void CondenseMatrix(Teuchos::RCP<LINALG::BlockSparseMatrixBase>& mat);

  /// Perform condensation in the right-hand side
  void CondenseRhs(Teuchos::RCP<Epetra_Vector>& mat);

  /// recover condensed primal slave-sided dofs
  void RecoverIncr(Teuchos::RCP<Epetra_Vector>& incr);

  private:
  std::vector<Teuchos::RCP<LINALG::SparseMatrix> > p_;
};
}




#endif
