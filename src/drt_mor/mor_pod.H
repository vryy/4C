/*!----------------------------------------------------------------------
\file mor_pod.H

\brief Model Order Reduction (MOR) using Proper Orthogonal Decomposition (POD)

\level 2

<pre>
\maintainer Martin Pfaller
</pre>

 *----------------------------------------------------------------------*/

#ifndef MOR_POD_H
#define MOR_POD_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Teuchos_ParameterList.hpp>

namespace DRT
{
class Discretization;
}

namespace LINALG
{
  class SparseMatrix;
//  class DefaultBlockMatrixStrategy;
//  class BlockSparseMatrix;
}
namespace UTILS
{
class MOR
{
public:

  /*!
      \brief Constructor
   */
  MOR(Teuchos::RCP<DRT::Discretization> discr);

  /*!
      \brief Destructor
   */
  ~MOR()
  {
    return;
  };

  /*!
      \brief Read in binary matrix from file
   */
  /* Read matrix from specified binary file
   * The binary file has to be formatted like this:
   * Number of Rows: int
   * Number of Columns: int
   * Values (row-wise): float
   */
  void ReadMatrix(std::string filename, Teuchos::RCP<Epetra_MultiVector>& projmatrix);

  Teuchos::RCP<LINALG::SparseMatrix> ReduceDiagnoal(Teuchos::RCP<LINALG::SparseMatrix> M);

  Teuchos::RCP<LINALG::SparseMatrix> ReduceOffDiagonal(Teuchos::RCP<LINALG::SparseMatrix> M);

  Teuchos::RCP<Epetra_MultiVector> ReduceRHS(Teuchos::RCP<Epetra_MultiVector> v);

  Teuchos::RCP<Epetra_Vector> ReduceResidual(Teuchos::RCP<Epetra_Vector> v);

  Teuchos::RCP<Epetra_Vector> ExtendSolution(Teuchos::RCP<Epetra_Vector> v);

  bool HaveMOR()
  {
    return havemor_;
  };

  Teuchos::RCP<Epetra_MultiVector> GetPODMatrix()
  {
    return projmatrix_;
  };

  int GetRedDim()
  {
    return projmatrix_->NumVectors();
  };

private:

  void MultiplyEpetraMultiVectors(
       Teuchos::RCP<Epetra_MultiVector>,
       char,
       Teuchos::RCP<Epetra_MultiVector>,
       char,
       Teuchos::RCP<Epetra_Map>,
       Teuchos::RCP<Epetra_Import>,
       Teuchos::RCP<Epetra_MultiVector>);

  void EpetraMultiVectorToLINALGSparseMatrix(
      Teuchos::RCP<Epetra_MultiVector> multivect,
      Teuchos::RCP<Epetra_Map> rangemap,
      Teuchos::RCP<Epetra_Map> domainmap,
      Teuchos::RCP<LINALG::SparseMatrix> sparsemat);

  bool IsOrthogonal(Teuchos::RCP<Epetra_MultiVector> M);

  Teuchos::RCP<DRT::Discretization> actdisc_;
  int myrank_;  //!< ID of actual processor in parallel
  Teuchos::ParameterList morparams_;
  bool havemor_;
  Teuchos::RCP<Epetra_MultiVector> projmatrix_; //!< projection matrix for POD


  Teuchos::RCP<Epetra_Map> structmapr_; ///< unique map of structure dofs after POD-MOR
  Teuchos::RCP<Epetra_Map> redstructmapr_; ///< full redundant map of structure dofs after POD-MOR
  Teuchos::RCP<Epetra_Import> structrimpo_; ///< importer for fully redundant map of structure dofs after POD-MOR into distributed one
  Teuchos::RCP<Epetra_Import> structrinvimpo_; ///< importer for distributed map of structure dofs after POD-MOR into fully redundant one

}; //class
}
#endif
