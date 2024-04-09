/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_FEM_GENERAL_EXTRACT_VALUES_HPP
#define FOUR_C_DISCRETIZATION_FEM_GENERAL_EXTRACT_VALUES_HPP

#include "baci_config.hpp"

#include "baci_lib_element.hpp"
#include "baci_lib_node.hpp"

BACI_NAMESPACE_OPEN

namespace CORE::FE
{
  /*!
  \brief Locally extract a subset of values from an Epetra_Vector

  Extracts lm.size() values from a distributed epetra vector and stores them into local.
  this is NOT a parallel method, meaning that all values to be extracted on a processor
  must be present in global on that specific processor. This usually means that global
  has to be in column map style.

  \param global (in): global distributed vector with values to be extracted
  \param local (out): vector or matrix holding values extracted from global
  \param lm     (in): vector containing global ids to be extracted. Size of lm
                      determines number of values to be extracted.
  */
  void ExtractMyValues(
      const Epetra_Vector& global, std::vector<double>& local, const std::vector<int>& lm);

  void ExtractMyValues(const Epetra_Vector& global, CORE::LINALG::SerialDenseVector& local,
      const std::vector<int>& lm);

  void ExtractMyValues(
      const Epetra_MultiVector& global, std::vector<double>& local, const std::vector<int>& lm);

  template <class matrix>
  void ExtractMyValues(
      const Epetra_Vector& global, std::vector<matrix>& local, const std::vector<int>& lm)
  {
    // safety check
    if (local[0].N() != 1 or local.size() * (unsigned)local[0].M() != lm.size())
      dserror("Received matrix vector of wrong size!");

    // loop over all nodes of current element
    for (unsigned inode = 0; inode < local[0].M(); ++inode)
    {
      // loop over all dofs of current node
      for (unsigned idof = 0; idof < local.size(); ++idof)
      {
        // extract local ID of current dof
        const int lid = global.Map().LID(lm[inode * local.size() + idof]);

        // safety check
        if (lid < 0)
          dserror("Proc %d: Cannot find gid=%d in Epetra_Vector", global.Comm().MyPID(),
              lm[inode * local.size() + idof]);

        // store current dof in local matrix vector consisting of ndof matrices of size nnode x 1,
        // where nnode denotes the number of element nodes and ndof denotes the number of degrees
        // of freedom per element node.
        local[idof](inode, 0) = global[lid];
      }
    }
  }

  template <class matrix>
  void ExtractMyValues(const Epetra_Vector& global, matrix& local, const std::vector<int>& lm)
  {
    // safety check
    if ((unsigned)(local.numRows() * local.numCols()) != lm.size())
      dserror("Received matrix of wrong size!");

    // loop over all columns of cal matrix
    for (unsigned icol = 0; icol < local.numCols(); ++icol)
    {
      // loop over all rows of local matrix
      for (unsigned irow = 0; irow < local.numRows(); ++irow)
      {
        // extract local ID of current dof
        const unsigned index = icol * local.numRows() + irow;
        const int lid = global.Map().LID(lm[index]);

        // safety check
        if (lid < 0)
          dserror("Proc %d: Cannot find gid=%d in Epetra_Vector", global.Comm().MyPID(), lm[index]);

        // store current dof in local matrix, which is filled column-wise with the dofs listed in
        // the lm vector
        local(irow, icol) = global[lid];
      }
    }
  }

  /// Locally extract a subset of values from a (column)-nodemap-based Epetra_MultiVector
  /*  \author henke
   *  \date 06/09
   */
  void ExtractMyNodeBasedValues(const DRT::Element* ele,  ///< pointer to current element
      std::vector<double>& local,                         ///< local vector on element-level
      const Epetra_MultiVector& global                    ///< global (multi) vector
  );


  /// Locally extract a subset of values from a (column)-nodemap-based Epetra_MultiVector
  /*  \author g.bau
   *  \date 08/08
   */
  void ExtractMyNodeBasedValues(const DRT::Element* ele,  ///< pointer to current element
      CORE::LINALG::SerialDenseVector& local,             ///< local vector on element-level
      const Teuchos::RCP<Epetra_MultiVector>& global,     ///< global vector
      const int nsd                                       ///< number of space dimensions
  );

  /// Locally extract a subset of values from a (column)-nodemap-based Epetra_MultiVector
  /*  \author schott
   *  \date 12/16
   */
  void ExtractMyNodeBasedValues(const DRT::Node* node,  ///< pointer to current element
      CORE::LINALG::SerialDenseVector& local,           ///< local vector on node-level
      const Teuchos::RCP<Epetra_MultiVector>& global,   ///< global vector
      const int nsd                                     ///< number of space dimensions
  );


  /// Locally extract a subset of values from a (column)-nodemap-based Epetra_MultiVector
  /// and fill a local matrix that has implemented the (.,.) operator
  /*  \author g.bau
   *  \date 04/09
   */
  template <class M>
  void ExtractMyNodeBasedValues(const DRT::Element* ele,  ///< pointer to current element
      M& localmatrix,                                     ///< local matrix on element-level
      const Teuchos::RCP<Epetra_MultiVector>& global,     ///< global vector
      const int nsd                                       ///< number of space dimensions
  )
  {
    if (global == Teuchos::null) dserror("received a TEUCHOS::null pointer");
    if (nsd > global->NumVectors())
      dserror("Requested %d of %d available columns", nsd, global->NumVectors());
    const int iel = ele->NumNode();  // number of nodes
    if (((int)localmatrix.numCols()) != iel) dserror("local matrix has wrong number of columns");
    if (((int)localmatrix.numRows()) != nsd) dserror("local matrix has wrong number of rows");

    for (int i = 0; i < nsd; i++)
    {
      // access actual component column of multi-vector
      double* globalcolumn = (*global)[i];
      // loop over the element nodes
      for (int j = 0; j < iel; j++)
      {
        const int nodegid = (ele->Nodes()[j])->Id();
        const int lid = global->Map().LID(nodegid);
        if (lid < 0)
          dserror(
              "Proc %d: Cannot find gid=%d in Epetra_Vector", (*global).Comm().MyPID(), nodegid);
        localmatrix(i, j) = globalcolumn[lid];
      }
    }
  }

  /*!
  \brief extract local values from global node-based (multi) vector

  This function returns a column vector!

  \author henke
 */
  template <class M>
  void ExtractMyNodeBasedValues(const DRT::Element* ele, M& local, const Epetra_MultiVector& global)
  {
    const int numnode = ele->NumNode();
    const int numcol = global.NumVectors();
    if (((int)local.N()) != 1) dserror("local matrix must have one column");
    if (((int)local.M()) != numnode * numcol) dserror("local matrix has wrong number of rows");

    // loop over element nodes
    for (int i = 0; i < numnode; ++i)
    {
      const int nodegid = (ele->Nodes()[i])->Id();
      const int lid = global.Map().LID(nodegid);
      if (lid < 0)
        dserror("Proc %d: Cannot find gid=%d in Epetra_Vector", global.Comm().MyPID(), nodegid);

      // loop over multi vector columns (numcol=1 for Epetra_Vector)
      for (int col = 0; col < numcol; col++)
      {
        double* globalcolumn = (global)[col];
        local((col + (numcol * i)), 0) = globalcolumn[lid];
      }
    }
  }
}  // namespace CORE::FE


BACI_NAMESPACE_CLOSE

#endif
