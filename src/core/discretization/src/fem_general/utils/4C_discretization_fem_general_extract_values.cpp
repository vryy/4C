/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 0


*/
/*---------------------------------------------------------------------*/

#include "4C_discretization_fem_general_extract_values.hpp"

#include "4C_io_control.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_FEVector.h>

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::FE::ExtractMyValues(
    const Epetra_Vector& global, std::vector<double>& local, const std::vector<int>& lm)
{
  const size_t ldim = lm.size();
  local.resize(ldim);
  for (size_t i = 0; i < ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid < 0)
      FOUR_C_THROW("Proc %d: Cannot find gid=%d in Epetra_Vector", global.Comm().MyPID(), lm[i]);
    local[i] = global[lid];
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::FE::ExtractMyValues(
    const Epetra_Vector& global, CORE::LINALG::SerialDenseVector& local, const std::vector<int>& lm)
{
  const size_t ldim = lm.size();
  local.size(ldim);
  for (size_t i = 0; i < ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid < 0)
      FOUR_C_THROW("Proc %d: Cannot find gid=%d in Epetra_Vector", global.Comm().MyPID(), lm[i]);
    local[i] = global[lid];
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::FE::ExtractMyValues(
    const Epetra_MultiVector& global, std::vector<double>& local, const std::vector<int>& lm)
{
  const int numcol = global.NumVectors();
  const size_t ldim = lm.size();

  local.resize(ldim * numcol);

  // loop over element nodes
  for (size_t i = 0; i < ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid < 0)
      FOUR_C_THROW(
          "Proc %d: Cannot find gid=%d in Epetra_MultiVector", global.Comm().MyPID(), lm[i]);

    // loop over multi vector columns (numcol=1 for Epetra_Vector)
    for (int col = 0; col < numcol; col++)
    {
      double* globalcolumn = (global)[col];
      local[col + (numcol * i)] = globalcolumn[lid];
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::FE::ExtractMyNodeBasedValues(
    const DRT::Element* ele, std::vector<double>& local, const Epetra_MultiVector& global)
{
  const int numnode = ele->NumNode();
  const int numcol = global.NumVectors();
  local.resize(numnode * numcol);

  // loop over element nodes
  for (int i = 0; i < numnode; ++i)
  {
    const int nodegid = (ele->Nodes()[i])->Id();
    const int lid = global.Map().LID(nodegid);
    if (lid < 0)
      FOUR_C_THROW("Proc %d: Cannot find gid=%d in Epetra_Vector", global.Comm().MyPID(), nodegid);

    // loop over multi vector columns (numcol=1 for Epetra_Vector)
    for (int col = 0; col < numcol; col++)
    {
      double* globalcolumn = (global)[col];
      local[col + (numcol * i)] = globalcolumn[lid];
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::FE::ExtractMyNodeBasedValues(const DRT::Element* ele,
    CORE::LINALG::SerialDenseVector& local, const Teuchos::RCP<Epetra_MultiVector>& global,
    const int nsd)
{
  if (global == Teuchos::null) FOUR_C_THROW("received a TEUCHOS::null pointer");
  if (nsd > global->NumVectors())
    FOUR_C_THROW("Requested %d of %d available columns", nsd, global->NumVectors());
  const int iel = ele->NumNode();  // number of nodes
  if (local.length() != (iel * nsd)) FOUR_C_THROW("vector size mismatch.");

  // TODO: might we do change the loops?
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
        FOUR_C_THROW(
            "Proc %d: Cannot find gid=%d in Epetra_MultiVector", global->Comm().MyPID(), nodegid);
      local(i + (nsd * j)) = globalcolumn[lid];
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::FE::ExtractMyNodeBasedValues(const DRT::Node* node,
    CORE::LINALG::SerialDenseVector& local, const Teuchos::RCP<Epetra_MultiVector>& global,
    const int nsd)
{
  if (global == Teuchos::null) FOUR_C_THROW("received a TEUCHOS::null pointer");
  if (nsd > global->NumVectors())
    FOUR_C_THROW("Requested %d of %d available columns", nsd, global->NumVectors());
  if (local.length() != nsd) FOUR_C_THROW("vector size mismatch.");

  const int nodegid = node->Id();
  const int lid = global->Map().LID(nodegid);

  for (int i = 0; i < nsd; i++)
  {
    // access actual component column of multi-vector
    double* globalcolumn = (*global)[i];

    local(i + nsd) = globalcolumn[lid];
  }
  return;
}

FOUR_C_NAMESPACE_CLOSE
