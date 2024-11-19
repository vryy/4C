// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_general_extract_values.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_FEVector.h>

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::FE::extract_my_values(const Core::LinAlg::Vector<double>& global,
    std::vector<double>& local, const std::vector<int>& lm)
{
  const size_t ldim = lm.size();
  local.resize(ldim);
  for (size_t i = 0; i < ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid < 0)
      FOUR_C_THROW("Proc %d: Cannot find gid=%d in Core::LinAlg::Vector<double>",
          Core::Communication::my_mpi_rank(global.Comm()), lm[i]);
    local[i] = global[lid];
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::FE::extract_my_values(const Core::LinAlg::Vector<double>& global,
    Core::LinAlg::SerialDenseVector& local, const std::vector<int>& lm)
{
  const size_t ldim = lm.size();
  local.size(ldim);
  for (size_t i = 0; i < ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid < 0)
      FOUR_C_THROW("Proc %d: Cannot find gid=%d in Core::LinAlg::Vector<double>",
          Core::Communication::my_mpi_rank(global.Comm()), lm[i]);
    local[i] = global[lid];
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::FE::extract_my_values(const Core::LinAlg::MultiVector<double>& global,
    std::vector<double>& local, const std::vector<int>& lm)
{
  const int numcol = global.NumVectors();
  const size_t ldim = lm.size();

  local.resize(ldim * numcol);

  // loop over element nodes
  for (size_t i = 0; i < ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid < 0)
      FOUR_C_THROW("Proc %d: Cannot find gid=%d in Core::LinAlg::MultiVector<double>",
          Core::Communication::my_mpi_rank(global.Comm()), lm[i]);

    // loop over multi vector columns (numcol=1 for Core::LinAlg::Vector<double>)
    for (int col = 0; col < numcol; col++)
    {
      local[col + (numcol * i)] = global(col)[lid];
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::FE::extract_my_node_based_values(const Core::Elements::Element* ele,
    std::vector<double>& local, const Core::LinAlg::MultiVector<double>& global)
{
  const int numnode = ele->num_node();
  const int numcol = global.NumVectors();
  local.resize(numnode * numcol);

  // loop over element nodes
  for (int i = 0; i < numnode; ++i)
  {
    const int nodegid = (ele->nodes()[i])->id();
    const int lid = global.Map().LID(nodegid);
    if (lid < 0)
      FOUR_C_THROW("Proc %d: Cannot find gid=%d in Core::LinAlg::Vector<double>",
          Core::Communication::my_mpi_rank(global.Comm()), nodegid);

    // loop over multi vector columns (numcol=1 for Core::LinAlg::Vector<double>)
    for (int col = 0; col < numcol; col++)
    {
      local[col + (numcol * i)] = global(col)[lid];
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::FE::extract_my_node_based_values(const Core::Elements::Element* ele,
    Core::LinAlg::SerialDenseVector& local, Core::LinAlg::MultiVector<double>& global,
    const int nsd)
{
  if (nsd > global.NumVectors())
    FOUR_C_THROW("Requested %d of %d available columns", nsd, global.NumVectors());
  const int iel = ele->num_node();  // number of nodes
  if (local.length() != (iel * nsd)) FOUR_C_THROW("vector size mismatch.");

  // TODO: might we do change the loops?
  for (int i = 0; i < nsd; i++)
  {
    // loop over the element nodes
    for (int j = 0; j < iel; j++)
    {
      const int nodegid = (ele->nodes()[j])->id();
      const int lid = global.Map().LID(nodegid);
      if (lid < 0)
        FOUR_C_THROW("Proc %d: Cannot find gid=%d in Core::LinAlg::MultiVector<double>",
            Core::Communication::my_mpi_rank(global.Comm()), nodegid);
      local(i + (nsd * j)) = global(i)[lid];
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::FE::extract_my_node_based_values(const Core::Nodes::Node* node,
    Core::LinAlg::SerialDenseVector& local, Core::LinAlg::MultiVector<double>& global,
    const int nsd)
{
  if (nsd > global.NumVectors())
    FOUR_C_THROW("Requested %d of %d available columns", nsd, global.NumVectors());
  if (local.length() != nsd) FOUR_C_THROW("vector size mismatch.");

  const int nodegid = node->id();
  const int lid = global.Map().LID(nodegid);

  for (int i = 0; i < nsd; i++)
  {
    // access actual component column of multi-vector
    local(i + nsd) = global(i)[lid];
  }
  return;
}

FOUR_C_NAMESPACE_CLOSE
