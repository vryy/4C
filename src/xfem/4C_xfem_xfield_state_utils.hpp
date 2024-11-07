// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_XFEM_XFIELD_STATE_UTILS_HPP
#define FOUR_C_XFEM_XFIELD_STATE_UTILS_HPP


#include "4C_config.hpp"

#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_solver_nonlin_nox_aux.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN


namespace XFEM
{
  /** \brief Destroy the Core::LinAlg::SparseOperator object and it's date
   *
   *  \author hiermeier
   *  \date 07/16 */
  inline void destroy_matrix(
      std::shared_ptr<Core::LinAlg::SparseOperator>& mat, bool throw_exception = true)
  {
    mat = nullptr;
  }

  /** \brief Destroy the Core::LinAlg::SparseMatrix object and it's date
   *
   *  \author schott
   *  \date 01/15 */
  inline void destroy_matrix(
      std::shared_ptr<Core::LinAlg::SparseMatrix>& mat, bool throw_exception = true)
  {
    mat = nullptr;
  }


  /** \brief Destroy the reference counted object and the reference counter
   *
   *  \author schott
   *  \date 01/15 */
  template <class OBJECT>
  inline void destroy_rcp_object(std::shared_ptr<OBJECT>& obj_rcp, bool throw_exception = true)
  {
    obj_rcp = nullptr;
  }


  /** \brief More efficient and memory safe Zero routine for system matrix
   *
   *  \author schott
   *  \date 01/15 */
  inline void zero_matrix(Core::LinAlg::SparseMatrix& mat)
  {
    if (mat.explicit_dirichlet())
    {
      mat.zero();  // matrix could have been changed due to Dirichlet conditions, go back to
                   // original Graph if savegraph == true
    }
    else
    {
      // do not create a new matrix via Zero() but zero entries
      mat.put_scalar(0.0);
    }
  }

}  // namespace XFEM

FOUR_C_NAMESPACE_CLOSE

#endif
