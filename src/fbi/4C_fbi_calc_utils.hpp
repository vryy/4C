// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FBI_CALC_UTILS_HPP
#define FOUR_C_FBI_CALC_UTILS_HPP

#include "4C_config.hpp"

#include <memory>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Elements
{
  class Element;
}

namespace Core::LinAlg
{
  class SerialDenseVector;
  class SerialDenseMatrix;
  class SparseOperator;
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace FBI
{
  namespace Utils
  {
    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/

    /**
     * \brief Get the local indices of the centerline DOFs of an element.
     * @param discret (in) Pointer to the fluid or structure discretization.
     * @param ele (in) Pointer to the fluid or beam element.
     * @param ele_centerline_dof_indices (out) Vector with local indices of centerline DOFs in the
     * element.
     * @param num_dof (out) Number total DOFs on the element.
     *
     */
    void get_fbi_element_centerline_dof_indices(Core::FE::Discretization const& discret,
        const Core::Elements::Element* ele, std::vector<unsigned int>& ele_centerline_dof_indices,
        unsigned int& num_dof);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/

    /**
     * \brief Assembles centerline DOF coupling contributions into element coupling matrices and
     * force vectors
     *
     * The function expects a vector object of length 2 containing discretizations
     *
     * \param[in] discretizations vector to the structure and then fluid discretization
     * \param[in] elegid global id for the beam and fluid element
     * \param[in] eleforce_centerlineDOFs force vectors containing contributions from the centerline
     * dofs in the case of beams, otherwise just the classical element DOFs
     * \param[in]
     * elestiffcenterlineDOFs stiffness matrices containing contributions from the centerline dofs
     * in the case of beams, otherwise just the classical element DOFs
     * \param[inout] eleforce
     * element force vector
     * \param[inout] elestiff element stiffness matrix
     *
     */

    void assemble_centerline_dof_force_stiff_into_fbi_element_force_stiff(
        const Core::FE::Discretization& discretization1,
        const Core::FE::Discretization& discretization2, std::vector<int> const& elegid,
        std::vector<Core::LinAlg::SerialDenseVector> const& eleforce_centerlineDOFs,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> const& elestiff_centerlineDOFs,
        std::vector<Core::LinAlg::SerialDenseVector>* eleforce,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>>* elestiff);

  }  // namespace Utils
}  // namespace FBI



FOUR_C_NAMESPACE_CLOSE

#endif
