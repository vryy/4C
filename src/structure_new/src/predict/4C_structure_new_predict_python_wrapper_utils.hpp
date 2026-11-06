// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_PREDICT_PYTHON_WRAPPER_UTILS_HPP
#define FOUR_C_STRUCTURE_NEW_PREDICT_PYTHON_WRAPPER_UTILS_HPP

#ifdef FOUR_C_WITH_PYBIND11

#include "4C_config.hpp"

#include <pybind11/embed.h>
#include <pybind11/numpy.h>

FOUR_C_NAMESPACE_OPEN

namespace Solid::Predict::PythonWrapperUtils
{

  //! Convenience struct for encapsulating the mesh information in one place.
  struct FOUR_C_API(HIDDEN) PythonMesh
  {
    pybind11::array_t<int> node_ids;
    pybind11::array_t<double> coordinates;
    pybind11::array_t<int> disp_dof_ids;
    pybind11::array_t<int> element_ids;
    pybind11::list connectivity;
  };

  /*!
   * \brief Convenience struct for encapsulating a lightweight mesh and some additional information
   * in one place.
   */
  struct FOUR_C_API(HIDDEN) PythonContext
  {
    std::string input_file;
    unsigned int problem_dim;
    int step_max;
    double time_max;
    PythonMesh mesh;
  };

  /*!
   * \brief Convenience struct for encapsulating the global state vectors and some additional
   * information in one place.
   */
  struct FOUR_C_API(HIDDEN) PythonState
  {
    int step;
    double time;
    double dt;
    int num_global_dofs;

    pybind11::array_t<const double> dis_n;
    pybind11::array_t<double> dis_np;
    pybind11::array_t<const double> vel_n;
    pybind11::array_t<double> vel_np;
    pybind11::array_t<const double> acc_n;
    pybind11::array_t<double> acc_np;
  };

  /*!
   * \brief Ensures that the Python wrapper module providing the helper structs PythonMesh,
   * PythonContext and PythonState is initialized.
   * \return A reference to the Python wrapper module.
   */
  pybind11::module_ ensure_python_wrapper_module();

}  // namespace Solid::Predict::PythonWrapperUtils

FOUR_C_NAMESPACE_CLOSE

#endif

#endif