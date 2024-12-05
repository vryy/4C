// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_EHL_UTILS_HPP
#define FOUR_C_EHL_UTILS_HPP


#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

#include <mpi.h>


FOUR_C_NAMESPACE_OPEN

namespace EHL
{
  namespace Utils
  {
    /// Function for checking that the different time steps are a
    /// multiplicative of each other
    int check_time_stepping(double dt1, double dt2);

    // Modification of time parameter list for problem with different time step size
    void change_time_parameter(MPI_Comm comm, Teuchos::ParameterList& ehlparams,
        Teuchos::ParameterList& lubricationdyn, Teuchos::ParameterList& sdyn);

  };  // namespace Utils

  //! prints the 4C EHL-logo on the screen
  void printlogo();

}  // namespace EHL


FOUR_C_NAMESPACE_CLOSE

#endif
