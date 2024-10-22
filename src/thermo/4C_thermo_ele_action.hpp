// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_THERMO_ELE_ACTION_HPP
#define FOUR_C_THERMO_ELE_ACTION_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

#include <iostream>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace Thermo
{
  /*--------------------------------------------------------------------------
   | enum that provides all possible thermo actions
   *--------------------------------------------------------------------------*/
  enum Action
  {
    none,
    calc_thermo_fint,
    calc_thermo_fintcapa,
    calc_thermo_finttang,
    calc_thermo_heatflux,
    postproc_thermo_heatflux,
    integrate_shape_functions,
    calc_thermo_update_istep,
    calc_thermo_reset_istep,
    calc_thermo_energy,
    calc_thermo_coupltang,
    calc_thermo_fintcond,
    calc_thermo_fext,
    calc_thermo_error,
  };  // enum Action

  /*--------------------------------------------------------------------------
   | enum that provides all possible thermo actions on a boundary
   *--------------------------------------------------------------------------*/
  enum BoundaryAction
  {
    ba_none,
    calc_thermo_fextconvection,
    calc_thermo_fextconvection_coupltang,
  };

  /*!
   * \brief translate to string for screen output
   */
  inline std::string action_to_string(const Action action)
  {
    switch (action)
    {
      case none:
        return "none";
      case calc_thermo_fint:
        return "calc_thermo_fint";
      case calc_thermo_fintcapa:
        return "calc_thermo_fintcapa";
      case calc_thermo_finttang:
        return "calc_thermo_finttang";
      case calc_thermo_heatflux:
        return "calc_thermo_heatflux";
      case postproc_thermo_heatflux:
        return "postproc_thermo_heatflux";
      case integrate_shape_functions:
        return "integrate_shape_functions";
      case calc_thermo_update_istep:
        return "calc_thermo_update_istep";
      case calc_thermo_reset_istep:
        return "calc_thermo_reset_istep";
      case calc_thermo_energy:
        return "calc_thermo_energy";
      case calc_thermo_coupltang:
        return "calc_thermo_coupltang";
      case calc_thermo_fintcond:
        return "calc_thermo_fintcond";
      default:
        FOUR_C_THROW("no string for action %d defined!", action);
    };
  }

  inline std::string boundary_action_to_string(const BoundaryAction baction)
  {
    switch (baction)
    {
      case ba_none:
        return "ba_none";
      case calc_thermo_fextconvection:
        return "calc_thermo_fextconvection";
      case calc_thermo_fextconvection_coupltang:
        return "calc_thermo_fextconvection_coupltang";
      default:
        FOUR_C_THROW("no string for the boundary action %d defined!", baction);
    };
  }

}  // namespace Thermo

FOUR_C_NAMESPACE_CLOSE

#endif
