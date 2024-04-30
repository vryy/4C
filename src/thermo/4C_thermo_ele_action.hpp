/*----------------------------------------------------------------------*/
/*! \file

 \brief enums defining action of the thermo elements and related helper functions

\level 1

 *----------------------------------------------------------------------*/


#ifndef FOUR_C_THERMO_ELE_ACTION_HPP
#define FOUR_C_THERMO_ELE_ACTION_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

#include <iostream>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace THR
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
    calc_normal_vectors,
    ba_integrate_shape_functions
  };

  /*!
   * \brief translate to string for screen output
   */
  inline std::string ActionToString(const Action action)
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

  inline std::string BoundaryActionToString(const BoundaryAction baction)
  {
    switch (baction)
    {
      case ba_none:
        return "ba_none";
      case calc_thermo_fextconvection:
        return "calc_thermo_fextconvection";
      case calc_thermo_fextconvection_coupltang:
        return "calc_thermo_fextconvection_coupltang";
      case calc_normal_vectors:
        return "calc_normal_vectors";
      case ba_integrate_shape_functions:
        return "ba_integrate_shape_functions";
      default:
        FOUR_C_THROW("no string for the boundary action %d defined!", baction);
    };
  }

}  // namespace THR

FOUR_C_NAMESPACE_CLOSE

#endif
