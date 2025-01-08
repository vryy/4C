// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_PARAMETER_TIMINT_HPP
#define FOUR_C_SCATRA_ELE_PARAMETER_TIMINT_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_parameter_base.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace Elements
  {
    /// Evaluation of general parameters (constant over time)
    class ScaTraEleParameterTimInt : public ScaTraEleParameterBase
    {
     public:
      //! singleton access method
      static ScaTraEleParameterTimInt* instance(const std::string& disname);

      //! set parameters
      void set_parameters(Teuchos::ParameterList& parameters) override;

      bool is_gen_alpha() const { return is_genalpha_; };
      bool is_stationary() const { return is_stationary_; };
      bool is_incremental() const { return is_incremental_; };
      double time() const { return time_; };
      double time_derivative_fac() const { return timederivativefac_; }
      double dt() const { return dt_; };
      double time_fac() const { return timefac_; };
      double time_fac_rhs() const { return timefacrhs_; };
      double time_fac_rhs_tau() const { return timefacrhstau_; };
      double alpha_f() const { return alpha_f_; };

     private:
      //! private constructor for singletons
      ScaTraEleParameterTimInt();

      bool is_genalpha_;
      bool is_stationary_;
      bool is_incremental_;
      double time_;
      double timederivativefac_;
      double dt_;
      double timefac_;
      double timefacrhs_;
      double timefacrhstau_;
      double alpha_f_;
    };  // class ScaTraEleParameterTimInt
  }  // namespace Elements
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
