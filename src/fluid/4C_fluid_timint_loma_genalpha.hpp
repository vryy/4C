// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_TIMINT_LOMA_GENALPHA_HPP
#define FOUR_C_FLUID_TIMINT_LOMA_GENALPHA_HPP


#include "4C_config.hpp"

#include "4C_fluid_timint_genalpha.hpp"
#include "4C_fluid_timint_loma.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TimIntLomaGenAlpha : public TimIntGenAlpha, public TimIntLoma
  {
   public:
    /// Standard Constructor
    TimIntLomaGenAlpha(const std::shared_ptr<Core::FE::Discretization>& actdis,
        const std::shared_ptr<Core::LinAlg::Solver>& solver,
        const std::shared_ptr<Teuchos::ParameterList>& params,
        const std::shared_ptr<Core::IO::DiscretizationWriter>& output, bool alefluid = false);


    /*!
    \brief initialization

    */
    void init() override;


   protected:
   private:
  };  // class TimIntLomaGenAlpha

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
