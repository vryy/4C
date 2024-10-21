// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_TIMINT_LOMA_BDF2_HPP
#define FOUR_C_FLUID_TIMINT_LOMA_BDF2_HPP


#include "4C_config.hpp"

#include "4C_fluid_timint_bdf2.hpp"
#include "4C_fluid_timint_loma.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TimIntLomaBDF2 : public TimIntBDF2, public TimIntLoma
  {
   public:
    /// Standard Constructor
    TimIntLomaBDF2(const Teuchos::RCP<Core::FE::Discretization>& actdis,
        const Teuchos::RCP<Core::LinAlg::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool alefluid = false);


    /*!
    \brief initialization

    */
    void init() override;


   protected:
   private:
  };  // class TimIntLomaBDF2

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
