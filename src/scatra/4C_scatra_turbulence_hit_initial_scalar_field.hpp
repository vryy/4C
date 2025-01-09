// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_TURBULENCE_HIT_INITIAL_SCALAR_FIELD_HPP
#define FOUR_C_SCATRA_TURBULENCE_HIT_INITIAL_SCALAR_FIELD_HPP

#include "4C_config.hpp"

#include "4C_inpar_scatra.hpp"
#include "4C_linalg_vector.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace ScaTra
{
  // forward declarations
  class ScaTraTimIntImpl;

  // initial condition for homogeneous isotropic turbulence
  // based on the Comte-Bellot - Corrsion experiment
  class HomoIsoTurbInitialScalarField
  {
   public:
    //! constructor
    HomoIsoTurbInitialScalarField(
        ScaTraTimIntImpl& timeint, const Inpar::ScaTra::InitialField initfield);

    //! calculate initial field
    void calculate_initial_field();

   protected:
    //! sort criterium for double values up to a tolerance of 10-9
    class LineSortCriterion
    {
     public:
      bool operator()(const double& p1, const double& p2) const { return (p1 < p2 - 1E-9); }

     protected:
     private:
    };

   private:
    //! estimate energy form given scalar variance spectrum (function for E_phi)
    double calculate_energy_from_spectrum(double k);

    //! scatra discretization
    std::shared_ptr<Core::FE::Discretization> discret_;

    //! state vectors to be initialized
    std::shared_ptr<Core::LinAlg::Vector<double>> phinp_;
    std::shared_ptr<Core::LinAlg::Vector<double>> phin_;

    //! type of energy spectrum for initialization
    Inpar::ScaTra::InitialField type_;

    //! number of resolved mode
    int nummodes_;

    //! vector of coordinates in one spatial direction (same for the other two directions)
    std::shared_ptr<std::vector<double>> coordinates_;
  };

}  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
