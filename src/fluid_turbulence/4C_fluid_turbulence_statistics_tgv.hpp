// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_TURBULENCE_STATISTICS_TGV_HPP
#define FOUR_C_FLUID_TURBULENCE_STATISTICS_TGV_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_ParameterList.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TurbulenceStatisticsTgv
  {
   public:
    /*!
    \brief Standard Constructor (public)
    */
    TurbulenceStatisticsTgv(std::shared_ptr<Core::FE::Discretization> actdis,
        Teuchos::ParameterList& params, const std::string& statistics_outfilename);

    /*!
    \brief Destructor

    */
    virtual ~TurbulenceStatisticsTgv() = default;

    //! @name functions for (spatial) averaging

    /*!
    \brief evaluate dissipation
    */
    void evaluate_residuals(
        std::map<std::string, std::shared_ptr<Core::LinAlg::Vector<double>>> statevecs,
        std::map<std::string, std::shared_ptr<Core::LinAlg::MultiVector<double>>> statetenss,
        const double thermpressaf, const double thermpressam, const double thermpressdtaf,
        const double thermpressdtam);

    //@}

    //! @name Miscellaneous

    /*!
    \brief Dump the result to file.
    */

    void dump_statistics(const int step);

    /*!
    \brief reset values
    */

    void clear_statistics();


   private:
    //! time step size
    double dt_;

    //! number of elements in sample plane
    int numele_;

    //! The discretization (required for nodes, dofs etc;)
    std::shared_ptr<Core::FE::Discretization> discret_;

    //! contains plane normal direction etc --- this is the original
    //! fluid dynamic parameterlist
    Teuchos::ParameterList& params_;
    //! name of statistics output file, despite the ending
    const std::string statistics_outfilename_;

    //! parameterlist for the element call when averages of residuals
    //! are calculated --- used for communication between element
    //! and averaging methods --- for fluid field
    Teuchos::ParameterList eleparams_;

    //! the dim_-coordinates of the homogeneous planes containing nodes
    std::shared_ptr<std::vector<double>> nodeplanes_;

    //!--------------------------------------------------
    //!  averaged resdiuals and subscale quantities
    //!--------------------------------------------------

    //! sum over all in plane element sizes
    std::shared_ptr<std::vector<double>> sumhk_;
    //! sum over all in plane element sizes for the Bazilevs
    //! parameter in the viscous regime
    std::shared_ptr<std::vector<double>> sumhbazilevs_;
    //! sum over all in plane stream lengths
    std::shared_ptr<std::vector<double>> sumstrle_;
    //! sum over all in plane stream lengths
    std::shared_ptr<std::vector<double>> sumgradle_;

    //! sum over all in plane residuals
    std::shared_ptr<std::vector<double>> sumtau_m_;
    //! sum over all in plane squared residuals
    std::shared_ptr<std::vector<double>> sumtau_c_;

    //! sum over all in plane mk (parameter for stabilization parameter, 1/3 for lin ele)
    std::shared_ptr<std::vector<double>> summk_;

    //! sum over all in plane residuals
    std::shared_ptr<std::vector<double>> sumres_;
    //! sum over all in plane squared residuals
    std::shared_ptr<std::vector<double>> sumres_sq_;
    //! sum over all in plane residuals norms
    std::shared_ptr<std::vector<double>> sumabsres_;
    //! sum over all in plane values of svel/tau
    std::shared_ptr<std::vector<double>> sumtauinvsvel_;
    //! sum over all in plane subscale velocities
    std::shared_ptr<std::vector<double>> sumsvelaf_;
    //! sum over all in plane squared subscale velocities
    std::shared_ptr<std::vector<double>> sumsvelaf_sq_;
    //! sum over all in plane subscale velocities norms
    std::shared_ptr<std::vector<double>> sumabssvelaf_;

    //! sum over all in plane residuals of the continuity equation
    std::shared_ptr<std::vector<double>> sumres_c_;
    //! sum over all in plane squared residuals of the continuity equation
    std::shared_ptr<std::vector<double>> sumres_c_sq_;
    //! sum over all in plane subscale pressure values at current timestep
    std::shared_ptr<std::vector<double>> sumspressnp_;
    //! sum over all in plane squared subscale pressure values at current timestep
    std::shared_ptr<std::vector<double>> sumspressnp_sq_;

    //! sum over all in plane averaged dissipation rates from pspg stabilisation
    std::shared_ptr<std::vector<double>> sum_eps_pspg_;
    //! sum over all in plane averaged dissipation rates from supg stabilisation
    std::shared_ptr<std::vector<double>> sum_eps_supg_;
    //! sum over all in plane averaged dissipation rates from cross term
    std::shared_ptr<std::vector<double>> sum_eps_cross_;
    //! sum over all in plane averaged dissipation rates from reynolds term
    std::shared_ptr<std::vector<double>> sum_eps_rey_;
    //! sum over all in plane averaged dissipation rates from least squares continuity term
    std::shared_ptr<std::vector<double>> sum_eps_graddiv_;
    //! sum over all in plane averaged dissipation rates from eddy viscosity model (Smagorinsky)
    std::shared_ptr<std::vector<double>> sum_eps_eddyvisc_;
    //! sum over all in plane averaged dissipation rates from eddy viscosity model (AVM3)
    std::shared_ptr<std::vector<double>> sum_eps_avm3_;
    //! sum over all in plane averaged dissipation rates from mfs model
    std::shared_ptr<std::vector<double>> sum_eps_mfs_;
    //! sum over all in plane averaged dissipation rates from ,fs model (forwardscatter)
    std::shared_ptr<std::vector<double>> sum_eps_mfscross_;
    //! sum over all in plane averaged dissipation rates from mfs model (backscatter)
    std::shared_ptr<std::vector<double>> sum_eps_mfsrey_;
    //! sum over all in plane averaged dissipation rates from Galerkin viscous term
    std::shared_ptr<std::vector<double>> sum_eps_visc_;
    //! sum over all in plane averaged dissipation rates from Galerkin convective term
    std::shared_ptr<std::vector<double>> sum_eps_conv_;

    //! sum over all in plane averaged supg+cross stress
    std::shared_ptr<std::vector<double>> sum_crossstress_;
    //! sum over all in plane averaged reynolds stress
    std::shared_ptr<std::vector<double>> sum_reystress_;
  };

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
