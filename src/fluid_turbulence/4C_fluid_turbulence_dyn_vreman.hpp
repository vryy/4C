// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_TURBULENCE_DYN_VREMAN_HPP
#define FOUR_C_FLUID_TURBULENCE_DYN_VREMAN_HPP


#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Teuchos_TimeMonitor.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN



namespace FLD
{
  class Boxfilter;

  class Vreman
  {
   public:
    /*!
    \brief Standard Constructor (public)

    */
    Vreman(std::shared_ptr<Core::FE::Discretization> actdis, Teuchos::ParameterList& params);

    /*!
    \brief Destructor

    */
    virtual ~Vreman() = default;



    void apply_filter_for_dynamic_computation_of_cv(
        const std::shared_ptr<const Core::LinAlg::Vector<double>> velocity,
        const std::shared_ptr<const Core::LinAlg::Vector<double>> scalar, const double thermpress,
        const std::shared_ptr<const Core::LinAlg::Vector<double>> dirichtoggle);

    void apply_filter_for_dynamic_computation_of_dt(
        const std::shared_ptr<const Core::LinAlg::Vector<double>> scalar, const double thermpress,
        const std::shared_ptr<const Core::LinAlg::Vector<double>> dirichtoggle,
        Teuchos::ParameterList& extraparams, const int ndsvel);

    void add_scatra(std::shared_ptr<Core::FE::Discretization> scatradis);

    void get_cv(double& Cv)
    {
      Cv = Cv_;
      return;
    }

    double Cv_;

   private:
    /// provide access to the box filter
    std::shared_ptr<FLD::Boxfilter> boxfilter();
    // Boxfilter
    std::shared_ptr<FLD::Boxfilter> boxf_;
    std::shared_ptr<FLD::Boxfilter> boxfsc_;



    //! the discretization
    std::shared_ptr<Core::FE::Discretization> discret_;
    //! parameterlist including time params, stabilization params and turbulence sublist
    Teuchos::ParameterList& params_;
    //! flag for physical type of fluid flow
    Inpar::FLUID::PhysicalType physicaltype_;
    // scatra specific
    std::shared_ptr<Core::FE::Discretization> scatradiscret_;


    double dyn_vreman_compute_cv();

    void dyn_vreman_compute_dt(Teuchos::ParameterList& extraparams);

    //@}

    //! @name vectors used for filtering (for dynamic Smagorinsky model)
    //        --------------------------

    //! the filtered reystress exported to column map
    std::shared_ptr<Core::LinAlg::MultiVector<double>> col_filtered_strainrate_;
    std::shared_ptr<Core::LinAlg::Vector<double>> col_filtered_expression_;
    std::shared_ptr<Core::LinAlg::MultiVector<double>> col_filtered_alphaij_;
    std::shared_ptr<Core::LinAlg::Vector<double>> col_filtered_alpha2_;
    std::shared_ptr<Core::LinAlg::MultiVector<double>> col_filtered_phi_;
    std::shared_ptr<Core::LinAlg::Vector<double>> col_filtered_phi2_;
    std::shared_ptr<Core::LinAlg::Vector<double>> col_filtered_phiexpression_;
    std::shared_ptr<Core::LinAlg::MultiVector<double>> col_filtered_alphaijsc_;
    //@}

    //! @name homogeneous flow specials
    //        -------------------------------

    //@}

  };  // end class Vreman

}  // end namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
