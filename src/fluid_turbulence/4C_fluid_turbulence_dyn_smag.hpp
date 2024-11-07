// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_TURBULENCE_DYN_SMAG_HPP
#define FOUR_C_FLUID_TURBULENCE_DYN_SMAG_HPP


#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Teuchos_TimeMonitor.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN



namespace FLD
{
  class Boxfilter;

  class DynSmagFilter
  {
   public:
    /*!
    \brief Standard Constructor (public)

    */
    DynSmagFilter(std::shared_ptr<Core::FE::Discretization> actdis, Teuchos::ParameterList& params);

    /*!
    \brief Destructor

    */
    virtual ~DynSmagFilter() = default;

    void add_scatra(std::shared_ptr<Core::FE::Discretization> scatradis);

    /*!
    \brief Perform box filter operation, compare filtered quantities
    to solution to get an estimate for Cs (using clpping), average
    over element layers in turbulent channel flows.

    This method initialises element quantities (standard case) or
    provides information for the element via the parameter list
    (in plane averaging for channel flow)

    \param solution     (in) velocity field to filter and to
                             determine Cs from
    \param dirichtoggle (in) information on dirichlet dofs to be
                             able to exclude boundary nodes from
                             filtering

    */
    void apply_filter_for_dynamic_computation_of_cs(
        const std::shared_ptr<const Core::LinAlg::Vector<double>> velocity,
        const std::shared_ptr<const Core::LinAlg::Vector<double>> scalar, const double thermpress,
        const std::shared_ptr<const Core::LinAlg::Vector<double>> dirichtoggle);

    void apply_filter_for_dynamic_computation_of_prt(
        const std::shared_ptr<const Core::LinAlg::Vector<double>> scalar, const double thermpress,
        const std::shared_ptr<const Core::LinAlg::Vector<double>> dirichtoggle,
        Teuchos::ParameterList& extraparams, const int ndsvel);


    /*!
    \brief Output of averaged velocity vector for paraview IO

    \param outvec  (in/out) vector in dofrowmap-format to use for
                            output of averaged solution


    */



   private:
    /*!
    \brief perform box filtering in five steps

    1) Integrate element Heaviside functions against the quantities
       which are filtered. Add the result to the nodevectors
       (we get a contribution for every node of the element)
       This is an element call!
    2) send/add values from slaves to masters
    3) zero out dirichlet nodes
    4) do normalization by division by the patchvolume
       (Heaviside function -> box filter function)
    5) Communication part: Export filtered quantities from row to
       column map

       \param velocity     (i) the velocity defining the
                               unfiltered quantities
       \param dirichtoggle (i) specifying which nodes have to be
                               set to zero

    */

    /// provide access to the box filter
    std::shared_ptr<FLD::Boxfilter> boxfilter();
    /*!
    \brief Compute Cs using the filtered quantities.
    This is an element call!

    For a turbulent channel flow, the averaging of the Smagorinsky
    constant is done in here.
    */
    void dyn_smag_compute_cs();

    /*!
    \brief Compute Prt using the filtered quantities.
    This is an element call!

    For a turbulent channel flow, the averaging of the turbulent
    Prandtl number is done in here.
    */
    void dyn_smag_compute_prt(Teuchos::ParameterList& extraparams, int& numele_layer);


    //! @name input arguments of the constructor
    //

    // Boxfilter
    std::shared_ptr<FLD::Boxfilter> boxf_;
    std::shared_ptr<FLD::Boxfilter> boxfsc_;

    //! the discretization
    std::shared_ptr<Core::FE::Discretization> discret_;
    //! parameterlist including time params, stabilization params and turbulence sublist
    Teuchos::ParameterList& params_;
    //! flag for physical type of fluid flow
    Inpar::FLUID::PhysicalType physicaltype_;
    //@}

    //! @name control parameters
    bool homdir_;
    std::string special_flow_homdir_;
    bool apply_dynamic_smagorinsky_;
    bool calc_ci_;
    //@}

    //! @name special scatra variables
    //! the discretization
    std::shared_ptr<Core::FE::Discretization> scatradiscret_;
    //@}

    //! @name vectors used for filtering (for dynamic Smagorinsky model)
    //        --------------------------

    //! the filtered vel exported to column map
    std::shared_ptr<Core::LinAlg::MultiVector<double>> col_filtered_vel_;
    //! the filtered reystress exported to column map
    std::shared_ptr<Core::LinAlg::MultiVector<double>> col_filtered_reynoldsstress_;
    //! the modeled subgrid stresses exported to column map
    std::shared_ptr<Core::LinAlg::MultiVector<double>> col_filtered_modeled_subgrid_stress_;
    //! the filtered velocities times rho exported to column map
    std::shared_ptr<Core::LinAlg::MultiVector<double>> col_filtered_dens_vel_;
    //! the filtered density exported to column map
    std::shared_ptr<Core::LinAlg::Vector<double>> col_filtered_dens_;
    //! the filtered strainrate times rho exported to column map
    std::shared_ptr<Core::LinAlg::Vector<double>> col_filtered_dens_strainrate_;
    //! the modeled fine scale velocities exported to column map
    std::shared_ptr<Core::LinAlg::MultiVector<double>> col_fs_vel_;
    //! the filtered density times temperature times velocity exported to column map (scalar)
    std::shared_ptr<Core::LinAlg::MultiVector<double>> col_filtered_dens_vel_temp_;
    //! the filtered density times temperature gradient times rate of strain exported to column map
    //! (scalar)
    std::shared_ptr<Core::LinAlg::MultiVector<double>> col_filtered_dens_rateofstrain_temp_;
    //  //! the filtered temperature gradient exported to column map (scalar)
    //  std::shared_ptr<Core::LinAlg::MultiVector<double>>      col_filtered_gradtemp_;
    //! the filtered temperature exported to column map (scalar)
    std::shared_ptr<Core::LinAlg::Vector<double>> col_filtered_temp_;
    //! the filtered density times temperature exported to column map (scalar)
    std::shared_ptr<Core::LinAlg::Vector<double>> col_filtered_dens_temp_;
    //@}

    //! @name homogeneous flow specials
    //        -------------------------------

    //! the direction coordinates for the above mentioned averaging procedure
    std::shared_ptr<std::vector<double>> dir1coords_;
    std::shared_ptr<std::vector<double>> dir2coords_;
    //@}

  };  // end class DynSmagFilter

}  // end namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
