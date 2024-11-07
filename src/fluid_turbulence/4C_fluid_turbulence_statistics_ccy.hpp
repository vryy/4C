// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_TURBULENCE_STATISTICS_CCY_HPP
#define FOUR_C_FLUID_TURBULENCE_STATISTICS_CCY_HPP


#include "4C_config.hpp"

#include "4C_linalg_serialdensematrix.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_ParameterList.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  template <typename T>
  class Vector;
}

namespace FLD
{
  class TurbulenceStatisticsCcy
  {
   public:
    /*!
    \brief Standard Constructor (public)

        o Create vector of radial coordinates of homogeneous shells

    o allocate all sum_something vectors

    o initialise the output (open/clear files, print header)


    */
    TurbulenceStatisticsCcy(std::shared_ptr<Core::FE::Discretization> actdis, bool alefluid,
        std::shared_ptr<Core::LinAlg::Vector<double>> dispnp, Teuchos::ParameterList& params,
        const std::string& statistics_outfilename, const bool withscatra);

    /*!
    \brief Destructor

    */
    virtual ~TurbulenceStatisticsCcy() = default;


    //! @name functions for (spatial) averaging

    /*!
    \brief Compute the in-shell mean values of first and second order
    moments for velocities, pressure (and transported scalar fields).
    */
    void do_time_sample(Core::LinAlg::Vector<double>& velnp,
        std::shared_ptr<Core::LinAlg::Vector<double>> scanp,
        std::shared_ptr<Core::LinAlg::Vector<double>> fullphinp);


    /*!
    \brief Compute in plane means of u,u^2 etc. (nodal quantities)

    The averages here are calculated as the arithmetic mean of
    point values (computed by interpolation)

    The calculated values are added to the pointsum**,pointsumsq** variables
    in the component corresponding to the plane.

    velnp is the solution vector provided by the time integration
    algorithm
    */
    void evaluate_pointwise_mean_values_in_planes();

    //@}

    //! @name Miscellaneous

    /*!
    \brief Compute a time average of the mean values over all steps
    since the last output. Dump the result to file.

    step on input is used to print the timesteps which belong to the
    statistic to the file

    */

    void time_average_means_and_output_of_statistics(int step);

    /*!
    \brief Reset sums and number of samples to 0

    */

    void clear_statistics();

    /*!
    \brief Provide the radius of the homogeneous shell for a
    flow around a rotating circular cylinder

    */
    std::vector<double> return_shell_plane_radius() { return (*nodeshells_); };

    //@}

    // Add results from scalar transport field solver to statistics
    void add_scatra_results(
        std::shared_ptr<Core::FE::Discretization> scatradis, Core::LinAlg::Vector<double>& phinp);

   protected:
    /*!
    \brief sort criterium for double values up to a tolerance of 10-6

    This is used to create sets of doubles (e.g. coordinates)

    */
    class PlaneSortCriterion
    {
     public:
      bool operator()(const double& p1, const double& p2) const { return (p1 < p2 - 1E-6); }

     protected:
     private:
    };

   private:
    //! direction normal to homogenous plane
    int dim_;

    //! number of samples taken
    int numsamp_;

    //! number of records written
    int countrecord_;

    //! The discretisation (required for nodes, dofs etc;)
    std::shared_ptr<Core::FE::Discretization> discret_;

    //! the scatra discretization
    std::shared_ptr<Core::FE::Discretization> scatradis_;

    //! node displacements due to mesh motion
    std::shared_ptr<Core::LinAlg::Vector<double>> dispnp_;

    //! contains plane normal direction etc --- this is the original
    //! fluid dynamic parameterlist
    Teuchos::ParameterList& params_;

    //! name of statistics output file, despite the ending
    const std::string statistics_outfilename_;


    //! parameterlist for the element call when averages of residuals
    //! are calculated --- used for communication between element
    //! and averaging methods
    Teuchos::ParameterList eleparams_;

    //! pointer to mean vel/pres field
    std::shared_ptr<Core::LinAlg::Vector<double>> meanvelnp_;

    //! the dim_-coordinates of the homogeneous planes containing nodes
    std::shared_ptr<std::vector<double>> nodeshells_;

    //! the dim_-coordinates of the homogeneous planes --- including
    // additional sampling planes
    std::shared_ptr<std::vector<double>> shellcoordinates_;

    //!--------------------------------------------------
    //!       the pointwise averaged stuff
    //!--------------------------------------------------
    //

    //! sum over u (over one plane in each component)
    std::shared_ptr<std::vector<double>> pointsumu_;
    //! sum over v (over one plane in each component)
    std::shared_ptr<std::vector<double>> pointsumv_;
    //! sum over w (over one plane in each component)
    std::shared_ptr<std::vector<double>> pointsumw_;
    //! sum over p (over one plane in each component)
    std::shared_ptr<std::vector<double>> pointsump_;

    //! sum over u^2 (over one plane in each component)
    std::shared_ptr<std::vector<double>> pointsumuu_;
    //! sum over v^2 (over one plane in each component)
    std::shared_ptr<std::vector<double>> pointsumvv_;
    //! sum over w^2 (over one plane in each component)
    std::shared_ptr<std::vector<double>> pointsumww_;
    //! sum over p^2 (over one plane in each component)
    std::shared_ptr<std::vector<double>> pointsumpp_;

    //! sum over uv (over one plane in each component)
    std::shared_ptr<std::vector<double>> pointsumuv_;
    //! sum over uw (over one plane in each component)
    std::shared_ptr<std::vector<double>> pointsumuw_;
    //! sum over vw (over one plane in each component)
    std::shared_ptr<std::vector<double>> pointsumvw_;

    // for additional scalar transport

    //! flag for mass transport statistics
    const bool withscatra_;

    //! number of scatra dofs per node
    int numscatradofpernode_;

    //! pointer to mean scalar field
    std::shared_ptr<Core::LinAlg::Vector<double>> meanscanp_;

    //! pointer to mean field of all scatra results
    std::shared_ptr<Core::LinAlg::Vector<double>> meanfullphinp_;

    //! sum over c (over one plane in each component)
    std::shared_ptr<std::vector<double>> pointsumc_;
    //! sum over c^2 (over one plane in each component)
    std::shared_ptr<std::vector<double>> pointsumcc_;

    //! sum over c (over one plane in each component)
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> pointsumphi_;
    //! sum over c^2 (over one plane in each component)
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> pointsumphiphi_;
  };

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
