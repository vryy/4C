// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_TURBULENCE_STATISTICS_SQC_HPP
#define FOUR_C_FLUID_TURBULENCE_STATISTICS_SQC_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_MpiComm.h>

#include <memory>

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TurbulenceStatisticsSqc
  {
   public:
    /*!
    \brief Standard Constructor (public)

        o Create sets for lines in x1- and x2-direction

    o Allocate distributed vector for squares

    */
    TurbulenceStatisticsSqc(std::shared_ptr<Core::FE::Discretization> actdis,
        Teuchos::ParameterList& params, const std::string& statistics_outfilename);

    /*!
    \brief Destructor

    */
    virtual ~TurbulenceStatisticsSqc() = default;


    //! @name functions for averaging

    /*!
    \brief The values of lift and drag and its squared values are added.
     This method allows to do the time average after a certain amount of
     timesteps.
    */
    void do_lift_drag_time_sample(double dragforce, double liftforce);

    /*!
    \brief The values of velocity and its squared values are added to
    global vectors. This method allows to do the time average of the
    nodal values after a certain amount of timesteps.
    */
    void do_time_sample(Core::LinAlg::Vector<double>& velnp);

    /*!
    \brief Dump the result to file.

    step on input is used to print the timesteps which belong to the
    statistic to the file
    */

    void dump_statistics(int step);

    /*!
    \brief Reset sums and number of samples to 0
    */

    void clear_statistics();


   protected:
    /*!
    \brief sort criterium for double values up to a tolerance of 10-9

    This is used to create sets of doubles (e.g. coordinates)

    */
    class LineSortCriterion
    {
     public:
      bool operator()(const double& p1, const double& p2) const { return (p1 < p2 - 1E-9); }

     protected:
     private:
    };

   private:
    //! number of samples taken
    int numsamp_;
    //! homogeneous direction for sampling
    std::string homdir_;

    //! bounds for extension of cavity in x3-direction
    double x3min_;
    double x3max_;

    //! sums over lift and drag values
    double lift_;
    double drag_;
    double liftsq_;
    double dragsq_;

    //! The discretisation (required for nodes, dofs etc;)
    std::shared_ptr<Core::FE::Discretization> discret_;

    //! parameter list
    Teuchos::ParameterList& params_;

    //! name of statistics output file, despite the ending
    const std::string statistics_outfilename_;


    //! pointer to vel/pres^2 field (space allocated in constructor)
    std::shared_ptr<Core::LinAlg::Vector<double>> squaredvelnp_;

    //! toogle vectors: sums are computed by scalarproducts
    std::shared_ptr<Core::LinAlg::Vector<double>> toggleu_;
    std::shared_ptr<Core::LinAlg::Vector<double>> togglev_;
    std::shared_ptr<Core::LinAlg::Vector<double>> togglew_;
    std::shared_ptr<Core::LinAlg::Vector<double>> togglep_;

    //! the coordinates of various lines
    std::shared_ptr<std::vector<double>> x1ccoordinates_;
    std::shared_ptr<std::vector<double>> x2ccoordinates_;
    std::shared_ptr<std::vector<double>> x2wcoordinates_;
    std::shared_ptr<std::vector<double>> clrcoordinates_;
    std::shared_ptr<std::vector<double>> ctbcoordinates_;
    //! all coordinates in x1- and x2-direction (required for averaging of Smagorinsky constant)
    std::shared_ptr<std::vector<double>> x1coordinates_;
    std::shared_ptr<std::vector<double>> x2coordinates_;

    //! sum over u
    std::shared_ptr<std::vector<double>> x1csumu_;
    std::shared_ptr<std::vector<double>> x2csumu_;
    std::shared_ptr<std::vector<double>> x2w1sumu_;
    std::shared_ptr<std::vector<double>> x2w2sumu_;
    std::shared_ptr<std::vector<double>> cyllsumu_;
    std::shared_ptr<std::vector<double>> cyltsumu_;
    std::shared_ptr<std::vector<double>> cylrsumu_;
    std::shared_ptr<std::vector<double>> cylbsumu_;
    //! sum over v
    std::shared_ptr<std::vector<double>> x1csumv_;
    std::shared_ptr<std::vector<double>> x2csumv_;
    std::shared_ptr<std::vector<double>> x2w1sumv_;
    std::shared_ptr<std::vector<double>> x2w2sumv_;
    std::shared_ptr<std::vector<double>> cyllsumv_;
    std::shared_ptr<std::vector<double>> cyltsumv_;
    std::shared_ptr<std::vector<double>> cylrsumv_;
    std::shared_ptr<std::vector<double>> cylbsumv_;
    //! sum over w
    std::shared_ptr<std::vector<double>> x1csumw_;
    std::shared_ptr<std::vector<double>> x2csumw_;
    std::shared_ptr<std::vector<double>> x2w1sumw_;
    std::shared_ptr<std::vector<double>> x2w2sumw_;
    std::shared_ptr<std::vector<double>> cyllsumw_;
    std::shared_ptr<std::vector<double>> cyltsumw_;
    std::shared_ptr<std::vector<double>> cylrsumw_;
    std::shared_ptr<std::vector<double>> cylbsumw_;
    //! sum over p
    std::shared_ptr<std::vector<double>> x1csump_;
    std::shared_ptr<std::vector<double>> x2csump_;
    std::shared_ptr<std::vector<double>> x2w1sump_;
    std::shared_ptr<std::vector<double>> x2w2sump_;
    std::shared_ptr<std::vector<double>> cyllsump_;
    std::shared_ptr<std::vector<double>> cyltsump_;
    std::shared_ptr<std::vector<double>> cylrsump_;
    std::shared_ptr<std::vector<double>> cylbsump_;

    //! sum over u^2
    std::shared_ptr<std::vector<double>> x1csumsqu_;
    std::shared_ptr<std::vector<double>> x2csumsqu_;
    std::shared_ptr<std::vector<double>> x2w1sumsqu_;
    std::shared_ptr<std::vector<double>> x2w2sumsqu_;
    std::shared_ptr<std::vector<double>> cyllsumsqu_;
    std::shared_ptr<std::vector<double>> cyltsumsqu_;
    std::shared_ptr<std::vector<double>> cylrsumsqu_;
    std::shared_ptr<std::vector<double>> cylbsumsqu_;
    //! sum over v^2
    std::shared_ptr<std::vector<double>> x1csumsqv_;
    std::shared_ptr<std::vector<double>> x2csumsqv_;
    std::shared_ptr<std::vector<double>> x2w1sumsqv_;
    std::shared_ptr<std::vector<double>> x2w2sumsqv_;
    std::shared_ptr<std::vector<double>> cyllsumsqv_;
    std::shared_ptr<std::vector<double>> cyltsumsqv_;
    std::shared_ptr<std::vector<double>> cylrsumsqv_;
    std::shared_ptr<std::vector<double>> cylbsumsqv_;
    //! sum over w^2
    std::shared_ptr<std::vector<double>> x1csumsqw_;
    std::shared_ptr<std::vector<double>> x2csumsqw_;
    std::shared_ptr<std::vector<double>> x2w1sumsqw_;
    std::shared_ptr<std::vector<double>> x2w2sumsqw_;
    std::shared_ptr<std::vector<double>> cyllsumsqw_;
    std::shared_ptr<std::vector<double>> cyltsumsqw_;
    std::shared_ptr<std::vector<double>> cylrsumsqw_;
    std::shared_ptr<std::vector<double>> cylbsumsqw_;
    //! sum over uv
    std::shared_ptr<std::vector<double>> x1csumuv_;
    std::shared_ptr<std::vector<double>> x2csumuv_;
    std::shared_ptr<std::vector<double>> x2w1sumuv_;
    std::shared_ptr<std::vector<double>> x2w2sumuv_;
    std::shared_ptr<std::vector<double>> cyllsumuv_;
    std::shared_ptr<std::vector<double>> cyltsumuv_;
    std::shared_ptr<std::vector<double>> cylrsumuv_;
    std::shared_ptr<std::vector<double>> cylbsumuv_;
    //! sum over uw
    std::shared_ptr<std::vector<double>> x1csumuw_;
    std::shared_ptr<std::vector<double>> x2csumuw_;
    std::shared_ptr<std::vector<double>> x2w1sumuw_;
    std::shared_ptr<std::vector<double>> x2w2sumuw_;
    std::shared_ptr<std::vector<double>> cyllsumuw_;
    std::shared_ptr<std::vector<double>> cyltsumuw_;
    std::shared_ptr<std::vector<double>> cylrsumuw_;
    std::shared_ptr<std::vector<double>> cylbsumuw_;
    //! sum over vw
    std::shared_ptr<std::vector<double>> x1csumvw_;
    std::shared_ptr<std::vector<double>> x2csumvw_;
    std::shared_ptr<std::vector<double>> x2w1sumvw_;
    std::shared_ptr<std::vector<double>> x2w2sumvw_;
    std::shared_ptr<std::vector<double>> cyllsumvw_;
    std::shared_ptr<std::vector<double>> cyltsumvw_;
    std::shared_ptr<std::vector<double>> cylrsumvw_;
    std::shared_ptr<std::vector<double>> cylbsumvw_;
    //! sum over p^2
    std::shared_ptr<std::vector<double>> x1csumsqp_;
    std::shared_ptr<std::vector<double>> x2csumsqp_;
    std::shared_ptr<std::vector<double>> x2w1sumsqp_;
    std::shared_ptr<std::vector<double>> x2w2sumsqp_;
    std::shared_ptr<std::vector<double>> cyllsumsqp_;
    std::shared_ptr<std::vector<double>> cyltsumsqp_;
    std::shared_ptr<std::vector<double>> cylrsumsqp_;
    std::shared_ptr<std::vector<double>> cylbsumsqp_;
  };

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
