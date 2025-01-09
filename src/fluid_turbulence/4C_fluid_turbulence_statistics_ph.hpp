// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_TURBULENCE_STATISTICS_PH_HPP
#define FOUR_C_FLUID_TURBULENCE_STATISTICS_PH_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_MpiComm.h>

#include <memory>

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TurbulenceStatisticsPh
  {
   public:
    /*!
    \brief Standard Constructor (public)

        o Create sets for lines

    o Allocate distributed vector for squares

    */
    TurbulenceStatisticsPh(std::shared_ptr<Core::FE::Discretization> actdis,
        Teuchos::ParameterList& params, const std::string& statistics_outfilename);

    /*!
    \brief Destructor

    */
    virtual ~TurbulenceStatisticsPh() = default;


    //! @name functions for averaging

    /*!
    \brief The values of velocity and its squared values are added to
    global vectors. This method allows to do the time average of the
    nodal values after a certain amount of timesteps.
    */
    void do_time_sample(
        Core::LinAlg::Vector<double>& velnp, Core::LinAlg::Vector<double>& stresses);

    /*!
    \brief Dump the result to file.

    step on input is used to print the timesteps which belong to the
    statistic to the file
    */

    void dump_statistics(int step);


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

    //! number of coordinates in x1- and x2-direction
    int numx1coor_;
    int numx2coor_;

    //! number of locations in x1- and x2-direction for statistical evaluation
    int numx1statlocations_;

    //! The discretisation (required for nodes, dofs etc;)
    std::shared_ptr<Core::FE::Discretization> discret_;

    //! parameter list
    Teuchos::ParameterList& params_;

    //! name of statistics output file, despite the ending
    const std::string statistics_outfilename_;

    //! pointer to vel/pres^2 field (space allocated in constructor)
    std::shared_ptr<Core::LinAlg::Vector<double>> squaredvelnp_;

    //! toggle vectors: sums are computed by scalarproducts
    std::shared_ptr<Core::LinAlg::Vector<double>> toggleu_;
    std::shared_ptr<Core::LinAlg::Vector<double>> togglev_;
    std::shared_ptr<Core::LinAlg::Vector<double>> togglew_;
    std::shared_ptr<Core::LinAlg::Vector<double>> togglep_;

    //! available x1- and x2-coordinates
    std::shared_ptr<std::vector<double>> x1coordinates_;
    std::shared_ptr<std::vector<double>> x2coordinates_;

    //! coordinates of locations in x1- and x2-direction for statistical evaluation
    Core::LinAlg::SerialDenseMatrix x1statlocations_;

    //! matrix for r-coordinates (columns are evaluation planes
    Core::LinAlg::SerialDenseMatrix x2statlocations_;

    //! set coordinates of locations in x1-direction for statistical evaluation
    std::shared_ptr<std::vector<double>> x1setstatlocations_;

    //! coordinates in x1-direction for sampling velocity gradient at the middle bottom

    //! matrices containing values
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> x1sumu_;
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> x1sump_;
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> x1sumf_;

    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> x2sumu_;
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> x2sumv_;
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> x2sumw_;
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> x2sump_;


    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> x2sumsqu_;
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> x2sumsqv_;
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> x2sumsqw_;
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> x2sumsqp_;

    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> x2sumuv_;
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> x2sumuw_;
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> x2sumvw_;
  };

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
