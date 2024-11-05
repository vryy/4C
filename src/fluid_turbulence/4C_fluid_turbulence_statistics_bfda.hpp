// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_TURBULENCE_STATISTICS_BFDA_HPP
#define FOUR_C_FLUID_TURBULENCE_STATISTICS_BFDA_HPP

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
  class TurbulenceStatisticsBfda
  {
   public:
    /*!
    \brief Standard Constructor (public)

        o Create vector "zcoordinates" with coordinates along z-axis
    o Create matrix "rcoordinates_" with coordinates of radial coordinates of evaluation planes
      columns are evaluation planes corresponding to the positions in "posEvaluation_"

    */
    TurbulenceStatisticsBfda(std::shared_ptr<Core::FE::Discretization> actdis,
        Teuchos::ParameterList& params, const std::string& statistics_outfilename);

    /*!
    \brief Destructor

    */
    virtual ~TurbulenceStatisticsBfda() = default;


    //! @name functions for averaging


    void do_time_sample(Core::LinAlg::Vector<double>& velnp);

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


    //! number of coordinates in z-direction
    int numzcoor_;

    //! number in radial direction for statistical evaluation
    int numrstatlocations_;

    //! The discretisation (required for nodes, dofs etc;)
    std::shared_ptr<Core::FE::Discretization> discret_;

    //! parameter list
    Teuchos::ParameterList& params_;

    //! name of statistics output file, despite the ending
    const std::string statistics_outfilename_;

    //! positions for evaluation
    std::vector<double> pos_evaluation_;

    //! positions for evaluation of actual mesh
    std::vector<double> act_pos_evaluation_;

    //! Toggle vectors
    std::shared_ptr<Core::LinAlg::Vector<double>> togglew_;
    std::shared_ptr<Core::LinAlg::Vector<double>> togglep_;

    //! vector for z-coordinates
    std::shared_ptr<std::vector<double>> zcoordinates_;
    //! matrix for r-coordinates (columns are evaluation planes
    Core::LinAlg::SerialDenseMatrix rcoordinates_;

    //! matrices for mean values along z-axis
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> zsumw_;
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> zsump_;

    //! matrices for mean values of evaluation planes
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> rsumw_;
    std::shared_ptr<Core::LinAlg::SerialDenseMatrix> rsump_;
  };

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
