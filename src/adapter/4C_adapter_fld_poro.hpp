// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_FLD_PORO_HPP
#define FOUR_C_ADAPTER_FLD_PORO_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_fluid_fpsi.hpp"
#include "4C_fem_condition.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_poroelast_utils.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class FluidPoro : public FluidFPSI
  {
   public:
    //! Constructor
    FluidPoro(std::shared_ptr<Fluid> fluid, std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Core::IO::DiscretizationWriter> output, bool isale, bool dirichletcond);

    //! Evaluate no penetration constraint
    /*!
     \param Cond_RHS                  (o) condition part of rhs
     \param ConstraintMatrix          (o) static part of Fluid matrix associated with constraints
     \param struct_vel_constraint_matrix (o) transient part of Fluid matrix associated with
     constraints \param condIDs                   (o) vector containing constraint dofs \param
     coupltype                 (i) coupling type, determines which matrix is to be evaluated (0==
     fluid-fluid, 1== fluid -structure)
     */
    void evaluate_no_penetration_cond(std::shared_ptr<Core::LinAlg::Vector<double>> Cond_RHS,
        std::shared_ptr<Core::LinAlg::SparseMatrix> ConstraintMatrix,
        std::shared_ptr<Core::LinAlg::SparseMatrix> struct_vel_constraint_matrix,
        std::shared_ptr<Core::LinAlg::Vector<double>> condVector, std::set<int>& condIDs,
        PoroElast::Coupltype coupltype = PoroElast::fluidfluid);

    //! calls the VelPresSplitter on the time integrator
    virtual std::shared_ptr<Core::LinAlg::MapExtractor> vel_pres_splitter();

    /*!
      \brief Write extra output for specified step and time.
             Useful if you want to write output every iteration in partitioned schemes.
             If no step and time is provided, standard Output of fluid field is invoked.

      \param step (in) : Pseudo-step for which extra output is written
      \param time (in) : Pseudo-time for which extra output is written

      \note This is a pure DEBUG functionality. Originally used in immersed method development.

      \warning This method partly re-implements redundantly few lines of the common fluid output()
      routine. \return void
    */
    virtual void output(const int step = -1, const double time = -1);

   private:
    /// fluid field
    const std::shared_ptr<Adapter::Fluid>& fluid_field() { return fluid_; }

    std::vector<Core::Conditions::Condition*>
        nopencond_;  ///< vector containing no penetration conditions
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
