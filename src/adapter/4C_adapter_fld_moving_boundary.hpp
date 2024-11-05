// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_FLD_MOVING_BOUNDARY_HPP
#define FOUR_C_ADAPTER_FLD_MOVING_BOUNDARY_HPP


#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_utils_result_test.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Discret
{
  class ResultTest;
}

namespace FLD
{
  namespace Utils
  {
    class MapExtractor;
  }
}  // namespace FLD

namespace Adapter
{
  class Fluid;

  /// generalized fluid base: fluid with moving boundaries
  /*!

    The generalized fluid is a fluid on a variable domain. The domain could be
    deforming (thus the fluid is solved on an ale mesh) or could be cut using
    xfem. Nevermind, the outside world sees a general fluid interface.

    \author u.kue
    \date 03/08
   */
  class FluidMovingBoundary
  {
   public:
    /// virtual destructor to get polymorph destruction
    virtual ~FluidMovingBoundary() = default;
    //! @name Misc

    /// direct access to discretization
    virtual std::shared_ptr<Core::FE::Discretization> discretization() = 0;

    virtual const std::shared_ptr<Adapter::Fluid>& fluid_field() = 0;

    /// communication object at the interface
    virtual std::shared_ptr<FLD::Utils::MapExtractor> const& interface() const = 0;

    //@}

    //! @name Time step helpers

    /// start new time step
    virtual void prepare_time_step() = 0;

    /// evaluate elements with given displacement
    // virtual void evaluate(std::shared_ptr<const Core::LinAlg::Vector<double>> vel) = 0;

    /// update at time step end
    virtual void update() = 0;

    /// output results
    virtual void output() = 0;

    /// read restart information for given time step
    virtual double read_restart(int step) = 0;

    //@}

    //! @name Solver calls

    /// nonlinear solve
    virtual void nonlinear_solve(std::shared_ptr<Core::LinAlg::Vector<double>> idisp = nullptr,
        std::shared_ptr<Core::LinAlg::Vector<double>> ivel = nullptr) = 0;

    /// nonlinear solve
    virtual void apply_interface_values(
        std::shared_ptr<Core::LinAlg::Vector<double>> idisp = nullptr,
        std::shared_ptr<Core::LinAlg::Vector<double>> ivel = nullptr)
    {
      FOUR_C_THROW("Not implemented in base class");
    }

    /// linear fluid solve with just a interface load
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> relaxation_solve(
        std::shared_ptr<Core::LinAlg::Vector<double>> idisp, double dt) = 0;

    //@}

    //! @name Extract interface forces

    /// After the fluid solve we need the forces at the FSI interface.
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> extract_interface_forces() = 0;

    //@}

    virtual std::shared_ptr<Core::LinAlg::Vector<double>> extract_interface_velnp() = 0;

    /// extract old velocities
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> extract_interface_veln() = 0;

    //! @name Number of Newton iterations
    //! For simplified FD MFNK solve we want to temporally limit the
    /// number of Newton steps inside the fluid solver

    virtual int itemax() const = 0;
    virtual void set_itemax(int itemax) = 0;

    //@}

    /// integrate FSI interface shape functions
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> integrate_interface_shape() = 0;

    /// create result test for encapulated fluid algorithm
    virtual std::shared_ptr<Core::Utils::ResultTest> create_field_test() = 0;
  };


  /// base of all algorithms that use a fluid on a variable domain
  class FluidMovingBoundaryBaseAlgorithm
  {
   public:
    /// constructor
    explicit FluidMovingBoundaryBaseAlgorithm(
        const Teuchos::ParameterList& prbdyn, std::string condname);

    /// virtual destructor to support polymorph destruction
    virtual ~FluidMovingBoundaryBaseAlgorithm() = default;
    /// fluid field solver
    const std::shared_ptr<FluidMovingBoundary>& mb_fluid_field() { return fluid_; }

   private:
    /// fluid field solver
    std::shared_ptr<FluidMovingBoundary> fluid_;
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
