// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_ALGORITHMBASE_HPP
#define FOUR_C_ADAPTER_ALGORITHMBASE_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

#include <mpi.h>

#include <string>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  /// base class for algorithms
  /*!

    Here we just keep the time variables as it is always the same.

    \author u.kue
    \date 07/09
   */
  class AlgorithmBase
  {
   public:
    explicit AlgorithmBase(MPI_Comm comm, const Teuchos::ParameterList& timeparams);

    /// virtual destruction
    virtual ~AlgorithmBase() = default;
    /// read restart data
    virtual void read_restart(int step) = 0;

    /// communicator
    MPI_Comm get_comm() const { return comm_; }

    /// tests if there are more time steps to do
    bool not_finished() const { return step_ < nstep_ and time_ + 1e-8 * dt_ < maxtime_; }

    /// time step size
    double dt() const { return dt_; }

    /// set time step size
    void set_dt(const double stepsize) { dt_ = stepsize; }

    /// current time
    double time() const { return time_; }

    /// current time step number
    int step() const { return step_; }

    /// set new time step in read_restart() or in a potential outer control
    void set_time_step(const double time,  ///< physical time to set
        const int step                     ///< time step number to set
    );

   protected:
    /// total number of time steps
    int n_step() const { return nstep_; }

    /// maximum simulation time
    double max_time() const { return maxtime_; }

    //! @name Time loop building blocks

    /// increment time and step value
    void increment_time_and_step()
    {
      step_ += 1;
      time_ += dt_;
    }

    /// start a new time step
    virtual void prepare_time_step() {}

    /// take current results for converged and save for next time step
    virtual void update() {}

    /// write output
    virtual void output() {}

    /// set method name for screen output
    void set_method(std::string method) { method_ = method; }

    //@}

    /// print time step header
    void print_header();

    /// return printscreen_
    int print_screen_every() { return printscreen_; }

   private:
    /// communication (mainly for screen output)
    MPI_Comm comm_;

    /// method name prepared for output
    std::string method_;

    /// print infos to standard out every printscreen_ steps
    int printscreen_;

    //! @name Time stepping variables

    /// current time step number
    int step_;

    /// total number of time steps
    int nstep_;

    /// current physical time at end of time step, i.e. \f$t_{n+1}\f$
    double time_;

    /// maximum simulation time
    double maxtime_;

    /// current time step size
    double dt_;

    //@}
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
