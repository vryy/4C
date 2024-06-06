/*----------------------------------------------------------------------*/
/*! \file

\brief ...

\level 2

*/

#ifndef FOUR_C_FSI_XFEM_FLUID_HPP
#define FOUR_C_FSI_XFEM_FLUID_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_fluid_xfem.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FSI
{
  /// Fluid on XFEM test algorithm
  class FluidXFEMAlgorithm : public Adapter::FluidMovingBoundaryBaseAlgorithm
  {
   public:
    explicit FluidXFEMAlgorithm(const Epetra_Comm& comm);

    /// time loop
    void Timeloop();

    /// communicator
    const Epetra_Comm& Comm() const { return comm_; }

    /// read restart data
    virtual void read_restart(int step);

   protected:
    /// time step size
    double dt() const { return dt_; }

    /// step number
    int step() const { return step_; }

    //! @name Time loop building blocks

    /// tests if there are more time steps to do
    bool not_finished() { return step_ < nstep_ and time_ <= maxtime_; }

    /// start a new time step
    void prepare_time_step();

    /// solve ale and fluid fields
    void solve();

    /// take current results for converged and save for next time step
    void update();

    /// write output
    void output();

    //@}

   private:
    /// comunication (mainly for screen output)
    const Epetra_Comm& comm_;

    //! @name Time stepping variables
    int step_;
    int nstep_;
    double time_;
    double maxtime_;
    double dt_;
    //@}
  };

}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
