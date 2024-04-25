/*----------------------------------------------------------------------*/
/*! \file

\brief Solve fluid problems on ALE mesh

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FSI_FLUID_ALE_HPP
#define FOUR_C_FSI_FLUID_ALE_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_fluid_ale.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FSI
{
  /// Fluid on Ale test algorithm
  class FluidAleAlgorithm : public ADAPTER::FluidMovingBoundaryBaseAlgorithm
  {
   public:
    explicit FluidAleAlgorithm(const Epetra_Comm& comm);

    /// time loop
    void Timeloop();

    /// communicator
    const Epetra_Comm& Comm() const { return comm_; }

    /// read restart data
    virtual void ReadRestart(int step);

   protected:
    /// time step size
    double Dt() const { return dt_; }

    /// step number
    int Step() const { return step_; }

    //! @name Time loop building blocks

    /// tests if there are more time steps to do
    bool NotFinished() { return step_ < nstep_ and time_ <= maxtime_; }

    /// start a new time step
    void PrepareTimeStep();

    /// solve ale and fluid fields
    void Solve();

    /// take current results for converged and save for next time step
    void Update();

    /// write output
    void Output();

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
