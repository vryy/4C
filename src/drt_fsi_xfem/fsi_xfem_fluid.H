/*----------------------------------------------------------------------*/
/*!
\file fsi_xfem_fluid.H

\brief ...

\level 2

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
*/

#ifndef FSI_FLUID_XFEM_H
#define FSI_FLUID_XFEM_H

#include "../drt_adapter/ad_fld_fluid_xfem.H"

namespace FSI
{
  /// Fluid on XFEM test algorithm
  class FluidXFEMAlgorithm : public ::ADAPTER::FluidMovingBoundaryBaseAlgorithm
  {
   public:
    explicit FluidXFEMAlgorithm(const Epetra_Comm& comm);

    virtual ~FluidXFEMAlgorithm();

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

#endif
