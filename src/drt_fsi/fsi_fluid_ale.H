#ifdef CCADISCRET

#ifndef FSI_FLUID_ALE_H
#define FSI_FLUID_ALE_H

#include <Epetra_Comm.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "fsi_fluid.H"
#include "fsi_ale.H"
#include "fsi_coupling.H"

namespace FSI
{

  /// Fluid on Ale test algorithm
  class FluidAleAlgorithm : public FluidBaseAlgorithm,
                            public AleBaseAlgorithm
  {
  public:
    explicit FluidAleAlgorithm(Epetra_Comm& comm);

    virtual ~FluidAleAlgorithm();

    /// time loop
    void Timeloop();

    /// communicator
    const Epetra_Comm& Comm() const { return comm_; }

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

    //! @name Transfer helpers

    virtual Teuchos::RCP<Epetra_Vector> AleToFluid(Teuchos::RCP<Epetra_Vector> iv) const;

    virtual Teuchos::RCP<Epetra_Vector> AleToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;

    //@}

    Coupling& FluidAleCoupling()       { return coupfa_; }

    const Coupling& FluidAleCoupling()       const { return coupfa_; }

  private:

    /// coupling of fluid and ale
    Coupling               coupfa_;

    /// comunication (mainly for screen output)
    Epetra_Comm& comm_;

    //! @name Time stepping variables
    int step_;
    int nstep_;
    double time_;
    double maxtime_;
    double dt_;
    //@}
  };
}

#endif
#endif
