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

  /// generalized fluid adapter base
  /*!

    The generalized fluid is a fluid on a variable domain. The domain could be
    deforming (thus the fluid is solved on an ale mesh) or could be cut using
    xfem. Nevermind, the outside world sees a general fluid interface.

    \author u.kue
    \date 03/08
   */
  class GeneralFluid
  {
  public:

    /// virtual destructor to get polymorph destruction
    virtual ~GeneralFluid() {}

    //! @name Misc

    /// direct access to discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;

    /// communication object at the interface
    virtual const LINALG::MapExtractor& Interface() const = 0;

    //@}

    //! @name Time step helpers

    /// start new time step
    virtual void PrepareTimeStep() = 0;

    /// evaluate elements with given displacement
    //virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> vel) const = 0;

    /// update at time step end
    virtual void Update() = 0;

    /// output results
    virtual void Output() = 0;

    /// read restart information for given time step
    virtual double ReadRestart(int step) = 0;

    /// calculation (and output) of lift&drag
    //virtual void LiftDrag() = 0;

    //@}

    //! @name Solver calls

    /// nonlinear solve
    virtual void NonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp=Teuchos::null,
                                Teuchos::RCP<Epetra_Vector> ivel=Teuchos::null) = 0;

    /// linear fluid solve with just a interface load
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                                        double dt) = 0;

    /// get the linear solver object used for this field
    //virtual Teuchos::RCP<LINALG::Solver> LinearSolver() = 0;

    //@}

    //! @name Extract interface forces

    /// After the fluid solve we need the forces at the FSI interface.
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces() = 0;

    //@}

    //! @name Number of Newton iterations
    //! For simplified FD MFNK solve we want to temporally limit the
    /// number of Newton steps inside the fluid solver

    virtual int Itemax() const = 0;
    virtual void SetItemax(int itemax) = 0;

    //@}

    /// integrate FSI interface shape functions
    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape() = 0;

    /// create result test for encapulated fluid algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest() = 0;
  };


  /// fluid on ale
  class FluidAleAdapter : public GeneralFluid
  {
  public:

    FluidAleAdapter(const Teuchos::ParameterList& prbdyn);

    /// fluid field
    Fluid& FluidField() { return fluid_.FluidField(); }
    const Fluid& FluidField() const { return fluid_.FluidField(); }

    /// ale field
    Ale& AleField() { return ale_.AleField(); }
    const Ale& AleField() const { return ale_.AleField(); }

    virtual Teuchos::RCP<DRT::Discretization> Discretization();
    virtual const LINALG::MapExtractor& Interface() const;

    virtual void PrepareTimeStep();
    virtual void Update();
    virtual void Output();
    virtual double ReadRestart(int step);

    virtual void NonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,Teuchos::RCP<Epetra_Vector> ivel);
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                                        double dt);

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();

    virtual int Itemax() const { return FluidField().Itemax(); }
    virtual void SetItemax(int itemax) { FluidField().SetItemax(itemax); }

    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape();

    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  protected:

    //! @name Transfer helpers

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const;

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> AleToFluidField(Teuchos::RCP<const Epetra_Vector> iv) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<const Epetra_Vector> iv) const;

    //@}

    //Coupling& FluidAleFieldCoupling()             { return coupfa_; }
    //const Coupling& FluidAleFieldCoupling() const { return coupfa_; }

  private:

    FluidBaseAlgorithm fluid_;
    AleBaseAlgorithm ale_;

    /// coupling of fluid and ale (whole field)
    Coupling               coupfa_;

    /// coupling of fluid and ale (interface only)
    Coupling               icoupfa_;

  };


  /// base of all algorithms that use a fluid on a variable domain
  class GeneralFluidBaseAlgorithm
  {
  public:

    /// constructor
    explicit GeneralFluidBaseAlgorithm(const Teuchos::ParameterList& prbdyn);

    /// virtual destructor to support polymorph destruction
    virtual ~GeneralFluidBaseAlgorithm() {}

    /// fluid field solver
    GeneralFluid& FluidField() { return *fluid_; }
    const GeneralFluid& FluidField() const { return *fluid_; }

  private:

    /// fluid field solver
    Teuchos::RCP<GeneralFluid> fluid_;
  };


  /// Fluid on Ale test algorithm
  class FluidAleAlgorithm : public GeneralFluidBaseAlgorithm
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

  private:

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
