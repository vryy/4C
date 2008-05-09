/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_moving_boundary.H

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_GENERAL_FLUID_H
#define ADAPTER_GENERAL_FLUID_H

#include <Epetra_Comm.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "adapter_fluid.H"
#include "adapter_ale.H"

namespace ADAPTER
{

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
    virtual ~FluidMovingBoundary() {}

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
    //virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> vel) = 0;

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

    virtual void RobinNonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,
				     Teuchos::RCP<Epetra_Vector> ivel,
				     Teuchos::RCP<Epetra_Vector> iforce)=0;

    /// linear fluid solve with just a interface load
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                                        double dt) = 0;

    /// get the linear solver object used for this field
    //virtual Teuchos::RCP<LINALG::Solver> LinearSolver() = 0;

    //@}

    //! @name Extract interface forces

    /// After the fluid solve we need the forces at the FSI interface.
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces() = 0;

    /// After the fluid solve we need the forces at the FSI interface.
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForcesRobin() = 0;

    //@}

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceFluidVelocity() = 0;

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


  /// base of all algorithms that use a fluid on a variable domain
  class FluidMovingBoundaryBaseAlgorithm
  {
  public:

    /// constructor
    explicit FluidMovingBoundaryBaseAlgorithm(const Teuchos::ParameterList& prbdyn,
                                       std::string condname);

    /// virtual destructor to support polymorph destruction
    virtual ~FluidMovingBoundaryBaseAlgorithm() {}

    /// fluid field solver
    FluidMovingBoundary& MBFluidField() { return *fluid_; }
    const FluidMovingBoundary& MBFluidField() const { return *fluid_; }

  private:

    /// fluid field solver
    Teuchos::RCP<FluidMovingBoundary> fluid_;
  };
}

#endif
#endif

