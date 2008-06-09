/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid.H

\brief Fluid field adapter

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_FLUID_H
#define ADAPTER_FLUID_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_lib/linalg_systemmatrix.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/drt_resulttest.H"

#include "adapter_utils.H"


/// Algorithm adapter classes for use in multi-field problems
namespace ADAPTER {

  /// general fluid field interface for multiphysics (FSI, ELCH, ...)
  /*!

  This is the FSI algorithm's view on a fluid algorithm. This pure virtual
  interface contains all the methods any FSI algorithm might want to
  call. The idea is to implement this interface with a concrete adapter class
  for each fluid algorithm we want to use for FSI.

  FSI is quite demanding when it comes to knowledge about the internal details
  of a fluid algorithm. Furthermore there are different coupling versions, all
  of them requiring a slightly different view. Yet these views have a lot in
  common, so a common adapter clas seems appropriate.

  We use this adapter interface instead of deriving from the fluid algorithm
  class. The good thing is that this way we keep control over the exported
  entities. The down side is that the fluid algorithm still has to grant
  access to a lot of internal state variables.

  Currently supported FSI couplings:

  - Dirichlet-Neumann coupling (fixed-point, Newton-Krylov, vector
    extrapolation)

  - Monolithic (Newton with overlapping blocks and block preconditioning)

  \warning Further cleanup is still needed.

  \sa Structure, Ale
  \author u.kue
  \date 11/07
  */
  class Fluid
  {
  public:

    /// virtual destructor to get polymorph destruction
    virtual ~Fluid();

    //! @name Vector access

    /// initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() = 0;

    /// rhs of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> RHS() = 0;

    /// unknown velocities (and pressures) at t(n+1)
    virtual Teuchos::RCP<const Epetra_Vector> Velnp() = 0;

    /// unknown velocities (and pressures) at t(n)
    virtual Teuchos::RCP<const Epetra_Vector> Veln() = 0;

    /// mesh displacements at t(n+1)
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp() = 0;

    //@}

    //! @name Misc

    /// dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() = 0;

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() = 0;

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() = 0;

    /// linearization of Navier-Stokes with respect to mesh movement
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> MeshMoveMatrix() = 0;

    /// direct access to discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;

    //virtual Teuchos::RCP<Epetra_Vector> StructCondRHS() const = 0;

    //@}

    //! @name Time step helpers

    /// start new time step
    virtual void PrepareTimeStep() = 0;

    /// evaluate elements with given displacement
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> vel) = 0;

    /// update at time step end
    virtual void Update() = 0;

    /// output results
    virtual void Output() = 0;

    /// read restart information for given time step
    virtual void ReadRestart(int step) = 0;

    /// current time value
    virtual double Time() = 0;

    /// current time step
    virtual int Step() = 0;

    /// calculation (and output) of lift&drag
    virtual void LiftDrag() = 0;

    //@}

    //! @name Solver calls

    /// nonlinear solve
    /*!
      Do the nonlinear solve for the time step. All boundary conditions have
      been set.
     */
    virtual void NonlinearSolve() = 0;

    /// linear fluid solve with just a interface load
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel) = 0;

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() = 0;

    //@}

    /// Map of all velocity dofs that are not Dirichlet-constrained
    virtual Teuchos::RCP<const Epetra_Map> InnerVelocityRowMap() = 0;

    /// Map of all velocity dofs
    virtual Teuchos::RCP<const Epetra_Map>    VelocityRowMap() = 0;

    /// Map of all pressure dofs
    virtual Teuchos::RCP<const Epetra_Map> PressureRowMap() = 0;

    /// the mesh map contains all velocity dofs that are covered by an ALE node
    virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm) = 0;

    /// scaling factor between velocity and displacement
    virtual double ResidualScaling() const = 0;

    /// time scaling factor to get back to real forces
    virtual double TimeScaling() const = 0;

    /// communication object at the interface
    virtual const LINALG::MapExtractor& Interface() const = 0;

    /// communication object at the free surface
    virtual const LINALG::MapExtractor& FreeSurface() const = 0;

    //! @name Extract the velocity-related part of a fluid vector (e.g. velnp, veln, residual)
    /// The idea is to have one function that does the extraction and call it
    /// with different vectors.

    /// Some applications need only access to velocity-related values of an fluid result vector.
    virtual Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres) = 0;

    //@}

    //! @name Apply interface values

    /// at the interface the velocity is prescribed as a Dirichlet condition
    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel) = 0;

    /// apply a generalized robin force to the interface
    virtual void ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> ivel,
                                          Teuchos::RCP<Epetra_Vector> iforce) = 0;

    //@}

    //! @name Extract interface values
    /// Maybe we do not need all of them?

    /// extract fluid velocity at the interface
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceFluidVelocity() = 0;

    /// extract fluid forces at the interface
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces() = 0;

    /// After the fluid solve we need the forces at the FSI interface.
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForcesRobin() = 0;

    //@}


    //! @name Extract mesh values

    /// tell the mesh displacement to the fluid solver
    virtual void ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp) = 0;

    /// tell the mesh velocity to the fluid solver
    virtual void ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel) = 0;

    //@}

    /// convert Delta d(n+1,i+1) to the fluid unknown at the interface
    virtual void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx) = 0;

    /// convert the fluid unknown to Delta d(n+1,i+1) at the interface
    virtual void VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx) = 0;

    //! @name Number of Newton iterations
    //! For simplified FD MFNK solve we want to temporally limit the
    /// number of Newton steps inside the fluid solver

    virtual int Itemax() const = 0;
    virtual void SetItemax(int itemax) = 0;

    //@}

    /// integrate FSI interface shape functions
    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape() = 0;

    /// switch fluid field to block matrix
    virtual void UseBlockMatrix(const LINALG::MultiMapExtractor& domainmaps,
                                const LINALG::MultiMapExtractor& rangemaps,
                                bool splitmatrix) = 0;

    /// create result test for encapulated fluid algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest() = 0;
  };
}

#endif
#endif
