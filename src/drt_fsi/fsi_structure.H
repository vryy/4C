/*----------------------------------------------------------------------*/
/*!
\file fsi_structure.H

\brief Structure field adapter

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef MFSI_STRUCTURE_H
#define MFSI_STRUCTURE_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_structure/strugenalpha.H"

#include "fsi_utils.H"


namespace FSI {

  /// general structural field interface
  /*!

  The point is to keep FSI as far apart from our field solvers as
  possible. Each structure field solver we want to use should get its own
  subclass of this. The FSI algorithm should be able to extract all the
  information from the structure field it needs using this interface.

  All FSI algorithms use this adapter to communicate with the structural
  field. There are different ways to use this adapter.

  In all cases you need to tell the structural algorithm about your time
  step. Therefore PrepareTimeStep(), Update() and Output() must be called at
  the appropriate position in the FSI algorithm.

  <h3>Dirichlet-Neumann coupled FSI</h3>

  A good starting displacement can be guessed with PredictInterfaceDispnp().

  Dirichlet-Neumann coupled FSI will need to Solve() the nonlinear
  structural problem for each time step after the fluid forces have been
  applied (ApplyInterfaceForces()). Solve() will be called many times for each
  time step until the interface equilibrium is reached. The structural
  algorithm has to preserve its state until Update() is called.

  After each Solve() you get the interface forces by ExtractInterfaceDispnp().

  A Dirichlet-Neumann FSI with steepest descent relaxation or matrix free
  Newton Krylov will want to solve the structural problem linearly without
  history and prescribed interface forces: RelaxationSolve().

  <h3>Monolithic FSI</h3>

  Monolithic FSI is based on Evaluate() of elements. This results in a new
  RHS() and a new SysMat(). Together with the InitialGuess() these form the
  building blocks for a block based Newton's method.

  \warning Further cleanup is still needed.

  \sa Fluid, Ale
  \author u.kue
  \date 11/07
  */
  class Structure
  {
  public:

    /// virtual to get polymorph destruction
    virtual ~Structure();

    //! @name Vector access

    /// initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() const = 0;

    /// rhs of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> RHS() const = 0;

    /// unknown displacements at t(n+1)
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp() const = 0;

    /// unknown displacements at t(n)
    virtual Teuchos::RCP<const Epetra_Vector> Dispn() const = 0;

    /// unknown displacements at t(n-1)
    virtual Teuchos::RCP<const Epetra_Vector> Dispnm() const = 0;

    //@}

    //! @name Misc

    /// dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() const = 0;

    /// domain map of system matrix (do we really need this?)
    virtual const Epetra_Map& DomainMap() = 0;

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() const = 0;

    /// direct access to discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;

    /// scaling factor for Lagrangian MFSI coupling (unused)
    virtual double DispIncrFactor() = 0;

    /// structure part of the rhs to enforce coupling condition at the interface
    /*!
      structure part of
      \f$ \mathbf{u}^{n+1}_{\Gamma} \Delta t = \mathbf{d}^{n+1}_{\Gamma} - \mathbf{d}^{n}_{\Gamma} \f$
    */
    //virtual Teuchos::RCP<Epetra_Vector> FluidCondRHS() const = 0;


    /// structure part of the rhs to enforce structure -- mesh coupling at the interface
    /*!
      \f$ \mathbf{d}^{G,n+1}_{\Gamma} = \mathbf{d}^{n+1}_{\Gamma} \f$
    */
    //virtual Teuchos::RCP<Epetra_Vector> MeshCondRHS() const = 0;

    /// communication object at the interface
    virtual const LINALG::MapExtractor& Interface() const = 0;

    //@}

    //! @name Time step helpers

    /// start new time step
    virtual void PrepareTimeStep() = 0;

    /// evaluate elements with given displacement
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disp) const = 0;

    /// update at time step end
    virtual void Update() = 0;

    /// output results
    virtual void Output() = 0;

    /// read restart information for given time step
    virtual void ReadRestart(int step) = 0;

    //@}

    //! @name Solver calls

    /// nonlinear solve
    /*!
      Do the nonlinear solve for the time step. All boundary conditions have
      been set.
     */
    virtual void Solve() = 0;

    /// linear structure solve with just a interface load
    /*!
      The very special solve done in steepest descent relaxation
      calculation (and matrix free Newton Krylov).

      \note Can only be called after a valid structural solve.
    */
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce) = 0;

    //@}

    //! @name Extract interface values

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn() = 0;

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispnp() = 0;

    /// Predictor for interface displacements
    virtual Teuchos::RCP<Epetra_Vector> PredictInterfaceDispnp() = 0;

    //@}

    //! @name Apply interface forces

    /// apply interface forces to structural solver
    /*!
      This prepares a new solve of the structural field within one time
      step. The middle values are newly created.

      \note This is not yet the most efficient implementation.
    */
    virtual void ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce) = 0;

    //@}
  };


  /// adapter to generalized alpha algorithm
  class StructureAdapter : public Structure
  {
  public:

    StructureAdapter(Teuchos::RCP<Teuchos::ParameterList> params,
                     Teuchos::RCP<DRT::Discretization> dis,
                     Teuchos::RCP<LINALG::Solver> solver,
                     Teuchos::RCP<IO::DiscretizationWriter> output);

    Teuchos::RCP<const Epetra_Vector> InitialGuess() const;
    Teuchos::RCP<const Epetra_Vector> RHS() const;
    Teuchos::RCP<const Epetra_Vector> Dispnp() const;
    Teuchos::RCP<const Epetra_Vector> Dispn() const;
    Teuchos::RCP<const Epetra_Vector> Dispnm() const;
    Teuchos::RCP<const Epetra_Map> DofRowMap() const;
    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() const;
    Teuchos::RCP<DRT::Discretization> Discretization();
    double DispIncrFactor();

    //virtual Teuchos::RCP<Epetra_Vector> FluidCondRHS() const;
    //virtual Teuchos::RCP<Epetra_Vector> MeshCondRHS() const;

    //Teuchos::RCP<Epetra_Vector> InterfaceDispn();

    void PrepareTimeStep();
    void Evaluate(Teuchos::RCP<const Epetra_Vector> disp) const;
    void Update();
    void Output();

    /// communication object at the interface
    virtual const LINALG::MapExtractor& Interface() const { return interface_; }

    virtual const Epetra_Map& DomainMap();

    virtual void ReadRestart(int step);

    virtual void Solve();

    /// linear structure solve with just a interface load
    /*!
      The very special solve done in steepest descent relaxation
      calculation (and matrix free Newton Krylov).

      \note Can only be called after a valid structural solve.
    */
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce);

    //! @name Extract interface values

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn();

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispnp();

    /// Predictor for interface displacements
    virtual Teuchos::RCP<Epetra_Vector> PredictInterfaceDispnp();

    //@}

    //! @name Apply interface forces

    /// apply interface forces to structural solver
    /*!
      This prepares a new solve of the structural field within one time
      step. The middle values are newly created.

      \note This is not yet the most efficient implementation.
    */
    virtual void ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce);

    //@}

  private:

    /// the actual structure algorithm
    mutable StruGenAlpha structure_;

    /// the interface map setup for interface <-> full translation
    LINALG::MapExtractor interface_;

    //! @name local copies of input parameters

    Teuchos::RCP<DRT::Discretization> dis_;
    Teuchos::RCP<Teuchos::ParameterList> params_;
    Teuchos::RCP<LINALG::Solver> solver_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    //@}
  };
}

#endif
#endif
