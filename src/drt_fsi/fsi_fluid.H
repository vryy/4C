/*----------------------------------------------------------------------*/
/*!
\file fsi_fluid.H

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

#ifndef FSI_FLUID_H
#define FSI_FLUID_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluidresulttest.H"

#include "fsi_utils.H"


namespace FSI {

  /// general fluid field interface for FSI
  /*!

  This is the FSI algorithm's view on a fluid algorithm. This pure virtual
  interface contains all the methods anny FSI algorithm might want to
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

  \warning A cleanup is needed badly.

  \author u.kue
  \date 11/07
  */
  class Fluid
  {
  public:

    /// virtual destructor to get polymorph destruction
    virtual ~Fluid();

    virtual Teuchos::RCP<Epetra_Vector> InitialGuess() = 0;
    virtual Teuchos::RCP<Epetra_Vector> RHS() = 0;
    virtual Teuchos::RCP<Epetra_Vector> Velnp() = 0;
    virtual Teuchos::RCP<Epetra_Vector> Veln() = 0;
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() = 0;
    virtual Teuchos::RCP<Epetra_Operator> SysMat() = 0;
    virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;

    virtual Teuchos::RCP<Epetra_Vector> StructCondRHS() = 0;

    virtual void PrepareTimeStep() = 0;
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> vel) = 0;
    virtual void Update() = 0;
    virtual void Output() = 0;

    /// solve the nonlinear problem
    virtual void NonlinearSolve() = 0;

    /// the interface map contains all velocity dofs at the interface
    virtual Teuchos::RCP<Epetra_Map> InterfaceMap() = 0;
    virtual Teuchos::RCP<Epetra_Map> InnerVelocityRowMap() = 0;
    virtual Teuchos::RCP<Epetra_Map> PressureRowMap() = 0;

    /// the mesh map contains all velocity dofs that are covered by an ALE node
    virtual void SetMeshMap(Teuchos::RCP<Epetra_Map> mm) = 0;

    virtual double ResidualScaling() = 0;

    /// communication object at the interface
    virtual const LINALG::MapExtractor& Interface() const = 0;

    virtual void ReadRestart(int step) = 0;
    virtual double Time() = 0;
    virtual int Step() = 0;
    virtual void LiftDrag() = 0;

    //! @name Extract interface forces

    /// After the fluid solve we need the forces at the FSI interface.
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces() = 0;

    //@}

    //! @name Extract interface values

    /// at the interface the velocity is prescribed as a Dirichlet condition
    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel) = 0;

    //@}

    //! @name Extract mesh values

    /// tell the mesh displacement to the fluid solver
    virtual void ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp) = 0;

    /// tell the mesh velocity to the fluid solver
    virtual void ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel) = 0;

    //@}

    //! @name Number of Newton iterations
    //! For simplified FD MFNK solve we want to temporally limit the
    /// number of Newton steps inside the fluid solver

    virtual int Itemax() const = 0;
    virtual void SetItemax(int itemax) = 0;

    //@}

    /// integrate FSI interface shape functions
    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape() = 0;

    /// linear fluid solve with just a interface load
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel) = 0;

    /// create result test for encapulated fluid algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest() = 0;
  };


  /// adapter to normal fluid algorithm
  class FluidAdapter : public Fluid
  {
  public:
    FluidAdapter(Teuchos::RCP<DRT::Discretization> dis,
                 Teuchos::RCP<LINALG::Solver> solver,
                 Teuchos::RCP<ParameterList> params,
                 Teuchos::RCP<IO::DiscretizationWriter> output);

    virtual Teuchos::RCP<Epetra_Vector> InitialGuess();
    virtual Teuchos::RCP<Epetra_Vector> RHS();
    virtual Teuchos::RCP<Epetra_Vector> Velnp();
    virtual Teuchos::RCP<Epetra_Vector> Veln();
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap();
    virtual Teuchos::RCP<Epetra_Operator> SysMat();
    virtual Teuchos::RCP<DRT::Discretization> Discretization();

    virtual Teuchos::RCP<Epetra_Vector> StructCondRHS();

    virtual void PrepareTimeStep();
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> vel);
    virtual void Update();
    virtual void Output();

    virtual void NonlinearSolve();

    /// the interface map contains all velocity dofs at the interface
    virtual Teuchos::RCP<Epetra_Map> InterfaceMap();
    virtual Teuchos::RCP<Epetra_Map> InnerVelocityRowMap();
    virtual Teuchos::RCP<Epetra_Map> PressureRowMap();

    /// the mesh map contains all velocity dofs that are covered by an ALE node
    virtual void SetMeshMap(Teuchos::RCP<Epetra_Map> mm);

    /// scaling factor needed to convert the residual to real forces
    virtual double ResidualScaling();

    /// communication object at the interface
    virtual const LINALG::MapExtractor& Interface() const { return interface_; }

    virtual void ReadRestart(int step);
    virtual double Time();
    virtual int Step();
    virtual void LiftDrag();

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();
    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel);
    virtual void ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp);
    virtual void ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel);

    virtual int Itemax() const;
    virtual void SetItemax(int itemax);

    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape();

    /// linear fluid solve with just a interface load
    /*!
      The very special solve done in steepest descent relaxation
      calculation (and matrix free Newton Krylov).

      \note Can only be called after a valid fluid solve.
    */
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel);

    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  private:

    /// the actual fluid algorithm
    FluidImplicitTimeInt fluid_;

    /// the interface map setup for interface <-> full translation
    LINALG::MapExtractor interface_;

    /// ALE dof map
    LINALG::MapExtractor meshmap_;

    //! @name local copies of input parameters
    Teuchos::RCP<DRT::Discretization> dis_;
    Teuchos::RCP<LINALG::Solver> solver_;
    Teuchos::RCP<Teuchos::ParameterList> params_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;
    //@}

    /// all velocity dofs not at the interface
    Teuchos::RefCountPtr<Epetra_Map> innervelmap_;
  };

}

#endif
#endif
