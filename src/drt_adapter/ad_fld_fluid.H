/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid.H

\brief Fluid field adapter

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
#ifndef ADAPTER_FLUID_H
#define ADAPTER_FLUID_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsGraph.h>

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_combust.H"
#include "../drt_inpar/inpar_poroelast.H"
#include "../drt_lib/drt_dserror.H"


namespace LINALG
{
  class SparseMatrix;
  class BlockSparseMatrixBase;
  class MapExtractor;
  class Solver;
}

namespace DRT
{
  class Discretization;
  class DofSet;
  class ResultTest;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace FLD
{
  class TurbulenceStatisticManager;
  class DynSmagFilter;
  namespace UTILS
  {
    class MapExtractor;
    class FluidXFluidMapExtractor;
  }
}

namespace COMBUST
{
  class InterfaceHandleCombust;
  class FlameFront;
}

namespace TOPOPT
{
  class Optimizer;
}

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

    /// true residual
    virtual Teuchos::RCP<const Epetra_Vector> TrueResidual() = 0;

    /// velocities (and pressures) at t(n+1)
    virtual Teuchos::RCP<const Epetra_Vector> Velnp() = 0;

    /// velocities (and pressures) at t(n+alpha_F)
    virtual Teuchos::RCP<const Epetra_Vector> Velaf() = 0;

    /// velocities (and pressures) at t(n)
    virtual Teuchos::RCP<const Epetra_Vector> Veln() = 0;

    /// velocities (and pressures) at t(n-1)
    virtual Teuchos::RCP<const Epetra_Vector> Velnm() = 0;

    /// accelerations at t(n+1)
    virtual Teuchos::RCP<const Epetra_Vector> Accnp() = 0;

    /// accelerations at t(n)
    virtual Teuchos::RCP<const Epetra_Vector> Accn() = 0;

    /// accelerations at t(n+alpha_M)
    virtual Teuchos::RCP<const Epetra_Vector> Accam() = 0;

    /// scalars at t(n+alpha_F)
    virtual Teuchos::RCP<const Epetra_Vector> Scaaf() = 0;

    /// scalars at t(n+alpha_M)
    virtual Teuchos::RCP<const Epetra_Vector> Scaam() = 0;

    /// history vector
    virtual Teuchos::RCP<const Epetra_Vector> Hist() = 0;

    /// mesh displacements at t(n+1)
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp() = 0;

    /// mesh displacements at t(n+1)
    virtual Teuchos::RCP<const Epetra_Vector> Dispn() = 0;

    /// convective velocity (= velnp - grid velocity)
    virtual Teuchos::RCP<const Epetra_Vector> ConvectiveVel() = 0;

    virtual Teuchos::RCP<const Epetra_Vector> GridVel() = 0;

    /// fine-scale velocity
    virtual Teuchos::RCP<const Epetra_Vector> FsVel() = 0;

    /// velocities (and pressures) at t(n) w/out enriched dofs
    virtual Teuchos::RCP<Epetra_Vector> StdVeln() = 0;

    /// velocities (and pressures) at t(n+1) w/out enriched dofs
    virtual Teuchos::RCP<Epetra_Vector> StdVelnp() = 0;

    /// velocities (and pressures) at t(n+alpha_F) w/out enriched dofs
    virtual Teuchos::RCP<Epetra_Vector> StdVelaf() = 0;

    //@}

    //! @name Misc

    /// dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() = 0;

    /// dof map of vector of unknowns for multiple dofsets
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds) = 0;

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() = 0;

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() = 0;

    /// linearization of Navier-Stokes with respect to mesh movement
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives() = 0;

    /// direct access to discretization
    virtual const Teuchos::RCP<DRT::Discretization>& Discretization() = 0;

    /// direct access to dofset
    virtual Teuchos::RCP<const DRT::DofSet> DofSet() = 0;

    /// Return MapExtractor for Dirichlet boundary conditions
    virtual Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() = 0;

    // should it be pure virtual?!
    /// Return dirichlet boundray condition (used for fluidfluid)
    virtual const Teuchos::RCP<const Epetra_Map>& FluidDirichMaps() = 0;

    /// Shorten the Dirichlet DOF set
    virtual void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove) = 0;

    //virtual Teuchos::RCP<Epetra_Vector> StructCondRHS() const = 0;

    /// set initial flow field
    virtual void SetInitialFlowField(const INPAR::FLUID::InitialField initfield,const int startfuncno) = 0;

    /// set initial flow field
    virtual void SetInitialPorosityField(const INPAR::POROELAST::InitialField initfield,const int startfuncno) = 0;

    ///  set fields for low-Mach-number flow within outer iteration loop
    virtual void SetIterLomaFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
                                   Teuchos::RCP<const Epetra_Vector> scalaram,
                                   Teuchos::RCP<const Epetra_Vector> scalardtam,
                                   Teuchos::RCP<const Epetra_Vector> fsscalaraf,
                                   const double             thermpressaf,
                                   const double             thermpressam,
                                   const double             thermpressdtaf,
                                   const double             thermpressdtam,
                                   Teuchos::RCP<DRT::Discretization> scatradis) = 0;

    /// set fields for low-Mach-number flow at end of time step
    virtual void SetTimeLomaFields(Teuchos::RCP<const Epetra_Vector> scalarnp,
                                   const double             thermpressnp,
                                   Teuchos::RCP<const Epetra_Vector> scatraresidual,
                                   Teuchos::RCP<DRT::Discretization> scatradis,
                                   const int                whichscalar = -1) = 0;

    /// set velocity field (separate computation)
    virtual void SetVelocityField(Teuchos::RCP<const Epetra_Vector> velnp) = 0;

    /// provide access to the turbulence statistic manager
    virtual Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager() = 0;
    /// provide access to the box filter class for dynamic Smaorinsky model
    virtual Teuchos::RCP<FLD::DynSmagFilter> DynSmagFilter() = 0;

    /// reset state vectors (needed for biofilm simulations)
    virtual void Reset() { dserror("reset function not implemented for this fluid adapter"); };

    //@}

    //! @name Time step helpers

    /// run a complete simulation (for fluid stand-alone simulations)
//    virtual void TimeLoop() = 0;
    virtual void Integrate() = 0;

    /// start new time step
    virtual void PrepareTimeStep() = 0;
    virtual void PrepareTimeandStep() = 0;

    /// preparatives for solve
    virtual void PrepareSolve() = 0;

    /// update fluid unknowns and evaluate elements
    ///
    /// there are two increments possible
    /// x^n+1_i+1 = x^n+1_i + iterinc, and
    ///
    /// x^n+1_i+1 = x^n     + stepinc
    ///
    /// with n and i being time and Newton iteration step
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> stepinc ///< increment between time step n and n+1
        ) = 0;

    /// convergence check
    virtual bool ConvergenceCheck(int itnum, int itmax, const double ittol) = 0;

    /// update at end of iteration step
    virtual void IterUpdate(const Teuchos::RCP<const Epetra_Vector> increment) = 0;

    /// update at end of time step
    virtual void Update() = 0;

    /// update velocity increment after Newton step
    virtual void UpdateNewton(Teuchos::RCP<const Epetra_Vector> vel) = 0;

    /// lift'n'drag forces, statistics time sample and
    /// output of solution and statistics
    virtual void StatisticsAndOutput() = 0;

    /// output results
    virtual void Output() = 0;

    /// output statistics
    virtual void StatisticsOutput() = 0;

    /// access to output
    virtual const Teuchos::RCP<IO::DiscretizationWriter>& DiscWriter() = 0;

    /// read restart information for given time step
    virtual void ReadRestart(int step) = 0;
    /// set restart
    virtual void SetRestart(const int step,
                    const double time,
                    Teuchos::RCP<const Epetra_Vector> readvelnp,
                    Teuchos::RCP<const Epetra_Vector> readveln,
                    Teuchos::RCP<const Epetra_Vector> readvelnm,
                    Teuchos::RCP<const Epetra_Vector> readaccnp,
                    Teuchos::RCP<const Epetra_Vector> readaccn) = 0;

    /// current time value
    virtual double Time() const = 0;

    /// current time step
    virtual int Step() const = 0;

    /// time step size
    virtual double Dt() const = 0;

    /// calculation (and output) of lift&drag
    virtual void LiftDrag() const = 0;

    /// this procs element evaluate time
    virtual const double EvalTime() const = 0;

    /// redistribute the fluid discretization and vectors according to nodegraph in std. mode
    virtual void Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph) = 0;

    /// rewrite vectors to the new dofmaps after redistributing the discretization
    virtual void TransferVectorsToNewDistribution(Teuchos::RCP<COMBUST::FlameFront> flamefront) = 0;

    //@}

    //! @name Solver calls

    /// nonlinear solve
    /*!
      Do the nonlinear solve for the time step. All boundary conditions have
      been set.
     */
    virtual void NonlinearSolve() = 0;

    /// predictor
    virtual void Predictor() = 0;

    /// (multiple) corrector
    virtual void MultiCorrector() = 0;

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

    /// Use ResidualScaling() to convert the implemented fluid residual to an actual force with unit Newton [N]
    virtual double ResidualScaling() const = 0;

    /// Velocity-displacement conversion at the fsi interface
    /*! Time integration of the fsi interface reads:
     *  \f$\mathbf{d}^{n+1} = \mathbf{d}^{n} + \tau(\mathbf{u}^{n+1}-\mathbf{u}^{n}) + \Delta t \mathbf{u}^{n}\f$
     *
     *  Currently, two time integration schemes for the fsi interface
     *  are implemented:
     *  - Backward-Euler: \f$\tau = \Delta t\f$
     *  - Trapezoidal rule: \f$\tau = \frac{\Delta t}{2}\f$
     *
     *  Use TimeScaling() to get \f$\tau=\frac{1}{\text{TimeScaling()}}\f$
     */
    virtual double TimeScaling() const = 0;

    /// return time integration factor
    virtual double TimIntParam() const = 0;

    /// communication object at the interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor>const& Interface() const = 0;

    /// return type of time integration scheme
    virtual const INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const = 0;

    //! @name Extract the velocity-related part of a fluid vector (e.g. velnp, veln, residual)
    /// The idea is to have one function that does the extraction and call it
    /// with different vectors.

    /// Some applications need only access to velocity-related values of an fluid result vector.
    virtual Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres) = 0;

    //@}

    //! @name Transfer the porosity field used in the topology optimization to the fluid

    /// set pseudo-porosity field in the fluid
    virtual void SetTopOptData(Teuchos::RCP<const Epetra_Vector> porosity,Teuchos::RCP<TOPOPT::Optimizer>& optimizer) = 0;
    //@}

    //! @name Import geometrical coupling information about an implicitely defined interface from a different field

    /// A fluid in multi field problems which involve an implicitly defined interface, e.g. via level set
    /// methods, might need to know the interface geometry
    virtual void ImportFlameFront(const Teuchos::RCP<COMBUST::FlameFront> flamefront, bool UpdateDofSet) = 0;

    /// set initial flow field for combustion problem
    virtual void SetInitialFlowField(const INPAR::COMBUST::InitialField initfield, const int initfuncno) = 0;

    /// clear time integration data for XFEM problems
    virtual void ClearTimeInt() = 0;

    //@}

    //! @name Apply interface values

    /// at the interface the velocity is prescribed as a Dirichlet condition
    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel) = 0;

    //@}

    //! @name Extract interface values
    /// Maybe we do not need all of them?

    /// extract fluid velocity at the interface from time step n+1
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp() = 0;

    /// extract fluid velocity at the interface from time step n
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln() = 0;

    /// extract fluid velocity at the free surface from time step n
    virtual Teuchos::RCP<Epetra_Vector> ExtractFreeSurfaceVeln() = 0;

    /// extract fluid forces at the interface
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces() = 0;

    virtual Teuchos::RCP<FLD::UTILS::FluidXFluidMapExtractor>const& XFluidFluidMapExtractor() = 0;

    virtual Teuchos::RCP<const Epetra_Vector> Stepinc() = 0;

    //@}


    //! @name Extract mesh values

    /// tell the mesh displacement to the fluid solver
    virtual void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp) = 0;

    /// tell the mesh displacement step increment to the fluid solver
    virtual void ApplyMeshDisplacementIncrement(Teuchos::RCP<const Epetra_Vector> dispstepinc) = 0;

    /// tell the mesh velocity to the fluid solver
    virtual void ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel) = 0;

    //@}

    /// convert Delta d(n+1,i+1) to the fluid unknown at the interface
    virtual void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx) = 0;

    /// Conversion of displacement to velocity at the interface considering predictors or inhomogeneous DBCs
    virtual void DisplacementToVelocity(
        Teuchos::RCP<Epetra_Vector> fcx,    ///< fluid displacement step increment at interface
        Teuchos::RCP<Epetra_Vector> ddgpre, ///< displacement increment due to structural predictor
        Teuchos::RCP<Epetra_Vector> dugpre  ///< velocity increment due to fluid predictor
    ){ dserror("Not implemented in base class."); }

    /// convert the fluid unknown to Delta d(n+1,i+1) at the interface
    virtual void VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx) = 0;

    /// Conversion of velocity to displacement at the interface considering predictors or inhomogeneous DBCs
    virtual void VelocityToDisplacement(
        Teuchos::RCP<Epetra_Vector> fcx,    ///< fluid velocity step increment at interface
        Teuchos::RCP<Epetra_Vector> ddgpre, ///< displacement increment due to structural predictor
        Teuchos::RCP<Epetra_Vector> dugpre  ///< velocity increment due to fluid predictor
    ){ dserror("Not implemented in base class."); }

    /// convert Delta d(n+1,i+1) to the fluid unknown at the free surface
    virtual void FreeSurfDisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx) = 0;

    /// convert the fluid unknown to Delta d(n+1,i+1) at the free surface
    virtual void FreeSurfVelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx) = 0;

    //! @name Number of Newton iterations
    //! For simplified FD MFNK solve we want to temporally limit the
    /// number of Newton steps inside the fluid solver

    /// return maximum for iteration steps
    virtual int Itemax() const = 0;
    /// set maximum for iteration steps
    virtual void SetItemax(int itemax) = 0;

    //@}

    /// integrate FSI interface shape functions
    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape() = 0;

    /// switch fluid field to block matrix
    virtual void UseBlockMatrix(bool splitmatrix) = 0;

    /// create result test for encapulated fluid algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest() = 0;
  };
}

#endif
