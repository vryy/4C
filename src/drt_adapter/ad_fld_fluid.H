/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter

\level 1

\maintainer Martin Kronbichler
*/
#ifndef ADAPTER_FLUID_H
#define ADAPTER_FLUID_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsGraph.h>

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_poroelast.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_fixedsizematrix.H"

// forward declarations
namespace LINALG
{
  class SparseMatrix;
  class BlockSparseMatrixBase;
  class MapExtractor;
  class Solver;
}  // namespace LINALG

namespace DRT
{
  class Discretization;
  class DofSet;
  class ResultTest;
  class Element;
}  // namespace DRT

namespace IO
{
  class DiscretizationWriter;
}

namespace FLD
{
  class TurbulenceStatisticManager;
  class DynSmagFilter;
  class Vreman;
  namespace UTILS
  {
    class MapExtractor;
  }
}  // namespace FLD

namespace ADAPTER
{
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

    /// initialize time integration
    virtual void Init() = 0;

    //! @name Vector access

    /// initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() = 0;

    /// rhs of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> RHS() = 0;

    /// true residual
    virtual Teuchos::RCP<const Epetra_Vector> TrueResidual() = 0;

    /// velocities (and pressures) at \f$t^{n+1}\f$ for write access
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelnp() = 0;

    /// velocities (and pressures) at \f$t^{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Velnp() = 0;

    /// velocities (and pressures) at \f$t^{n+\alpha_F}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Velaf() = 0;

    /// velocities (and pressures) at \f$t^n\f$
    virtual Teuchos::RCP<const Epetra_Vector> Veln() = 0;

    /// velocities (and pressures) at \f$t^{n-1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Velnm() = 0;

    /// accelerations at \f$t^{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Accnp() = 0;

    /// accelerations at \f$t^n\f$
    virtual Teuchos::RCP<const Epetra_Vector> Accn() = 0;

    /// accelerations at \f$t^{n-1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Accnm() = 0;

    /// accelerations at \f$t^{n+\alpha_M}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Accam() = 0;

    /// scalars at \f$t^{n+\alpha_F}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Scaaf() = 0;

    /// scalars at \f$t^{n+\alpha_M}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Scaam() = 0;

    /// history vector
    virtual Teuchos::RCP<const Epetra_Vector> Hist() = 0;

    /// mesh displacements at \f$t^{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp() = 0;

    /// mesh displacements at \f$t^{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispn() = 0;

    /// convective velocity (= velnp - grid velocity)
    virtual Teuchos::RCP<const Epetra_Vector> ConvectiveVel() = 0;

    /// grid velocity at \f$t^{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> GridVel() = 0;

    /// grid velocity at \f$t^{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> GridVeln() = 0;

    /// fine-scale velocity
    virtual Teuchos::RCP<const Epetra_Vector> FsVel() = 0;

    /// velocities (and pressures) at \f$t^{n}\f$ w/out enriched dofs
    virtual Teuchos::RCP<Epetra_Vector> StdVeln() = 0;

    /// velocities (and pressures) at \f$t^{n+1}\f$ w/out enriched dofs
    virtual Teuchos::RCP<Epetra_Vector> StdVelnp() = 0;

    /// velocities (and pressures) at \f$t^{n+\alpha_F}\f$ w/out enriched dofs
    virtual Teuchos::RCP<Epetra_Vector> StdVelaf() = 0;

    //@}

    //! @name Misc

    /// dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() = 0;

    /// dof map of vector of unknowns for multiple dofsets
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds) = 0;

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() = 0;

    /// direct access to merged system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemSparseMatrix() = 0;

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

    /// set initial flow field
    virtual void SetInitialFlowField(
        const INPAR::FLUID::InitialField initfield, const int startfuncno) = 0;

    /// set initial flow field
    virtual void SetInitialPorosityField(
        const INPAR::POROELAST::InitialField initfield, const int startfuncno) = 0;

    /// apply external forces to the fluid
    virtual void ApplyExternalForces(Teuchos::RCP<Epetra_MultiVector> fext) = 0;

    /// apply contribution to neumann loads of the fluid (similar to ApplyExternalForces but without
    /// residual scaling)
    virtual void AddContributionToExternalLoads(
        const Teuchos::RCP<const Epetra_Vector> contributing_vector) = 0;

    /// expand dirichlet dbc set by provided map containing dofs to add
    virtual void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd) = 0;

    /// contract dirichlet set by provided map containing dofs to remove
    virtual void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove) = 0;

    ///  set scalar fields within outer iteration loop
    virtual void SetIterScalarFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
        Teuchos::RCP<const Epetra_Vector> scalaram, Teuchos::RCP<const Epetra_Vector> scalardtam,
        Teuchos::RCP<DRT::Discretization> scatradis, int dofset = 0) = 0;

    /// set scalar fields within outer iteration loop for low-Mach-number flow
    virtual void SetLomaIterScalarFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
        Teuchos::RCP<const Epetra_Vector> scalaram, Teuchos::RCP<const Epetra_Vector> scalardtam,
        Teuchos::RCP<const Epetra_Vector> fsscalaraf, const double thermpressaf,
        const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
        Teuchos::RCP<DRT::Discretization> scatradis) = 0;

    /// set scalar fields
    virtual void SetScalarFields(Teuchos::RCP<const Epetra_Vector> scalarnp,
        const double thermpressnp, Teuchos::RCP<const Epetra_Vector> scatraresidual,
        Teuchos::RCP<DRT::Discretization> scatradis, const int whichscalar = -1) = 0;

    /// set velocity field (separate computation)
    virtual void SetVelocityField(Teuchos::RCP<const Epetra_Vector> velnp) = 0;

    /// provide access to the turbulence statistic manager
    virtual Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager() = 0;
    /// provide access to the box filter class for dynamic Smaorinsky model
    virtual Teuchos::RCP<FLD::DynSmagFilter> DynSmagFilter() = 0;
    virtual Teuchos::RCP<FLD::Vreman> Vreman() = 0;

    /// reset state vectors (needed for biofilm simulations)
    virtual void Reset(bool completeReset = false, int numsteps = 1, int iter = -1) = 0;

    /// set fluid displacement vector due to biofilm growth
    virtual void SetFldGrDisp(Teuchos::RCP<Epetra_Vector> fluid_growth_disp) = 0;
    //@}

    //! @name Time step helpers

    /// run a complete simulation (for fluid stand-alone simulations)
    //    virtual void TimeLoop() = 0;
    virtual void Integrate() = 0;

    /// start new time step
    virtual void PrepareTimeStep() = 0;

    /// increment the time and the step
    virtual void IncrementTimeAndStep() = 0;

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
        Teuchos::RCP<const Epetra_Vector> stepinc  ///< increment between time step n and n+1
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

    /// access to map extractor for velocity and pressure
    virtual const Teuchos::RCP<LINALG::MapExtractor> GetVelPressSplitter() = 0;

    /// read restart information for given time step
    virtual void ReadRestart(int step) = 0;

    /// set restart
    virtual void SetRestart(const int step, const double time,
        Teuchos::RCP<const Epetra_Vector> readvelnp, Teuchos::RCP<const Epetra_Vector> readveln,
        Teuchos::RCP<const Epetra_Vector> readvelnm, Teuchos::RCP<const Epetra_Vector> readaccnp,
        Teuchos::RCP<const Epetra_Vector> readaccn) = 0;

    /// current time value
    virtual double Time() const = 0;

    /// current time step
    virtual int Step() const = 0;

    /// time step size
    virtual double Dt() const = 0;

    //! @name Time step size adaptivity in monolithic FSI
    //@{

    /*! Do one step with auxiliary time integration scheme
     *
     *  Do a single time step with the user given auxiliary time integration
     *  scheme. Result is stored in \p locerrvelnp_ and is used later to
     *  estimate the local discretization error of the marching time integration
     *  scheme.
     *
     *  \author mayr.mt \date 12/2013
     */
    virtual void TimeStepAuxiliar() = 0;

    /*! Indicate norms of local discretization error
     *
     *  \author mayr.mt \date 12/2013
     */
    virtual void IndicateErrorNorms(
        double& err,       ///< L2-norm of temporal discretization error based on all DOFs
        double& errcond,   ///< L2-norm of temporal discretization error based on interface DOFs
        double& errother,  ///< L2-norm of temporal discretization error based on interior DOFs
        double& errinf,    ///< L-inf-norm of temporal discretization error based on all DOFs
        double&
            errinfcond,  ///< L-inf-norm of temporal discretization error based on interface DOFs
        double& errinfother  ///< L-inf-norm of temporal discretization error based on interior DOFs
        ) = 0;

    /// set time step size
    virtual void SetDt(const double dtold) = 0;

    /// set time and step
    virtual void SetTimeStep(const double time,  ///< time to set
        const int step                           ///< time step number to set
        ) = 0;

    //@}

    /*!
    \brief Reset time step

    In case of time step size adaptivity, time steps might have to be repeated.
    Therefore, we need to reset the solution back to the initial solution of the
    time step.

    \author mayr.mt
    \date 08/2013
    */
    virtual void ResetStep() = 0;

    /*!
    \brief Reset time and step in case that a time step has to be repeated

    Fluid field increments time and step at the beginning of a time step. If a time
    step has to be repeated, we need to take this into account and decrease time and
    step beforehand. They will be incremented right at the beginning of the repetition
    and, thus, everything will be fine. Currently, this is needed for time step size
    adaptivity in FSI.

    \author mayr.mt
    \date 08/2013
     */
    virtual void ResetTime(const double dtold) = 0;

    /// this procs element evaluate time
    virtual double EvalTime() const = 0;

    /// redistribute the fluid discretization and vectors according to nodegraph in std. mode
    virtual void Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph) = 0;


    //@}

    //! @name Solver calls

    /// nonlinear solve
    /*!
      Do the nonlinear solve for the time step. All boundary conditions have
      been set.
     */
    virtual void Solve() = 0;

    /// linear fluid solve with just a interface load
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel) = 0;

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() = 0;

    /// do an intermediate solution step
    virtual void CalcIntermediateSolution() = 0;

    //@}

    /// Map of all velocity dofs that are not Dirichlet-constrained
    virtual Teuchos::RCP<const Epetra_Map> InnerVelocityRowMap() = 0;

    /// Map of all velocity dofs
    virtual Teuchos::RCP<const Epetra_Map> VelocityRowMap() = 0;

    /// Map of all pressure dofs
    virtual Teuchos::RCP<const Epetra_Map> PressureRowMap() = 0;

    /// the mesh map contains all velocity dofs that are covered by an ALE node
    virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm) = 0;

    /// Use ResidualScaling() to convert the implemented fluid residual to an actual force with unit
    /// Newton [N]
    virtual double ResidualScaling() const = 0;

    /// Velocity-displacement conversion at the fsi interface
    /*! Time integration of the fsi interface reads:
     *  \f$\mathbf{d}^{n+1} = \mathbf{d}^{n} + \tau(\mathbf{u}^{n+1}-\mathbf{u}^{n}) + \Delta t
     * \mathbf{u}^{n}\f$
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

    /// communication object at the interface (neglecting pressure dofs)
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor> const& Interface() const = 0;

    /// communication object at the interface (including pressure dofs)
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor> const& FPSIInterface() const = 0;

    /// return type of time integration scheme
    virtual INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const = 0;

    //! @name Extract the velocity-related part of a fluid vector (e.g. velnp, veln, residual)
    /// The idea is to have one function that does the extraction and call it
    /// with different vectors.

    /// Some applications need only access to velocity-related values of an fluid result vector.
    virtual Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(
        Teuchos::RCP<const Epetra_Vector> velpres) = 0;

    /// Some applications need only access to pressure-related values of an fluid result vector.
    virtual Teuchos::RCP<const Epetra_Vector> ExtractPressurePart(
        Teuchos::RCP<const Epetra_Vector> velpres) = 0;

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

    //! @name Conversion between displacement and velocity at interface

    /// convert Delta d(n+1,i+1) to the fluid unknown at the interface
    virtual void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx) = 0;

    /// convert the fluid unknown to Delta d(n+1,i+1) at the interface
    virtual void VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx) = 0;

    /// convert Delta d(n+1,i+1) to the fluid unknown at the free surface
    virtual void FreeSurfDisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx) = 0;

    /// convert the fluid unknown to Delta d(n+1,i+1) at the free surface
    virtual void FreeSurfVelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx) = 0;

    //@}

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

    /// calculate error in comparison to analytical solution
    virtual void CalculateError() = 0;

    /// return physical type of fluid algorithm
    virtual INPAR::FLUID::PhysicalType PhysicalType() const = 0;
  };
}  // namespace ADAPTER

#endif
