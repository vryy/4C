/*-----------------------------------------------------------*/
/*! \file

\brief Base class for all fluid time integrations

\maintainer Martin Kronbichler

\level 1

*/
/*-----------------------------------------------------------*/

#ifndef FLUID_TIMINT_H
#define FLUID_TIMINT_H

#include "../drt_adapter/ad_fld_fluid.H"

#include "../drt_lib/drt_dserror.H"

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

// forward declarations
namespace LINALG
{
  class Sparsematrix;
  class BlockSparseMatrixBase;
  class MapExtractor;
  class KrylovProjector;

}  // namespace LINALG

namespace DRT
{
  class Discretization;
  class DofSet;
  class Condition;
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
    class KSPMapExtractor;
  }

  class TimInt : public ADAPTER::Fluid
  {
   public:
    TimInt(const Teuchos::RCP<DRT::Discretization>& discret,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output);

    /*!
    \brief Destructor

    */
    virtual ~TimInt();

    virtual void Init() = 0;

    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> RHS() = 0;
    virtual Teuchos::RCP<const Epetra_Vector> TrueResidual()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelnp()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> Velnp() = 0;
    virtual Teuchos::RCP<const Epetra_Vector> Velaf()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> Veln() = 0;
    virtual Teuchos::RCP<const Epetra_Vector> Velnm()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> Accnp()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> Accn()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> Accnm()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> Accam()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> Scaaf()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> Scaam()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> Hist()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> GridVel()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> GridVeln()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> Dispn()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> ConvectiveVel()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> FsVel()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<Epetra_Vector> StdVeln()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<Epetra_Vector> StdVelnp()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<Epetra_Vector> StdVelaf()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() { return DofRowMap(0); }
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds);
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemSparseMatrix()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual const Teuchos::RCP<DRT::Discretization>& Discretization() { return discret_; }
    virtual Teuchos::RCP<const DRT::DofSet> DofSet()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Vector> Stepinc()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual void Integrate() = 0;
    virtual void PrepareTimeStep()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    }
    virtual void PrepareSolve()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    }
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> vel)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    }
    virtual bool ConvergenceCheck(int itnum, int itmax, const double ittol)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return false;
    }
    virtual void IterUpdate(const Teuchos::RCP<const Epetra_Vector> increment)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    }
    virtual void Update() = 0;
    virtual void StatisticsAndOutput()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    }
    virtual void Output()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    }
    virtual void StatisticsOutput()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    }
    virtual const Teuchos::RCP<IO::DiscretizationWriter>& DiscWriter() { return output_; }
    virtual const Teuchos::RCP<LINALG::MapExtractor> GetVelPressSplitter()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }


    virtual void Solve() = 0;

    virtual void CalcIntermediateSolution() { dserror("Not implemented in the base class"); }
    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Map> InnerVelocityRowMap()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }  // only used for FSI
    virtual Teuchos::RCP<const Epetra_Map> VelocityRowMap()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<const Epetra_Map> PressureRowMap()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }

    /// preparations for Krylov space projection
    virtual void SetupKrylovSpaceProjection(DRT::Condition* kspcond)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    }
    virtual void UpdateKrylovSpaceProjection()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    }
    virtual void CheckMatrixNullspace()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    }


    /// the mesh map contains all velocity dofs that are covered by an ALE node
    virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm, const int nds_master = 0)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /// scaling factor needed to convert the residual to real forces
    virtual double ResidualScaling() const = 0;

    virtual double TimeScaling() const
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return 1337.0;
    }

    /// return time integration factor
    virtual double TimIntParam() const = 0;

    /// communication object at the interface (neglecting pressure dofs)
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor> const& Interface() const
    {
      dserror("Implemented in the fluid wrapper and derived classes");
      static Teuchos::RCP<FLD::UTILS::MapExtractor> ret = Teuchos::null;
      return ret;
    }

    /// communication object at the interface needed for fpsi problems (including pressure dofs)
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor> const& FPSIInterface() const
    {
      dserror("Implemented in the fluid wrapper and derived classes");
      static Teuchos::RCP<FLD::UTILS::MapExtractor> ret = Teuchos::null;
      return ret;
    }

    virtual void ReadRestart(int step)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    virtual void SetRestart(const int step, const double time,
        Teuchos::RCP<const Epetra_Vector> readvelnp, Teuchos::RCP<const Epetra_Vector> readveln,
        Teuchos::RCP<const Epetra_Vector> readvelnm, Teuchos::RCP<const Epetra_Vector> readaccnp,
        Teuchos::RCP<const Epetra_Vector> readaccn)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    virtual double Time() const { return time_; }
    virtual int Step() const { return step_; }
    virtual double Dt() const { return dta_; }

    //! increment time and step value
    virtual void IncrementTimeAndStep();

    //! @name Time step size adaptivity in monolithic FSI
    //@{

    /*! Do one step with auxiliary time integration scheme
     *
     *  Do a single time step with the user given auxiliary time integration
     *  scheme. Result is stored in #locerrvelnp_ and is used later to estimate
     *  the local discretization error of the marching time integration scheme.
     *
     *  \author mayr.mt \date 12/2013
     */
    virtual void TimeStepAuxiliar()
    {
      dserror("We do this in the Adapter until time adaptivity is available in the fluid field.");
    }

    /*! Indicate norms of temporal discretization error
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
    )
    {
      dserror("We do this in the Adapter until time adaptivity is available in the fluid field.");
    }

    //@}

    //! set time step size
    virtual void SetDt(const double dtnew)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    //! set time and step
    virtual void SetTimeStep(const double time,  ///< time to set
        const int step                           ///< time step number to set
    )
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /*!
    \brief Reset time step

    In case of time step size adaptivity, time steps might have to be repeated.
    Therefore, we need to reset the solution back to the initial solution of the
    time step.

    \author mayr.mt
    \date 08/2013
    */
    virtual void ResetStep()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

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
    virtual void ResetTime(const double dtold)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    virtual void LiftDrag() const = 0;
    virtual double EvalTime() const
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return 0.0;
    }
    virtual void Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForcesRobin()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<Epetra_Vector> ExtractFreeSurfaceVeln()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }
    virtual void ApplyInterfaceRobinValue(
        Teuchos::RCP<Epetra_Vector> ivel, Teuchos::RCP<Epetra_Vector> iforce)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }
    virtual void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }
    virtual void ApplyMeshDisplacementIncrement(Teuchos::RCP<const Epetra_Vector> dispstepinc)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }
    virtual void ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    virtual void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }
    virtual void VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    virtual void FreeSurfDisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }
    virtual void FreeSurfVelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    int Itemax() const { return itemax_; }
    virtual void SetItemax(int itemax)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /*!
    \brief return type of time integration scheme

    */
    INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const { return timealgo_; }

    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }

    virtual void UseBlockMatrix(bool splitmatrix)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /// linear fluid solve with just a interface load
    /*!
      The very special solve done in steepest descent relaxation
      calculation (and matrix free Newton Krylov).

      \note Can only be called after a valid fluid solve.
    */
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }

    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }

    virtual Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(
        Teuchos::RCP<const Epetra_Vector> velpres)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }

    virtual Teuchos::RCP<const Epetra_Vector> ExtractPressurePart(
        Teuchos::RCP<const Epetra_Vector> velpres)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }

    /// set initial flow field
    virtual void SetInitialFlowField(
        const INPAR::FLUID::InitialField initfield, const int startfuncno)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    }
    /// set initial porosity field
    virtual void SetInitialPorosityField(
        const INPAR::POROELAST::InitialField initfield, const int startfuncno)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /// apply external forces to the fluid
    virtual void ApplyExternalForces(Teuchos::RCP<Epetra_MultiVector> fext)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /// apply external forces to the fluid
    virtual void AddContributionToExternalLoads(
        const Teuchos::RCP<const Epetra_Vector> contributing_vector)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /// expand dirichlet set
    virtual void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    };

    /// contract dirichlet set
    virtual void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    ///  set scalar fields within outer iteration loop
    virtual void SetIterScalarFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
        Teuchos::RCP<const Epetra_Vector> scalaram, Teuchos::RCP<const Epetra_Vector> scalardtam,
        Teuchos::RCP<DRT::Discretization> scatradis, int dofset)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    }

    virtual void SetLomaIterScalarFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
        Teuchos::RCP<const Epetra_Vector> scalaram, Teuchos::RCP<const Epetra_Vector> scalardtam,
        Teuchos::RCP<const Epetra_Vector> fsscalaraf, const double thermpressaf,
        const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
        Teuchos::RCP<DRT::Discretization> scatradis)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    }

    /// set scalar fields
    virtual void SetScalarFields(Teuchos::RCP<const Epetra_Vector> scalarnp,
        const double thermpressnp, Teuchos::RCP<const Epetra_Vector> scatraresidual,
        Teuchos::RCP<DRT::Discretization> scatradis, const int whichscalar = -1)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    }

    /// set velocity field (separate computation)
    virtual void SetVelocityField(Teuchos::RCP<const Epetra_Vector> velnp)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /// provide access to the turbulence statistic manager
    virtual Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    };
    /// provide access to the box filter for dynamic Smagorinsky model
    virtual Teuchos::RCP<FLD::DynSmagFilter> DynSmagFilter()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    };
    virtual Teuchos::RCP<FLD::Vreman> Vreman() { return Teuchos::null; };

    /// update velocity increment after Newton step
    virtual void UpdateNewton(Teuchos::RCP<const Epetra_Vector> vel)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /// reset data for restart of simulation at beginning
    virtual void Reset(bool completeReset = false, int numsteps = 1, int iter = -1)
    {
      dserror("reset function not implemented for this fluid adapter");
    };

    // set fluid displacement vector due to biofilm growth
    virtual void SetFldGrDisp(Teuchos::RCP<Epetra_Vector> fluid_growth_disp)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    virtual void CalculateError()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    INPAR::FLUID::PhysicalType PhysicalType() const { return physicaltype_; }

   protected:
    //! fluid discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! linear solver
    Teuchos::RCP<LINALG::Solver> solver_;

    //! parameter list
    Teuchos::RCP<Teuchos::ParameterList> params_;

    //! output writer
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    //! @name Time loop stuff
    //@{

    double time_;  /// physical time
    int step_;     /// timestep number
    double dta_;   /// time step size of current time step

    int stepmax_;     ///< maximal number of timesteps
    double maxtime_;  ///< maximal physical computation time
    int itemax_;      /// maximum number of nonlinear iterations

    //@}

    int uprestart_;  ///< write restart data every uprestart_ steps
    int upres_;      ///< write result every upres_ steps

    INPAR::FLUID::TimeIntegrationScheme timealgo_;  ///< time algorithm flag
    INPAR::FLUID::PhysicalType
        physicaltype_;  ///< flag for physical type of fluid flow (standard: incompressible)

    int myrank_;  ///< the processor ID from the communicator

    //! @name variables for Krylov Space projection
    //@{

    //! flag setting whether Krylov projection needs to be updated or not
    bool updateprojection_;

    //! Krylov projector himself
    Teuchos::RCP<LINALG::KrylovProjector> projector_;

    /// Krylov space projection map extractor
    Teuchos::RCP<FLD::UTILS::KSPMapExtractor> kspsplitter_;

    //@}

  };  // class TimInt

}  // namespace FLD

#endif  // #ifndef FLUIDIMPLICITINTEGRATION_H
