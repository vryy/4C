/*-----------------------------------------------------------*/
/*! \file

\brief Base class for all fluid time integrations


\level 1

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TIMINT_HPP
#define FOUR_C_FLUID_TIMINT_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_fluid.hpp"
#include "4C_fluid_discretization_runtime_output_params.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class Sparsematrix;
  class BlockSparseMatrixBase;
  class MapExtractor;
  class KrylovProjector;
}  // namespace CORE::LINALG

namespace DRT
{
  class Discretization;
  class DofSet;
  class Condition;
}  // namespace DRT

namespace IO
{
  class DiscretizationWriter;
  class DiscretizationVisualizationWriterMesh;
}  // namespace IO

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
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output);


    void Init() override = 0;

    Teuchos::RCP<const Epetra_Vector> InitialGuess() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> RHS() override = 0;
    Teuchos::RCP<const Epetra_Vector> TrueResidual() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<Epetra_Vector> WriteAccessVelnp() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> Velnp() override = 0;
    Teuchos::RCP<const Epetra_Vector> Velaf() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> Veln() override = 0;
    Teuchos::RCP<const Epetra_Vector> Velnm() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> Accnp() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> Accn() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> Accnm() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> Accam() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> Scaaf() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> Scaam() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> Hist() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> GridVel() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> GridVeln() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> Dispnp() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> Dispn() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> ConvectiveVel() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> FsVel() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<Epetra_Vector> StdVeln() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<Epetra_Vector> StdVelnp() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<Epetra_Vector> StdVelaf() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Map> DofRowMap() override { return DofRowMap(0); }
    Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds) override;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemMatrix() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemSparseMatrix() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> BlockSystemMatrix() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> ShapeDerivatives() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    const Teuchos::RCP<DRT::Discretization>& Discretization() override { return discret_; }
    Teuchos::RCP<const DRT::DofSet> DofSet() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> Stepinc() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const CORE::LINALG::MapExtractor> GetDBCMapExtractor() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    void Integrate() override = 0;
    void PrepareTimeStep() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    }
    void PrepareSolve() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    }
    void Evaluate(Teuchos::RCP<const Epetra_Vector> vel) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    }
    virtual bool ConvergenceCheck(int itnum, int itmax, const double ittol)
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return false;
    }
    void IterUpdate(const Teuchos::RCP<const Epetra_Vector> increment) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    }
    void Update() override = 0;
    void StatisticsAndOutput() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    }
    void Output() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    }
    void StatisticsOutput() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    }
    const Teuchos::RCP<IO::DiscretizationWriter>& DiscWriter() override { return output_; }
    Teuchos::RCP<CORE::LINALG::MapExtractor> GetVelPressSplitter() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }


    void Solve() override = 0;

    void CalcIntermediateSolution() override { FOUR_C_THROW("Not implemented in the base class"); }
    /// get the linear solver object used for this field
    Teuchos::RCP<CORE::LINALG::Solver> LinearSolver() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Map> InnerVelocityRowMap() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }  // only used for FSI
    Teuchos::RCP<const Epetra_Map> VelocityRowMap() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Map> PressureRowMap() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }

    /// preparations for Krylov space projection
    virtual void SetupKrylovSpaceProjection(DRT::Condition* kspcond)
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    }
    virtual void UpdateKrylovSpaceProjection()
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    }
    virtual void CheckMatrixNullspace()
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    }


    /// the mesh map contains all velocity dofs that are covered by an ALE node
    void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm, const int nds_master = 0) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /// scaling factor needed to convert the residual to real forces
    double ResidualScaling() const override = 0;

    double TimeScaling() const override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return 1337.0;
    }

    /// return time integration factor
    double TimIntParam() const override = 0;

    /// communication object at the interface (neglecting pressure dofs)
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& Interface() const override
    {
      FOUR_C_THROW("Implemented in the fluid wrapper and derived classes");
      static Teuchos::RCP<FLD::UTILS::MapExtractor> ret = Teuchos::null;
      return ret;
    }

    /// communication object at the interface needed for fpsi problems (including pressure dofs)
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& FPSIInterface() const override
    {
      FOUR_C_THROW("Implemented in the fluid wrapper and derived classes");
      static Teuchos::RCP<FLD::UTILS::MapExtractor> ret = Teuchos::null;
      return ret;
    }

    void ReadRestart(int step) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    void SetRestart(const int step, const double time, Teuchos::RCP<const Epetra_Vector> readvelnp,
        Teuchos::RCP<const Epetra_Vector> readveln, Teuchos::RCP<const Epetra_Vector> readvelnm,
        Teuchos::RCP<const Epetra_Vector> readaccnp,
        Teuchos::RCP<const Epetra_Vector> readaccn) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    double Time() const override { return time_; }
    int Step() const override { return step_; }
    double Dt() const override { return dta_; }

    //! increment time and step value
    void IncrementTimeAndStep() override;

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
    void TimeStepAuxiliar() override
    {
      FOUR_C_THROW(
          "We do this in the Adapter until time adaptivity is available in the fluid field.");
    }

    /*! Indicate norms of temporal discretization error
     *
     *  \author mayr.mt \date 12/2013
     */
    void IndicateErrorNorms(
        double& err,       ///< L2-norm of temporal discretization error based on all DOFs
        double& errcond,   ///< L2-norm of temporal discretization error based on interface DOFs
        double& errother,  ///< L2-norm of temporal discretization error based on interior DOFs
        double& errinf,    ///< L-inf-norm of temporal discretization error based on all DOFs
        double&
            errinfcond,  ///< L-inf-norm of temporal discretization error based on interface DOFs
        double& errinfother  ///< L-inf-norm of temporal discretization error based on interior DOFs
        ) override
    {
      FOUR_C_THROW(
          "We do this in the Adapter until time adaptivity is available in the fluid field.");
    }

    //@}

    //! set time step size
    void SetDt(const double dtnew) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    //! set time and step
    void SetTimeStep(const double time,  ///< time to set
        const int step                   ///< time step number to set
        ) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
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
    void ResetStep() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
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
    void ResetTime(const double dtold) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    virtual void LiftDrag() const = 0;
    double EvalTime() const override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return 0.0;
    }
    void Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForcesRobin()
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<Epetra_Vector> ExtractFreeSurfaceVeln() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }
    virtual void ApplyInterfaceRobinValue(
        Teuchos::RCP<Epetra_Vector> ivel, Teuchos::RCP<Epetra_Vector> iforce)
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }
    /// Apply initial mesh displacement
    void ApplyInitialMeshDisplacement(Teuchos::RCP<const Epetra_Vector> initfluiddisp) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }
    void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }
    void ApplyMeshDisplacementIncrement(Teuchos::RCP<const Epetra_Vector> dispstepinc) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }
    void ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }
    void VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    void FreeSurfDisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }
    void FreeSurfVelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    int Itemax() const override { return itemax_; }
    void SetItemax(int itemax) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /*!
    \brief return type of time integration scheme

    */
    INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const override { return timealgo_; }

    Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }

    void UseBlockMatrix(bool splitmatrix) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /// linear fluid solve with just a interface load
    /*!
      The very special solve done in steepest descent relaxation
      calculation (and matrix free Newton Krylov).

      \note Can only be called after a valid fluid solve.
    */
    Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }

    Teuchos::RCP<DRT::ResultTest> CreateFieldTest() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }

    Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(
        Teuchos::RCP<const Epetra_Vector> velpres) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }

    Teuchos::RCP<const Epetra_Vector> ExtractPressurePart(
        Teuchos::RCP<const Epetra_Vector> velpres) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }

    /// set initial flow field
    void SetInitialFlowField(
        const INPAR::FLUID::InitialField initfield, const int startfuncno) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    }
    /// set initial porosity field
    void SetInitialPorosityField(
        const INPAR::POROELAST::InitialField initfield, const int startfuncno) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /// apply external forces to the fluid
    void ApplyExternalForces(Teuchos::RCP<Epetra_MultiVector> fext) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /// apply external forces to the fluid
    void AddContributionToExternalLoads(
        const Teuchos::RCP<const Epetra_Vector> contributing_vector) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /// expand dirichlet set
    void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    };

    /// contract dirichlet set
    void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    ///  set scalar fields within outer iteration loop
    void SetIterScalarFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
        Teuchos::RCP<const Epetra_Vector> scalaram, Teuchos::RCP<const Epetra_Vector> scalardtam,
        Teuchos::RCP<DRT::Discretization> scatradis, int dofset) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    }

    void SetLomaIterScalarFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
        Teuchos::RCP<const Epetra_Vector> scalaram, Teuchos::RCP<const Epetra_Vector> scalardtam,
        Teuchos::RCP<const Epetra_Vector> fsscalaraf, const double thermpressaf,
        const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
        Teuchos::RCP<DRT::Discretization> scatradis) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    }

    /// set scalar fields
    void SetScalarFields(Teuchos::RCP<const Epetra_Vector> scalarnp, const double thermpressnp,
        Teuchos::RCP<const Epetra_Vector> scatraresidual,
        Teuchos::RCP<DRT::Discretization> scatradis, const int whichscalar = -1) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    }

    /// set velocity field (separate computation)
    void SetVelocityField(Teuchos::RCP<const Epetra_Vector> velnp) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /// provide access to the turbulence statistic manager
    Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    };
    /// provide access to the box filter for dynamic Smagorinsky model
    Teuchos::RCP<FLD::DynSmagFilter> DynSmagFilter() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    };
    Teuchos::RCP<FLD::Vreman> Vreman() override { return Teuchos::null; };

    /// update velocity increment after Newton step
    void UpdateNewton(Teuchos::RCP<const Epetra_Vector> vel) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    /// reset data for restart of simulation at beginning
    void Reset(bool completeReset = false, int numsteps = 1, int iter = -1) override
    {
      FOUR_C_THROW("reset function not implemented for this fluid adapter");
    };

    // set fluid displacement vector due to biofilm growth
    void SetFldGrDisp(Teuchos::RCP<Epetra_Vector> fluid_growth_disp) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    void CalculateError() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return;
    }

    INPAR::FLUID::PhysicalType PhysicalType() const override { return physicaltype_; }

   protected:
    //! fluid discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! linear solver
    Teuchos::RCP<CORE::LINALG::Solver> solver_;

    //! parameter list
    Teuchos::RCP<Teuchos::ParameterList> params_;

    //! output writer
    Teuchos::RCP<IO::DiscretizationWriter> output_;


    /// runtime output writer
    Teuchos::RCP<IO::DiscretizationVisualizationWriterMesh> runtime_output_writer_;

    /// runtime output parameter
    DRT::ELEMENTS::FluidRuntimeOutputParams runtime_output_params_;

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
    Teuchos::RCP<CORE::LINALG::KrylovProjector> projector_;

    /// Krylov space projection map extractor
    Teuchos::RCP<FLD::UTILS::KSPMapExtractor> kspsplitter_;

    //@}

  };  // class TimInt

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
