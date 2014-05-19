/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_fluid_fsi.H

\brief
Fluid field adapter for embedded (ALE-)fluid-fluid problems using XFEM

<pre>
Maintainer: Raffaela Kruse
            kruse@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ADAPTER_FLUID_FLUID_FSI_H
#define ADAPTER_FLUID_FLUID_FSI_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include "ad_fld_wrapper.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_inpar/inpar_xfem.H"

namespace DRT
{
  class Condition;
}

namespace LINALG
{
  class Solver;
  class MapExtractor;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace FLD
{
  class XFluidFluid;
  namespace UTILS
  {
    class MapExtractor;
  }
}

namespace ADAPTER {

  class FluidFluidFSI : public FluidWrapper
  {

  public:

    /// constructor
    FluidFluidFSI(Teuchos::RCP<Fluid> fluid,
        Teuchos::RCP<DRT::Discretization> embfluiddis,
        Teuchos::RCP<DRT::Discretization> bgfluiddis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        bool isale,
        bool dirichletcond,
        bool monolithicfluidfluidfsi);

    /// initialize and prepare maps
    virtual void Init();

    /// save results of current time step, do XFEM cut and refresh the
    /// merged fluid map extractor
    virtual void Update();

    /// @name Accessors

    /// get time scaling from FSI interface interpolation
    virtual double TimeScaling() const;

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

    /// communication object at the interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor>const& Interface() const { return mergedfluidinterface_; }

    /// communication object at the interface (dummy)
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor>const& FPSIInterface() const
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      static Teuchos::RCP<FLD::UTILS::MapExtractor> ret = Teuchos::null;
      return ret;
    }

    /// get Dirichlet-free map of fluid velocity DOF
    virtual Teuchos::RCP<const Epetra_Map> InnerVelocityRowMap();

    /// get interface forces
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();

    /// get interface velocity at time step n+1
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp();

    /// get interface velocity at time step n
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();

    /// get Dirichlet map extractor for embedded fluid
    virtual Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor();

    /// get block system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

    /// get discretization writer for output
    virtual const Teuchos::RCP<IO::DiscretizationWriter>& DiscWriter();

    /// get the velocity row maps for merged background & embedded fluid;
    /// please note: in the context of monolithic fluid-fluid-fsi with fluid split,
    /// only the embedded velocity row maps are obtained
    virtual Teuchos::RCP<const Epetra_Map> VelocityRowMap();

    /// get map extractor for background/embedded fluid
    virtual Teuchos::RCP<FLD::UTILS::FluidXFluidMapExtractor>const& XFluidFluidMapExtractor();

    //@}

    /// set old ALE displacement state before relaxation
    virtual void ApplyEmbFixedMeshDisplacement(Teuchos::RCP<const Epetra_Vector> disp);

    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
    { dserror("ADAPTER::FluidFluidFSI::RelaxationSolve not implemented."); return Teuchos::null; }

    /// set ALE-mesh map
    virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm);

    /// apply FSI interface velocities
    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel);

    /// apply embedded fluid mesh displacement
    virtual void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp);

    virtual void ApplyMeshDisplacementIncrement(Teuchos::RCP<const Epetra_Vector> dispstepinc)
    {dserror("ApplyMeshDisplacementIncrement for XFluid-Fluid problems.");}

    /// apply embedded fluid mesh velocity
    virtual void ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel);

    /// conversion of displacement to velocity at the interface without predictors or inhomogeneous DBCs
    virtual void DisplacementToVelocity(
      Teuchos::RCP<Epetra_Vector> fcx ///< fluid displacement step increment
      );

    /// evaluate the fluid and update the merged fluid/FSI DOF-map extractor in case of a change in the DOF-maps
    virtual void Evaluate(
      Teuchos::RCP<const Epetra_Vector> stepinc ///< solution increment between time step n and n+1
      );

    /// conversion of velocity to displacement at the interface without predictors or inhomogeneous DBCs
    virtual void VelocityToDisplacement(
      Teuchos::RCP<Epetra_Vector> fcx     ///< fluid velocity step increment at interface
      );

    /// request fluid system matrix & shapederivatives as blockmatrices when called with true
    /// (required in fluidsplit xffsi)
    virtual void UseBlockMatrix(bool splitmatrix = false);

    /// remove Dirichlet conditions: legacy method for the case of coincident fluid/structure FSI nodes
    virtual void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove);

  private:

    /// setup of map extractor to distinguish between FSI DOF-map and
    /// merged inner embedded fluid and background fluid DOF-map
    void SetupInterfaceExtractor();

    /// type cast pointer to XFluidFluid
    Teuchos::RCP<FLD::XFluidFluid>         xfluidfluid_;

    /// map extractor for the embedded fluid discretization
    Teuchos::RCP<FLD::UTILS::MapExtractor> interface_;

    /// map extractor for merged fluid maps (to keep fsi interface-DOF apart from
    /// merged inner DOF (inner embedded fluid together with background fluid)
    Teuchos::RCP<FLD::UTILS::MapExtractor> mergedfluidinterface_;

    /// interface force at old time level t_n
    Teuchos::RCP<Epetra_Vector>            interfaceforcen_;

    /// ALE dof map
    Teuchos::RCP<LINALG::MapExtractor>     meshmap_;

    //! @name local copies of input parameters
    // fluid discretization
    Teuchos::RCP<DRT::Discretization>      embfluiddis_;

    // xfluid discretization
    Teuchos::RCP<DRT::Discretization>      bgfluiddis_;

    /// solver
    Teuchos::RCP<LINALG::Solver>           solver_;

    /// parameters
    Teuchos::RCP<Teuchos::ParameterList>   params_;

    /// output writer
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    //@}

    /// Dirichlet-free, merged velocity DOF
    Teuchos::RCP<Epetra_Map>       innervelmap_;

    /// flag that indicates, whether fluid interface velocities shall
    /// be marked as Dirichlet values
    bool dirichletcond_;

    /// type of monolithic XFluid-Fluid approach (decides whether ALE-mesh is fixed during
    /// Newton iteration)
    enum INPAR::XFEM::Monolithic_xffsi_Approach monolithic_approach_;

    /// flag, that indicates an XFFSI problem (in contrast to pure embedded fluid problem)
    bool monolithicfluidfluidfsi_;

    /*!
     * \brief flag, that indicates a fluidsplit context
     * in communication with an FSI class, there are some occasions
     * (namely: output of Lagrangian multiplier and access to inner
     * embedded fluid velocities) where we need the information about
     * the FSI approach in here
     */
    bool isfluidsplit_;

  };
}

#endif
