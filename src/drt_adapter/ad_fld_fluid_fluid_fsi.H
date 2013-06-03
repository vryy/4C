/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_fsi.H

\brief Fluid field adapter for fsi

Can only be used in conjunction with fluidimplicittimeint!

<pre>
Maintainer:
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ADAPTER_FLUID_FLUID_FSI_H
#define ADAPTER_FLUID_FLUID_FSI_H


#include "ad_fld_wrapper.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_inpar/inpar_xfem.H"


#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

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

    /// Constructor
    FluidFluidFSI(Teuchos::RCP<Fluid> fluid,
        Teuchos::RCP<DRT::Discretization> embfluiddis,
        Teuchos::RCP<DRT::Discretization> bgfluiddis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        bool isale,
        bool dirichletcond,
        bool monolithicfluidfluidfsi);

    virtual double TimeScaling() const;

    virtual void Update();

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

    /// communication object at the interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor>const& Interface() const { return interface_; }

    virtual Teuchos::RCP<const Epetra_Map> InnerVelocityRowMap();

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp();

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();

    virtual Teuchos::RCP<FLD::UTILS::FluidXFluidMapExtractor>const& XFluidFluidMapExtractor();

    virtual void ApplyEmbFixedMeshDisplacement(Teuchos::RCP<const Epetra_Vector> disp);

    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel);

    virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm);

    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel);

    virtual void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp);

    virtual void ApplyMeshDisplacementIncrement(Teuchos::RCP<const Epetra_Vector> dispstepinc){dserror("not implemented!");}

    virtual void ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel);

    /// Conversion of displacement to velocity at the interface without predictors or inhomogeneous DBCs
    virtual void DisplacementToVelocity(
      Teuchos::RCP<Epetra_Vector> fcx      ///< fluid displacement step increment
      );

    /// Conversion of velocity to displacement at the interface without predictors or inhomogeneous DBCs
    virtual void VelocityToDisplacement(
      Teuchos::RCP<Epetra_Vector> fcx     ///< fluid velocity step increment at interface
      );

    virtual void UseBlockMatrix(bool splitmatrix);

    virtual void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove);

    virtual Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor();

    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

    //modified to assign fsilambda to the correct (embedded fluid)
    //discretization in case of fluid split
    virtual const Teuchos::RCP<IO::DiscretizationWriter>& DiscWriter();

  protected:
    /// type cast pointer to XFluidFluid
    Teuchos::RCP<FLD::XFluidFluid>         xfluidfluid_;

    /// the interface map setup for fsi interface, free surface, interior translation
    Teuchos::RCP<FLD::UTILS::MapExtractor> interface_;

    /// interface force at old time level t_n
    Teuchos::RCP<Epetra_Vector>            interfaceforcen_;

    /// ALE dof map
    Teuchos::RCP<LINALG::MapExtractor>     meshmap_;

    // FSI-Cond Map of ALE (used for monolithic-fluidfluid-fsi)
    Teuchos::RCP<LINALG::MapExtractor>     fsiale_;

    //! @name local copies of input parameters
    // fluid discretization
    Teuchos::RCP<DRT::Discretization>      embfluiddis_;
    // xfluid discretization
    Teuchos::RCP<DRT::Discretization>      bgfluiddis_;

    Teuchos::RCP<LINALG::Solver>           solver_;
    Teuchos::RCP<Teuchos::ParameterList>   params_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;
    //@}

    /// all velocity dofs not at the interface
    Teuchos::RCP<Epetra_Map>       innervelmap_;

    enum INPAR::XFEM::Monolithic_xffsi_Approach monolithic_approach_; // xfluidfluid monolithic approach

    bool monolithicfluidfluidfsi_;

    bool isembdiscwriter_;

  };
}

#endif
