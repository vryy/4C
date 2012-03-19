/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_fsi.H

\brief Fluid field adapter for fsi

Can only be used in conjunction with fluidimplicittimeint!

<pre>
Maintainer:
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ADAPTER_FLUID_FLUID_FSI_H
#define ADAPTER_FLUID_FLUID_FSI_H


#include "adapter_fluid_wrapper.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"


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
        Teuchos::RCP<ParameterList> params,
        bool isale,
        bool dirichletcond,
        bool monolithicfluidfluidfsi);

    virtual double TimeScaling() const;

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

    /// communication object at the interface
    virtual const FLD::UTILS::MapExtractor& Interface() const { return *interface_; }

    virtual Teuchos::RCP<const Epetra_Map> InnerVelocityRowMap();

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceFluidVelocity();

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();

    virtual void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp);

    virtual void ApplyMeshDisplacementIncrement(Teuchos::RCP<const Epetra_Vector> dispstepinc){dserror("not implemented!");};

    virtual void ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel);

    virtual void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx);

    virtual void VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx);

    virtual void UseBlockMatrix(bool splitmatrix);

  protected:
    /// type cast pointer to XFluidFluid
    Teuchos::RCP<FLD::XFluidFluid>         xfluidfluid_;

    /// the interface map setup for fsi interface, free surface, interior translation
    Teuchos::RCP<FLD::UTILS::MapExtractor> interface_;

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
    Teuchos::RefCountPtr<Epetra_Map>       innervelmap_;

    enum INPAR::XFEM::Monolithic_xffsi_Approach monolithic_approach_; // xfluidfluid monolithic approach

    bool monolithicfluidfluidfsi_;
  };
}

#endif
