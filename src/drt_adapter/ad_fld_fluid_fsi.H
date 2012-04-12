/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fsi.H

\brief Fluid field adapter for fsi

Can only be used in conjunction with fluidimplicittimeint!

<pre>
Maintainer:
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ADAPTER_FLUID_FSI_H
#define ADAPTER_FLUID_FSI_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include "ad_fld_wrapper.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

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
  class FluidImplicitTimeInt;
  namespace UTILS
  {
    class MapExtractor;
  }
}

namespace ADAPTER {

  class FluidFSI : public FluidWrapper
  {
  public:

    /// Constructor
    FluidFSI(Teuchos::RCP<Fluid> fluid,
        Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output,
        bool isale,
        bool dirichletcond);

    virtual Teuchos::RCP<const Epetra_Map> DofRowMap();

    virtual Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds);

    virtual double TimeScaling() const;

    virtual void Update();

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel);

    /// communication object at the interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor>const& Interface() const { return interface_; }

    virtual Teuchos::RCP<const Epetra_Map> InnerVelocityRowMap();

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceFluidVelocity();

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();

    virtual Teuchos::RCP<Epetra_Vector> ExtractFreeSurfaceVeln();

    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel);

    virtual void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp);

    virtual void ApplyMeshDisplacementIncrement(Teuchos::RCP<const Epetra_Vector> dispstepinc){dserror("not implemented!");};

    virtual void ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel);

    virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm);

    virtual void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx);

    virtual void VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx);

    virtual void FreeSurfDisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx);

    virtual void FreeSurfVelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx);

    virtual Teuchos::RCP<Epetra_Vector>   IntegrateInterfaceShape();

    virtual void UseBlockMatrix(bool splitmatrix);

  protected:
    /// A castet pointer to the fluid itself
    Teuchos::RCP<FLD::FluidImplicitTimeInt> fluidimpl_;

    /// the interface map setup for fsi interface, free surface, interior translation
    Teuchos::RCP<FLD::UTILS::MapExtractor> interface_;
    Teuchos::RCP<Epetra_Vector>            interfaceforcen_;  /// interface force at old time level t_n

    /// ALE dof map
    Teuchos::RCP<LINALG::MapExtractor>     meshmap_;

    //! @name local copies of input parameters
    Teuchos::RCP<DRT::Discretization>      dis_;
    Teuchos::RCP<LINALG::Solver>           solver_;
    Teuchos::RCP<Teuchos::ParameterList>   params_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;
    //@}

    /// all velocity dofs not at the interface
    Teuchos::RCP<Epetra_Map>               innervelmap_;

  };
}

#endif
