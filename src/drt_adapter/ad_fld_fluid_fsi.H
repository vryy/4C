/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_fsi.H

\brief Fluid field adapter for fsi

Can only be used in conjunction with #FluidImplicitTimeInt

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10362
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

  /*! \brief Fluid field adapter for fsi
   *
   *
   *  Can only be used in conjunction with #FLD::FluidImplicitTimeInt
   */
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

    /// Velocity-displacement conversion at the fsi interface
    virtual double TimeScaling() const;

    /// take current results for converged and save for next time step
    virtual void Update();

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver();

    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel);

    /// communication object at the interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor>const& Interface() const { return interface_; }

    virtual Teuchos::RCP<const Epetra_Map> InnerVelocityRowMap();

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();

    /// Return interface velocity at new time level n+1
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp();

    /// Return interface velocity at old time level n
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();

    virtual Teuchos::RCP<Epetra_Vector> ExtractFreeSurfaceVeln();

    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel);

    virtual void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp);

    virtual void ApplyMeshDisplacementIncrement(Teuchos::RCP<const Epetra_Vector> dispstepinc){dserror("not implemented!");};

    virtual void ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel);

    virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm);

    //! @name Conversion between displacement and velocity at interface

    //! Conversion of displacement to velocity at the interface without predictors or inhomogeneous DBCs
    //!
    //! All input vectors have to live on the fluid field map.
    virtual void DisplacementToVelocity(
        Teuchos::RCP<Epetra_Vector> fcx     ///< interface displacement step increment
    );

    //! Conversion of velocity to displacement at the interface without predictors or inhomogeneous DBCs
    //!
    //! All input vectors have to live on the fluid field map.
    virtual void VelocityToDisplacement(
        Teuchos::RCP<Epetra_Vector> fcx     ///< interface velocity step increment at interface
    );

    virtual void FreeSurfDisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx);

    virtual void FreeSurfVelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx);

    //@}

    virtual Teuchos::RCP<Epetra_Vector>   IntegrateInterfaceShape();

    virtual void UseBlockMatrix(bool splitmatrix);

    virtual void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove);

    //! Project the velocity field into a divergence free subspace
    //!
    //! Project the velocity field into a divergence free subspace
    //! while interface and Dirichlet DOFS are not affected.
    //! The projection is done by the following operation:
    //!
    //! \$f v_{divfree} = \left(I - B(B^TB)^{-1}B^T)\right) v + B(B^TB)^{-1} R\$f
    //!
    //! The vector \$f R \$f ensures that interface and Dirichlet DOFs are not modified.
    //!
    //! \autor: mayr.mt
    //! \date  06/2012
    void ProjVelToDivZero();

    /// reset state vectors
    virtual void Reset(
        bool completeReset = false,
        bool newFiles = false,
        int iter = -1);

    /// calculate error in comparison to analytical solution
    virtual void CalculateError();

  protected:
    /// A casted pointer to the fluid itself
    Teuchos::RCP<FLD::FluidImplicitTimeInt> fluidimpl_;

    //! @name local copies of input parameters
    Teuchos::RCP<DRT::Discretization>      dis_;
    Teuchos::RCP<Teuchos::ParameterList>   params_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;
    //@}

    //! \brief interface map setup for fsi interface, free surface, interior translation
    //!
    //! Note: full map contains velocity AND pressure DOFs
    Teuchos::RCP<FLD::UTILS::MapExtractor> interface_;

    /// interface force at old time level t_n
    Teuchos::RCP<Epetra_Vector>            interfaceforcen_;

    /// ALE dof map
    Teuchos::RCP<LINALG::MapExtractor>     meshmap_;

    /// all velocity dofs not at the interface
    Teuchos::RCP<Epetra_Map>               innervelmap_;

  };
}

#endif
