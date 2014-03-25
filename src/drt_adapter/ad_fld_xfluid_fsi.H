/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_fsi.H

\brief Fluid field adapter for fsi

Can only be used in conjunction with XFluid!

<pre>
Maintainer:  Benedikt Schott
             schott@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ADAPTER_FLUID_XFLUID_FSI_H
#define ADAPTER_FLUID_XFLUID_FSI_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include "ad_fld_wrapper.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

/// forward declarations

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
  class XFluid;
  namespace UTILS
  {
    class MapExtractor;
  }
}

namespace ADAPTER {

  class XFluidFSI : public FluidWrapper
  {
  public:

    /// Constructor
    XFluidFSI(Teuchos::RCP<Fluid>              fluid,
        Teuchos::RCP<DRT::Discretization>      xfluiddis,
        Teuchos::RCP<DRT::Discretization>      soliddis,
        Teuchos::RCP<LINALG::Solver>           solver,
        Teuchos::RCP<Teuchos::ParameterList>   params,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    /// communication object at the interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor>const& Interface() const { return interface_; }

    /// communication object at the interface (dummy)
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor>const& FPSIInterface() const
    { dserror("Not implemented in the base class, may be overridden by a subclass.");
      static Teuchos::RCP<FLD::UTILS::MapExtractor> ret = Teuchos::null;
      return ret;
    }

    /// Velocity-displacement conversion at the fsi interface
    virtual double TimeScaling() const;

    /// Return interface forces
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();

    /// Return interface velocity at old time level n
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();

    /// Return interface velocity at new time level n+1
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp() { dserror("Not implemented, yet!");
                                                                  return Teuchos::null; }

    /// apply the interface velocities to the fluid
    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel);

    /// apply the interface displacements to the fluid
    virtual void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp);

    /// convert increment of displacement to increment in velocity
    virtual void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx);



    /// return coupling matrix between fluid and structure as sparse matrices
    Teuchos::RCP<LINALG::SparseMatrix> C_Struct_Fluid_Matrix();
    Teuchos::RCP<LINALG::SparseMatrix> C_Fluid_Struct_Matrix();
    Teuchos::RCP<LINALG::SparseMatrix> C_Struct_Struct_Matrix();

    Teuchos::RCP<Epetra_Vector> RHS_Struct_Vec();

  protected:
    /// A casted pointer to the fluid itself
    Teuchos::RCP<FLD::XFluid>              xfluid_;

    /// the interface map setup for fsi interface, free surface, interior translation
    Teuchos::RCP<FLD::UTILS::MapExtractor> interface_;

    /// ALE dof map
    Teuchos::RCP<LINALG::MapExtractor>     meshmap_;

    //! @name local copies of input parameters
    Teuchos::RCP<DRT::Discretization>      xfluiddis_;
    Teuchos::RCP<DRT::Discretization>      soliddis_;
    Teuchos::RCP<DRT::Discretization>      boundarydis_;
    Teuchos::RCP<LINALG::Solver>           solver_;
    Teuchos::RCP<Teuchos::ParameterList>   params_;
    //@}

  };
}

#endif
