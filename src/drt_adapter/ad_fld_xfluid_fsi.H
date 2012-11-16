/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_fsi.H

\brief Fluid field adapter for fsi

Can only be used in conjunction with XFluid!

<pre>
Maintainer:
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
        Teuchos::RCP<Teuchos::ParameterList>            params,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    /// communication object at the interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor>const& Interface() const { return interface_; }

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp() { dserror("Not implemented, yet!");
                                                                  return Teuchos::null; }

    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel);
    virtual void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp);



  protected:
    /// A castet pointer to the fluid itself
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
