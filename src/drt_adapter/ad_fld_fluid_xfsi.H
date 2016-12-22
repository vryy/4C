/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_xfsi.H

\brief Fluid field adapter for fsi

Can only be used in conjunction with XFluid!

<pre>
\level 1
\maintainer  Benedikt Schott
             schott@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef AD_FLD_FLUID_XFSI_H
#define AD_FLD_FLUID_XFSI_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include "ad_fld_wrapper.H"
#include "../drt_fluid_xfluid/xfluid.H"
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

namespace XFEM
{
  class MeshCouplingFSI;
}

namespace ADAPTER {

  class XFluidFSI : public FluidWrapper
  {
  public:

    /// Constructor
    XFluidFSI(Teuchos::RCP<Fluid>              fluid,
        const std::string                      coupling_name, // name of the FSI coupling condition
        Teuchos::RCP<LINALG::Solver>           solver,
        Teuchos::RCP<Teuchos::ParameterList>   params,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    /// initialize algorithm
    virtual void Init();

    /// communication object at the interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor>const& StructInterface() const { return structinterface_; }

    /// communication object at the interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor>const& Interface() const { return interface_; }

    /// communication object at the interface without pressure dofs for FPSI problems
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor>const& FPSIInterface() const { return fpsiinterface_; }

    /// Velocity-displacement conversion at the fsi interface
    virtual double TimeScaling() const;

    /// Return interface forces
    virtual Teuchos::RCP<Epetra_Vector> ExtractStructInterfaceForces();

    /// Return interface velocity at old time level n
    virtual Teuchos::RCP<Epetra_Vector> ExtractStructInterfaceVeln();

    /// Return interface velocity at new time level n+1
    virtual Teuchos::RCP<Epetra_Vector> ExtractStructInterfaceVelnp() { dserror("Not implemented, yet!");
                                                                  return Teuchos::null; }

    /// apply the interface velocities to the fluid
    virtual void ApplyStructInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel);

    /// apply the interface displacements to the fluid
    virtual void ApplyStructMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp);

    /// convert increment of displacement to increment in velocity
    virtual void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx);

    /// apply the interface displacements to the fluid
    virtual void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp);

    virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm);

    /// return coupling matrix between fluid and structure as sparse matrices
    Teuchos::RCP<LINALG::SparseMatrix> C_Struct_Fluid_Matrix();
    Teuchos::RCP<LINALG::SparseMatrix> C_Fluid_Struct_Matrix();
    Teuchos::RCP<LINALG::SparseMatrix> C_Struct_Struct_Matrix();

    Teuchos::RCP<const Epetra_Vector> RHS_Struct_Vec();

    Teuchos::RCP<FLD::XFluid> MyFluid() {return xfluid_; }

    /// return boundary discretization
    const Teuchos::RCP<DRT::Discretization> BoundaryDiscretization();

    bool NewtonRestartMonolithic(){ return xfluid_->NewtonRestartMonolithic(); }

    Teuchos::RCP<std::map<int, int> > GetPermutationMap(){ return xfluid_->GetPermutationMap(); }

    /// GmshOutput for background mesh and cut mesh
    void GmshOutput(
        const std::string & name,                       ///< name for output file
        const int step,                                 ///< step number
        const int count,                                ///< counter for iterations within a global time step
        Teuchos::RCP<Epetra_Vector> vel,                ///< vector holding velocity and pressure dofs
        Teuchos::RCP<Epetra_Vector> acc = Teuchos::null ///< vector holding accelerations
    );

  protected:
    /// A casted pointer to the fluid itself
    Teuchos::RCP<FLD::XFluid>              xfluid_;

    /// the interface map setup for fsi interface, free surface, interior translation
    Teuchos::RCP<FLD::UTILS::MapExtractor> interface_;

    /// the interface map setup for fsi interface, free surface, interior translation
    Teuchos::RCP<FLD::UTILS::MapExtractor> structinterface_;

    /// the interface map setup for fpsi interface
    Teuchos::RCP<FLD::UTILS::MapExtractor> fpsiinterface_;

    /// ALE dof map
    Teuchos::RCP<LINALG::MapExtractor>     meshmap_;
    Teuchos::RCP<Epetra_Map> permfluidmap_;
    Teuchos::RCP<Epetra_Map> fullfluidmap_;

    //! @name local copies of input parameters
    std::string coupling_name_;              /// the name of the XFEM::MeshCoupling object
    Teuchos::RCP<XFEM::MeshCouplingFSI>    mesh_coupling_fsi_;
    Teuchos::RCP<LINALG::Solver>           solver_;
    Teuchos::RCP<Teuchos::ParameterList>   params_;
    //@}

  };
}

#endif
