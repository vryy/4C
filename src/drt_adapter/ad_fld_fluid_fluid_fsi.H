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

#include "ad_fld_fluid_fsi.H"
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
    class FluidXFluidMapExtractor;
  }
}

namespace ADAPTER {

  class FluidFluidFSI : public FluidFSI
  {

  public:

    /// constructor
    FluidFluidFSI(Teuchos::RCP<Fluid> fluid,
        Teuchos::RCP<DRT::Discretization> embfluiddis,
        Teuchos::RCP<DRT::Discretization> bgfluiddis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output,
        bool isale,
        bool dirichletcond);

    /// initialize and prepare maps
    virtual void Init();

    /// save results of current time step, do XFEM cut and refresh the
    /// merged fluid map extractor
    virtual void Update();

    /// @name Accessors
    //@{

    // get merged xfluid-fluid dof row map
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap();

    /// communication object at the interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor>const& Interface() const { return mergedfluidinterface_; }

    /// communication object at the interface
    virtual Teuchos::RCP<FLD::UTILS::FsiMapExtractor>const& FsiInterface() const { return fsiinterface_; }

    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
    { dserror("Do not call RexationSolve for XFFSI."); return Teuchos::null; }

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces()
    { dserror("Do not call ExtractInterfaceForces for XFFSI."); return Teuchos::null; }

    /// get block system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

    /// get the velocity row map of the embedded fluid
    virtual Teuchos::RCP<const Epetra_Map> VelocityRowMap();

    // access to embedded discretization
    const Teuchos::RCP<DRT::Discretization>& Discretization();

    // return discretization writer of embedded fluid discretization (for special purpose output)
    const Teuchos::RCP<IO::DiscretizationWriter>& DiscWriter();

    /// get map extractor for background/embedded fluid
    virtual Teuchos::RCP<FLD::UTILS::FluidXFluidMapExtractor>const& XFluidFluidMapExtractor();

    //@}

    /// set old ALE displacement state before relaxation
    virtual void ApplyEmbFixedMeshDisplacement(Teuchos::RCP<const Epetra_Vector> disp);

    /// apply embedded fluid mesh displacement
    virtual void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp);

    /// apply embedded fluid mesh velocity
    virtual void ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel);

    /// evaluate the fluid and update the merged fluid/FSI DOF-map extractor in case of a change in the DOF-maps
    virtual void Evaluate(
      Teuchos::RCP<const Epetra_Vector> stepinc ///< solution increment between time step n and n+1
      );

    /// request fluid system matrix & shapederivatives as blockmatrices when called with true
    /// (indicated monolithic XFFSI with fluidsplit)
    virtual void UseBlockMatrix(bool split_fluidsysmat = false);

    void UpdateSlaveDOF(Teuchos::RCP<Epetra_Vector>& f){};

  private:

    /// setup of map extractor to distinguish between FSI DOF-map and
    /// merged inner embedded fluid and background fluid DOF-map
    void SetupInterface();

    /// prepare underlying extended shape derivatives matrix, that is based
    /// on the merged fluid dof-map (with background fluid dof set to zero),
    /// as it may change
    void PrepareShapeDerivatives();

    /// type cast pointer to XFluidFluid
    Teuchos::RCP<FLD::XFluidFluid>         xfluidfluid_;

    /// fsi map extractor for embedded fluid
    Teuchos::RCP<FLD::UTILS::FsiMapExtractor> embfluidfsiinterface_;

    /// fsi map extractor for merged fluid maps (to keep fsi interface-DOF apart from
    /// merged inner DOF (inner embedded fluid together with background fluid)
    Teuchos::RCP<FLD::UTILS::MapExtractor> mergedfluidinterface_;

    Teuchos::RCP<FLD::UTILS::FsiMapExtractor> fsiinterface_;

    /// type of monolithic XFluid-Fluid approach (decides whether ALE-mesh is fixed during
    /// Newton iteration)
    enum INPAR::XFEM::Monolithic_xffsi_Approach monolithic_approach_;
  };
}

#endif
