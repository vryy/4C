/*!----------------------------------------------------------------------
\file xfluid_state.H
\brief State class for (in)stationary XFEM fluid problems

Attention:
These classes are still prototypes and have to be completed.

<pre>
Maintainer:  Raffaela Kruse
             kruse@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15240
</pre>

*----------------------------------------------------------------------*/

#ifndef XFLUID_STATE_H_
#define XFLUID_STATE_H_

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_cut.H"

namespace DRT
{
  class Discretization;
}
namespace LINALG
{
  class SparseMatrix;
  class MultiMapExtractor;
  class MapExtractor;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace XFEM
{
  class FluidWizard;
  class FluidWizardMesh;
  class FluidWizardLevelSet;
  class FluidDofSet;
}

namespace FLD
{

namespace UTILS
{
  class KSPMapExtractor;
}
/**
 * Container class for the state vectors and maps of the intersected background
 * fluid - tied to a specific intersection state (interface position,
 * independent from the fact, in which form the interface is given (boundary mesh or level-set field)).
 */
class XFluidState
{
  friend class XFluid;
  friend class XFluidFluid;
  friend class XFluidResultTest2;

public:

  /*!
   ctor for one-sided problems
   @param xfluiddofrowmap dof-rowmap of intersected fluid
   */
  explicit XFluidState(Teuchos::RCP<const Epetra_Map> & xfluiddofrowmap);

  /*!
   ctor for two-sided problems
   @param xfluiddofrowmap dof-rowmap of intersected fluid
   @param slavedofrowmap dof-row map of discretization to couple with
   @param slavedofcolmap dof-column map of discretization to couple with
   */
  explicit XFluidState(Teuchos::RCP<const Epetra_Map> & xfluiddofrowmap,
                       const Epetra_Map &               slavedofrowmap,
                       const Epetra_Map &               slavedofcolmap);

  /// dtor
  virtual ~XFluidState(){}

  /// setup map extractors for dirichlet maps & velocity/pressure maps
  void SetupMapExtractors(
    const Teuchos::RCP<DRT::Discretization> & xfluiddiscret,
    const double & time
    );

  /// setup map extractors for krylov space projection
  void SetupKSPMapExtractor(
    const Teuchos::RCP<DRT::Discretization> & xfluiddiscret
    );

protected:

  /// initialize all state members based on the xfluid dof-rowmap
  virtual void InitStateVectors();

  /// initialize the system matrix of the intersected fluid
  void InitSystemMatrix();

  /*!
   \brief initialize coupling matrices & rhs
   @param slavedofrowmap dof-row map of discretization to couple with
   @param slavedofcolmap dof-column map of discretization to couple with
   */
  void InitCouplingMatricesAndRhs(const Epetra_Map & slavedofrowmap,
                                  const Epetra_Map & slavedofcolmap);

  /// dof-rowmap of intersected fluid
  Teuchos::RCP<const Epetra_Map>      xfluiddofrowmap_;

  /// system matrix (internally EpetraFECrs)
  Teuchos::RCP<LINALG::SparseMatrix>  sysmat_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  Teuchos::RCP<Epetra_Vector>         zeros_;

  /// the vector containing body and surface forces
  Teuchos::RCP<Epetra_Vector>         neumann_loads_;

  /// (standard) residual vector (rhs for the incremental form),
  Teuchos::RCP<Epetra_Vector>         residual_;

  /// true (rescaled) residual vector without zeros at dirichlet positions
  Teuchos::RCP<Epetra_Vector>         trueresidual_;

  /// rhs: right hand side vector for linearized flow solution
  Teuchos::RCP<Epetra_Vector>         rhs_;

  /// nonlinear iteration increment vector
  Teuchos::RCP<Epetra_Vector>         incvel_;

  //! @name acceleration/(scalar time derivative) at time n+1, n and n+alpha_M/(n+alpha_M/n)
  //@{
  Teuchos::RCP<Epetra_Vector>         accnp_;
  Teuchos::RCP<Epetra_Vector>         accn_;
  Teuchos::RCP<Epetra_Vector>         accam_;
  //@}

  //! @name velocity and pressure at time n+1, n, n-1 and n+alpha_F
  //@{
  Teuchos::RCP<Epetra_Vector>         velnp_;
  Teuchos::RCP<Epetra_Vector>         veln_;
  Teuchos::RCP<Epetra_Vector>         velnm_;
  Teuchos::RCP<Epetra_Vector>         velaf_;
  //@}

  //! @name scalar at time n+alpha_F/n+1 and n+alpha_M/n
  //@{
  Teuchos::RCP<Epetra_Vector>         scaaf_;
  Teuchos::RCP<Epetra_Vector>         scaam_;
  //@}

  //! @name displacemets at time n+1, n and n-1 (if we have an XFEM-ALE-fluid)
  //@{
  Teuchos::RCP<Epetra_Vector>         dispnp_;
  Teuchos::RCP<Epetra_Vector>         dispn_;
  Teuchos::RCP<Epetra_Vector>         dispnm_;
  //@}

  /// grid velocity (if we have an XFEM-ALE-fluid)
  Teuchos::RCP<Epetra_Vector>         gridv_;

  /// histvector --- a linear combination of velnm, veln (BDF)
  ///                or veln, accn (One-Step-Theta)
  Teuchos::RCP<Epetra_Vector>         hist_;

  //! @name coupling matrices x: xfluid, s: coupling slave (structure, ALE-fluid, xfluid-element with other active dofset, etc.)
  //@{
  Teuchos::RCP<LINALG::SparseMatrix>  C_sx_;             ///< slave - xfluid coupling block
  Teuchos::RCP<LINALG::SparseMatrix>  C_xs_;             ///< xfluid - slave coupling block
  Teuchos::RCP<LINALG::SparseMatrix>  C_ss_;             ///< slave - slave coupling block

  Teuchos::RCP<Epetra_Vector>         rhC_s_;            ///< slave rhs block
  Teuchos::RCP<Epetra_Vector>         rhC_s_col_;        ///< slave rhs block based on col map
  //@}

  //! @name map extractors
  //@{
  /// maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<LINALG::MapExtractor>        dbcmaps_;

  /// velocity/pressure map extractor, used for convergence check
  Teuchos::RCP<LINALG::MapExtractor>        velpressplitter_;

  /// Krylov space projection map extractor
  Teuchos::RCP<FLD::UTILS::KSPMapExtractor> kspsplitter_;
  //@}
};

/**
 * Builder class for XFluidState.
 * Creates the appropriate wizard & handles the cut state (level-set field or boundary discretization).
 */
class XFluidStateCreator
{

public:

  /// ctor
  XFluidStateCreator(
    INPAR::CUT::VCellGaussPts VolumeCellGaussPointBy,
    INPAR::CUT::BCellGaussPts BoundCellGaussPointBy,
    bool gmsh_cut_out,
    int maxnumdofsets,
    int minnumdofsets ) :
      VolumeCellGaussPointBy_(VolumeCellGaussPointBy),
      BoundCellGaussPointBy_(BoundCellGaussPointBy),
      gmsh_cut_out_(gmsh_cut_out),
      maxnumdofsets_(maxnumdofsets),
      minnumdofsets_(minnumdofsets)
  {
  }

  /// create a xfem fluid wizard (mesh)
  void SetWizard(
      const Teuchos::RCP<DRT::Discretization> & discret,
      const Teuchos::RCP<DRT::Discretization> & boundarydiscret);

  /// create a xfem fluid wizard (level-set)
  void SetWizard(
      const Teuchos::RCP<DRT::Discretization> & discret);

  /// create a state-object after a cut
  Teuchos::RCP<XFluidState> Create(
    const Teuchos::RCP<DRT::Discretization> & discret,
    const Teuchos::RCP<DRT::Discretization> & boundarydiscret,
    const Epetra_Vector &                     idispcol,
    Teuchos::ParameterList &                  solver_params,
    const int                                 step,
    const double &                            time,
    bool                                      coupling = true);

  /// create a state-object after a cut
  Teuchos::RCP<XFluidState> Create(
    const Teuchos::RCP<DRT::Discretization> & discret
    )
  {
    dserror("Not implemented yet!");
    return Teuchos::null;
  }

  /// access to the fluid wizard (boundary mesh)
  Teuchos::RCP<XFEM::FluidWizard> Wizard()
  {
    if (wizard_ == Teuchos::null)
      dserror("Xfem fluid wizard is uninitialized!");
    return wizard_;
  }

  //! access to the xfem-dofset
  Teuchos::RCP<XFEM::FluidDofSet> DofSet()
  {
    return dofset_;
  }

private:

  const INPAR::CUT::VCellGaussPts VolumeCellGaussPointBy_;
  const INPAR::CUT::BCellGaussPts BoundCellGaussPointBy_;

  /// is gmsh-output active?
  const bool gmsh_cut_out_;

  //! @name size limits for dofsets with variable size
  //@{
  const int maxnumdofsets_;
  int minnumdofsets_;
  //@}

  /// xfem-fluid dofset
  Teuchos::RCP<XFEM::FluidDofSet>  dofset_;

  /// xfem fluid wizard (boundary mesh)
  Teuchos::RCP<XFEM::FluidWizard>  wizard_;

};

//class XFluidMeshState : public XFluidState
//{
//
//public:
//
//  /*!
//   ctor
//   @param xfluiddofrowmap dof-rowmap of intersected fluid
//   */
//  explicit XFluidMeshState( Teuchos::RCP<const Epetra_Map> &  xfluiddofrowmap ) :
//    XFluidState(xfluiddofrowmap)
//  {
//  }
//
//  /*!
//   ctor
//   @param xfluiddofrowmap dof-rowmap of intersected fluid
//   @param slavedofrowmap  dof-rowmap from elements to couple with
//   @param slavedofrowmap  dof-colmap from elements to couple with
//   */
//  explicit XFluidMeshState( Teuchos::RCP<const Epetra_Map> &  xfluiddofrowmap,
//                            const Epetra_Map &                slavedofrowmap,
//                            const Epetra_Map &                slavedofcolmap);
//
//private:
//
//  /// initialize all state members, that are tied to the fact, that the interface is described by a mesh
//  void InitStateVectors();
//
//};
//
//class XFluidLevelSetState : public XFluidState
//{
//
//public:
//
//  /*!
//   \brief Create a container class based on the background fluid dof-rowmap
//   @param xfluiddofrowmap dof-rowmap of intersected fluid
//   */
//  explicit XFluidLevelSetState( Teuchos::RCP<const Epetra_Map> &  xfluiddofrowmap ) :
//    XFluidState(xfluiddofrowmap)
//  {
//  }
//
//private:
//
//  /// initialize all state members, that are tied to the fact, that the interface is described by a level
//  void InitStateVectors();
//
//};
}

#endif
