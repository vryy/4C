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

#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include <Teuchos_RCP.hpp>

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_cut.H"

namespace DRT
{
  class Discretization;
  class DiscretizationXFEM;
}

namespace GEO
{
  class CutWizard;
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
  class XFEMDofSet;
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
  friend class XFluidResultTest;
  friend class XFluidFluidResultTest;

  friend class XFluidStateCreator;

public:

  /*!
   ctor for one-sided problems
   @param xfluiddofrowmap dof-rowmap of intersected fluid
   */
  explicit XFluidState(Teuchos::RCP<GEO::CutWizard>& wizard, Teuchos::RCP<XFEM::XFEMDofSet>& dofset, Teuchos::RCP<const Epetra_Map> & xfluiddofrowmap);

  /// dtor
  virtual ~XFluidState(){}

  /// setup map extractors for dirichlet maps & velocity/pressure maps
  void SetupMapExtractors(
    const Teuchos::RCP<DRT::Discretization> & xfluiddiscret,
    const double & time
    );

  /// zero all coupling matrices and rhs vectors
  void ZeroCouplingMatricesAndRhs();

  /// access to the cut wizard
  Teuchos::RCP<GEO::CutWizard> Wizard()
  {
    if (wizard_ == Teuchos::null)
      dserror("Cut wizard is uninitialized!");
    return wizard_;
  }


  //! access to the xfem-dofset
  Teuchos::RCP<XFEM::XFEMDofSet> DofSet()
  {
    return dofset_;
  }

protected:

  /// initialize all state members based on the xfluid dof-rowmap
  virtual void InitStateVectors();

  /// initialize the system matrix of the intersected fluid
  void InitSystemMatrix();


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
//  Teuchos::RCP<Epetra_Vector>         dispn_;
//  Teuchos::RCP<Epetra_Vector>         dispnm_;
  //@}

  /// grid velocity (if we have an XFEM-ALE-fluid) (set from the adapter!)
  Teuchos::RCP<Epetra_Vector>         gridvnp_;

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

  //@}


private:

  /// XFEM dofset
  Teuchos::RCP<XFEM::XFEMDofSet> dofset_;

  /// cut wizard
  Teuchos::RCP<GEO::CutWizard> wizard_;

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

  /// create a state-object after a cut
  Teuchos::RCP<XFluidState> Create(
      const Teuchos::RCP<DRT::DiscretizationXFEM>&      xdiscret,           //!< xfluid background discretization
      const Teuchos::RCP<DRT::Discretization> &         boundarydiscret,    //!< boundary discretization, a discretization whose surface elements cut the background mesh
      Teuchos::RCP<const Epetra_Vector>                 back_disp_col,      //!< col vector holding background ALE displacements for backdis
      Teuchos::RCP<const Epetra_Vector>                 cutter_disp_col,    //!< col vector holding interface displacements for cutterdis
      Teuchos::RCP<const Epetra_Vector>                 back_levelset_col,  //!< col vector holding nodal level-set values based on backdis
      std::map<int, LINALG::Matrix<3,1> > &             tip_nodes,          //!< tip nodes for crack application
      Teuchos::ParameterList &                          solver_params,
      const int                                         step,
      const double &                                    time
      );

  /*!
   \brief initialize coupling matrices
   @param slave discretization for matrices to couple with
   */
  void InitCouplingMatrices(const Teuchos::RCP<DRT::Discretization>&  slavediscret);

  /*!
   \brief initialize coupling rhs vectors
   @param slave discretization for rhs-vector to couple with
   */
  void InitCouplingRhs(const Teuchos::RCP<DRT::Discretization>&  slavediscret);

  /*!
   \brief initialize ALE state vectors
   @param dispnp and grivnp vectors w.r.t initial full dofrowmap
   */
  void InitALEStateVectors(
      const Teuchos::RCP<DRT::DiscretizationXFEM>& xdiscret,
      Teuchos::RCP<const Epetra_Vector> dispnp_initmap,
      Teuchos::RCP<const Epetra_Vector> gridvnp_initmap
  );

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

  /// state class
  Teuchos::RCP<FLD::XFluidState> state_;

};

}

#endif
