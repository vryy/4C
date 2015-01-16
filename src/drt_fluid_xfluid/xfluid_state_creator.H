/*!----------------------------------------------------------------------
\file xfluid_state_creator.H
\brief Creates a state object for (in)stationary XFEM fluid problems

<pre>
Maintainer:  Raffaela Kruse and Benedikt Schott
             [kruse,schott]@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15240
</pre>

*----------------------------------------------------------------------*/

#ifndef XFLUID_STATE_CREATOR_H_
#define XFLUID_STATE_CREATOR_H_

#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include <Teuchos_RCP.hpp>

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_cut.H"

#include "../linalg/linalg_utils.H"

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

class XFluidState;
class XFluidFluidState;

namespace UTILS
{
  class KSPMapExtractor;
}

/**
 * Builder class for XFluid(Fluid)State.
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
    int minnumdofsets,
    bool include_inner) :
      VolumeCellGaussPointBy_(VolumeCellGaussPointBy),
      BoundCellGaussPointBy_(BoundCellGaussPointBy),
      gmsh_cut_out_(gmsh_cut_out),
      maxnumdofsets_(maxnumdofsets),
      minnumdofsets_(minnumdofsets),
      include_inner_(include_inner)
  {
  }

  /// create a state-object after a cut (pure XFEM fluid)
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

  /// create a state-object after a cut (XFEM fluid with embedded fluid mesh)
    Teuchos::RCP<XFluidFluidState> Create(
      const Teuchos::RCP<DRT::DiscretizationXFEM>&      xdiscret,           //!< xfluid background discretization
      const Teuchos::RCP<DRT::Discretization> &         boundarydiscret,    //!< boundary discretization, a discretization whose surface elements cut the background mesh
      const Teuchos::RCP<DRT::Discretization> &         embfluiddiscret,    //!< embedded fluid discretization
      Teuchos::RCP<const Epetra_Vector>                 back_disp_col,      //!< col vector holding background ALE displacements for backdis
      Teuchos::RCP<const Epetra_Vector>                 cutter_disp_col,    //!< col vector holding interface displacements for cutterdis
      Teuchos::RCP<const Epetra_Vector>                 back_levelset_col,  //!< col vector holding nodal level-set values based on backdis
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

  /// create wizard, perform cut, create new dofset and update xfem discretization
  void CreateNewCutState(
      Teuchos::RCP<XFEM::XFEMDofSet> &                  dofset,             //!< xfem dofset obtained from the new wizard
      Teuchos::RCP<GEO::CutWizard> &                    wizard,             //!< cut wizard associated with current intersection state
      const Teuchos::RCP<DRT::DiscretizationXFEM>&      xdiscret,           //!< xfluid background discretization
      const Teuchos::RCP<DRT::Discretization> &         boundarydiscret,    //!< boundary discretization, a discretization whose surface elements cut the background mesh
      Teuchos::RCP<const Epetra_Vector>                 back_disp_col,      //!< col vector holding background ALE displacements for backdis
      Teuchos::RCP<const Epetra_Vector>                 cutter_disp_col,    //!< col vector holding interface displacements for cutterdis
      Teuchos::RCP<const Epetra_Vector>                 back_levelset_col,  //!< col vector holding nodal level-set values based on backdis
      Teuchos::ParameterList &                          solver_params,
      const int                                         step,
      std::map<int, LINALG::Matrix<3,1> >               tip_nodes = std::map<int, LINALG::Matrix<3,1> >()  //!< tip nodes for crack application
  );


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

  bool include_inner_;

};

}

#endif
