/*----------------------------------------------------------------------*/
/*! \file

\brief Creates a state object for (in)stationary XFEM fluid problems

\level 0

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_XFLUID_STATE_CREATOR_HPP
#define FOUR_C_FLUID_XFLUID_STATE_CREATOR_HPP

#include "4C_config.hpp"

#include "4C_inpar_cut.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_inpar_xfem.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Discretization;
}  // namespace DRT

namespace CORE::GEO
{
  class CutWizard;
}

namespace CORE::LINALG
{
  class SparseMatrix;
  class MultiMapExtractor;
  class MapExtractor;
}  // namespace CORE::LINALG

namespace IO
{
  class DiscretizationWriter;
}

namespace XFEM
{
  class ConditionManager;
  class DiscretizationXFEM;
  class XFEMDofSet;
}  // namespace XFEM

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
   * Creates the appropriate wizard & handles the cut state (level-set field or boundary
   * discretization).
   */
  class XFluidStateCreator
  {
   public:
    /// ctor
    XFluidStateCreator(Teuchos::RCP<XFEM::ConditionManager> condition_manager,
        Teuchos::ParameterList& params_xfem, int maxnumdofsets, int minnumdofsets,
        bool include_inner)
        : condition_manager_(condition_manager),
          nodal_dofset_strategy_(CORE::UTILS::IntegralValue<INPAR::CUT::NodalDofSetStrategy>(
              params_xfem, "NODAL_DOFSET_STRATEGY")),
          volume_cell_gauss_point_by_(CORE::UTILS::IntegralValue<INPAR::CUT::VCellGaussPts>(
              params_xfem, "VOLUME_GAUSS_POINTS_BY")),
          bound_cell_gauss_point_by_(CORE::UTILS::IntegralValue<INPAR::CUT::BCellGaussPts>(
              params_xfem, "BOUNDARY_GAUSS_POINTS_BY")),
          gmsh_cut_out_(CORE::UTILS::IntegralValue<int>(params_xfem, "GMSH_CUT_OUT")),
          maxnumdofsets_(maxnumdofsets),
          minnumdofsets_(minnumdofsets),
          include_inner_(include_inner)
    {
    }

    /// create a state-object after a cut (pure XFEM fluid)
    Teuchos::RCP<XFluidState> Create(const Teuchos::RCP<XFEM::DiscretizationXFEM>&
                                         xdiscret,  //!< xfluid background discretization
        Teuchos::RCP<const Epetra_Vector>
            back_disp_col,  //!< col vector holding background ALE displacements for backdis
        Teuchos::ParameterList& solver_params,  //!< solver parameters
        const int step,                         //!< current time step
        const double& time                      //!< current time
    );

    /// create a state-object after a cut (XFEM fluid with embedded fluid mesh)
    Teuchos::RCP<XFluidFluidState> Create(const Teuchos::RCP<XFEM::DiscretizationXFEM>&
                                              xdiscret,  //!< xfluid background discretization
        const Teuchos::RCP<DRT::Discretization>&
            embfluiddiscret,  //!< embedded fluid discretization
        Teuchos::RCP<const Epetra_Vector>
            back_disp_col,  //!< col vector holding background ALE displacements for backdis
        Teuchos::ParameterList& solver_params,  //!< solver parameters
        const int step,                         //!< current time step
        const double& time                      //!< current time
    );


   private:
    /// create wizard, perform cut, create new dofset and update xfem discretization
    void create_new_cut_state(
        Teuchos::RCP<XFEM::XFEMDofSet>& dofset,  //!< xfem dofset obtained from the new wizard
        Teuchos::RCP<CORE::GEO::CutWizard>&
            wizard,  //!< cut wizard associated with current intersection state
        const Teuchos::RCP<XFEM::DiscretizationXFEM>&
            xdiscret,  //!< xfluid background discretization
        Teuchos::RCP<const Epetra_Vector>
            back_disp_col,  //!< col vector holding background ALE displacements for backdis
        Teuchos::ParameterList& solver_params,  //!< solver parameters
        const int step                          //!< current time step
    );


    //! condition manager which handles all coupling objects and the coupling/boundary conditions
    Teuchos::RCP<XFEM::ConditionManager> condition_manager_;

    //! strategy for nodal dofset management
    const INPAR::CUT::NodalDofSetStrategy nodal_dofset_strategy_;

    const INPAR::CUT::VCellGaussPts volume_cell_gauss_point_by_;
    const INPAR::CUT::BCellGaussPts bound_cell_gauss_point_by_;

    /// is gmsh-output active?
    const bool gmsh_cut_out_;

    //! @name size limits for dofsets with variable size
    //@{
    const int maxnumdofsets_;
    int minnumdofsets_;
    //@}

    bool include_inner_;
  };

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
