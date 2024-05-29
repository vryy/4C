/*----------------------------------------------------------------------*/
/*! \file

\brief Service class for XFluid-related output.

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_XFLUID_OUTPUTSERVICE_HPP
#define FOUR_C_FLUID_XFLUID_OUTPUTSERVICE_HPP


/*! header inclusions */
#include "4C_config.hpp"

#include "4C_inpar_cut.hpp"
#include "4C_inpar_xfem.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Discretization;
  class DiscretizationXFEM;
}  // namespace DRT

namespace CORE::Elements
{
  class Element;
}

namespace CORE::Dofsets
{
  class IndependentDofSet;
}

namespace CORE::LINALG
{
  class MapExtractor;
}

namespace CORE::GEO
{
  class CutWizard;

  namespace CUT
  {
    class ElementHandle;
    class VolumeCell;
  }  // namespace CUT
}  // namespace CORE::GEO

namespace XFEM
{
  class ConditionManager;
  class XfemEdgeStab;
}  // namespace XFEM

namespace FLD
{
  class XFluidState;

  /*!
   * \brief Class handles output of XFluid and derived classes
   */
  class XFluidOutputService
  {
   public:
    XFluidOutputService(const Teuchos::RCP<DRT::DiscretizationXFEM>& discret,
        const Teuchos::RCP<XFEM::ConditionManager>& cond_manager);

    virtual ~XFluidOutputService() = default;

    /// prepare standard output
    void prepare_output();

    /// standard output routine
    void Output(int step, double time, bool write_restart_data,
        Teuchos::RCP<const FLD::XFluidState> state,
        Teuchos::RCP<Epetra_Vector> dispnp = Teuchos::null,
        Teuchos::RCP<Epetra_Vector> gridvnp = Teuchos::null);

    /// Gmsh solution output
    virtual void GmshSolutionOutput(const std::string& filename_base,  ///< name for output file
        int step,                                                      ///< step number
        const Teuchos::RCP<FLD::XFluidState>& state,                   ///< state
        int count = -1){};

    /// Gmsh solution output for previous time step
    virtual void gmsh_solution_output_previous(
        const std::string& filename_base,             ///< name for output file
        int step,                                     ///< step number
        const Teuchos::RCP<FLD::XFluidState>& state,  ///< state
        int count = -1){};

    /// Gmsh output of solution (debug)
    virtual void gmsh_solution_output_debug(
        const std::string& filename_base,  ///< name for output file
        int step,                          ///< step number
        int count,                         ///< counter for iterations within a global time step
        const Teuchos::RCP<FLD::XFluidState>& state  ///< state
    ){};

    /// Gmsh output of residual (debug)
    virtual void gmsh_residual_output_debug(
        const std::string& filename_base,  ///< name for output file
        int step,                          ///< step number
        int count,                         ///< counter for iterations within a global time step
        const Teuchos::RCP<FLD::XFluidState>& state  ///< state
    ){};

    /// Gmsh output of increment (debug)
    virtual void gmsh_increment_output_debug(
        const std::string& filename_base,  ///< name for output file
        int step,                          ///< step number
        int count,                         ///< counter for iterations within a global time step
        const Teuchos::RCP<FLD::XFluidState>& state  ///< state
    ){};

    /// Gmsh output of discretization
    virtual void gmsh_output_discretization(bool print_faces, int step,
        std::map<int, CORE::LINALG::Matrix<3, 1>>* curr_pos = nullptr){};

    /// Main output routine for gmsh output
    virtual void GmshOutput(const std::string& filename_base,  ///< name for output file
        const std::string& prefix,                             ///< data prefix
        int step,                                              ///< step number
        int count,  ///< counter for iterations within a global time step
        const Teuchos::RCP<CORE::GEO::CutWizard>& wizard,  ///< cut wizard
        Teuchos::RCP<const Epetra_Vector> vel,  ///< vector holding velocity and pressure dofs
        Teuchos::RCP<const Epetra_Vector> acc = Teuchos::null  ///< vector holding acceleration
    ){};

    /// Gmsh output for EOS
    virtual void GmshOutputEOS(int step,            ///< step number
        Teuchos::RCP<XFEM::XfemEdgeStab> edge_stab  ///< stabilization handler
    ){};

   protected:
    //! @name XFEM discretization
    //@{
    const Teuchos::RCP<DRT::DiscretizationXFEM> discret_;
    //@}

    //! XFEM condition manager
    const Teuchos::RCP<XFEM::ConditionManager> cond_manager_;

    //! dofset for fluid output
    Teuchos::RCP<CORE::Dofsets::IndependentDofSet> dofset_out_;

    //! output vector (mapped to initial fluid dofrowmap)
    Teuchos::RCP<Epetra_Vector> outvec_fluid_;

    //! vel-pres splitter for output purpose
    Teuchos::RCP<CORE::LINALG::MapExtractor> velpressplitter_out_;

    bool firstoutputofrun_;

    //! how many restart steps have already been written
    int restart_count_;
  };

  /*!
   * \brief Class handles output of XFluid and derived classes, capable of handling gmsh output
   */
  class XFluidOutputServiceGmsh : public XFluidOutputService
  {
   public:
    XFluidOutputServiceGmsh(Teuchos::ParameterList& params_xfem,
        const Teuchos::RCP<DRT::DiscretizationXFEM>& discret,
        const Teuchos::RCP<XFEM::ConditionManager>& cond_manager, const bool include_inner);

    /// Gmsh solution output
    void GmshSolutionOutput(const std::string& filename_base,  ///< name for output file
        int step,                                              ///< step number
        const Teuchos::RCP<FLD::XFluidState>& state,           ///< state
        int count = -1) override;

    /// Gmsh solution output for previous time step
    void gmsh_solution_output_previous(const std::string& filename_base,  ///< name for output file
        int step,                                                         ///< step number
        const Teuchos::RCP<FLD::XFluidState>& state,                      ///< state
        int count = -1) override;

    /// Gmsh output of solution (debug)
    void gmsh_solution_output_debug(const std::string& filename_base,  ///< name for output file
        int step,                                                      ///< step number
        int count,  ///< counter for iterations within a global time step
        const Teuchos::RCP<FLD::XFluidState>& state  ///< state
        ) override;

    /// Gmsh output of residual (debug)
    void gmsh_residual_output_debug(const std::string& filename_base,  ///< name for output file
        int step,                                                      ///< step number
        int count,  ///< counter for iterations within a global time step
        const Teuchos::RCP<FLD::XFluidState>& state  ///< state
        ) override;

    /// Gmsh output of increment (debug)
    void gmsh_increment_output_debug(const std::string& filename_base,  ///< name for output file
        int step,                                                       ///< step number
        int count,  ///< counter for iterations within a global time step
        const Teuchos::RCP<FLD::XFluidState>& state  ///< state
        ) override;

    /// Gmsh output of discretization
    void gmsh_output_discretization(bool print_faces, int step,
        std::map<int, CORE::LINALG::Matrix<3, 1>>* curr_pos = nullptr) override;

    /// Main output routine for gmsh output
    void GmshOutput(const std::string& filename_base,  ///< name for output file
        const std::string& prefix,                     ///< data prefix (e.g. "SOL")
        int step,                                      ///< step number
        int count,  ///< counter for iterations within a global time step
        const Teuchos::RCP<CORE::GEO::CutWizard>& wizard,  ///< cut wizard
        Teuchos::RCP<const Epetra_Vector> vel,  ///< vector holding velocity and pressure dofs
        Teuchos::RCP<const Epetra_Vector> acc = Teuchos::null,  ///< vector holding acceleration
        Teuchos::RCP<const Epetra_Vector> dispnp =
            Teuchos::null  ///< vector holding ale displacements
    );

    /// Gmsh output for EOS
    void GmshOutputEOS(int step,                    ///< step number
        Teuchos::RCP<XFEM::XfemEdgeStab> edge_stab  ///< stabilization handler
        ) override;

   private:
    /// Gmsh output function for elements without an CORE::GEO::CUT::ElementHandle
    void gmsh_output_element(DRT::Discretization& discret,  ///< background fluid discretization
        std::ofstream& vel_f,                               ///< output file stream for velocity
        std::ofstream& press_f,                             ///< output file stream for pressure
        std::ofstream& acc_f,                               ///< output file stream for acceleration
        CORE::Elements::Element* actele,                    ///< element
        std::vector<int>& nds,                              ///< vector holding the nodal dofsets
        Teuchos::RCP<const Epetra_Vector> vel,  ///< vector holding velocity and pressure dofs
        Teuchos::RCP<const Epetra_Vector> acc = Teuchos::null,  ///< vector holding acceleration
        Teuchos::RCP<const Epetra_Vector> dispnp =
            Teuchos::null  ///< vector holding ale displacements
    );

    /// Gmsh output function for volumecells
    void gmsh_output_volume_cell(DRT::Discretization& discret,  ///< background fluid discretization
        std::ofstream& vel_f,                                   ///< output file stream for velocity
        std::ofstream& press_f,                                 ///< output file stream for pressure
        std::ofstream& acc_f,                      ///< output file stream for acceleration
        CORE::Elements::Element* actele,           ///< element
        CORE::GEO::CUT::ElementHandle* e,          ///< elementhandle
        CORE::GEO::CUT::VolumeCell* vc,            ///< volumecell
        const std::vector<int>& nds,               ///< vector holding the nodal dofsets
        Teuchos::RCP<const Epetra_Vector> velvec,  ///< vector holding velocity and pressure dofs
        Teuchos::RCP<const Epetra_Vector> accvec = Teuchos::null  ///< vector holding acceleration
    );

    /// Gmsh output function for boundarycells
    void gmsh_output_boundary_cell(
        DRT::Discretization& discret,                     ///< background fluid discretization
        std::ofstream& bound_f,                           ///< output file stream for boundary mesh
        CORE::GEO::CUT::VolumeCell* vc,                   ///< volumecell
        const Teuchos::RCP<CORE::GEO::CutWizard>& wizard  ///< cut wizard
    );

    //! @name flags for detailed gmsh output
    const bool gmsh_sol_out_;           ///< Gmsh solution output for each timestep
    const bool gmsh_ref_sol_out_;       ///< Gmsh reference solution output
    const bool gmsh_debug_out_;         ///< Gmsh debug output (increment, residual, etc.)
    const bool gmsh_debug_out_screen_;  ///< print information about output to screen
    const bool
        gmsh_eos_out_;  ///< output for edge-oriented stabilization and ghost-penalty stabilization
    const bool gmsh_discret_out_;  ///< output of XFEM discretization
    const int gmsh_step_diff_;     ///< no. of kept steps
    //@}

    //! integration approach
    const INPAR::CUT::VCellGaussPts volume_cell_gauss_point_by_;

    //! include elements with inside position?
    const bool include_inner_;
  };

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
