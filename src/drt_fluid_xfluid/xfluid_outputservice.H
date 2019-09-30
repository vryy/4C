/*----------------------------------------------------------------------*/
/*! \file

\brief Service class for XFluid-related output.

\level 2

\maintainer  Martin Kronbichler
             kronbichler@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
*/
/*----------------------------------------------------------------------*/

#ifndef XFLUID_OUTPUTSERVICE_H_
#define XFLUID_OUTPUTSERVICE_H_


/*! header inclusions */
#include "Epetra_Vector.h"

#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_cut.H"

#include "../linalg/linalg_fixedsizematrix.H"

// forward declarations
namespace DRT
{
  class Discretization;
  class DiscretizationXFEM;
  class Element;
  class IndependentDofSet;
}  // namespace DRT

namespace LINALG
{
  class MapExtractor;
}

namespace GEO
{
  class CutWizard;

  namespace CUT
  {
    class ElementHandle;
    class VolumeCell;
  }  // namespace CUT
}  // namespace GEO

namespace XFEM
{
  class ConditionManager;
  class XFEM_EdgeStab;
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

    virtual ~XFluidOutputService(){};

    /// prepare standard output
    void PrepareOutput();

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
    virtual void GmshSolutionOutputPrevious(
        const std::string& filename_base,             ///< name for output file
        int step,                                     ///< step number
        const Teuchos::RCP<FLD::XFluidState>& state,  ///< state
        int count = -1){};

    /// Gmsh output of solution (debug)
    virtual void GmshSolutionOutputDebug(
        const std::string& filename_base,  ///< name for output file
        int step,                          ///< step number
        int count,                         ///< counter for iterations within a global time step
        const Teuchos::RCP<FLD::XFluidState>& state  ///< state
    ){};

    /// Gmsh output of residual (debug)
    virtual void GmshResidualOutputDebug(
        const std::string& filename_base,  ///< name for output file
        int step,                          ///< step number
        int count,                         ///< counter for iterations within a global time step
        const Teuchos::RCP<FLD::XFluidState>& state  ///< state
    ){};

    /// Gmsh output of increment (debug)
    virtual void GmshIncrementOutputDebug(
        const std::string& filename_base,  ///< name for output file
        int step,                          ///< step number
        int count,                         ///< counter for iterations within a global time step
        const Teuchos::RCP<FLD::XFluidState>& state  ///< state
    ){};

    /// Gmsh output of discretization
    virtual void GmshOutputDiscretization(
        bool print_faces, int step, std::map<int, LINALG::Matrix<3, 1>>* curr_pos = NULL){};

    /// Main output routine for gmsh output
    virtual void GmshOutput(const std::string& filename_base,  ///< name for output file
        const std::string& prefix,                             ///< data prefix
        int step,                                              ///< step number
        int count,  ///< counter for iterations within a global time step
        const Teuchos::RCP<GEO::CutWizard>& wizard,  ///< cut wizard
        Teuchos::RCP<const Epetra_Vector> vel,       ///< vector holding velocity and pressure dofs
        Teuchos::RCP<const Epetra_Vector> acc = Teuchos::null  ///< vector holding acceleration
    ){};

    /// Gmsh output for EOS
    virtual void GmshOutputEOS(int step,             ///< step number
        Teuchos::RCP<XFEM::XFEM_EdgeStab> edge_stab  ///< stabilization handler
    ){};

   protected:
    //! @name XFEM discretization
    //@{
    const Teuchos::RCP<DRT::DiscretizationXFEM> discret_;
    //@}

    //! XFEM condition manager
    const Teuchos::RCP<XFEM::ConditionManager> cond_manager_;

    //! dofset for fluid output
    Teuchos::RCP<DRT::IndependentDofSet> dofset_out_;

    //! output vector (mapped to initial fluid dofrowmap)
    Teuchos::RCP<Epetra_Vector> outvec_fluid_;

    //! vel-pres splitter for output purpose
    Teuchos::RCP<LINALG::MapExtractor> velpressplitter_out_;

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
        int count = -1);

    /// Gmsh solution output for previous time step
    void GmshSolutionOutputPrevious(const std::string& filename_base,  ///< name for output file
        int step,                                                      ///< step number
        const Teuchos::RCP<FLD::XFluidState>& state,                   ///< state
        int count = -1);

    /// Gmsh output of solution (debug)
    void GmshSolutionOutputDebug(const std::string& filename_base,  ///< name for output file
        int step,                                                   ///< step number
        int count,  ///< counter for iterations within a global time step
        const Teuchos::RCP<FLD::XFluidState>& state  ///< state
    );

    /// Gmsh output of residual (debug)
    void GmshResidualOutputDebug(const std::string& filename_base,  ///< name for output file
        int step,                                                   ///< step number
        int count,  ///< counter for iterations within a global time step
        const Teuchos::RCP<FLD::XFluidState>& state  ///< state
    );

    /// Gmsh output of increment (debug)
    void GmshIncrementOutputDebug(const std::string& filename_base,  ///< name for output file
        int step,                                                    ///< step number
        int count,  ///< counter for iterations within a global time step
        const Teuchos::RCP<FLD::XFluidState>& state  ///< state
    );

    /// Gmsh output of discretization
    void GmshOutputDiscretization(
        bool print_faces, int step, std::map<int, LINALG::Matrix<3, 1>>* curr_pos = NULL);

    /// Main output routine for gmsh output
    void GmshOutput(const std::string& filename_base,  ///< name for output file
        const std::string& prefix,                     ///< data prefix (e.g. "SOL")
        int step,                                      ///< step number
        int count,  ///< counter for iterations within a global time step
        const Teuchos::RCP<GEO::CutWizard>& wizard,  ///< cut wizard
        Teuchos::RCP<const Epetra_Vector> vel,       ///< vector holding velocity and pressure dofs
        Teuchos::RCP<const Epetra_Vector> acc = Teuchos::null,  ///< vector holding acceleration
        Teuchos::RCP<const Epetra_Vector> dispnp =
            Teuchos::null  ///< vector holding ale displacements
    );

    /// Gmsh output for EOS
    void GmshOutputEOS(int step,                     ///< step number
        Teuchos::RCP<XFEM::XFEM_EdgeStab> edge_stab  ///< stabilization handler
    );

   private:
    /// Gmsh output function for elements without an GEO::CUT::ElementHandle
    void GmshOutputElement(DRT::Discretization& discret,  ///< background fluid discretization
        std::ofstream& vel_f,                             ///< output file stream for velocity
        std::ofstream& press_f,                           ///< output file stream for pressure
        std::ofstream& acc_f,                             ///< output file stream for acceleration
        DRT::Element* actele,                             ///< element
        std::vector<int>& nds,                            ///< vector holding the nodal dofsets
        Teuchos::RCP<const Epetra_Vector> vel,  ///< vector holding velocity and pressure dofs
        Teuchos::RCP<const Epetra_Vector> acc = Teuchos::null,  ///< vector holding acceleration
        Teuchos::RCP<const Epetra_Vector> dispnp =
            Teuchos::null  ///< vector holding ale displacements
    );

    /// Gmsh output function for volumecells
    void GmshOutputVolumeCell(DRT::Discretization& discret,  ///< background fluid discretization
        std::ofstream& vel_f,                                ///< output file stream for velocity
        std::ofstream& press_f,                              ///< output file stream for pressure
        std::ofstream& acc_f,                      ///< output file stream for acceleration
        DRT::Element* actele,                      ///< element
        GEO::CUT::ElementHandle* e,                ///< elementhandle
        GEO::CUT::VolumeCell* vc,                  ///< volumecell
        const std::vector<int>& nds,               ///< vector holding the nodal dofsets
        Teuchos::RCP<const Epetra_Vector> velvec,  ///< vector holding velocity and pressure dofs
        Teuchos::RCP<const Epetra_Vector> accvec = Teuchos::null  ///< vector holding acceleration
    );

    /// Gmsh output function for boundarycells
    void GmshOutputBoundaryCell(DRT::Discretization& discret,  ///< background fluid discretization
        std::ofstream& bound_f,                     ///< output file stream for boundary mesh
        GEO::CUT::VolumeCell* vc,                   ///< volumecell
        const Teuchos::RCP<GEO::CutWizard>& wizard  ///< cut wizard
    );

    //! @name flags for detailed gmsh output
    const bool gmsh_sol_out_;           ///< Gmsh solution output for each timestep
    const bool gmsh_ref_sol_out_;       ///< Gmsh reference solution output
    const bool gmsh_debug_out_;         ///< Gmsh debug output (increment, residual, etc.)
    const bool gmsh_debug_out_screen_;  ///< print information about output to screen
    const bool
        gmsh_EOS_out_;  ///< output for edge-oriented stabilization and ghost-penalty stabilization
    const bool gmsh_discret_out_;  ///< output of XFEM discretization
    const int gmsh_step_diff_;     ///< no. of kept steps
    //@}

    //! integration approach
    const INPAR::CUT::VCellGaussPts VolumeCellGaussPointBy_;

    //! include elements with inside position?
    const bool include_inner_;
  };

}  // namespace FLD

#endif
