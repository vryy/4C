/*----------------------------------------------------------------------------*/
/**
\file xfem_state_creator.H

\brief Merge the mutual functionality of the XCONTACT state creators.

\maintainer Matthias Mayr

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_XFEM_XFEM_STATE_CREATOR_H_
#define SRC_DRT_XFEM_XFEM_STATE_CREATOR_H_


#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_cut.H"
#include <Teuchos_RCP.hpp>

// forward declarations
class Epetra_Vector;
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos
namespace DRT
{
  class Discretization;
  class DiscretizationInterface;
  class DiscretizationXFEM;
}  // namespace DRT
namespace GEO
{
  class CutWizard;
}  // namespace GEO
namespace XFEM
{
  class ConditionManager;
  class XFEMDofSet;
  class XFieldState;

  enum StateStatus
  {
    state_unchanged = 0,
    state_changed = 1,
    state_undefined = -1
  };

  static inline std::string StateStatus2String(enum StateStatus status)
  {
    switch (status)
    {
      case state_unchanged:
        return "state_unchanged";
      case state_changed:
        return "state_changed";
      case state_undefined:
        return "state_undefined";
      default:
        return "Unknown XFEM::StateStatus";
    }
    exit(EXIT_FAILURE);
  }

  class StateCreator
  {
   public:
    //! constructor
    StateCreator();

    //! destructor
    virtual ~StateCreator(){};

    //! initialize member variables
    void Init(const Teuchos::RCP<XFEM::ConditionManager>& condition_manager,
        const Teuchos::ParameterList& params_xfem, int numdof, int maxnumdofsets, int minnumdofsets,
        bool include_inner);

    //! setup member variables
    virtual void Setup();

    /*! \brief create a state-object after a cut (pure XFEM discretization)
     *
     *  \param xfieldiscret  : xFEM background discretization
     *  \param back_disp_col : col vector holding background ALE displacements for backdis
     *  \param solver_params : solver parameters
     *  \param step          : current time step
     *  \param time          : current time
     */
    Teuchos::RCP<XFieldState> Create(const Teuchos::RCP<DRT::DiscretizationXFEM>& xfielddiscret,
        Teuchos::RCP<const Epetra_Vector> back_disp_col, Teuchos::ParameterList& solver_params,
        int step, double time, bool dosetup);

    /*! create a state-object after a cut (XFEM with embedded standard mesh)
     *
     *  \param xfieldiscret  : xFEM (background) discretization
     *  \param fieldiscret   : FEM (embedded / additional) discretization
     *  \param back_disp_col : col vector holding background ALE displacements for backdis
     *  \param solver_params : solver parameters
     *  \param step          : current time step
     *  \param time          : current time
     *  \param dosetup       : call setup on the new state object
     */
    Teuchos::RCP<XFEM::XFieldState> Create(
        const Teuchos::RCP<DRT::DiscretizationInterface>& xfielddiscret,
        const Teuchos::RCP<DRT::DiscretizationInterface>& fielddiscret,
        Teuchos::RCP<const Epetra_Vector> back_disp_col, Teuchos::ParameterList& solver_params,
        int step, double time, bool dosetup);
    Teuchos::RCP<XFieldState> Create(const Teuchos::RCP<DRT::DiscretizationXFEM>& xfielddiscret,
        const Teuchos::RCP<DRT::Discretization>& fielddiscret,
        Teuchos::RCP<const Epetra_Vector> back_disp_col, Teuchos::ParameterList& solver_params,
        int step, double time, bool dosetup);

    /*! create a new cut wizard and complete the given state-object
     *
     *  \param state         : state object which will be completed
     *  \param xfieldiscret  : xFEM (background) discretization
     *  \param fieldiscret   : FEM (embedded/additional) discretization
     *  \param back_disp_col : col vector holding background ALE displacements for backdis
     *  \param solver_params : solver parameters
     *  \param step          : current time step
     *  \param time          : current time
     *  \param dosetup       : call setup on the new state object
     */
    void CompleteState(const Teuchos::RCP<XFEM::XFieldState>& state,
        const Teuchos::RCP<DRT::DiscretizationInterface>& xfielddiscret,
        const Teuchos::RCP<DRT::DiscretizationInterface>& fielddiscret,
        Teuchos::RCP<const Epetra_Vector> back_disp_col, Teuchos::ParameterList& solver_params,
        int step, double time, bool dosetup);
    void CompleteState(const Teuchos::RCP<XFEM::XFieldState>& state,
        const Teuchos::RCP<DRT::DiscretizationXFEM>& xfielddiscret,
        const Teuchos::RCP<DRT::Discretization>& fielddiscret,
        Teuchos::RCP<const Epetra_Vector> back_disp_col, Teuchos::ParameterList& solver_params,
        int step, double time, bool dosetup);

   protected:
    //! check the initialization indicator
    inline void CheckInit() const
    {
      if (not isinit_) dserror("Call Init() first!");
    }

    //! check the initialization and setup indicators
    inline void CheckInitSetup() const
    {
      if (not issetup_ or not isinit_) dserror("Call Init() and Setup() first!");
    }

    void CreateCutWizard(const Teuchos::RCP<DRT::DiscretizationXFEM>& xdiscret,
        const Teuchos::RCP<const Epetra_Vector>& back_disp_col,
        Teuchos::RCP<GEO::CutWizard>& wizard);


    //! factory routine to create the right xfield state corresponding to the problem type
    Teuchos::RCP<XFEM::XFieldState> CreateXFieldState() const;

    //! factory routine to create the right xfield/field state corresponding to the problem type
    Teuchos::RCP<XFEM::XFieldState> CreateXFieldFieldState() const;

   private:
    enum StateStatus FinishBackgroundDiscretization(const Teuchos::RCP<XFEM::XFEMDofSet>& xdofset,
        Teuchos::ParameterList& solver_params, DRT::DiscretizationXFEM& xdiscret);

    /*! \brief create wizard, perform cut, create new dofset and update xfem discretization
     *
     *  \param xdofset        : xfem dofset obtained from the new wizard
     *  \param wizard        : cut wizard associated with current intersection state
     *  \param xdiscret      : xflield (background) discretization
     *  \param back_disp_col : col vector holding background ALE displacements for backdis
     *  \param solver_params : solver parameters
     *  \param step          : current time step
     *
     *  \return \c state_changed if the discretization dof set has been modified */
    enum StateStatus CreateNewCutState(Teuchos::RCP<XFEM::XFEMDofSet>& xdofset,
        Teuchos::RCP<GEO::CutWizard>& wizard, const Teuchos::RCP<DRT::DiscretizationXFEM>& xdiscret,
        const Teuchos::RCP<const Epetra_Vector>& back_disp_col,
        Teuchos::ParameterList& solver_params, int step);

    /*! \brief add cutter states to the cut wizard if a condition manager was initialized
     *
     *  \param wizard: cut wizard associated with current intersection state
     */
    void AddCutterStates(Teuchos::RCP<GEO::CutWizard>& wizard);

   protected:
    //! init indicator
    bool isinit_;

    //! setup indicator
    bool issetup_;

    //! condition manager which handles all coupling objects and the coupling/boundary conditions
    Teuchos::RCP<XFEM::ConditionManager> condition_manager_;

    //! strategy for nodal dofset management
    INPAR::CUT::NodalDofSetStrategy nodal_dofset_strategy_;

    INPAR::CUT::VCellGaussPts volume_cell_gauss_point_by_;
    INPAR::CUT::BCellGaussPts bound_cell_gauss_point_by_;

    /// is gmsh-output active?
    bool gmsh_cut_out_;

    /// number of dofs of the underlying standard field discretization
    int numdof_;

    //! @name size limits for dofsets with variable size
    //@{
    int maxnumdofsets_;
    int minnumdofsets_;
    //@}

    bool include_inner_;

  };  // class StateCreator
}  // namespace XFEM


#endif /* SRC_DRT_XFEM_XFEM_STATE_CREATOR_H_ */
