/*----------------------------------------------------------------------------*/
/**
\brief State creator for the xcontact simulation

\maintainer Matthias Mayr

\level 3
*/
/*----------------------------------------------------------------------------*/


#ifndef SRC_DRT_CONTACT_XCONTACT_XCONTACT_STATE_CREATOR_H_
#define SRC_DRT_CONTACT_XCONTACT_XCONTACT_STATE_CREATOR_H_

#include "../drt_xfem/xfem_state_creator.H"

namespace ADAPTER
{
  class Field;
}  // namespace ADAPTER
namespace XSTR
{
  class MultiDiscretizationWrapper;
}  // namespace XSTR

namespace XCONTACT
{
  class StateCreator : public XFEM::StateCreator
  {
   public:
    StateCreator();

    virtual ~StateCreator(){};

    /*! \brief Recreate a state object if necessary
     *
     *  \param xstate        : Old state object, which will be destroyed and recreated,
     *                         if  necessary.
     *  \param xfieldiscret  : xFEM (background) discretization
     *  \param fieldiscret   : FEM (embedded / additional) discretization
     *  \param back_disp_col : col vector holding background ALE displacements for backdis
     *  \param solver_params : solver parameters
     *  \param step          : current time step
     *  \param time          : current time
     *  \param dosetup       : call setup on the new state object
     */
    void Recreate(Teuchos::RCP<XFEM::XFieldState>& xstate,
        const Teuchos::RCP<ADAPTER::Field>& xfield,
        const Teuchos::RCP<DRT::DiscretizationInterface>& full_discret,
        const Teuchos::RCP<const Epetra_Vector>& back_disp_col,
        const Teuchos::RCP<const Epetra_Vector>& levelset_field_row,
        Teuchos::ParameterList& solver_params, int step, double time, bool dosetup);
    void Recreate(Teuchos::RCP<XFEM::XFieldState>& xstate,
        const Teuchos::RCP<ADAPTER::Field>& xfield, XSTR::MultiDiscretizationWrapper& full_discret,
        const Teuchos::RCP<const Epetra_Vector>& back_disp_col,
        const Teuchos::RCP<const Epetra_Vector>& levelset_field_row,
        Teuchos::ParameterList& solver_params, int step, double time, bool dosetup);

   protected:
    void CreateCutWizard(const Teuchos::RCP<DRT::DiscretizationXFEM>& xdiscret,
        const Teuchos::RCP<const Epetra_Vector>& back_disp_col,
        const Teuchos::RCP<const Epetra_Vector>& levelset_field_row,
        Teuchos::RCP<GEO::CutWizard>& wizard);

   private:
    enum XFEM::StateStatus FinishBackgroundDiscretization(
        const Teuchos::RCP<XFEM::XFEMDofSet>& xdofset, Teuchos::ParameterList& solver_params,
        XSTR::MultiDiscretizationWrapper& full_discret);

    enum XFEM::StateStatus CreateNewCutState(Teuchos::RCP<XFEM::XFEMDofSet>& xdofset,
        Teuchos::RCP<GEO::CutWizard>& wizard, XSTR::MultiDiscretizationWrapper& full_discret,
        const Teuchos::RCP<const Epetra_Vector>& back_disp_col,
        const Teuchos::RCP<const Epetra_Vector>& levelset_field_row,
        Teuchos::ParameterList& solver_params, int step);

    Teuchos::RCP<XFEM::XFieldState> CreateNewXFieldState(const Teuchos::RCP<ADAPTER::Field>& xfield,
        const Teuchos::RCP<GEO::CutWizard>& wizard, const Teuchos::RCP<XFEM::XFEMDofSet>& xdofset,
        XSTR::MultiDiscretizationWrapper& full_discret) const;

    void DestroyXFieldState(const Teuchos::RCP<ADAPTER::Field>& xfield) const;

    void ResetXFieldNonStandardDofs(const Teuchos::RCP<ADAPTER::Field>& xfield) const;

  };  // class StateCreator
}  // namespace XCONTACT



#endif /* SRC_DRT_CONTACT_XCONTACT_XCONTACT_STATE_CREATOR_H_ */
