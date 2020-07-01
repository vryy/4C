/*----------------------------------------------------------------------------*/
/*! \file

\brief Interface between the xcontact algorithm and the structure algorithm



\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_ADAPTER_AD_STR_XCONTACT_H_
#define SRC_DRT_ADAPTER_AD_STR_XCONTACT_H_


#include "../drt_lib/drt_discret_interface.H"
#include "ad_str_wrapper.H"
#include "ad_str_structure_new.H"

namespace XFEM
{
  class XFieldState;
}  // namespace XFEM
namespace XCONTACT
{
  class MultiDiscretizationWrapper;
  class CutWizard;
  class StateCreator;
}  // namespace XCONTACT
namespace STR
{
  namespace MODELEVALUATOR
  {
    class XContact;
  }  // namespace MODELEVALUATOR
}  // namespace STR
namespace ADAPTER
{
  class StructureNew;

  class StructureXContact : public StructureWrapper
  {
   public:
    /// constructor
    explicit StructureXContact(const Teuchos::RCP<Structure>& structure);

    //! initialize member variables
    void Init(const Teuchos::ParameterList& p_xfem_general, const int& num_dof_per_node);

    //! setup member variables
    virtual void Setup();

    //! prepare the structural time step and detect overlap or adhesion regions
    virtual void PrepareTimeStep();

    bool IsComingIntoContact();

    /** \brief Solve the structural/contact part of the inequality level-set problem
     *
     *  ( 1 ) Solve the direct contact problem
     *  ( 2 ) Solve m-times the sensibility problem
     *  ( 3 ) Solve the shape optimization problem
     *
     *  \author hiermeier \date 06/16 */
    virtual INPAR::STR::ConvergenceStatus Solve();

    /** Set the solution values of the scalar field/level set of the step \f$n+1\f$
     *
     * ToDo Add additional input variables.
     *
     * \date 06/16 */
    void SetScaTraValuesInStructure_Np(const Epetra_Vector& phi_np);

    /// access the weighted gap vector with slave normal dof row map layout (read-only)
    const Epetra_Vector& GetWeightedGap() const;

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

    const Epetra_Comm& Comm() const;

    inline XCONTACT::StateCreator& StateCreator()
    {
      if (state_creator_ptr_.is_null()) dserror("The XFEM::StateCreator was not yet created!");
      return *state_creator_ptr_;
    }

    //! reference to the new underlying structure time integration strategy object
    inline StructureNew& Structure() { return *structure_ptr_; }
    inline const StructureNew& Structure() const { return *structure_ptr_; }

    //! reference to the xstr multi discretization wrapper
    inline XCONTACT::MultiDiscretizationWrapper& MultiDiscret() { return *multi_discret_ptr_; }
    inline const XCONTACT::MultiDiscretizationWrapper& MultiDiscret() const
    {
      return *multi_discret_ptr_;
    }

    //! get the xcontact model (read and write)
    STR::MODELEVALUATOR::XContact& XContactModel();

    //! get the xcontact model (read-only)
    const STR::MODELEVALUATOR::XContact& XContactModel() const;

    //! @name control and maintain the global state objects
    //! @{
    void CompleteInitialState();

    //! recreate the state object
    void RecreateState();
    //! @}

    void PreSolve();

   protected:
    //! initialization indicator
    bool isinit_;

    //! setup indicator
    bool issetup_;

   private:
    //! Pointer to the new structure object
    Teuchos::RCP<StructureNew> structure_ptr_;

    //! Pointer to the structural multi-discretization wrapper
    Teuchos::RCP<XCONTACT::MultiDiscretizationWrapper> multi_discret_ptr_;

    Teuchos::RCP<const Teuchos::ParameterList> p_xfem_general_ptr_;

    //!
    Teuchos::RCP<XCONTACT::StateCreator> state_creator_ptr_;

    Teuchos::RCP<XFEM::XFieldState> state_ptr_;

    int num_dof_per_node_;

    int max_num_dof_sets_;

    int state_count_;

    bool iscontact_;
  };
}  // namespace ADAPTER


#endif /* SRC_DRT_ADAPTER_AD_STR_XCONTACT_H_ */
