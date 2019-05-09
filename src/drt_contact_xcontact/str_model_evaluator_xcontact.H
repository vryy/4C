/*----------------------------------------------------------------------------*/
/**
\brief Model evaluator for the eXtended contact formulation

\maintainer Matthias Mayr

\level 3
*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_XCONTACT_H_
#define SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_XCONTACT_H_

#include "../drt_structure_new/str_model_evaluator_contact.H"

// forward declaration
namespace XCONTACT
{
  class Strategy;
}  // namespace XCONTACT
namespace STR
{
  namespace MODELEVALUATOR
  {
    class XContact : public Contact
    {
     public:
      //! constructor
      XContact();

      //! destructor
      virtual ~XContact(){};

      virtual void Setup();

      void SetContactStatus(const bool& is_in_contact);

      /// evaluate the weighted gap values at each slave node
      void EvaluateWeightedGap();

      /// access the weighted gap values (slave normal dof row map layout)
      const Epetra_Vector& GetWeightedGap() const;

      void SetLevelSetValuesPtr(const Teuchos::RCP<Epetra_Vector>& levelset_values_ptr);

      const Epetra_Vector& LevelSetValues() const;

     protected:
      virtual void CheckPseudo2D() const;

      //! @name Accessors
      //! @{
      const XCONTACT::Strategy& XStrategy() const;

      XCONTACT::Strategy& XStrategy();
      //! @}

     private:
     private:
      Teuchos::RCP<XCONTACT::Strategy> xstrategy_ptr_;

      Teuchos::RCP<const Epetra_Vector> levelset_values_ptr_;

    };  // class XContact
  }     // namespace MODELEVALUATOR
}  // namespace STR


#endif /* SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_XCONTACT_H_ */
