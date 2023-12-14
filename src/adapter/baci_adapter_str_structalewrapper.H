/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for Structure-ALE problems.



\level 3
*/
/*----------------------------------------------------------------------*/


#ifndef BACI_ADAPTER_STR_STRUCTALEWRAPPER_H
#define BACI_ADAPTER_STR_STRUCTALEWRAPPER_H


#include "baci_adapter_str_wrapper.H"
#include "baci_struct_ale_str_model_evaluator.H"

BACI_NAMESPACE_OPEN

// forward declarations
namespace STR
{
  namespace MODELEVALUATOR
  {
    class StructAle;
    class Generic;
  }  // namespace MODELEVALUATOR
}  // namespace STR


namespace ADAPTER
{
  // forward declaration
  class StructureNew;

  /*! \brief Wrapper for structural time integration in case of Struct-Ale.
   *
   *  This wrapper is constructed in case Struct-Ale is required. Struct-Ale is
   *  needed, e.g. in \ref SSI_Part2WC_PROTRUSIONFORMATION , \ref WEAR::Algorithm ,
   *  and biofilm.
   *
   * \date 12/2016
   */
  class StructAleWrapper : public StructureWrapper
  {
   public:
    /// constructor
    explicit StructAleWrapper(Teuchos::RCP<Structure> structure);

    /// return reference to material displacement RCP
    const Teuchos::RCP<Epetra_Vector>& GetMaterialDisplacementNpPtr();

    /// update the material displacements
    void UpdateMaterialDisplacements(Teuchos::RCP<Epetra_Vector> dispmat);

   protected:
    Teuchos::RCP<STR::MODELEVALUATOR::StructAle> GetStructAleModelEvaluatorPtr()
    {
      return Teuchos::rcp_dynamic_cast<STR::MODELEVALUATOR::StructAle>(struct_ale_model_evaluator_);
    };

   private:
    /// adapter for new structural time integration
    Teuchos::RCP<ADAPTER::StructureNew> structure_;

    /// structural model evaluator object
    Teuchos::RCP<STR::MODELEVALUATOR::Generic> struct_ale_model_evaluator_;

  };  // class StructAleWrapper
}  // namespace ADAPTER


BACI_NAMESPACE_CLOSE

#endif  // ADAPTER_STR_STRUCTALEWRAPPER_H
