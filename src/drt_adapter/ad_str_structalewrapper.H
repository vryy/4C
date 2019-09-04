/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for Structure-ALE problems.

\maintainer Anh-Tu Vuong


\level 3
*/
/*----------------------------------------------------------------------*/


#ifndef STRUCT_ALE_STRUCTURE_WRAPPER_H
#define STRUCT_ALE_STRUCTURE_WRAPPER_H


#include "ad_str_wrapper.H"
#include "../drt_struct_ale/struct_ale_str_model_evaluator.H"


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
   * \author rauch (rauch@lnm.mw.tum.de)          */
  class StructAleWrapper : public StructureWrapper
  {
   public:
    /// constructor
    explicit StructAleWrapper(Teuchos::RCP<Structure> structure);

    /// Destructor
    virtual ~StructAleWrapper() { ; }

    /// return reference to material displacement RCP
    const Teuchos::RCP<Epetra_Vector>& GetMaterialDisplacementNpPtr();

    /// update the material displacements
    void UpdateMaterialDisplacements(Teuchos::RCP<Epetra_Vector> dispmat);

   protected:
    Teuchos::RCP<::STR::MODELEVALUATOR::StructAle> GetStructAleModelEvaluatorPtr()
    {
      return Teuchos::rcp_dynamic_cast<::STR::MODELEVALUATOR::StructAle>(
          struct_ale_model_evaluator_);
    };

   private:
    /// adapter for new structural time integration
    Teuchos::RCP<ADAPTER::StructureNew> structure_;

    /// structural model evaluator object
    Teuchos::RCP<::STR::MODELEVALUATOR::Generic> struct_ale_model_evaluator_;

  };  // class StructAleWrapper
}  // namespace ADAPTER


#endif /* STRUCT_ALE_STRUCTURE_WRAPPER_H */
