/*-----------------------------------------------------------*/
/*! \file
\brief Model evaluator for Structure-ALE problems.


\level 3
*/
/*-----------------------------------------------------------*/


#ifndef BACI_STRUCT_ALE_STR_MODEL_EVALUATOR_HPP
#define BACI_STRUCT_ALE_STR_MODEL_EVALUATOR_HPP

#include "baci_config.hpp"

#include "baci_structure_new_model_evaluator_structure.hpp"

BACI_NAMESPACE_OPEN


namespace STR
{
  namespace MODELEVALUATOR
  {
    /*! \brief Specialized model evaluator in case of Struct-Ale.
     *
     *  This model evaluator is constructed in case Struct-Ale is required. Struct-Ale
     *  is needed, e.g. in \ref SSI_Part2WC_PROTRUSIONFORMATION , \ref WEAR::Algorithm ,
     *  and biofilm.
     *
     *  This model evaluator holds a vector containing the material displacements
     *  \ref material_displacements_np_ . We provide access to this vector via
     *  \ref GetMaterialDisplacementNpPtr . In the corresponding structural adapter
     *  \ref StructAleWrapper the link between your algorithm and the this
     *  model evaluator is established.
     *
     *  This model evaluator is constructed instead of the standard model evaluator
     *  \ref STR::MODELEVALUATOR::Structure() in \ref STR::MODELEVALUATOR::Factory::
        BuildStructureModelEvaluator() .
     *
     * \date 12/2016
     */
    class StructAle : public Structure
    {
     public:
      //! constructor
      StructAle();

      //! get pointer to mat. displ. vector at time level n+1 (full structural map).
      const Teuchos::RCP<Epetra_Vector>& GetMaterialDisplacementNpPtr()
      {
        return material_displacements_np_;
      };

      //! setup class variables [derived]
      void Setup() override;

     protected:
      //! @name Functions which are derived from the base Structure class
      //! @{

      //! pre-operator EvalutaInternal
      void PreEvaluateInternal() override;

      //! @}


     private:
      //! material displacements at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> material_displacements_np_;

    };  // class StructAle

  }  // namespace MODELEVALUATOR
}  // namespace STR


BACI_NAMESPACE_CLOSE

#endif  // STRUCT_ALE_STR_MODEL_EVALUATOR_H
