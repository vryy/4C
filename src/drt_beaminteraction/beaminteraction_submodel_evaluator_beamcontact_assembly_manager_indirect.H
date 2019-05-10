/*-----------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_beamcontact_assembly_manager_indirect.H

\brief Class to assemble pair based contributions into global matrices. The pairs in this class can
not be directly assembled into the global matrices. They have to be assembled into the global
coupling matrices M and D first.

\maintainer Ivo Steinbrecher

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BEAMINTERACTION_SUBMODEL_EVALUATOR_BEAMCONTACT_ASSEMBLY_INDIRECT_H_
#define BEAMINTERACTION_SUBMODEL_EVALUATOR_BEAMCONTACT_ASSEMBLY_INDIRECT_H_


#include "beaminteraction_submodel_evaluator_beamcontact_assembly_manager.H"


// Forward declaration.
namespace BEAMINTERACTION
{
  class BeamToSolidMortarManager;
  class BeamContactParams;
}  // namespace BEAMINTERACTION


namespace BEAMINTERACTION
{
  namespace SUBMODELEVALUATOR
  {
    /**
     * \brief This class collects local coupling terms of the pairs (D and M) and assembles them
     * into the global coupling matrices. Those global coupling matrices are then multiplied with
     * each other and added to the global force vector and stiffness matrix.
     */
    class BeamContactAssemblyManagerInDirect : public BeamContactAssemblyManager
    {
     public:
      /**
       * \brief Constructor.
       * @param assembly_contact_elepairs (in) Vector with element pairs to be evaluated by this
       * class.
       */
      BeamContactAssemblyManagerInDirect(
          std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> assembly_contact_elepairs,
          Teuchos::RCP<DRT::Discretization> discret,
          Teuchos::RCP<BEAMINTERACTION::BeamContactParams> beam_contact_params_ptr);

      /**
       * \brief Destructor.
       */
      virtual ~BeamContactAssemblyManagerInDirect(){};

      /**
       * \brief Evaluate all force and stiffness terms and add them to the global matrices.
       * @param discret (in) Pointer to the disretization.
       * @param fe_sysvec (out) Global force vector.
       * @param fe_sysmat (out) Global stiffness matrix.
       * @param disp (in) Current displacement vector.
       */
      void EvaluateForceStiff(Teuchos::RCP<DRT::Discretization> discret,
          Teuchos::RCP<Epetra_FEVector> fe_sysvec, Teuchos::RCP<LINALG::SparseMatrix> fe_sysmat,
          Teuchos::RCP<const Epetra_Vector> disp) override;

      /**
       * \brief Return a const pointer to the mortar manager.
       */
      inline Teuchos::RCP<const BEAMINTERACTION::BeamToSolidMortarManager> GetMortarManager() const
      {
        return mortar_manager_;
      }

     private:
      //! Pointer to the mortar manager. This object stores the relevant mortar matrices.
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidMortarManager> mortar_manager_;
    };

  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

#endif
