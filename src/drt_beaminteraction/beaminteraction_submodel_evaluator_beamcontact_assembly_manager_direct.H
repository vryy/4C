/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble pair based contributions into global matrices. The pairs in this class can
be directly assembled into the global matrices.

\maintainer Ivo Steinbrecher

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BEAMINTERACTION_SUBMODEL_EVALUATOR_BEAMCONTACT_ASSEMBLY_DIRECT_H_
#define BEAMINTERACTION_SUBMODEL_EVALUATOR_BEAMCONTACT_ASSEMBLY_DIRECT_H_


#include "beaminteraction_submodel_evaluator_beamcontact_assembly_manager.H"


namespace BEAMINTERACTION
{
  namespace SUBMODELEVALUATOR
  {
    /**
     * \brief This class collects local force and stiffness terms of the pairs and adds them
     * directly into the global force vector and stiffness matrix.
     */
    class BeamContactAssemblyManagerDirect : public BeamContactAssemblyManager
    {
     public:
      /**
       * \brief Constructor.
       * @param assembly_contact_elepairs (in) Vector with element pairs to be evaluated by this
       * class.
       */
      BeamContactAssemblyManagerDirect(
          const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>&
              assembly_contact_elepairs);

      /**
       * \brief Destructor.
       */
      virtual ~BeamContactAssemblyManagerDirect(){};

      /**
       * \brief Evaluate all force and stiffness terms and add them to the global matrices.
       * @param discret (in) Pointer to the disretization.
       * @param data_state (in) Beam interaction data state.
       * @param fe_sysvec (out) Global force vector.
       * @param fe_sysmat (out) Global stiffness matrix.
       */
      void EvaluateForceStiff(Teuchos::RCP<DRT::Discretization> discret,
          const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>& data_state,
          Teuchos::RCP<Epetra_FEVector> fe_sysvec,
          Teuchos::RCP<LINALG::SparseMatrix> fe_sysmat) override;

     protected:
      //! Vector of pairs to be evaluated by this class.
      std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> assembly_contact_elepairs_;
    };
  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

#endif
