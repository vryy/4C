/*-----------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_beamcontact_assembly_manager_direct.H

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
          std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> assembly_contact_elepairs);

      /**
       * \brief Destructor.
       */
      virtual ~BeamContactAssemblyManagerDirect(){};

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
    };

  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

#endif
