/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble pair based contributions into global matrices.

\maintainer Ivo Steinbrecher

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BEAMINTERACTION_SUBMODEL_EVALUATOR_BEAMCONTACT_ASSEMBLY_H_
#define BEAMINTERACTION_SUBMODEL_EVALUATOR_BEAMCONTACT_ASSEMBLY_H_


#include "../drt_lib/drt_dserror.H"
#include "Teuchos_RCP.hpp"


// Forward declarations.
class Epetra_FEVector;
class Epetra_Vector;
namespace LINALG
{
  class SparseMatrix;
}
namespace DRT
{
  class Discretization;
}
namespace BEAMINTERACTION
{
  class BeamContactPair;
}
namespace STR
{
  namespace MODELEVALUATOR
  {
    class BeamInteractionDataState;
  }
}  // namespace STR


namespace BEAMINTERACTION
{
  namespace SUBMODELEVALUATOR
  {
    /**
     * \brief This class assembles the contribution of beam contact pairs into the global force
     * vector and stiffness matrix. The method EvaluateForceStiff has to be overloaded in the
     * derived classes to implement the correct assembly method.
     */
    class BeamContactAssemblyManager
    {
     public:
      /**
       * \brief Constructor.
       * @param assembly_contact_elepairs (in) Vector with element pairs to be evaluated by this
       * class.
       */
      BeamContactAssemblyManager(
          std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> assembly_contact_elepairs);

      /**
       * \brief Destructor.
       */
      virtual ~BeamContactAssemblyManager(){};

      /**
       * \brief Evaluate all force and stiffness terms and add them to the global matrices.
       * @param discret (in) Pointer to the disretization.
       * @param data_state (in) Beam interaction data state.
       * @param fe_sysvec (out) Global force vector.
       * @param fe_sysmat (out) Global stiffness matrix.
       */
      virtual void EvaluateForceStiff(Teuchos::RCP<DRT::Discretization> discret,
          const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>& data_state,
          Teuchos::RCP<Epetra_FEVector> fe_sysvec, Teuchos::RCP<LINALG::SparseMatrix> fe_sysmat)
      {
        dserror("Not implemented!");
      }

     protected:
      //! Vector of pairs to be evaluated by this class.
      std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> assembly_contact_elepairs_;
    };

  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

#endif
