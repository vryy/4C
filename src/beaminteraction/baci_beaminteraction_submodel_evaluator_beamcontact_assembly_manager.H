/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble pair based contributions into global matrices.


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BACI_BEAMINTERACTION_SUBMODEL_EVALUATOR_BEAMCONTACT_ASSEMBLY_MANAGER_H
#define BACI_BEAMINTERACTION_SUBMODEL_EVALUATOR_BEAMCONTACT_ASSEMBLY_MANAGER_H


#include "baci_utils_exceptions.H"

#include <Teuchos_RCP.hpp>

#include <vector>

// Forward declarations.
class Epetra_FEVector;
class Epetra_Vector;

BACI_NAMESPACE_OPEN

namespace CORE::LINALG
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
       */
      BeamContactAssemblyManager();

      /**
       * \brief Destructor.
       */
      virtual ~BeamContactAssemblyManager() = default;

      /**
       * \brief Evaluate all force and stiffness terms and add them to the global matrices.
       * @param discret (in) Pointer to the disretization.
       * @param data_state (in) Beam interaction data state.
       * @param fe_sysvec (out) Global force vector.
       * @param fe_sysmat (out) Global stiffness matrix.
       */
      virtual void EvaluateForceStiff(Teuchos::RCP<DRT::Discretization> discret,
          const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>& data_state,
          Teuchos::RCP<Epetra_FEVector> fe_sysvec,
          Teuchos::RCP<CORE::LINALG::SparseMatrix> fe_sysmat)
      {
        dserror("Not implemented!");
      }

      virtual double GetEnergy(const Teuchos::RCP<const Epetra_Vector>& disp) const { return 0.0; }
    };

  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

BACI_NAMESPACE_CLOSE

#endif
