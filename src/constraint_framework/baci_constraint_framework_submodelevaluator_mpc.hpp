/*-----------------------------------------------------------*/
/*! \file

\brief  Manage linear multipoint constraint equations including
        periodic displacement boundary conditions.

\level 3
 */
/*-----------------------------------------------------------*/

#ifndef FOUR_C_CONSTRAINT_FRAMEWORK_SUBMODELEVALUATOR_MPC_HPP
#define FOUR_C_CONSTRAINT_FRAMEWORK_SUBMODELEVALUATOR_MPC_HPP

#include "baci_config.hpp"

#include "baci_constraint_framework_submodelevaluator_base.hpp"
#include "baci_inpar_mpc_rve.hpp"
#include "baci_io_pstream.hpp"
#include "baci_lib_discret.hpp"
#include "baci_structure_new_model_evaluator_generic.hpp"

#include <boost/algorithm/string.hpp>
#include <Epetra_CrsMatrix.h>
#include <Teuchos_RCPDecl.hpp>

#include <map>

FOUR_C_NAMESPACE_OPEN


namespace CONSTRAINTS::SUBMODELEVALUATOR
{
  class MultiPointConstraintEquationBase;


  class RveMultiPointConstraintManager : public ConstraintBase
  {
   public:
    /*!
    \brief Standard Constructor
    */
    RveMultiPointConstraintManager(
        Teuchos::RCP<const DRT::Discretization> disc_ptr, CORE::LINALG::SparseMatrix* st_ptr);

    //! @name Public evaluation methods

    /*!
      \brief Perform basic checks of the input conditions and parameters
    */
    void CheckInput();

    /*!
     * \brief Reset the constraint stiffness matrix and delete node pairs
     */
    void Reset() override;

    //@}

   private:
    //! @name member variables

    //! Map of the Corner Node IDs Ni
    std::map<std::string, DRT::Node*> rveRefNodeMap_;

    //! RVE reference length vectors
    std::array<double, 2> r_xmxp_, r_ymyp_;

    //! Vector with all Conditions
    std::vector<Teuchos::RCP<DRT::Condition>> pointLinearCoupledEquationConditions_,
        pointPeriodicRveRefConditions_, linePeriodicRveConditions_, surfacePeriodicRveConditions_;

    //! Tolerance for the opposing edge node search
    double nodeSearchToler_ = 0.25;  // #ToDo: Add .dat parameter

    //! Parameter List for the rveType
    Teuchos::ParameterList mpcParameterList_;

    //! Dimension of the rve boundary
    enum INPAR::RVE_MPC::rveDimension rveDim_;

    //! Type of reference vector definition
    enum INPAR::RVE_MPC::rveReferenceDeformationDefinition rveRefType_;

    //@}

   private:
    //! @name Private evaluation methods

    //! find the opposite edge node pairs of the periodic rve
    void BuildPeriodicMPCs(std::map<std::string, const std::vector<int>*>& rveBoundaryNodeIdMap_,
        std::map<std::string, int>& rveCornerNodeIdMap_);

    //! add linear mpcs to the mpcList
    int BuildLinearMPCs();

    //! find a node that is member of edge1 and edge2
    int FindPeriodicRveCornerNodes(const std::vector<int>* edge1, const std::vector<int>* edge2);

    //! find a node that is member of surf1 thru surf3
    int FindPeriodicRveCornerNodes(const std::vector<int>* surf1, const std::vector<int>* surf2,
        const std::vector<int>* surf3);

    //! find the nodes containted i a mpc for the pbcs
    int FindOppositeEdgeNode(const int nodeID, INPAR::RVE_MPC::rveEdgeIdentifiers edge,
        std::map<std::string, const std::vector<int>*>& rveBoundaryNodeIdMap_);

    //! find the corner nodes of the periodic rve
    void BuildPeriodicRveCornerNodeMap(
        std::map<std::string, const std::vector<int>*>& rveBoundaryNodeIdMap,
        std::map<std::string, int>& rveCornerNodeIdMap_);

    //! retrive design line periodic rve boundary condition
    void BuildPeriodicRveBoundaryNodeMap(
        std::map<std::string, const std::vector<int>*>& rveBoundaryNodeIdMap_);
    //@}
  };
}  // namespace CONSTRAINTS::SUBMODELEVALUATOR
#endif
FOUR_C_NAMESPACE_CLOSE