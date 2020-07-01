/*---------------------------------------------------------------------*/
/*! \file
\brief Steepest ascent solution strategy based on the augmented contact
       formulation.

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_STEEPEST_ASCENT_STRATEGY_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_STEEPEST_ASCENT_STRATEGY_H_

#include "contact_aug_steepest_ascent_sp_strategy.H"

namespace CONTACT
{
  namespace AUG
  {
    class LagrangeMultiplierFunction;
    class PenaltyUpdate;
    namespace STEEPESTASCENT
    {
      // forward declarations
      class Interface;

      /*--------------------------------------------------------------------------*/
      /** \brief Condensed variant of the modified Newton approach.
       *
       * \author hiermeier \date 03/17 */
      class Strategy : public ::CONTACT::AUG::STEEPESTASCENT_SP::Strategy
      {
        /** The combo_strategy is a wrapper class for a set of augmented Lagrangian
         *  strategies and needs access to all methods. */
        friend class CONTACT::AUG::ComboStrategy;

       public:
        /// constructor
        Strategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
            const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
            const Teuchos::ParameterList& params, const plain_interface_set& interfaces, int dim,
            const Teuchos::RCP<const Epetra_Comm>& comm, int maxdof);

        /// destructor
        virtual ~Strategy(){};

        INPAR::CONTACT::SolvingStrategy Type() const override
        {
          return INPAR::CONTACT::solution_steepest_ascent;
        }

       protected:
        // un-do changes from the base class
        void EvalStrContactRHS() override { AUG::Strategy::EvalStrContactRHS(); }

        /// derived
        Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtrForNormCheck(
            const enum DRT::UTILS::VecBlockType& bt) const override;

        /// derived
        void AddContributionsToConstrRHS(Epetra_Vector& augConstrRhs) const override;

        /// derived
        Teuchos::RCP<LINALG::SparseMatrix> GetMatrixBlockPtr(
            const enum DRT::UTILS::MatBlockType& bt,
            const CONTACT::ParamsInterface* cparams = NULL) const override;

        /// derived
        void AddContributionsToMatrixBlockDisplDispl(LINALG::SparseMatrix& kdd,
            const CONTACT::ParamsInterface* cparams = NULL) const override;

        /// derived
        void RunPostApplyJacobianInverse(const CONTACT::ParamsInterface& cparams,
            const Epetra_Vector& rhs, Epetra_Vector& result, const Epetra_Vector& xold,
            const NOX::NLN::Group& grp) override;

        /// derived
        void RemoveCondensedContributionsFromRhs(Epetra_Vector& str_rhs) const override;

       private:
        void AugmentDirection(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
            Epetra_Vector& dir_mutable);

        Teuchos::RCP<Epetra_Vector> ComputeActiveLagrangeIncrInNormalDirection(
            const Epetra_Vector& displ_incr);

        Teuchos::RCP<Epetra_Vector> ComputeInactiveLagrangeIncrInNormalDirection(
            const Epetra_Vector& displ_incr, const Epetra_Vector& zold);

        void PostAugmentDirection(
            const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, Epetra_Vector& dir);
      };  //  class Strategy

    }  // namespace STEEPESTASCENT
  }    // namespace AUG
}  // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_STEEPEST_ASCENT_STRATEGY_H_ */
