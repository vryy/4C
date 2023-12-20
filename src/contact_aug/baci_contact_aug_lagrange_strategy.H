/*----------------------------------------------------------------------*/
/*! \file
\brief Lagrange contact solving strategy with standard Lagrangian
       multipliers based on the augmented Lagrange formulation.

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_CONTACT_AUG_LAGRANGE_STRATEGY_H
#define BACI_CONTACT_AUG_LAGRANGE_STRATEGY_H

#include "baci_config.H"

#include "baci_contact_aug_strategy.H"

BACI_NAMESPACE_OPEN

namespace CONTACT
{
  namespace AUG
  {
    namespace LAGRANGE
    {
      /*--------------------------------------------------------------------------*/
      /** \brief Standard Lagrange strategy based on the augmented Lagrangian
       *  strategy.
       *
       *  \author hiermeier \date 03/17 */
      class Strategy : public CONTACT::AUG::Strategy
      {
       public:
        /// constructor
        Strategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
            const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
            const Teuchos::ParameterList& params, const plain_interface_set& interfaces, int dim,
            const Teuchos::RCP<const Epetra_Comm>& comm, int maxdof);

        /// return the current solving strategy type
        INPAR::CONTACT::SolvingStrategy Type() const override
        {
          return INPAR::CONTACT::solution_std_lagrange;
        }

       protected:
        double InactiveScaleFactor() const override { return 1.0; }

        /** \brief Assemble the structural contact rhs [derived]
         *
         *  In contradiction to the base class only the the Lagrange multiplier
         *  forces are considered.
         *
         *  \author hiermeier \date 03/17 */
        void EvalStrContactRHS() override;

      };  // class Strategy
    }     // namespace LAGRANGE
  }       // namespace AUG
}  // namespace CONTACT


BACI_NAMESPACE_CLOSE

#endif  // CONTACT_AUG_LAGRANGE_STRATEGY_H
