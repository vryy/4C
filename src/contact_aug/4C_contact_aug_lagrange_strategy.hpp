/*----------------------------------------------------------------------*/
/*! \file
\brief Lagrange contact solving strategy with standard Lagrangian
       multipliers based on the augmented Lagrange formulation.

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_LAGRANGE_STRATEGY_HPP
#define FOUR_C_CONTACT_AUG_LAGRANGE_STRATEGY_HPP

#include "4C_config.hpp"

#include "4C_contact_aug_strategy.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  namespace Aug
  {
    namespace Lagrange
    {
      /*--------------------------------------------------------------------------*/
      /** \brief Standard Lagrange strategy based on the augmented Lagrangian
       *  strategy.
       *
       *  \author hiermeier \date 03/17 */
      class Strategy : public CONTACT::Aug::Strategy
      {
       public:
        /// constructor
        Strategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
            const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
            const Teuchos::ParameterList& params, const plain_interface_set& interfaces, int dim,
            const Teuchos::RCP<const Epetra_Comm>& comm, int maxdof);

        /// return the current solving strategy type
        Inpar::CONTACT::SolvingStrategy Type() const override
        {
          return Inpar::CONTACT::solution_std_lagrange;
        }

       protected:
        double inactive_scale_factor() const override { return 1.0; }

        /** \brief Assemble the structural contact rhs [derived]
         *
         *  In contradiction to the base class only the the Lagrange multiplier
         *  forces are considered.
         *
         *  \author hiermeier \date 03/17 */
        void evaluate_str_contact_rhs() override;

      };  // class Strategy
    }     // namespace Lagrange
  }       // namespace Aug
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
