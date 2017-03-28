/*----------------------------------------------------------------------*/
/*!
\file contact_aug_lagrange_strategy.H

\brief Lagrange contact solving strategy with standard Lagrangian
       multipliers based on the augmented Lagrange formulation.

\level 3

\maintainer Michael Hiermeier

\date Mar 28, 2017

*/
/*----------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_LAGRANGE_STRATEGY_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_LAGRANGE_STRATEGY_H_

#include "contact_augmented_strategy.H"

namespace CONTACT {
namespace AUG {
namespace LAGRANGE {

  /*--------------------------------------------------------------------------*/
  /** \brief Standard Lagrange strategy based on the augmented Lagrangian
   *  strategy.
   *
   *  \author hiermeier \date 03/17 */
  class Strategy : public ::CONTACT::AUG::Strategy
  {
  public:
    /// constructor
    Strategy(
        const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap,
        const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params,
        const plain_interface_set& interfaces,
        int dim,
        const Teuchos::RCP<const Epetra_Comm>& comm,
        int maxdof );

    /// destructor
    virtual ~Strategy() {};

    /// return the current solving strategy type
    virtual INPAR::CONTACT::SolvingStrategy Type() const
    {
      return INPAR::CONTACT::solution_std_lagrange;
    }

  protected:
    /** \brief Assemble the structural contact rhs [derived]
     *
     *  In contradiction to the base class only the the Lagrange multiplier
     *  forces are considered.
     *
     *  \author hiermeier \date 03/17 */
    virtual void EvalStrContactRHS();

  };  // class Strategy
} // namespace LAGRANGE
} // namespace AUG
} // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_LAGRANGE_STRATEGY_H_ */
