/*-----------------------------------------------------------*/
/*! \file

\brief Structural single step solver for explicit dynamics

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_STRUCTURE_NEW_NLN_SOLVER_SINGLESTEP_H
#define BACI_STRUCTURE_NEW_NLN_SOLVER_SINGLESTEP_H

#include "baci_solver_nonlin_nox_statustest_factory.H"
#include "baci_structure_new_nln_solver_nox.H"

namespace Teuchos
{
  class ParameterList;
}

BACI_NAMESPACE_OPEN

namespace STR::NLN::SOLVER
{
  /*! \brief Full step with single iteration via NOX for explicit structural dynamics
   *
   */
  class SingleStep : public Nox
  {
   public:
    //! derived from the base class
    void Setup() override;

    //! derived from the base class
    INPAR::STR::ConvergenceStatus Solve() override;

   protected:
    //! Reset the non-linear solver parameters and variables
    void ResetParams() override;

   private:
    //! set the full newton parameters in the nox parameter list
    void SetSingleStepParams();

    void SetSingleStepParams(Teuchos::ParameterList& p);

  };  // class SingleStep
}  // namespace STR::NLN::SOLVER


BACI_NAMESPACE_CLOSE

#endif /* BACI_STRUCTURE_NEW_NLN_SOLVER_SINGLESTEP_H */
