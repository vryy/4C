/*-----------------------------------------------------------*/
/*! \file

\brief Factory for nonlinear solvers in structural dynamics


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_STRUCTURE_NEW_NLN_SOLVER_FACTORY_H
#define BACI_STRUCTURE_NEW_NLN_SOLVER_FACTORY_H

#include "baci_config.H"

#include "baci_inpar_structure.H"

BACI_NAMESPACE_OPEN

namespace STR
{
  namespace NLN
  {
    namespace SOLVER
    {
      class Generic;
      /*! \brief Factory to build the desired non-linear structural solver
       *
       *  \author Michael Hiermeier */
      class Factory
      {
       public:
        //! constructor
        Factory();

        //! destructor
        virtual ~Factory() = default;

        //! build the specific nonlinear solver
        Teuchos::RCP<STR::NLN::SOLVER::Generic> BuildNlnSolver(
            const enum INPAR::STR::NonlinSolTech& nlnSolType) const;
      };

      /*! Non-member function, which relates to the STR::NLN::SOLVER::Factory class
       *  Please call this method from outside! */
      Teuchos::RCP<STR::NLN::SOLVER::Generic> BuildNlnSolver(
          const enum INPAR::STR::NonlinSolTech& nlnSolType);
    }  // namespace SOLVER
  }    // namespace NLN
}  // namespace STR


BACI_NAMESPACE_CLOSE

#endif  // STRUCTURE_NEW_NLN_SOLVER_FACTORY_H
