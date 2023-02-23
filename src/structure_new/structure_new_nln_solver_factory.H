/*-----------------------------------------------------------*/
/*! \file

\brief Factory for nonlinear solvers in structural dynamics


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STRUCTURE_NEW_NLN_SOLVER_FACTORY_H
#define STRUCTURE_NEW_NLN_SOLVER_FACTORY_H

#include "inpar_structure.H"

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
        virtual ~Factory(){};

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


#endif  // STRUCTURE_NEW_NLN_SOLVER_FACTORY_H
