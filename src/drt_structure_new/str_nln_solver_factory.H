/*-----------------------------------------------------------*/
/*!
\file str_nln_solver_factory.H

\maintainer Michael Hiermeier

\date Aug 13, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_NLN_SOLVER_FACTORY_H_
#define STR_NLN_SOLVER_FACTORY_H_

#include "../drt_inpar/inpar_structure.H"

// forward declaration
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


#endif /* STR_NLN_SOLVER_FACTORY_H_ */
