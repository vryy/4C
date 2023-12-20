/*-----------------------------------------------------------*/
/*! \file

\brief NOX's Newton with full step


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_STRUCTURE_NEW_NLN_SOLVER_FULLNEWTON_H
#define BACI_STRUCTURE_NEW_NLN_SOLVER_FULLNEWTON_H

#include "baci_config.H"

#include "baci_solver_nonlin_nox_statustest_factory.H"
#include "baci_structure_new_nln_solver_nox.H"

namespace Teuchos
{
  class ParameterList;
}

BACI_NAMESPACE_OPEN

namespace STR
{
  namespace NLN
  {
    namespace SOLVER
    {
      /*! \brief Newton's method with full step via NOX for structural dynamics
       *
       */
      class FullNewton : public Nox
      {
       public:
        //! constructor
        FullNewton();


        //! derived from the base class
        void Setup() override;

       private:
        //! set the full newton parameters in the nox parameter list
        void SetFullNewtonParams();

      };  // class FullNewton
    }     // namespace SOLVER
  }       // namespace NLN
}  // namespace STR


BACI_NAMESPACE_CLOSE

#endif  // STRUCTURE_NEW_NLN_SOLVER_FULLNEWTON_H
