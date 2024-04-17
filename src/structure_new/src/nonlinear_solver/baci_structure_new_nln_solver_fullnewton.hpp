/*-----------------------------------------------------------*/
/*! \file

\brief NOX's Newton with full step


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_NLN_SOLVER_FULLNEWTON_HPP
#define FOUR_C_STRUCTURE_NEW_NLN_SOLVER_FULLNEWTON_HPP

#include "baci_config.hpp"

#include "baci_solver_nonlin_nox_statustest_factory.hpp"
#include "baci_structure_new_nln_solver_nox.hpp"

namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

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


FOUR_C_NAMESPACE_CLOSE

#endif
