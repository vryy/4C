/*-----------------------------------------------------------*/
/*! \file

\brief NOX's Newton with full step


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_NLN_SOLVER_FULLNEWTON_HPP
#define FOUR_C_STRUCTURE_NEW_NLN_SOLVER_FULLNEWTON_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_statustest_factory.hpp"
#include "4C_structure_new_nln_solver_nox.hpp"

namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

namespace Solid
{
  namespace Nln
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
        void setup() override;

       private:
        //! set the full newton parameters in the nox parameter list
        void set_full_newton_params();

      };  // class FullNewton
    }     // namespace SOLVER
  }       // namespace Nln
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
