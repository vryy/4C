/*-----------------------------------------------------------*/
/*! \file

\brief uzawa (DEPRECATED)


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_NLN_SOLVER_UZAWA_HPP
#define FOUR_C_STRUCTURE_NEW_NLN_SOLVER_UZAWA_HPP

#include "baci_config.hpp"

#include "baci_structure_new_nln_solver_generic.hpp"

FOUR_C_NAMESPACE_OPEN

namespace STR
{
  namespace NLN
  {
    namespace SOLVER
    {
      // Can be deleted in an upcoming commit, since unnecessary/deprecated   Hiermeier 01/12/2015
      class Uzawa : public Generic
      {
       public:
        Uzawa(){};
      };
    }  // namespace SOLVER
  }    // namespace NLN
}  // namespace STR

FOUR_C_NAMESPACE_CLOSE

#endif
