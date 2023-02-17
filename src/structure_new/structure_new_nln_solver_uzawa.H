/*-----------------------------------------------------------*/
/*! \file

\brief uzawa (DEPRECATED)


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STRUCTURE_NEW_NLN_SOLVER_UZAWA_H
#define STRUCTURE_NEW_NLN_SOLVER_UZAWA_H

#include "structure_new_nln_solver_generic.H"

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

        virtual ~Uzawa(){};
      };
    }  // namespace SOLVER
  }    // namespace NLN
}  // namespace STR

#endif  // STRUCTURE_NEW_NLN_SOLVER_UZAWA_H
