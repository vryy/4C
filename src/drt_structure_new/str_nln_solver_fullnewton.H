/*-----------------------------------------------------------*/
/*! \file

\brief NOX's Newton with full step


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_NLN_SOLVER_FULLNEWTON_H_
#define STR_NLN_SOLVER_FULLNEWTON_H_

#include "str_nln_solver_nox.H"
#include "../solver_nonlin_nox/nox_nln_statustest_factory.H"

namespace Teuchos
{
  class ParameterList;
}

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

        //! destructor
        virtual ~FullNewton(){};

        //! derived from the base class
        virtual void Setup();

       private:
        //! set the full newton parameters in the nox parameter list
        void SetFullNewtonParams();

      };  // class FullNewton
    }     // namespace SOLVER
  }       // namespace NLN
}  // namespace STR


#endif /* STR_NLN_SOLVER_FULLNEWTON_H_ */
