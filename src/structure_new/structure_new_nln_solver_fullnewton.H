/*-----------------------------------------------------------*/
/*! \file

\brief NOX's Newton with full step


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STRUCTURE_NEW_NLN_SOLVER_FULLNEWTON_H
#define STRUCTURE_NEW_NLN_SOLVER_FULLNEWTON_H

#include "structure_new_nln_solver_nox.H"
#include "solver_nonlin_nox_statustest_factory.H"

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
        void Setup() override;

       private:
        //! set the full newton parameters in the nox parameter list
        void SetFullNewtonParams();

      };  // class FullNewton
    }     // namespace SOLVER
  }       // namespace NLN
}  // namespace STR


#endif  // STRUCTURE_NEW_NLN_SOLVER_FULLNEWTON_H
