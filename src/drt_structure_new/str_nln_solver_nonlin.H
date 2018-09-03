/*-----------------------------------------------------------*/
/*!
\file str_nln_solver_nonlin.H

\maintainer Matthias Mayr

\date Aug 13, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_NLN_SOLVER_NONLIN_H_
#define STR_NLN_SOLVER_NONLIN_H_

#include "str_nln_solver_generic.H"

// forward declaration
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos
namespace NLNSOL
{
  class NlnProblem;
  class NlnOperatorBase;
}  // namespace NLNSOL
namespace NOX
{
  class Utils;
  namespace Epetra
  {
    class Group;
    class Vector;
    class LinearSystem;
  }  // namespace Epetra
}  // namespace NOX
namespace STR
{
  namespace NLN
  {
    namespace SOLVER
    {
      class Nonlin : public Generic
      {
       public:
        Nonlin();

        virtual ~Nonlin(){};

        //! derived from the base class
        virtual void Setup();

        //! derived from the base class
        virtual void Reset();

        //! derived from the base class
        virtual INPAR::STR::ConvergenceStatus Solve();

       private:
        Teuchos::RCP<Teuchos::ParameterList> NoxCreatePrintParameters(const bool verbose) const;

        Teuchos::RCP<NOX::Epetra::LinearSystem> NoxCreateLinearSystem(
            Teuchos::ParameterList& nlParams, NOX::Epetra::Vector& noxsoln);

       private:
        //! interface to the non-linear problem
        Teuchos::RCP<NLNSOL::NlnProblem> nlnproblem_;

        //! non-linear operator to solve the non-linear problem
        Teuchos::RCP<NLNSOL::NlnOperatorBase> nlnoperator_;

        //! nox solution vector
        Teuchos::RCP<NOX::Epetra::Vector> noxsoln_;

        //! nox utils
        Teuchos::RCP<NOX::Utils> noxutils_;

        //! nonlin parameter list
        Teuchos::RCP<Teuchos::ParameterList> params_;
      };  // class Nonlin
    }     // namespace SOLVER
  }       // namespace NLN
}  // namespace STR


#endif /* STR_NLN_SOLVER_NONLIN_H_ */
