/*-----------------------------------------------------------*/
/*!
\file nox_nln_prepostoperator.H

\maintainer Michael Hiermeier

\date Jul 29, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_PREPOSTOPERATOR_H_
#define NOX_NLN_PREPOSTOPERATOR_H_

#include <NOX_Abstract_PrePostOperator.H>   // base class

namespace NOX {
namespace NLN {
  class PrePostOperator : public NOX::Abstract::PrePostOperator
  {
  public:
    //! constructor
    PrePostOperator();

    //! destructor
    virtual ~PrePostOperator() {};

    /*! \brief derived

       User defined method that will be executed at the start of a call to NOX::Solver::Generic::step().
    */
    virtual void runPreIterate(const NOX::Solver::Generic& solver);

    //! derived
    virtual void runPreSolve(const NOX::Solver::Generic& nlnSolver);

  };
} // namespace NLN
} // namespace NOX


#endif /* NOX_NLN_PREPOSTOPERATOR_H_ */
