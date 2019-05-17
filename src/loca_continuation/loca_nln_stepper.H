/*-----------------------------------------------------------*/
/*!

\brief LOCA stepper interface

\maintainer Martin Kronbichler

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_LOCA_CONTINUATION_LOCA_NLN_STEPPER_H_
#define SRC_LOCA_CONTINUATION_LOCA_NLN_STEPPER_H_

#include <LOCA_Stepper.H>

// forward declaration
namespace NOX
{
  namespace NLN
  {
    class GlobalData;
    namespace INNER
    {
      namespace StatusTest
      {
        class Generic;
      }  // namespace StatusTest
    }    // namespace INNER
  }      // namespace NLN
}  // namespace NOX

namespace LOCA
{
  namespace NLN
  {
    class Stepper : public LOCA::Stepper
    {
     public:
      //! Constructor with LOCA::StatusTest
      Stepper(const Teuchos::RCP<LOCA::GlobalData>& loca_gdata,
          const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& initial_guess,
          const Teuchos::RCP<LOCA::StatusTest::Abstract>& loca_test,
          const Teuchos::RCP<NOX::StatusTest::Generic>& nox_otest,
          const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& nox_itest,
          const Teuchos::RCP<NOX::NLN::GlobalData>& nln_gdata,
          const Teuchos::RCP<Teuchos::ParameterList>& p);

      //! Destructor
      virtual ~Stepper(){};

      //! derived
      virtual bool reset(const Teuchos::RCP<LOCA::GlobalData>& loca_gdata,
          const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& initial_guess,
          const Teuchos::RCP<LOCA::StatusTest::Abstract>& loca_test,
          const Teuchos::RCP<NOX::StatusTest::Generic>& nox_otest,
          const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& nox_itest,
          const Teuchos::RCP<NOX::NLN::GlobalData>& nln_gdata,
          const Teuchos::RCP<Teuchos::ParameterList>& p);

      //! forbid this method, since it is deprecated
      virtual bool reset(const Teuchos::RCP<LOCA::GlobalData>& loca_gdata,
          const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& initial_guess,
          const Teuchos::RCP<NOX::StatusTest::Generic>& nox_test,
          const Teuchos::RCP<Teuchos::ParameterList>& p);

     protected:
      //! reset solver
      virtual bool resetSolver();

      //! derived
      virtual LOCA::Abstract::Iterator::IteratorStatus start();

      //! derived
      virtual LOCA::Abstract::Iterator::IteratorStatus finish(
          LOCA::Abstract::Iterator::IteratorStatus itStatus);

      //! derived
      virtual LOCA::Abstract::Iterator::StepStatus preprocess(
          LOCA::Abstract::Iterator::StepStatus stepStatus);

     private:
      //! Private to prohibit copying
      Stepper(const Stepper&);

      //! Private to prohibit copying
      Stepper& operator=(const Stepper&);

     protected:
      Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> noxInnerStatusTestPtr_;

      Teuchos::RCP<NOX::NLN::GlobalData> noxNlnGlobalDataPtr_;


    };  // class Stepper
  }     // namespace NLN
}  // namespace LOCA

#endif /* SRC_LOCA_CONTINUATION_LOCA_NLN_STEPPER_H_ */
