/*-----------------------------------------------------------*/
/*!
\file nox_nln_inner_statustest_factory.H

\brief factory for user defined NOX inner status tests

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_INNER_STATUSTEST_FACTORY_H_
#define NOX_NLN_INNER_STATUSTEST_FACTORY_H_

#include <Teuchos_RCP.hpp>

// forward declaration
namespace Teuchos
{
  class ParameterList;
}

namespace NOX
{
  class Utils;
  namespace NLN
  {
    namespace INNER
    {
      namespace StatusTest
      {
        class Generic;

        class Factory
        {
         public:
          //! Constructor.
          Factory();

          //! Destructor.
          virtual ~Factory(){};

          //! Returns a inner status test set from a parameter list.
          Teuchos::RCP<Generic> BuildInnerStatusTests(Teuchos::ParameterList& p,
              const NOX::Utils& utils,
              std::map<std::string, Teuchos::RCP<Generic>>* tagged_tests) const;

         protected:
          //! Build the Armijo sufficient decrease test
          Teuchos::RCP<Generic> BuildArmijoTest(
              Teuchos::ParameterList& p, const NOX::Utils& u) const;

          //! Build the filter test
          Teuchos::RCP<Generic> BuildFilterTest(
              Teuchos::ParameterList& p, const NOX::Utils& u) const;

          //! Build the upper bound test
          Teuchos::RCP<Generic> BuildUpperBoundTest(
              Teuchos::ParameterList& p, const NOX::Utils& u) const;

          /*! \brief Build the inner_statustest_combo object

          */
          Teuchos::RCP<Generic> BuildComboTest(Teuchos::ParameterList& p, const NOX::Utils& u,
              std::map<std::string, Teuchos::RCP<Generic>>* tagged_tests) const;

          /// \brief Build volume change test
          Teuchos::RCP<Generic> BuildVolumeChangeTest(
              Teuchos::ParameterList& p, const NOX::Utils& u) const;

          bool CheckAndTagTest(const Teuchos::ParameterList& p, const Teuchos::RCP<Generic>& test,
              std::map<std::string, Teuchos::RCP<Generic>>* tagged_tests) const;

         private:
          //! Throws formated error
          void throwError(const std::string& functionName, const std::string& errorMsg) const;
        };

        /*! \brief Nonmember helper function for the NOX::NLN::INNER::StatusTest::Factory.

        \relates NOX::NLN::INNER::StatusTest::Factory

        */
        Teuchos::RCP<Generic> BuildInnerStatusTests(Teuchos::ParameterList& p,
            const NOX::Utils& utils,
            std::map<std::string, Teuchos::RCP<Generic>>* tagged_tests = 0);
      }  // namespace StatusTest
    }    // namespace INNER
  }      // namespace NLN
}  // namespace NOX


#endif /* NOX_NLN_INNER_STATUSTEST_FACTORY_H_ */
