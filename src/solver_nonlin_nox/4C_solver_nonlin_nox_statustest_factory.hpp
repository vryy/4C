/*-----------------------------------------------------------*/
/*! \file

\brief Create the desired outer status test.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_STATUSTEST_FACTORY_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_STATUSTEST_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_forward_decl.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace NLN
  {
    namespace StatusTest
    {
      class Factory
      {
       public:
        //! Constructor.
        Factory();

        //! Destructor.
        virtual ~Factory() = default;

        //! Returns a outer status test set from a parameter list.
        Teuchos::RCP<::NOX::StatusTest::Generic> BuildOuterStatusTests(Teuchos::ParameterList& p,
            const ::NOX::Utils& utils,
            std::map<std::string, Teuchos::RCP<::NOX::StatusTest::Generic>>* tagged_tests) const;

       protected:
        /*! \brief New implementation of the \c BuildNormFTest function.
         *
         *  The underlying Status Test is capable of a variety of quantities at the same
         *  time and combines the RelativeNormF and NormF test.
         */
        Teuchos::RCP<::NOX::StatusTest::Generic> BuildNormFTest(
            Teuchos::ParameterList& p, const ::NOX::Utils& u, const bool& relativeNormF) const;
        Teuchos::RCP<::NOX::StatusTest::Generic> BuildNormFTest(
            Teuchos::ParameterList& p, const ::NOX::Utils& u) const;

        /*! \brief New implementation of the \c BuildNormUpdateTest function.
         *
         *  The underlying Status Test is capable of a variety of quantities at the same
         *  time.
         */
        Teuchos::RCP<::NOX::StatusTest::Generic> BuildNormUpdateTest(
            Teuchos::ParameterList& p, const ::NOX::Utils& u) const;

        /*! \brief New implementation of the \c BuildNormWRMSTest function.
         *
         *  The underlying Status Test is capable of a variety of quantities at the same
         *  time.
         */
        Teuchos::RCP<::NOX::StatusTest::Generic> BuildNormWRMSTest(
            Teuchos::ParameterList& p, const ::NOX::Utils& u) const;

        /*! \brief Simple active set test
         *
         *  This active set test checks the active set status of the corresponding
         *  quantity. Basically we check if the active set did not change.
         *
         *  OPTIONAL: Based on the implementation of the internal get function, also
         *  a cycling of the active set can be detected and printed to the screen.
         */
        Teuchos::RCP<::NOX::StatusTest::Generic> BuildActiveSetTest(
            Teuchos::ParameterList& p, const ::NOX::Utils& u) const;

        /*! \brief Derived buildComboTest.
         *
         *  Restatement of the base class function, because the internal recursive call
         *  has to be redefined. Note that the base class input parameter
         *  <tt>Number of Tests<\tt> is not used in this version.
         */
        Teuchos::RCP<::NOX::StatusTest::Generic> BuildComboTest(Teuchos::ParameterList& p,
            const ::NOX::Utils& u,
            std::map<std::string, Teuchos::RCP<::NOX::StatusTest::Generic>>* tagged_tests) const;

        //! Checks if a tag is present in the parameter list and adds the test to the tagged_test
        //! std::map if true.  Returns true if a tag was present.
        bool CheckAndTagTest(const Teuchos::ParameterList& p,
            const Teuchos::RCP<::NOX::StatusTest::Generic>& test,
            std::map<std::string, Teuchos::RCP<::NOX::StatusTest::Generic>>* tagged_tests) const;

       private:
        //! Throws formated error
        void throwError(const std::string& functionName, const std::string& errorMsg) const;

       private:
        /*! \brief Reference to the base class ::NOX::StatusTest::Factory
         *
         *  Used for a direct call instead of the member function to prevent too
         *  many calls of the factory constructor.
         */
        Teuchos::RCP<const ::NOX::StatusTest::Factory> noxfactory_;

      };  // class Factory
      /*! \brief Nonmember helper function for the NOX::NLN::StatusTest::Factory.
       *
       *  \relates NOX::NLN::StatusTest::Factory
       *
       */
      Teuchos::RCP<::NOX::StatusTest::Generic> BuildOuterStatusTests(Teuchos::ParameterList& p,
          const ::NOX::Utils& utils,
          std::map<std::string, Teuchos::RCP<::NOX::StatusTest::Generic>>* tagged_tests = nullptr);

    }  // namespace StatusTest
  }    // namespace NLN
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
