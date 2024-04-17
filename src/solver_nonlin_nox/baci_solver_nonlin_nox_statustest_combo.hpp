/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_STATUSTEST_COMBO_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_STATUSTEST_COMBO_HPP

#include "baci_config.hpp"

#include <NOX_StatusTest_Combo.H>  // base class

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace NLN
  {
    namespace StatusTest
    {
      /*! \brief Arbitrary combination of status tests.
       *
       * Please see "NOX_StatusTest_Combo.H" for more information. This
       * is a restatement with no changes to the functionality of the class.
       *
       * \author Michael Hiermeier
       */

      class Combo : public ::NOX::StatusTest::Combo
      {
       public:
        //! Constructor. Optional argument is the error stream for output.
        Combo(ComboType t, const ::NOX::Utils* u = nullptr);

        //! Constructor with a single test.
        Combo(ComboType t, const Teuchos::RCP<Generic>& a, const ::NOX::Utils* u = nullptr);

        //! Constructor with two tests.
        Combo(ComboType t, const Teuchos::RCP<Generic>& a, const Teuchos::RCP<Generic>& b,
            const ::NOX::Utils* u = nullptr);

        //! Add another test to this combination.
        /*!
          Calls isSafe() to determine if it is safe to add \c a to the combination.
        */
        Combo& addStatusTest(const Teuchos::RCP<::NOX::StatusTest::Generic>& a) override;

        virtual Combo& addStatusTest(
            const Teuchos::RCP<::NOX::StatusTest::Generic>& a, const bool& init);

        virtual const std::vector<Teuchos::RCP<::NOX::StatusTest::Generic>>& GetTestVector() const;

       protected:
        //! Check whether or not it is safe to add a to this list of tests.
        bool isSafe(Generic& a);

       protected:
        /*! \brief Vector of generic status tests
         *
         * This a copy of the base class vector. It generates almost no overhead,
         * because it is a vector of pointers and we need it, because we have no
         * direct access to the base class vector (private). */
        std::vector<Teuchos::RCP<::NOX::StatusTest::Generic>> tests_;

        //! Ostream used to print errors
        ::NOX::Utils utils_;

      };  // class Combo

    }  // namespace StatusTest
  }    // namespace NLN
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
