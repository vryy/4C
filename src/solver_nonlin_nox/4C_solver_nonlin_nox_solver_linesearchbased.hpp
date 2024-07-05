/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_SOLVER_LINESEARCHBASED_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_SOLVER_LINESEARCHBASED_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_forward_decl.hpp"

#include <NOX_Solver_LineSearchBased.H>  // base class

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    namespace StatusTest
    {
      enum QuantityType : int;
    }  // namespace StatusTest
    namespace Inner
    {
      namespace StatusTest
      {
        class Generic;
      }  // namespace StatusTest
    }    // namespace Inner
    namespace Solver
    {
      class LineSearchBased : public ::NOX::Solver::LineSearchBased
      {
       public:
        //! Constructor
        /*!
          See reset(::NOX::Abstract::Group&, ::NOX::StatusTest::Generic&, Teuchos::ParameterList&)
          for description
         */
        LineSearchBased(const Teuchos::RCP<::NOX::Abstract::Group>& grp,
            const Teuchos::RCP<::NOX::StatusTest::Generic>& outerTests,
            const Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>& innerTests,
            const Teuchos::RCP<Teuchos::ParameterList>& params);

        //! Necessary for the pre/post operator class
        //! Be careful: This version allows you to actually change the parameter list entries
        inline virtual const Teuchos::RCP<Teuchos::ParameterList>& get_list_ptr() const
        {
          return paramsPtr;
        };

        //! Returns the stopping test.
        virtual const ::NOX::StatusTest::Generic& get_outer_status_test() const;

        ::NOX::StatusTest::StatusType getStatus() const override;

        //! \brief Returns a pointer to one specific stopping test %T.
        /** If there is no outer test of type T, a nullptr pointer will be returned.
         *
         *  \author hiermeier \date 04/17 */
        template <class T>
        ::NOX::StatusTest::Generic* get_outer_status_test() const;

        //! \brief Returns a pointer to one specific stopping test %T.
        /** If there is no outer test of type T with quantity type qtype, a nullptr
         *  pointer will be returned.
         *
         *  \author hiermeier \date 08/18 */
        template <class T>
        ::NOX::StatusTest::Generic* get_outer_status_test_with_quantity(
            const NOX::Nln::StatusTest::QuantityType qtype) const;

        //! Returns the ::NOX::Utils object
        [[nodiscard]] const ::NOX::Utils& get_utils() const;

        //! Returns the global status of the given test name
        template <class T>
        [[nodiscard]] ::NOX::StatusTest::StatusType get_status() const;

        /// access the direction object
        ::NOX::Direction::Generic& get_direction() const;

       protected:
        //! initialize additional variables or overwrite base class initialization
        virtual void init(const Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>& innerTests);

        void printUpdate() override;
      };  // class LineSearchBased
    }     // namespace Solver
  }       // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
