/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_SOLVER_NONLIN_NOX_SOLVER_LINESEARCHBASED_H
#define BACI_SOLVER_NONLIN_NOX_SOLVER_LINESEARCHBASED_H

#include <NOX_Solver_LineSearchBased.H>  // base class

// forward declaration
namespace NOX
{
  namespace Abstract
  {
    class Group;
  }  // namespace Abstract

  namespace StatusTest
  {
    class Generic;
  }  // namespace StatusTest
}  // namespace NOX

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

namespace NOX
{
  namespace NLN
  {
    namespace StatusTest
    {
      enum QuantityType : int;
    }  // namespace StatusTest
    namespace INNER
    {
      namespace StatusTest
      {
        class Generic;
      }  // namespace StatusTest
    }    // namespace INNER
    namespace Solver
    {
      class LineSearchBased : public NOX::Solver::LineSearchBased
      {
       public:
        //! Constructor
        /*!
          See reset(NOX::Abstract::Group&, NOX::StatusTest::Generic&, Teuchos::ParameterList&) for
          description
         */
        LineSearchBased(const Teuchos::RCP<NOX::Abstract::Group>& grp,
            const Teuchos::RCP<NOX::StatusTest::Generic>& outerTests,
            const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests,
            const Teuchos::RCP<Teuchos::ParameterList>& params);


        //! Necessary for the pre/post operator class
        //! Be careful: This version allows you to actually change the parameter list entries
        inline virtual const Teuchos::RCP<Teuchos::ParameterList>& GetListPtr() const
        {
          return paramsPtr;
        };

        //! Returns the stopping test.
        virtual const NOX::StatusTest::Generic& GetOuterStatusTest() const;

        NOX::StatusTest::StatusType getStatus() const override;

        //! \brief Returns a pointer to one specific stopping test %T.
        /** If there is no outer test of type T, a nullptr pointer will be returned.
         *
         *  \author hiermeier \date 04/17 */
        template <class T>
        NOX::StatusTest::Generic* GetOuterStatusTest() const;

        //! \brief Returns a pointer to one specific stopping test %T.
        /** If there is no outer test of type T with quantity type qtype, a nullptr
         *  pointer will be returned.
         *
         *  \author hiermeier \date 08/18 */
        template <class T>
        NOX::StatusTest::Generic* GetOuterStatusTestWithQuantity(
            const NOX::NLN::StatusTest::QuantityType qtype) const;

        //! Returns the NOX::Utils object
        [[nodiscard]] const NOX::Utils& GetUtils() const;

        //! Returns the global status of the given test name
        template <class T>
        [[nodiscard]] NOX::StatusTest::StatusType GetStatus() const;

        /// access the direction object
        NOX::Direction::Generic& GetDirection() const;

       protected:
        //! initialize additional variables or overwrite base class initialization
        virtual void init(const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests);

        void printUpdate() override;
      };  // class LineSearchBased
    }     // namespace Solver
  }       // namespace NLN
}  // namespace NOX

#endif  // SOLVER_NONLIN_NOX_SOLVER_LINESEARCHBASED_H
