/*-----------------------------------------------------------*/
/*!
\file nox_nln_solver_linesearchbased.H

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_SOLVER_LINESEARCHBASED_H_
#define NOX_NLN_SOLVER_LINESEARCHBASED_H_

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

        //! Destructor
        virtual ~LineSearchBased(){};

        //! Necessary for the pre/post operator class
        //! Be careful: This version allows you to actually change the parameter list entries
        inline virtual const Teuchos::RCP<Teuchos::ParameterList>& GetMutableListPtr() const
        {
          return paramsPtr;
        };

        //! Returns the stopping test.
        virtual const NOX::StatusTest::Generic& GetOuterStatusTest() const;

        virtual NOX::StatusTest::StatusType getStatus() const;

        //! \brief Returns a pointer to one specific stopping test %T.
        /** If there is no outer test of type T, a NULL pointer will be returned.
         *
         *  \author hiermeier \date 04/17 */
        template <class T>
        NOX::StatusTest::Generic* GetOuterStatusTest() const;

        //! \brief Returns a pointer to one specific stopping test %T.
        /** If there is no outer test of type T with quantity type qtype, a NULL
         *  pointer will be returned.
         *
         *  \author hiermeier \date 08/18 */
        template <class T>
        NOX::StatusTest::Generic* GetOuterStatusTestWithQuantity(
            const NOX::NLN::StatusTest::QuantityType qtype) const;

        //! Returns the NOX::Utils object
        const NOX::Utils& GetUtils() const;

        //! Returns the global status of the given test name
        template <class T>
        NOX::StatusTest::StatusType GetStatus() const;

        /// access the direction object
        NOX::Direction::Generic& GetDirection() const;

       protected:
        //! initialize additional variables or overwrite base class initialization
        virtual void init(const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests);

        // derived
        virtual void printUpdate();
      };  // class LineSearchBased
    }     // namespace Solver
  }       // namespace NLN
}  // namespace NOX

#endif /* NOX_NLN_CONSTRAINT_SOLVER_LINESEARCHBASED_H_ */
