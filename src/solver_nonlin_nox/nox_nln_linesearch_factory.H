/*-----------------------------------------------------------*/
/*! \file

\brief Factory to create the desired Line Search object.

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_LINESEARCH_FACTORY_H_
#define NOX_NLN_LINESEARCH_FACTORY_H_

#include <Teuchos_RCP.hpp>
#include <NOX_Common.H>

// Forward declarations
namespace Teuchos
{
  class ParameterList;
}
namespace NOX
{
  class Utils;

  namespace LineSearch
  {
    class Generic;
  }
}  // namespace NOX

namespace NOX
{
  class GlobalData;
  namespace StatusTest
  {
    class Generic;
  }
  namespace NLN
  {
    namespace INNER
    {
      namespace StatusTest
      {
        class Generic;
      }  // namespace StatusTest
    }    // namespace INNER
    namespace LineSearch
    {
      class Factory
      {
       public:
        //! Constructor
        Factory();

        //! Destructor
        ~Factory(){};

        /*! \brief Factory to build a line search object.

            @param gd A global data pointer that contains the top level parameter list.  Without
           storing this inside the line searchobject, there is no guarantee that the second
           parameter \c params will still exist.  It can be deleted by the top level RCP.
            @param params General nln parameterlist.

        */
        Teuchos::RCP<NOX::LineSearch::Generic> BuildLineSearch(
            const Teuchos::RCP<NOX::GlobalData>& gd,
            const Teuchos::RCP<NOX::StatusTest::Generic> outerTests,
            const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> innerTests,
            Teuchos::ParameterList& params);

       private:
        // checks if the inner status test pointer is initialized
        void InnerStatusTestIsRequired(
            const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests) const;
      };
      /*! Nonmember function to build a line search object.

      \relates NOX::NLNSOL::Constraint::LineSearch::Factory

      */
      Teuchos::RCP<NOX::LineSearch::Generic> BuildLineSearch(
          const Teuchos::RCP<NOX::GlobalData>& gd,
          const Teuchos::RCP<NOX::StatusTest::Generic> outerTests,
          const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> innerTests,
          Teuchos::ParameterList& params);

    }  // namespace LineSearch
  }    // namespace NLN
}  // namespace NOX

#endif /* NOX_NLN_LINESEARCH_FACTORY_H_ */
