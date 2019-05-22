/*-----------------------------------------------------------*/
/*!
\file nox_nln_meritfunction_factory.H

\brief Factory to create a merit function evaluation object.

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_MERITFUNCTION_FACTORY_H_
#define NOX_NLN_MERITFUNCTION_FACTORY_H_

#include <Teuchos_RCP.hpp>
#include "nox_nln_globaldata.H"

// Forward declarations
namespace Teuchos
{
  class ParameterList;
}
namespace NOX
{
  class Utils;
  namespace MeritFunction
  {
    class Generic;
  }  // namespace MeritFunction
  namespace NLN
  {
    namespace MeritFunction
    {
      class Factory
      {
       public:
        //! constructor
        Factory();

        //! destructor
        virtual ~Factory(){};

        /** \brief get a valid merit function pointer
         *
         *  choose between constraint and unconstraint */
        Teuchos::RCP<NOX::MeritFunction::Generic> BuildMeritFunction(
            const NOX::NLN::GlobalData& noxNlnGlobalData) const;

       private:
        //! unconstraint factory
        Teuchos::RCP<NOX::MeritFunction::Generic> BuildUnconstrainedMeritFunction(
            const std::string& mftype, const NOX::NLN::GlobalData& noxNlnGlobalData) const;

        //! constraint factory
        Teuchos::RCP<NOX::MeritFunction::Generic> BuildConstrainedMeritFunction(
            const std::string& mftype, const NOX::NLN::GlobalData& noxNlnGlobalData) const;

      };  // class Factory

      /*! \brief Non-member function to build a merit function object.

      \relates NOX::NLNSOL::Constraint::MeritFunction::Factory

      */
      Teuchos::RCP<NOX::MeritFunction::Generic> BuildMeritFunction(
          const NOX::NLN::GlobalData& noxNlnGlobalData);
    }  // namespace MeritFunction
  }    // namespace NLN
}  // namespace NOX

#endif /* NOX_NLN_MERITFUNCTION_FACTORY_H_ */
