/*-----------------------------------------------------------*/
/*! \file

\brief Factory to create a merit function evaluation object.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_MERITFUNCTION_FACTORY_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_MERITFUNCTION_FACTORY_HPP

#include "baci_config.hpp"

#include "baci_solver_nonlin_nox_globaldata.hpp"

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace NOX
{
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
        virtual ~Factory() = default;

        /** \brief get a valid merit function pointer
         *
         *  choose between constraint and unconstraint */
        Teuchos::RCP<::NOX::MeritFunction::Generic> BuildMeritFunction(
            const NOX::NLN::GlobalData& noxNlnGlobalData) const;

       private:
        //! unconstraint factory
        Teuchos::RCP<::NOX::MeritFunction::Generic> BuildUnconstrainedMeritFunction(
            const std::string& mftype, const NOX::NLN::GlobalData& noxNlnGlobalData) const;

        //! constraint factory
        Teuchos::RCP<::NOX::MeritFunction::Generic> BuildConstrainedMeritFunction(
            const std::string& mftype, const NOX::NLN::GlobalData& noxNlnGlobalData) const;

      };  // class Factory

      /*! \brief Non-member function to build a merit function object.

      \relates NOX::NLNSOL::Constraint::MeritFunction::Factory

      */
      Teuchos::RCP<::NOX::MeritFunction::Generic> BuildMeritFunction(
          const NOX::NLN::GlobalData& noxNlnGlobalData);
    }  // namespace MeritFunction
  }    // namespace NLN
}  // namespace NOX

BACI_NAMESPACE_CLOSE

#endif  // SOLVER_NONLIN_NOX_MERITFUNCTION_FACTORY_H
