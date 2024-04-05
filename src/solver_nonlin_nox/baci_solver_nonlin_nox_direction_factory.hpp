/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_DIRECTION_FACTORY_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_DIRECTION_FACTORY_HPP

#include "baci_config.hpp"

#include "baci_solver_nonlin_nox_forward_decl.hpp"

#include <NOX_Direction_UserDefinedFactory.H>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace NOX
{
  namespace NLN
  {
    namespace Direction
    {
      /*!
       \brief Factory to build direction objects derived from ::NOX::Direction::Generic.

       This factory class is closely related to the NOX_Direction_Factory.H. The main
       difference is that it allows to use direction methods which differ from the
       default NOX package.

       \author Michael Hiermeier
       */
      class Factory : public ::NOX::Direction::UserDefinedFactory
      {
       public:
        //! Constructor
        Factory();

        /*! \brief Factory to build user-defined direction objects.

            \note All direction methods which differ from the default NOX direction
            routines are user-defined. That means the user is the BACI user and the BACI
            programmer.

            @param gd A global data pointer that contains the top level parameter list.  Without
           storing this inside the direction object, there is no guarantee that the second parameter
           \c params will still exist.  It can be deleted by the top level RCP.
            @param params Sublist with direction construction parameters.

        */
        Teuchos::RCP<::NOX::Direction::Generic> buildDirection(
            const Teuchos::RCP<::NOX::GlobalData>& gd,
            Teuchos::ParameterList& params) const override;
      };

      /*! Nonmember function to build a direction object.

      \relates NOX::NLN::Direction::Factory

      */
      Teuchos::RCP<::NOX::Direction::Generic> BuildDirection(
          const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params);

    }  // namespace Direction
  }    // namespace NLN
}  // namespace NOX

BACI_NAMESPACE_CLOSE

#endif  // SOLVER_NONLIN_NOX_DIRECTION_FACTORY_H
