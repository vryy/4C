/*-----------------------------------------------------------------------*/
/*! \file

\brief

\level 1

*/

/*----------------------------------------------------------------------*/

#ifndef FOUR_C_INPAR_VOLMORTAR_HPP
#define FOUR_C_INPAR_VOLMORTAR_HPP


/*----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCPDecl.hpp>

FOUR_C_NAMESPACE_OPEN

namespace INPAR
{
  namespace VOLMORTAR
  {
    /// set the volmortar parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace VOLMORTAR
}  // namespace INPAR

FOUR_C_NAMESPACE_CLOSE

#endif
