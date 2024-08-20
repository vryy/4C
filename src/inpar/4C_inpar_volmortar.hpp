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

namespace Inpar
{
  namespace VolMortar
  {
    /// set the volmortar parameters
    void set_valid_parameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace VolMortar
}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE

#endif
