/*----------------------------------------------------------------------*/
/*! \file
\file inpar_io.H

\brief Input parameters for global IO section


\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_INPAR_IO_HPP
#define FOUR_C_INPAR_IO_HPP

#include "baci_config.hpp"

#include "baci_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace IO
  {
    /*! \brief Define valid parameter for global IO control
     *
     * @param[in/out] list Parameter list to be filled with valid parameters and their defaults
     */
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace IO
}  // namespace INPAR

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
