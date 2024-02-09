/*----------------------------------------------------------------------*/
/*! \file
\file inpar_io.H

\brief Input parameters for global IO section


\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_INPAR_IO_HPP
#define BACI_INPAR_IO_HPP

#include "baci_config.hpp"

#include "baci_inpar_parameterlist_utils.hpp"

BACI_NAMESPACE_OPEN

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
BACI_NAMESPACE_CLOSE

#endif  // INPAR_IO_H
