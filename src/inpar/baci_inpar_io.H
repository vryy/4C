/*----------------------------------------------------------------------*/
/*! \file
\file inpar_io.H

\brief Input parameters for global IO section


\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_INPAR_IO_H
#define BACI_INPAR_IO_H

#include "baci_config.H"

#include "baci_inpar_parameterlist_utils.H"

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
