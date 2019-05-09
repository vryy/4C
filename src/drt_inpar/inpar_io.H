/*----------------------------------------------------------------------*/
/*!
\file inpar_io.H

\brief Input parameters for global IO section

\maintainer Matthias Mayr

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef INPAR_IO_H
#define INPAR_IO_H

#include "inpar_parameterlist_utils.H"

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
#endif  // #ifndef INPAR_IO_H
