/*----------------------------------------------------------------------------*/
/*!
\file nln_utils.H

\brief Collection of utilities for nonlinear solver objects

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_UTILS_H
#define NLN_UTILS_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------------*/
/* forward declarations */

namespace Teuchos
{
  class ParameterList;
} // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of methods */
namespace NLNSOL
{
namespace UTILS
{
  /*! \brief Create Teuchos::ParameterList from XML-File
   *
   *  \author mayr.mt \date 06/2014
   */
  Teuchos::RCP<Teuchos::ParameterList> CreateParamListFromXML();
} // namespace UTILS
} // namespace NLNSOL

#endif /* NLN_UTILS_H */
