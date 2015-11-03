/*--------------------------------------------------------------------------*/
/*!
\file inpar_reynolds.H

\brief Reynolds dynamic parameters

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_INPAR_INPAR_REYNOLDS_H_
#define SRC_DRT_INPAR_INPAR_REYNOLDS_H_

#include "inpar_parameterlist_utils.H"

namespace INPAR
{

  /// Input parameters for REYNOLDS
  namespace REYNOLDS
  {

    /// set the reynolds parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace REYNOLDS
}  // namespace INPAR

#endif /* SRC_DRT_INPAR_INPAR_REYNOLDS_H_ */
