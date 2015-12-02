/*--------------------------------------------------------------------------*/
/*!
\file inpar_lubrication.H

\brief Lubrication dynamic parameters

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_INPAR_INPAR_LUBRICATION_H_
#define SRC_DRT_INPAR_INPAR_LUBRICATION_H_

#include "inpar_parameterlist_utils.H"

namespace INPAR
{

  /// Input parameters for LUBRICATION
  namespace LUBRICATION
  {

  /// compute error compared to analytical solution
  enum CalcError
  {
    calcerror_no,
    calcerror_byfunction
  };

    /// set the lubrication parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace LUBRICATION
}  // namespace INPAR

#endif /* SRC_DRT_INPAR_INPAR_LUBRICATION_H_ */
