/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_std.H

\brief Setting of general scatra parameter for element evaluation

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_PARAMETER_STD_H
#define SCATRA_ELE_PARAMETER_STD_H

#include "scatra_ele_parameter.H"


namespace DRT
{
namespace ELEMENTS
{

class ScaTraEleParameterStd : public ScaTraEleParameter
{

public:

  virtual ~ScaTraEleParameterStd() {}

  /// Singleton access method
  static ScaTraEleParameterStd* Instance( bool create=true );

  /// called upon destruction
  virtual void Done();

private:

protected:

  /// protected Constructor since we are a Singleton.
  ScaTraEleParameterStd();

};

}
}

#endif
