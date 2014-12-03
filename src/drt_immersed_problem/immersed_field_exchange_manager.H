/*!----------------------------------------------------------------------
\file immersed_field_exchange_manager.H

\brief manage access to and provide data in immersed problems

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/

#ifndef IMMERSED_FIELD_EXCHANGE_MANAGER_H
#define IMMERSED_FIELD_EXCHANGE_MANAGER_H

#include "../drt_lib/drt_singletondestruction.H"
#include <cstdlib>

namespace IMMERSED
{

class ImmersedFieldExchangeManager
{

public:

  virtual ~ImmersedFieldExchangeManager() {}

  /// Singleton access method
  static ImmersedFieldExchangeManager* Instance( bool create=true );

  /// called upon destruction
  virtual void Done();

private:

  // none

protected:

  /// protected Constructor since we are a Singleton.
  ImmersedFieldExchangeManager();

};

}


#endif
