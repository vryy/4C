/*!----------------------------------------------------------------------
\file drt_singletondestruction.cpp
\brief Implementation

<pre>
\brief Implementation
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "drt_singletondestruction.H"
#include "drt_globalproblem.H"

DRT::SingletonDestruction::SingletonDestruction()
{
  DRT::Problem::Instance()->Register(this);
}
