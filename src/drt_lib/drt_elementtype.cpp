/*----------------------------------------------------------------------*/
/*!
\file drt_elementtype.cpp

\brief Type definitions for elements

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "drt_elementtype.H"

DRT::ElementType::ElementType()
  : ParObjectType()
{

}

int DRT::ElementType::Initialize(DRT::Discretization& dis)
{
  return 0;
}
