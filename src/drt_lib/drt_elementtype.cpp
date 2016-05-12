/*----------------------------------------------------------------------*/
/*!
\file drt_elementtype.cpp

\brief Implementation Type definitions for elements

<pre>
\brief Implementation
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

