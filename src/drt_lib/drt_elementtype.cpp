/*----------------------------------------------------------------------*/
/*!
\file drt_elementtype.cpp

\brief Type definitions for elements

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
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

