/*----------------------------------------------------------------------*/
/*!

\brief Type definitions for elements

\level 0

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#include "drt_elementtype.H"

DRT::ElementType::ElementType() : ParObjectType() {}

int DRT::ElementType::Initialize(DRT::Discretization& dis) { return 0; }
