/*----------------------------------------------------------------------*/
/*! \file

\brief Type definitions for elements

\level 0


*/
/*----------------------------------------------------------------------*/

#include "baci_lib_elementtype.H"

BACI_NAMESPACE_OPEN

DRT::ElementType::ElementType() : ParObjectType() {}

int DRT::ElementType::Initialize(DRT::Discretization& dis) { return 0; }

BACI_NAMESPACE_CLOSE
