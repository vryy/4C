/*----------------------------------------------------------------------*/
/*! \file

\brief Type definitions for elements

\level 0


*/
/*----------------------------------------------------------------------*/

#include "4C_lib_elementtype.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ElementType::ElementType() : ParObjectType() {}

int DRT::ElementType::Initialize(DRT::Discretization& dis) { return 0; }

FOUR_C_NAMESPACE_CLOSE
