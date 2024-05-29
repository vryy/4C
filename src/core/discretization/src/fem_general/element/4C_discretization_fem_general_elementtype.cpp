/*----------------------------------------------------------------------*/
/*! \file

\brief Type definitions for elements

\level 0


*/
/*----------------------------------------------------------------------*/

#include "4C_discretization_fem_general_elementtype.hpp"

FOUR_C_NAMESPACE_OPEN

CORE::Elements::ElementType::ElementType() : ParObjectType() {}

int CORE::Elements::ElementType::Initialize(DRT::Discretization& dis) { return 0; }

FOUR_C_NAMESPACE_CLOSE
