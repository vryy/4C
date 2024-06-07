/*----------------------------------------------------------------------*/
/*! \file

\brief Type definitions for elements

\level 0


*/
/*----------------------------------------------------------------------*/

#include "4C_fem_general_elementtype.hpp"

FOUR_C_NAMESPACE_OPEN

Core::Elements::ElementType::ElementType() : ParObjectType() {}

int Core::Elements::ElementType::Initialize(Discret::Discretization& dis) { return 0; }

FOUR_C_NAMESPACE_CLOSE
