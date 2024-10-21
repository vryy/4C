#include "4C_fem_general_elementtype.hpp"

FOUR_C_NAMESPACE_OPEN

Core::Elements::ElementType::ElementType() : ParObjectType() {}

int Core::Elements::ElementType::initialize(Core::FE::Discretization& dis) { return 0; }

FOUR_C_NAMESPACE_CLOSE
