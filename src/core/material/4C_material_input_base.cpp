/*----------------------------------------------------------------------*/
/*! \file
\brief A material parameter container

\level 1

*/
/*---------------------------------------------------------------------*/


#include "4C_material_input_base.hpp"

FOUR_C_NAMESPACE_OPEN


CORE::MAT::PAR::Material::Material(const int id, const CORE::Materials::MaterialType type)
    : InputParameterContainer(), id_(id), type_(type)
{
}

FOUR_C_NAMESPACE_CLOSE
