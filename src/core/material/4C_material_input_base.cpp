/*----------------------------------------------------------------------*/
/*! \file
\brief A material parameter container

\level 1

*/
/*---------------------------------------------------------------------*/


#include "4C_material_input_base.hpp"

FOUR_C_NAMESPACE_OPEN


CORE::MAT::PAR::Material::Material(
    const int id, const CORE::Materials::MaterialType type, const std::string name)
    : InputParameterContainer(), id_(id), type_(type), name_(name), params_(Teuchos::null)
{
}

FOUR_C_NAMESPACE_CLOSE
