/*----------------------------------------------------------------------*/
/*! \file
\brief A material parameter container

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_material_parameter_base.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


Core::Mat::PAR::Parameter::Parameter(Data data) : data_(std::move(data)) {}


FOUR_C_NAMESPACE_CLOSE
