/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for poro problems

\level 3

*------------------------------------------------------------------------------------------------*/

#ifndef BACI_INPAR_PORO_HPP
#define BACI_INPAR_PORO_HPP

#include "baci_config.hpp"

BACI_NAMESPACE_OPEN

namespace INPAR::PORO
{
  //! poro element implementation type
  enum class PoroType
  {
    undefined,
    pressure_based,
    pressure_velocity_based
  };

}  // namespace INPAR::PORO

BACI_NAMESPACE_CLOSE

#endif