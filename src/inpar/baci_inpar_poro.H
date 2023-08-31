/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for poro problems

\level 3

*------------------------------------------------------------------------------------------------*/

#ifndef BACI_INPAR_PORO_H
#define BACI_INPAR_PORO_H

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

#endif