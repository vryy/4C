/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for poro problems

\level 3

*------------------------------------------------------------------------------------------------*/

#ifndef INPAR_PORO_H
#define INPAR_PORO_H

namespace INPAR::PORO
{
  //! poro element implementation type
  enum class PoroType
  {
    porotype_undefined,
    porotype_pressure_based,
    porotype_pressure_velocity_based
  };

}  // namespace INPAR::PORO

#endif