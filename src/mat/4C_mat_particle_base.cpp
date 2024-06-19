/*---------------------------------------------------------------------------*/
/*! \file
\brief particle material base

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "4C_mat_particle_base.hpp"

#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
Mat::PAR::ParticleMaterialBase::ParticleMaterialBase(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      initRadius_(matdata.parameters.get<double>("INITRADIUS")),
      initDensity_(matdata.parameters.get<double>("INITDENSITY"))
{
  // empty constructor
}

FOUR_C_NAMESPACE_CLOSE
