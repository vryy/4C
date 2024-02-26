/*---------------------------------------------------------------------------*/
/*! \file
\brief particle material base

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "baci_mat_particle_base.hpp"

#include "baci_mat_par_bundle.hpp"

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
MAT::PAR::ParticleMaterialBase::ParticleMaterialBase(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      initRadius_(*matdata->Get<double>("INITRADIUS")),
      initDensity_(*matdata->Get<double>("INITDENSITY"))
{
  // empty constructor
}

BACI_NAMESPACE_CLOSE
