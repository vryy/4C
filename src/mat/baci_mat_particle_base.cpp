/*---------------------------------------------------------------------------*/
/*! \file
\brief particle material base

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "baci_mat_particle_base.H"

#include "baci_mat_par_bundle.H"

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
MAT::PAR::ParticleMaterialBase::ParticleMaterialBase(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      initRadius_(matdata->GetDouble("INITRADIUS")),
      initDensity_(matdata->GetDouble("INITDENSITY"))
{
  // empty constructor
}

BACI_NAMESPACE_CLOSE
