/*---------------------------------------------------------------------------*/
/*!
\file particle_material_base.cpp

\brief particle material base

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_material_base.H"

#include "../drt_mat/matpar_bundle.H"

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
