/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_phase_change.cpp

\brief phase change handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_phase_change.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"
#include "../drt_particle_engine/particle_object.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHPhaseChangeBase::SPHPhaseChangeBase()
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init phase change handler                                  sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHPhaseChangeBase::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup phase change handler                                 sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHPhaseChangeBase::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();
}

/*---------------------------------------------------------------------------*
 | write restart of phase change handler                      sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHPhaseChangeBase::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of phase change handler                       sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHPhaseChangeBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}
