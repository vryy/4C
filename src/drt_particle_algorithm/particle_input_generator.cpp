/*---------------------------------------------------------------------------*/
/*!
\file particle_input_generator.cpp

\brief particle input generator for particle simulations

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_input_generator.H"

#include "../drt_particle_engine/particle_object.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::InputGenerator::InputGenerator(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : myrank_(comm.MyPID()), params_(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init input generator                                       sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::InputGenerator::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | generate particles                                         sfuchs 09/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::InputGenerator::GenerateParticles(
    std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& particlesgenerated) const
{
  /*
   * generate initial particles in addition to particles read in from input file
   *
   * attention: either generate particles only on one processor or be sure that your particles are
   * not generated twice (or even more) on different processors
   *
   * note: think about reserving (not resizing!) the vector particlesgenerated in advance if you
   * have a rough estimate of the total number of particles being generated
   *
   * add a generated particle to the vector particlesgenerated using the function
   * AddGeneratedParticle()
   */
}

/*---------------------------------------------------------------------------*
 | add generated particle                                     sfuchs 09/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::InputGenerator::AddGeneratedParticle(const std::vector<double>& position,
    const PARTICLEENGINE::TypeEnum particletype,
    std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& particlesgenerated) const
{
  // safety check
  if (position.size() != 3)
    dserror("particle can not be generated since position vector needs three entries!");

  // fill position into particle states map
  std::map<PARTICLEENGINE::StateEnum, std::vector<double>> particlestates;
  particlestates.insert(std::make_pair(PARTICLEENGINE::Position, position));

  // create and init new particle object
  PARTICLEENGINE::ParticleObjShrdPtr particleobject =
      std::make_shared<PARTICLEENGINE::ParticleObject>();
  particleobject->Init(particletype, particlestates);

  // append generated particle
  particlesgenerated.push_back(particleobject);
}
