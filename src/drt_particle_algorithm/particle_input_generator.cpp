/*---------------------------------------------------------------------------*/
/*!

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
   * note: take care of setting a global id (that is unique and not used by particles potentially
   * read in from the input file)
   *
   * add a generated particle to the vector particlesgenerated using the function
   * AddGeneratedParticle()
   */
}

/*---------------------------------------------------------------------------*
 | add generated particle                                     sfuchs 09/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::InputGenerator::AddGeneratedParticle(const std::vector<double>& position,
    const PARTICLEENGINE::TypeEnum particletype, int globalid,
    std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& particlesgenerated) const
{
  // safety check
  if (position.size() != 3)
    dserror("particle can not be generated since position vector needs three entries!");

  // allocate memory to hold particle states
  PARTICLEENGINE::ParticleStates particlestates;
  particlestates.assign((PARTICLEENGINE::Position + 1), std::vector<double>(0));

  // set position state
  particlestates[PARTICLEENGINE::Position] = position;

  // construct and store generated particle object
  particlesgenerated.emplace_back(
      std::make_shared<PARTICLEENGINE::ParticleObject>(particletype, globalid, particlestates));
}
