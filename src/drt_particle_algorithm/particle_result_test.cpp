/*---------------------------------------------------------------------------*/
/*!
\file particle_result_test.cpp

\brief particle result test for particle simulations

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
#include "particle_result_test.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_container_bundle.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_linedefinition.H"

#include "../drt_io/io_pstream.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::ResultTest::ResultTest(const Epetra_Comm& comm)
    : DRT::ResultTest("PARTICLE"), comm_(comm)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init result test                                           sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ResultTest::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup result test                                          sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ResultTest::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;
}

/*---------------------------------------------------------------------------*
 | test special quantity                                      sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ResultTest::TestSpecial(
    DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // extract string of particle type
  std::string type;
  res.ExtractString("TYPE", type);

  // extract reference position of particle
  std::vector<double> refpos(3, 0);
  res.ExtractDoubleVector("POS", refpos);

  // extract position tolerance of particle
  double positiontolerance = 0.0;
  res.ExtractDouble("POSTOLERANCE", positiontolerance);

  // get enum of particle types
  PARTICLEENGINE::TypeEnum particleType = PARTICLEENGINE::EnumFromTypeName(type);

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // get container of owned particles of current particle type
  PARTICLEENGINE::ParticleContainerShrdPtr container =
      particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

  // get number of particles stored in container
  int particlestored = container->ParticlesStored();

  // get dimension of particle position
  int posstatedim = PARTICLEENGINE::EnumToStateDim(PARTICLEENGINE::Position);

  // init minimum distance
  double minimumdistance = 1.0;
  int indexofparticlewithminimumdistance = -1;

  // declare pointer variables for current particle
  const double* pos;

  // only if particles are stored on this processor
  if (particlestored > 0)
  {
    // get pointer to particle states
    pos = container->GetPtrToParticleState(PARTICLEENGINE::Position, 0);

    // find particle
    for (int i = 0; i < particlestored; ++i)
    {
      // distance vector between reference position and current particle
      std::vector<double> currentdistancevector(3);
      for (int dim = 0; dim < posstatedim; ++dim)
        currentdistancevector[dim] = (pos[posstatedim * i + dim] - refpos[dim]);

      // get distance between reference position and current particle
      double currentdistance =
          std::sqrt(std::pow(currentdistancevector[0], 2) + std::pow(currentdistancevector[1], 2) +
                    std::pow(currentdistancevector[2], 2));

      // store current particle index and the distance
      if (currentdistance < minimumdistance)
      {
        minimumdistance = currentdistance;
        indexofparticlewithminimumdistance = i;
      }
    }
  }

  // get minimum distance over all processors
  double minimumdistanceallprocs = 0.0;
  comm_.MinAll(&minimumdistance, &minimumdistanceallprocs, 1);

  // decide which processor is the owner of the particle
  if (minimumdistance == minimumdistanceallprocs)
  {
    // prepare stringstream containing general information on particle position
    std::stringstream msghead;
    msghead << std::left << std::setw(9) << "PARTICLE"
            << ": " << std::left << std::setw(8) << "position"
            << "\t";

    if (minimumdistance < positiontolerance)
    {
      // position found
      IO::cout << msghead.str() << "\t is CORRECT"
               << ", abs(diff)=" << std::setw(24) << std::setprecision(17) << std::scientific
               << minimumdistance << " <" << std::setw(24) << std::setprecision(17)
               << std::scientific << positiontolerance << "\n";

      // get result
      std::string quantity;
      res.ExtractString("QUANTITY", quantity);

      // init actual result
      double actresult = 0.0;

      // component of result
      int dim = 0;

      // declare enum of particle state
      PARTICLEENGINE::StateEnum particleState;

      // velocity
      if (quantity == "velx" or quantity == "vely" or quantity == "velz")
      {
        // get enum of particle state
        particleState = PARTICLEENGINE::Velocity;

        // get component of result
        if (quantity == "velx")
          dim = 0;
        else if (quantity == "vely")
          dim = 1;
        else if (quantity == "velz")
          dim = 2;
      }
      // acceleration
      else if (quantity == "accx" or quantity == "accy" or quantity == "accz")
      {
        // get enum of particle state
        particleState = PARTICLEENGINE::Acceleration;

        // get component of result
        if (quantity == "accx")
          dim = 0;
        else if (quantity == "accy")
          dim = 1;
        else if (quantity == "accz")
          dim = 2;
      }
      // radius
      else if (quantity == "radius")
      {
        // get enum of particle state
        particleState = PARTICLEENGINE::Radius;

        // get component of result
        dim = 0;
      }
      // density
      else if (quantity == "density")
      {
        // get enum of particle state
        particleState = PARTICLEENGINE::Density;

        // get component of result
        dim = 0;
      }
      // pressure
      else if (quantity == "pressure")
      {
        // get enum of particle state
        particleState = PARTICLEENGINE::Pressure;

        // get component of result
        dim = 0;
      }
      else
        dserror("result check failed with unknown quantity '%s'!", quantity.c_str());

      // get pointer to particle state
      const double* state = container->GetPtrToParticleState(particleState, 0);

      // get dimension of particle state
      int statedim = PARTICLEENGINE::EnumToStateDim(particleState);

      // get actual result
      actresult = state[statedim * indexofparticlewithminimumdistance + dim];

      // compare values
      const int err = CompareValues(actresult, "SPECIAL", res);
      nerr += err;
      test_count++;
    }
    else
    {
      // position not found
      IO::cout << msghead.str() << "\t is WRONG"
               << ", abs(diff)=" << std::setw(24) << std::setprecision(17) << std::scientific
               << minimumdistance << " >" << std::setw(24) << std::setprecision(17)
               << std::scientific << positiontolerance << "\n"
               << "--> actposition: x=" << std::setw(24) << std::setprecision(17) << std::scientific
               << pos[posstatedim * indexofparticlewithminimumdistance + 0]
               << " y=" << std::setw(24) << std::setprecision(17) << std::scientific
               << pos[posstatedim * indexofparticlewithminimumdistance + 1]
               << " z=" << std::setw(24) << std::setprecision(17) << std::scientific
               << pos[posstatedim * indexofparticlewithminimumdistance + 2] << "\n"
               << "--> refposition: x=" << std::setw(24) << std::setprecision(17) << std::scientific
               << refpos[0] << " y=" << std::setw(24) << std::setprecision(17) << std::scientific
               << refpos[1] << " z=" << std::setw(24) << std::setprecision(17) << std::scientific
               << refpos[2] << "\n";

      // increase counter
      nerr++;
      test_count++;
    }
  }
}
