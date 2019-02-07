/*---------------------------------------------------------------------------*/
/*!
\file drt_particlereader.cpp

\brief functionality to read particles from file

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
#include "drt_particlereader.H"

#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"
#include "../drt_particle_engine/particle_object.H"

#include "../drt_io/io_pstream.H"

#include <Epetra_Time.h>


/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
DRT::INPUT::ParticleReader::ParticleReader(
    const DRT::INPUT::DatFileReader& reader, std::string sectionname)
    : reader_(reader), comm_(reader.Comm()), sectionname_(sectionname)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | do the actual reading of particles                         sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void DRT::INPUT::ParticleReader::Read(std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& particles)
{
  const int myrank = comm_->MyPID();
  const int numproc = comm_->NumProc();
  std::string inputfile_name = reader_.MyInputfileName();

  int numparticles = reader_.ExcludedSectionLength(sectionname_);

  // proceed only if particles are given in .dat-file
  if (numparticles > 0)
  {
    Epetra_Time time(*comm_);

    if (!myrank && !reader_.MyOutputFlag()) IO::cout << "Read and create particles\n" << IO::flush;

    // read in the particles block-wise:
    // EITHER one block per processor so that the number of blocks is numproc
    // OR number of blocks is numparticles if less particles than procs are read in
    // determine a rough blocksize
    int nblock = std::min(numproc, numparticles);
    int bsize = std::max(numparticles / nblock, 1);

    // an upper limit for bsize
    int maxblocksize = 200000;

    if (bsize > maxblocksize)
    {
      // without an additional increase of nblock by 1 the last block size
      // could reach a maximum value of (2*maxblocksize)-1, potentially
      // violating the intended upper limit!
      nblock = 1 + static_cast<int>(numparticles / maxblocksize);
      bsize = maxblocksize;
    }

    // open input file at the right position
    // note that stream is valid on proc 0 only!
    std::ifstream file;
    if (myrank == 0)
    {
      file.open(inputfile_name.c_str());
      file.seekg(reader_.ExcludedSectionPosition(sectionname_));
    }

    if (!myrank && !reader_.MyOutputFlag())
    {
      printf("numparticle %d nblock %d bsize %d\n", numparticles, nblock, bsize);
      fflush(stdout);
    }

    std::string tmp;
    std::string tmp2;
    int filecount = 0;

    int particlecounter = 0;

    PARTICLEENGINE::TypeEnum particleType;
    int globalid(0);
    PARTICLEENGINE::ParticleStates particlestates;

    // note that the last block is special....
    for (int block = 0; block < nblock; ++block)
    {
      double t1 = time.ElapsedTime();

      if (0 == myrank)
      {
#if defined(HAVE_PARMETIS)
        if (!reader_.MyOutputFlag()) printf("block %d ", block);
#endif

        int bcount = 0;
        for (; file; ++filecount)
        {
          file >> tmp;

          if (tmp == "TYPE")
          {
            std::string type;
            std::vector<double> pos(3);

            // read in particle type and position
            file >> type >> tmp2 >> pos[0] >> pos[1] >> pos[2];

            if (tmp2 != "POS") dserror("expected particle position!");

            // get enum of particle types
            particleType = PARTICLEENGINE::EnumFromTypeName(type);

            // insert position state
            particlestates.clear();
            particlestates.insert(std::make_pair(PARTICLEENGINE::Position, pos));

            // set global id
            globalid = particlecounter;

            // construct and init particleobject
            PARTICLEENGINE::ParticleObjShrdPtr particleobject =
                std::make_shared<PARTICLEENGINE::ParticleObject>();
            particleobject->Init(particleType, globalid, particlestates);

            // store read in particles
            particles.push_back(particleobject);

            ++particlecounter;

            ++bcount;
            if (block != nblock - 1)  // last block takes all the rest
              if (bcount == bsize)    // block is full
              {
                ++filecount;
                break;
              }
          }
          else if (tmp.find("--") == 0)
            break;
          else
            dserror("unexpected word '%s'", tmp.c_str());
        }
      }

      double t2 = time.ElapsedTime();
      if (!myrank && !reader_.MyOutputFlag())
      {
        printf("reading %10.5e secs\n", t2 - t1);
        fflush(stdout);
      }
    }

    if (!myrank && !reader_.MyOutputFlag())
      printf("in............................................. %10.5e secs\n", time.ElapsedTime());
  }
}
