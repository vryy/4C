/*---------------------------------------------------------------------------*/
/*! \file

\brief functionality to read particles from file

\level 3

\maintainer  Sebastian Fuchs

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
    if (!myrank)
    {
      file.open(inputfile_name.c_str());
      file.seekg(reader_.ExcludedSectionPosition(sectionname_));
    }

    std::string line;
    int filecount = 0;
    bool endofsection = false;

    if (!myrank && !reader_.MyOutputFlag())
    {
      printf("numparticle %d nblock %d bsize %d\n", numparticles, nblock, bsize);
      fflush(stdout);
    }

    // note that the last block is special....
    for (int block = 0; block < nblock; ++block)
    {
      double t1 = time.ElapsedTime();

      if (!myrank and !endofsection)
      {
        int bcount = 0;

        for (; getline(file, line); ++filecount)
        {
          if (line.find("--") == 0)
          {
            endofsection = true;
            break;
          }
          else
          {
            std::istringstream linestream;
            linestream.str(line);

            PARTICLEENGINE::TypeEnum particletype;
            int globalid(-1);
            PARTICLEENGINE::ParticleStates particlestates;

            std::string typelabel;
            std::string type;

            std::string poslabel;
            std::vector<double> pos(3);

            // read in particle type and position
            linestream >> typelabel >> type >> poslabel >> pos[0] >> pos[1] >> pos[2];

            if (typelabel != "TYPE") dserror("expected particle type label 'TYPE'!");

            if (poslabel != "POS") dserror("expected particle position label 'POS'!");

            // get enum of particle type
            particletype = PARTICLEENGINE::EnumFromTypeName(type);

            // allocate memory to hold particle position state
            particlestates.resize(PARTICLEENGINE::Position + 1);

            // set position state
            particlestates[PARTICLEENGINE::Position] = pos;

            // optional particle radius
            {
              std::string radlabel;
              std::vector<double> rad(1);

              // read in optional particle radius
              linestream >> radlabel;
              if (linestream)
              {
                if (radlabel != "RAD") dserror("expected particle radius label 'RAD'!");

                linestream >> rad[0];

                if (not linestream) dserror("expected radius if radius label 'RAD' is set!");
              }

              // allocate memory to hold particle radius state
              if (particlestates.size() < (PARTICLEENGINE::Radius + 1))
                particlestates.resize(PARTICLEENGINE::Radius + 1);

              // set radius state
              particlestates[PARTICLEENGINE::Radius] = rad;
            }

            // construct and store read in particle object
            particles.emplace_back(std::make_shared<PARTICLEENGINE::ParticleObject>(
                particletype, globalid, particlestates));

            ++bcount;
            if (block != nblock - 1)  // last block takes all the rest
              if (bcount == bsize)    // block is full
              {
                ++filecount;
                break;
              }
          }
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
