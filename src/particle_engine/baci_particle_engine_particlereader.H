/*---------------------------------------------------------------------------*/
/*! \file

\brief functionality to read particles from file

\level 3


*/
/*---------------------------------------------------------------------------*/

#ifndef BACI_PARTICLE_ENGINE_PARTICLEREADER_H
#define BACI_PARTICLE_ENGINE_PARTICLEREADER_H

#include "baci_io_inputreader.H"
#include "baci_particle_engine_typedefs.H"

namespace DRT::INPUT
{
  class ParticleReader
  {
   public:
    //! construct a reader that reads a given section
    ParticleReader(const DatFileReader& reader, std::string sectionname);

    //! do the actual reading of particles
    void Read(std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& particles);

   private:
    //! the main dat file reader
    const DatFileReader& reader_;

    //! my comm
    Teuchos::RCP<Epetra_Comm> comm_;

    //! my section to read
    std::string sectionname_;
  };

}  // namespace DRT::INPUT

#endif
