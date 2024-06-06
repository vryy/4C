/*---------------------------------------------------------------------------*/
/*! \file

\brief functionality to read particles from file

\level 3


*/
/*---------------------------------------------------------------------------*/

#ifndef FOUR_C_PARTICLE_ENGINE_PARTICLEREADER_HPP
#define FOUR_C_PARTICLE_ENGINE_PARTICLEREADER_HPP

#include "4C_config.hpp"

#include "4C_io_inputreader.hpp"
#include "4C_particle_engine_typedefs.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Input
{
  class ParticleReader
  {
   public:
    //! construct a reader that reads a given section
    ParticleReader(const Core::IO::DatFileReader& reader, std::string sectionname);

    //! do the actual reading of particles
    void Read(std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& particles);

   private:
    //! the main dat file reader
    const Core::IO::DatFileReader& reader_;

    //! my comm
    Teuchos::RCP<Epetra_Comm> comm_;

    //! my section to read
    std::string sectionname_;
  };

}  // namespace Input

FOUR_C_NAMESPACE_CLOSE

#endif
