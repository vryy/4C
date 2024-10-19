// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_PARTICLE_ENGINE_PARTICLEREADER_HPP
#define FOUR_C_PARTICLE_ENGINE_PARTICLEREADER_HPP

#include "4C_config.hpp"

#include "4C_particle_engine_typedefs.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  class DatFileReader;
}

namespace Input
{
  class ParticleReader
  {
   public:
    //! construct a reader that reads a given section
    ParticleReader(Core::IO::DatFileReader& reader, std::string sectionname);

    //! do the actual reading of particles
    void read(std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& particles);

   private:
    //! the main dat file reader
    Core::IO::DatFileReader& reader_;

    //! my comm
    Teuchos::RCP<Epetra_Comm> comm_;

    //! my section to read
    std::string sectionname_;
  };

}  // namespace Input

FOUR_C_NAMESPACE_CLOSE

#endif
