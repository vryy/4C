/*---------------------------------------------------------------------------*/
/*! \file

\brief functionality to read particles from file

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
#ifndef BACI_LIB_PARTICLEREADER_H
#define BACI_LIB_PARTICLEREADER_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
#include "baci_inpar.H"
#include "baci_lib_inputreader.H"
#include "baci_particle_engine_typedefs.H"

/*---------------------------------------------------------------------------*
 | helper class to read all particles of different types      sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
namespace DRT
{
  namespace INPUT
  {
    class ParticleReader
    {
     public:
      //! construct a reader that reads a given section
      ParticleReader(const DatFileReader& reader, std::string sectionname);

      //! destructor
      virtual ~ParticleReader() { ; };

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

  }  // namespace INPUT

}  // namespace DRT

/*---------------------------------------------------------------------------*/
#endif  // LIB_PARTICLEREADER_H
