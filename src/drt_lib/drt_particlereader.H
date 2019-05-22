/*---------------------------------------------------------------------------*/
/*!

\brief functionality to read particles from file

\level 3

\maintainer  Sebastian Fuchs

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
#ifndef DRT_PARTICLEREADER_H
#define DRT_PARTICLEREADER_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
#include "drt_inputreader.H"
#include "../drt_inpar/inpar.H"

#include "../drt_particle_engine/particle_typedefs.H"

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
#endif  // DRT_PARTICLEREADER_H
