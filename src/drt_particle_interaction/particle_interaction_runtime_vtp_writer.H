/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for particle interaction in vtk/vtp format at runtime

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_RUNTIME_VTP_WRITER_H
#define PARTICLE_INTERACTION_RUNTIME_VTP_WRITER_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
#include "../drt_lib/drt_dserror.H"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>
#include <unordered_map>
#include <memory>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
class RuntimeVtpWriter;

namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class InteractionWriter
  {
   public:
    //! constructor
    explicit InteractionWriter(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~InteractionWriter();

    //! init interaction writer
    void Init();

    //! setup interaction writer
    void Setup();

    //! write restart of interaction writer
    void WriteRestart(const int step, const double time) const;

    //! read restart of interaction writer
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! register specific runtime vtp writer
    void RegisterSpecificRuntimeVtpWriter(const std::string& fieldname);

    //! set current write result flag
    void SetCurrentWriteResultFlag(bool writeresultsthisstep)
    {
      writeresultsthisstep_ = writeresultsthisstep;
    };

    //! get current write result flag
    inline bool GetCurrentWriteResultFlag() const { return writeresultsthisstep_; }

    //! get specific runtime vtp writer
    inline RuntimeVtpWriter* GetSpecificRuntimeVtpWriter(const std::string& fieldname)
    {
#ifdef DEBUG
      if (not runtime_vtpwriters_.count(fieldname))
        dserror("no runtime vtp writer for field '%s' stored!", fieldname.c_str());
#endif

      return runtime_vtpwriters_[fieldname].get();
    }

    // write particle interaction runtime output
    void WriteParticleInteractionRuntimeOutput(const int step, const double time) const;

   private:
    //! communication
    const Epetra_Comm& comm_;

    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! setup time of runtime vtp writer
    double setuptime_;

    //! result control flag
    bool writeresultsthisstep_;

    //! holds all vtp writer objects
    std::unordered_map<std::string, std::shared_ptr<RuntimeVtpWriter>> runtime_vtpwriters_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
