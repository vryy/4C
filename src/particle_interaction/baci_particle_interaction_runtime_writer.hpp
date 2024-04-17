/*---------------------------------------------------------------------------*/
/*! \file
\brief write output for particle interaction at runtime
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_INTERACTION_RUNTIME_WRITER_HPP
#define FOUR_C_PARTICLE_INTERACTION_RUNTIME_WRITER_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_io_visualization_manager.hpp"
#include "baci_utils_exceptions.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

#include <memory>
#include <unordered_map>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
  class RuntimeCsvWriter;
}  // namespace IO

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class InteractionWriter final
  {
   public:
    //! constructor
    explicit InteractionWriter(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! init interaction writer
    void Init();

    //! setup interaction writer
    void Setup();

    //! read restart of interaction writer
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! register specific runtime output writer
    void RegisterSpecificRuntimeOutputWriter(const std::string& fieldname);

    //! register specific runtime csv writer
    void RegisterSpecificRuntimeCsvWriter(const std::string& fieldname);

    //! set current write result flag
    void SetCurrentWriteResultFlag(bool writeresultsthisstep)
    {
      writeresultsthisstep_ = writeresultsthisstep;
    };

    //! get current write result flag
    inline bool GetCurrentWriteResultFlag() const { return writeresultsthisstep_; }

    //! get specific runtime output writer
    inline IO::VisualizationManager* GetSpecificRuntimeOutputWriter(const std::string& fieldname)
    {
#ifdef BACI_DEBUG
      if (not runtime_visualization_managers_.count(fieldname))
        dserror("no runtime output writer for field '%s' stored!", fieldname.c_str());
#endif

      return runtime_visualization_managers_[fieldname].get();
    }

    //! get specific runtime csv writer
    inline IO::RuntimeCsvWriter* GetSpecificRuntimeCsvWriter(const std::string& fieldname)
    {
#ifdef BACI_DEBUG
      if (not runtime_csvwriters_.count(fieldname))
        dserror("no runtime csv writer for field '%s' stored!", fieldname.c_str());
#endif

      return runtime_csvwriters_[fieldname].get();
    }

    // write particle interaction runtime output
    void WriteParticleInteractionRuntimeOutput(const int step, const double time) const;

   private:
    //! communication
    const Epetra_Comm& comm_;

    //! setup time of runtime output writer
    double setuptime_;

    //! result control flag
    bool writeresultsthisstep_;

    //! holds all visualization output writer objects
    std::unordered_map<std::string, std::shared_ptr<IO::VisualizationManager>>
        runtime_visualization_managers_;

    //! holds all csv writer objects
    std::unordered_map<std::string, std::shared_ptr<IO::RuntimeCsvWriter>> runtime_csvwriters_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
