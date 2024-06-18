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
#include "4C_config.hpp"

#include "4C_io_visualization_manager.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

#include <memory>
#include <unordered_map>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace Core::IO
{
  class DiscretizationReader;
  class RuntimeCsvWriter;
}  // namespace Core::IO

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace ParticleInteraction
{
  class InteractionWriter final
  {
   public:
    //! constructor
    explicit InteractionWriter(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! init interaction writer
    void init();

    //! setup interaction writer
    void setup();

    //! read restart of interaction writer
    void read_restart(const std::shared_ptr<Core::IO::DiscretizationReader> reader);

    //! register specific runtime output writer
    void register_specific_runtime_output_writer(const std::string& fieldname);

    //! register specific runtime csv writer
    void register_specific_runtime_csv_writer(const std::string& fieldname);

    //! set current write result flag
    void set_current_write_result_flag(bool writeresultsthisstep)
    {
      writeresultsthisstep_ = writeresultsthisstep;
    };

    //! get current write result flag
    inline bool get_current_write_result_flag() const { return writeresultsthisstep_; }

    //! get specific runtime output writer
    inline Core::IO::VisualizationManager* get_specific_runtime_output_writer(
        const std::string& fieldname)
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (not runtime_visualization_managers_.count(fieldname))
        FOUR_C_THROW("no runtime output writer for field '%s' stored!", fieldname.c_str());
#endif

      return runtime_visualization_managers_[fieldname].get();
    }

    //! get specific runtime csv writer
    inline Core::IO::RuntimeCsvWriter* get_specific_runtime_csv_writer(const std::string& fieldname)
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (not runtime_csvwriters_.count(fieldname))
        FOUR_C_THROW("no runtime csv writer for field '%s' stored!", fieldname.c_str());
#endif

      return runtime_csvwriters_[fieldname].get();
    }

    // write particle interaction runtime output
    void write_particle_interaction_runtime_output(const int step, const double time) const;

   private:
    //! communication
    const Epetra_Comm& comm_;

    //! setup time of runtime output writer
    double setuptime_;

    //! result control flag
    bool writeresultsthisstep_;

    //! holds all visualization output writer objects
    std::unordered_map<std::string, std::shared_ptr<Core::IO::VisualizationManager>>
        runtime_visualization_managers_;

    //! holds all csv writer objects
    std::unordered_map<std::string, std::shared_ptr<Core::IO::RuntimeCsvWriter>>
        runtime_csvwriters_;
  };

}  // namespace ParticleInteraction

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
