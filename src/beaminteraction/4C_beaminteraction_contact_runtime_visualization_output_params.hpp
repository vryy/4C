/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container holding all input parameters for vtk-based visualization of beam contact

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_BEAMINTERACTION_CONTACT_RUNTIME_VISUALIZATION_OUTPUT_PARAMS_HPP
#define FOUR_C_BEAMINTERACTION_CONTACT_RUNTIME_VISUALIZATION_OUTPUT_PARAMS_HPP

#include "4C_config.hpp"

#include "4C_inpar_beamcontact.hpp"
#include "4C_io_visualization_parameters.hpp"

FOUR_C_NAMESPACE_OPEN

namespace BEAMINTERACTION
{
  class BeamContactRuntimeVisualizationOutputParams
  {
   public:
    //! constructor
    explicit BeamContactRuntimeVisualizationOutputParams(double restart_time);

    //! destructor
    virtual ~BeamContactRuntimeVisualizationOutputParams() = default;

    //! initialize with the stuff coming from input file
    void Init();

    //! setup member variables
    void Setup();

    /**
     * \brief Return the container holding the general output parameters
     */
    const Core::IO::VisualizationParameters& get_visualization_parameters() const
    {
      return visualization_parameters_;
    }

    /// output interval regarding steps: write output every INTERVAL_STEPS steps
    int output_interval_in_steps() const
    {
      throw_error_if_not_init_and_setup();
      return output_interval_steps_;
    };

    /// whether to write output in every iteration of the nonlinear solver
    bool output_every_iteration() const
    {
      throw_error_if_not_init_and_setup();
      return output_every_iteration_;
    };

    /// whether to write output for contact forces
    bool is_write_contact_forces() const
    {
      throw_error_if_not_init_and_setup();
      return output_forces_;
    };

    /// whether to write output for gaps
    bool IsWriteGaps() const
    {
      throw_error_if_not_init_and_setup();
      return output_gaps_;
    };


   private:
    //! returns the isinit_ flag
    inline const bool& is_init() const { return isinit_; };

    //! returns the issetup_ flag
    inline const bool& is_setup() const { return issetup_; };

    //! asserts the init and setup status
    void throw_error_if_not_init_and_setup() const;

    //! asserts the init status
    void throw_error_if_not_init() const;


   private:
    bool isinit_;

    bool issetup_;

    //! General visualization parameters
    Core::IO::VisualizationParameters visualization_parameters_;

    /// output interval regarding steps: write output every INTERVAL_STEPS steps
    int output_interval_steps_;

    /// whether to write output in every iteration of the nonlinear solver
    bool output_every_iteration_;

    /// whether to write forces
    bool output_forces_;

    /// whether to write gaps
    bool output_gaps_;
  };

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
