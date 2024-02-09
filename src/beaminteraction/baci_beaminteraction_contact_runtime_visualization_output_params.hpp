/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container holding all input parameters for vtk-based visualization of beam contact

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BACI_BEAMINTERACTION_CONTACT_RUNTIME_VISUALIZATION_OUTPUT_PARAMS_HPP
#define BACI_BEAMINTERACTION_CONTACT_RUNTIME_VISUALIZATION_OUTPUT_PARAMS_HPP

#include "baci_config.hpp"

#include "baci_inpar_beamcontact.hpp"
#include "baci_io_visualization_parameters.hpp"

BACI_NAMESPACE_OPEN

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
    const IO::VisualizationParameters& GetVisualizationParameters() const
    {
      return visualization_parameters_;
    }

    /// output interval regarding steps: write output every INTERVAL_STEPS steps
    int OutputIntervalInSteps() const
    {
      ThrowErrorIfNotInitAndSetup();
      return output_interval_steps_;
    };

    /// whether to write output in every iteration of the nonlinear solver
    bool OutputEveryIteration() const
    {
      ThrowErrorIfNotInitAndSetup();
      return output_every_iteration_;
    };

    /// whether to write output for contact forces
    bool IsWriteContactForces() const
    {
      ThrowErrorIfNotInitAndSetup();
      return output_forces_;
    };

    /// whether to write output for gaps
    bool IsWriteGaps() const
    {
      ThrowErrorIfNotInitAndSetup();
      return output_gaps_;
    };


   private:
    //! returns the isinit_ flag
    inline const bool& IsInit() const { return isinit_; };

    //! returns the issetup_ flag
    inline const bool& IsSetup() const { return issetup_; };

    //! asserts the init and setup status
    void ThrowErrorIfNotInitAndSetup() const;

    //! asserts the init status
    void ThrowErrorIfNotInit() const;


   private:
    bool isinit_;

    bool issetup_;

    //! General visualization parameters
    IO::VisualizationParameters visualization_parameters_;

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

BACI_NAMESPACE_CLOSE

#endif
