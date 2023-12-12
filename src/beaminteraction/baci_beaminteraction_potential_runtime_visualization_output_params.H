/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container for input parameters for visualization of potential-based beam interactions

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BACI_BEAMINTERACTION_POTENTIAL_RUNTIME_VISUALIZATION_OUTPUT_PARAMS_H
#define BACI_BEAMINTERACTION_POTENTIAL_RUNTIME_VISUALIZATION_OUTPUT_PARAMS_H

#include "baci_inpar_beampotential.H"
#include "baci_io_visualization_parameters.H"

BACI_NAMESPACE_OPEN

namespace BEAMINTERACTION
{
  /*!
   *  */
  class BeamToBeamPotentialRuntimeOutputParams
  {
   public:
    //! constructor
    BeamToBeamPotentialRuntimeOutputParams();

    //! destructor
    virtual ~BeamToBeamPotentialRuntimeOutputParams() = default;

    //! initialize with the stuff coming from input file
    void Init(const Teuchos::ParameterList& beam_contact_visualization_output_paramslist);

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

    /// whether to write output for forces
    bool IsWriteForces() const
    {
      ThrowErrorIfNotInitAndSetup();
      return output_forces_;
    };

    /// whether to write output for moments
    bool IsWriteMoments() const
    {
      ThrowErrorIfNotInitAndSetup();
      return output_moments_;
    };

    /// whether to write forces/moments separately for each element pair
    bool IsWriteForcesMomentsPerElementPair() const
    {
      ThrowErrorIfNotInitAndSetup();
      return write_force_moment_per_elepair_;
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

    /// whether to write moments
    bool output_moments_;

    /// whether to write forces/moments separately for each element pair
    bool write_force_moment_per_elepair_;
  };

}  // namespace BEAMINTERACTION

BACI_NAMESPACE_CLOSE

#endif
