/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container for input parameters for visualization of potential-based beam interactions

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_BEAMINTERACTION_POTENTIAL_RUNTIME_VISUALIZATION_OUTPUT_PARAMS_HPP
#define FOUR_C_BEAMINTERACTION_POTENTIAL_RUNTIME_VISUALIZATION_OUTPUT_PARAMS_HPP

#include "4C_config.hpp"

#include "4C_inpar_beampotential.hpp"
#include "4C_io_visualization_parameters.hpp"

FOUR_C_NAMESPACE_OPEN

namespace BEAMINTERACTION
{
  /*!
   *  */
  class BeamToBeamPotentialRuntimeOutputParams
  {
   public:
    //! constructor
    explicit BeamToBeamPotentialRuntimeOutputParams(double restart_time);

    //! destructor
    virtual ~BeamToBeamPotentialRuntimeOutputParams() = default;

    //! initialize with the stuff coming from input file
    void Init(const Teuchos::ParameterList& beam_contact_visualization_output_paramslist);

    //! setup member variables
    void Setup();

    /**
     * \brief Return the container holding the general output parameters
     */
    const CORE::IO::VisualizationParameters& get_visualization_parameters() const
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

    /// whether to write output for forces
    bool IsWriteForces() const
    {
      throw_error_if_not_init_and_setup();
      return output_forces_;
    };

    /// whether to write output for moments
    bool IsWriteMoments() const
    {
      throw_error_if_not_init_and_setup();
      return output_moments_;
    };

    /// whether to write forces/moments separately for each element pair
    bool is_write_forces_moments_per_element_pair() const
    {
      throw_error_if_not_init_and_setup();
      return write_force_moment_per_elepair_;
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
    CORE::IO::VisualizationParameters visualization_parameters_;

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

FOUR_C_NAMESPACE_CLOSE

#endif
