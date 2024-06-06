/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container holding all input parameters relevant for potential based beam interactions

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_BEAMINTERACTION_POTENTIAL_PARAMS_HPP
#define FOUR_C_BEAMINTERACTION_POTENTIAL_PARAMS_HPP

#include "4C_config.hpp"

#include "4C_inpar_beampotential.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace BEAMINTERACTION
{
  class BeamToBeamPotentialRuntimeOutputParams;

  /*!
   *  */
  class BeamPotentialParams
  {
   public:
    //! constructor
    BeamPotentialParams();

    //! destructor
    virtual ~BeamPotentialParams() = default;

    //! initialize with the stuff coming from input file
    void Init(double restart_time);

    //! setup member variables
    void Setup();

    //! returns the isinit_ flag
    inline bool is_init() const { return isinit_; }

    //! returns the issetup_ flag
    inline bool is_setup() const { return issetup_; }

    //! asserts the init and setup status
    void throw_error_if_not_init_and_setup() const;

    //! asserts the init status
    void throw_error_if_not_init() const;

    inline std::vector<double> const& potential_law_exponents() const
    {
      throw_error_if_not_init_and_setup();
      return *pot_law_exponents_;
    }

    inline std::vector<double> const& potential_law_prefactors() const
    {
      throw_error_if_not_init_and_setup();
      return *pot_law_prefactors_;
    }

    inline enum Inpar::BEAMPOTENTIAL::BeamPotentialType PotentialType() const
    {
      throw_error_if_not_init_and_setup();
      return potential_type_;
    }

    inline enum Inpar::BEAMPOTENTIAL::BeamPotentialStrategy Strategy() const
    {
      throw_error_if_not_init_and_setup();
      return strategy_;
    }

    inline double CutoffRadius() const
    {
      throw_error_if_not_init_and_setup();
      return cutoff_radius_;
    }

    inline enum Inpar::BEAMPOTENTIAL::BeamPotentialRegularizationType RegularizationType() const
    {
      throw_error_if_not_init_and_setup();
      return regularization_type_;
    }

    inline double regularization_separation() const
    {
      throw_error_if_not_init_and_setup();
      return regularization_separation_;
    }

    inline int number_integration_segments() const
    {
      throw_error_if_not_init_and_setup();
      return num_integration_segments_;
    }

    inline int NumberGaussPoints() const
    {
      throw_error_if_not_init_and_setup();
      return num_gp_s_;
    }

    inline bool UseFAD() const
    {
      throw_error_if_not_init_and_setup();
      return use_fad_;
    }

    inline enum Inpar::BEAMPOTENTIAL::MasterSlaveChoice ChoiceMasterSlave() const
    {
      throw_error_if_not_init_and_setup();
      return choice_master_slave_;
    }

    //! whether to write visualization output for beam contact
    inline bool RuntimeOutput() const
    {
      throw_error_if_not_init_and_setup();
      return visualization_output_;
    }

    //! get the data container for parameters regarding visualization output
    inline Teuchos::RCP<const BEAMINTERACTION::BeamToBeamPotentialRuntimeOutputParams>
    get_beam_potential_visualization_output_params() const
    {
      throw_error_if_not_init_and_setup();
      return params_runtime_visualization_output_btb_potential_;
    }

   private:
    bool isinit_;

    bool issetup_;

    //! exponents of the summands of a potential law in form of a power law
    // Todo maybe change to integer?
    Teuchos::RCP<std::vector<double>> pot_law_exponents_;

    //! prefactors of the summands of a potential law in form of a power law
    Teuchos::RCP<std::vector<double>> pot_law_prefactors_;

    //! type of applied potential (volume, surface)
    enum Inpar::BEAMPOTENTIAL::BeamPotentialType potential_type_;

    //! strategy to evaluate interaction potential
    enum Inpar::BEAMPOTENTIAL::BeamPotentialStrategy strategy_;

    //! neglect all contributions at separation larger than this cutoff radius
    double cutoff_radius_;

    //! type of regularization to use for force law at separations below specified separation
    enum Inpar::BEAMPOTENTIAL::BeamPotentialRegularizationType regularization_type_;

    //! use specified regularization type for separations smaller than this value
    double regularization_separation_;

    //! number of integration segments to be used per beam element
    int num_integration_segments_;

    //! number of Gauss points to be used per integration segment
    int num_gp_s_;

    //! use automatic differentiation via FAD
    bool use_fad_;

    //! rule how to assign the role of master and slave to beam elements (if applicable)
    enum Inpar::BEAMPOTENTIAL::MasterSlaveChoice choice_master_slave_;

    //! whether to write visualization output at runtime
    bool visualization_output_;

    //! data container for input parameters related to visualization output of beam contact at
    //! runtime
    Teuchos::RCP<BEAMINTERACTION::BeamToBeamPotentialRuntimeOutputParams>
        params_runtime_visualization_output_btb_potential_;
  };

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
