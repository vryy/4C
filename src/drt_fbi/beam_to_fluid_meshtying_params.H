/*----------------------------------------------------------------------*/
/*! \file
\file beam_to_fluid_meshtying_params.H

\brief Data container holding all beam to fluid volume meshtying input parameters.

\level 3
\maintainer Nora Hagmeyer
*/
// End doxygen header.


#ifndef BEAM_TO_FLUID_MESHTYING_PARAMS_H
#define BEAM_TO_FLUID_MESHTYING_PARAMS_H


#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_inpar/inpar_fbi.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_beaminteraction/beam_contact_params.H"
namespace FBI
{
  /**
   * \brief Class for beam to fluid meshtying parameters.
   */
  class BeamToFluidMeshtyingParams : public BEAMINTERACTION::BeamContactParams  // todo fix this
  {
   public:
    /**
     * \brief Constructor.
     */
    BeamToFluidMeshtyingParams();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToFluidMeshtyingParams(){};

    /**
     * \brief Initialize with the stuff coming from input file.
     */
    void Init();

    /**
     * \brief Setup member variables.
     */
    void Setup();

    /**
     * \brief Returns the isinit_ flag.
     */
    inline const bool& IsInit() const { return isinit_; };

    /**
     * \brief Returns the issetup_ flag.
     */
    inline const bool& IsSetup() const { return issetup_; };

    /**
     * \brief Checks the init and setup status.
     */
    inline void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup()) dserror("Call Init() and Setup() first!");
    }

    /**
     * \brief Checks the init status.
     */
    inline void CheckInit() const
    {
      if (!IsInit()) dserror("Init() has not been called, yet!");
    }

    /**
     * \brief Returns the contact discretization method.
     */
    inline INPAR::FBI::BeamToFluidConstraintEnforcement GetConstraintEnforcement() const
    {
      return constraint_enforcement_;
    }

    /**
     * \brief Returns constraints enforcement strategy.
     */
    inline INPAR::FBI::BeamToFluidDiscretization GetContactDiscretization() const
    {
      return meshtying_discretization_;
    }

    /**
     * \brief Returns the penalty parameter.
     * @return penalty parameter.
     */
    inline double GetPenaltyParameter() const { return penalty_parameter_; }

    /**
     * \brief Returns the Gauss rule.
     * @return gauss rule.
     */
    inline DRT::UTILS::GaussRule1D GetGaussRule() const { return gauss_rule_; }

   private:
    //! Flag if Init was called.
    bool isinit_;

    //! Flag if Setup was called.
    bool issetup_;

    //! Enforcement strategy for constraints.
    INPAR::FBI::BeamToFluidConstraintEnforcement constraint_enforcement_;

    //! Discretization used for the contact.
    INPAR::FBI::BeamToFluidDiscretization meshtying_discretization_;

    //! Penalty parameter.
    double penalty_parameter_;

    //! Gauss rule to be used.
    DRT::UTILS::GaussRule1D gauss_rule_;
  };

}  // namespace FBI

#endif
