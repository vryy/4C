/*!
\file beam_to_solid_volume_meshtying_params.H

\brief Data container holding all beam to solid volume meshtying input parameters.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VOLUME_MESHTYING_PARAMS_H
#define BEAM_TO_SOLID_VOLUME_MESHTYING_PARAMS_H


#include "../drt_fem_general/drt_utils_integration.H"

#include "../drt_lib/drt_dserror.H"


namespace BEAMINTERACTION
{
  /**
   * \brief Class for beam to solid meshtying parameters.
   */
  class BeamToSolidVolumeMeshtyingParams
  {
   public:
    /**
     * \brief Constructor.
     */
    BeamToSolidVolumeMeshtyingParams();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidVolumeMeshtyingParams(){};

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
     * \brief Returns the penalty parameter.
     * @return penalty parameter.
     */
    inline double PenaltyParameter() const { return penalty_parameter_; }

    /**
     * \brief Returns the Gauss rule.
     * @return gauss rule.
     */
    inline DRT::UTILS::GaussRule1D GaussRule() const { return gauss_rule_; }

   private:
    //! Flag if Init was called.
    bool isinit_;

    //! Flag if Setup was called.
    bool issetup_;

    //! Penalty parameter.
    double penalty_parameter_;

    //! Gauss rule to be used.
    DRT::UTILS::GaussRule1D gauss_rule_;
  };

}  // namespace BEAMINTERACTION

#endif
