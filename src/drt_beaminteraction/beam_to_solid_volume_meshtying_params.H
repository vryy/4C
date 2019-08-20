/*----------------------------------------------------------------------*/
/*! \file

\brief Data container holding all beam to solid volume meshtying input parameters.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VOLUME_MESHTYING_PARAMS_H
#define BEAM_TO_SOLID_VOLUME_MESHTYING_PARAMS_H


#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_inpar/inpar_beaminteraction.H"

#include "../drt_lib/drt_dserror.H"


// Forward declaration.
namespace BEAMINTERACTION
{
  class BeamToSolidVolumeMeshtyingVtkOutputParams;
}


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
     * \brief Returns the contact discretization method.
     */
    inline INPAR::BEAMINTERACTION::BeamToSolidVolumeConstraintEnforcement GetConstraintEnforcement()
        const
    {
      return constraint_enforcement_;
    }

    /**
     * \brief Returns constraints enforcement strategy.
     */
    inline INPAR::BEAMINTERACTION::BeamToSolidVolumeContactDiscretization GetContactDiscretization()
        const
    {
      return contact_discretization_;
    }

    /**
     * \brief Returns the shape function for the mortar Lagrange-multiplicators.
     */
    inline INPAR::BEAMINTERACTION::BeamToSolidVolumeMortarShapefunctions
    GetMortarShapeFunctionType() const
    {
      return mortar_shape_function_;
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

    /**
     * \brief Returns the number of integration points along the circumfence of the beams cross
     * section.
     * @return Number of points.
     */
    inline unsigned int GetNumberOfIntegrationPointsCircumfence() const
    {
      return integration_points_circumfence_;
    }

    /**
     * \brief Returns a pointer to the vtk output parameters.
     * @return Pointer to vtk output parameters.
     */
    Teuchos::RCP<BeamToSolidVolumeMeshtyingVtkOutputParams> GetVtkOuputParamsPtr();

   private:
    //! Flag if Init was called.
    bool isinit_;

    //! Flag if Setup was called.
    bool issetup_;

    //! Enforcement strategy for constraints.
    INPAR::BEAMINTERACTION::BeamToSolidVolumeConstraintEnforcement constraint_enforcement_;

    //! Discretization used for the contact.
    INPAR::BEAMINTERACTION::BeamToSolidVolumeContactDiscretization contact_discretization_;

    //! Shape function for the mortar Lagrange-multiplicators
    INPAR::BEAMINTERACTION::BeamToSolidVolumeMortarShapefunctions mortar_shape_function_;

    //! Penalty parameter.
    double penalty_parameter_;

    //! Gauss rule to be used.
    DRT::UTILS::GaussRule1D gauss_rule_;

    //! Number of integration points along the circumfencial direction in the cross section
    //! projection.
    unsigned int integration_points_circumfence_;

    //! Pointer to the vtk output parameters for beam to solid volume meshtying.
    Teuchos::RCP<BeamToSolidVolumeMeshtyingVtkOutputParams> output_params_ptr_;
  };

}  // namespace BEAMINTERACTION

#endif
