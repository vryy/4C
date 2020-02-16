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
#include "../drt_inpar/inpar_beam_to_solid.H"

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
    inline INPAR::BEAMTOSOLID::BeamToSolidConstraintEnforcement GetConstraintEnforcement() const
    {
      return constraint_enforcement_;
    }

    /**
     * \brief Returns constraints enforcement strategy.
     */
    inline INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization GetContactDiscretization() const
    {
      return contact_discretization_;
    }

    /**
     * \brief Returns the shape function for the mortar Lagrange-multiplicators.
     */
    inline INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions GetMortarShapeFunctionType() const
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
     * \brief Returns the number of integration points along the circumference of the beams cross
     * section.
     * @return Number of points.
     */
    inline unsigned int GetNumberOfIntegrationPointsCircumference() const
    {
      return integration_points_circumference_;
    }

    /**
     * \brief Returns a pointer to the vtk output parameters.
     * @return Pointer to vtk output parameters.
     */
    Teuchos::RCP<BeamToSolidVolumeMeshtyingVtkOutputParams> GetVtkOuputParamsPtr();

    /**
     * \brief Returns if the restart configuration should be coupled.
     */
    inline bool GetCoupleRestartState() const { return couple_restart_state_; }

   private:
    //! Flag if Init was called.
    bool isinit_;

    //! Flag if Setup was called.
    bool issetup_;

    //! Enforcement strategy for constraints.
    INPAR::BEAMTOSOLID::BeamToSolidConstraintEnforcement constraint_enforcement_;

    //! Discretization used for the contact.
    INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization contact_discretization_;

    //! Shape function for the mortar Lagrange-multiplicators
    INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions mortar_shape_function_;

    //! Penalty parameter.
    double penalty_parameter_;

    //! Gauss rule to be used.
    DRT::UTILS::GaussRule1D gauss_rule_;

    //! Number of integration points along the circumferencial direction in the cross section
    //! projection.
    unsigned int integration_points_circumference_;

    //! Pointer to the vtk output parameters for beam to solid volume meshtying.
    Teuchos::RCP<BeamToSolidVolumeMeshtyingVtkOutputParams> output_params_ptr_;

    //! If the coupling terms should be evaluated with the restart configuration.
    bool couple_restart_state_;
  };

}  // namespace BEAMINTERACTION

#endif
