/*----------------------------------------------------------------------*/
/*! \file

\brief Data container holding all beam to solid volume meshtying input parameters.

\level 3
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VOLUME_MESHTYING_PARAMS_H
#define BEAM_TO_SOLID_VOLUME_MESHTYING_PARAMS_H


#include "beaminteraction_beam_to_solid_params_base.H"


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
  class BeamToSolidVolumeMeshtyingParams : public BeamToSolidParamsBase
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
    void Init() override;

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
     * \brief Returns the type of rotational coupling to be used.
     */
    inline INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling GetRotationalCouplingType() const
    {
      return rotation_coupling_;
    }

    /**
     * \brief Returns the shape function for the mortar Lagrange-multiplicators to interpolate the
     * rotational coupling.
     */
    inline INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions GetMortarShapeFunctionRotationType()
        const
    {
      return mortar_shape_function_rotation_;
    }

    /**
     * \brief Returns the penalty parameter for rotational coupling.
     */
    inline double GetRotationalCouplingPenaltyParameter() const
    {
      return rotational_coupling_penalty_parameter_;
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
    //! Number of integration points along the circumferencial direction in the cross section
    //! projection.
    unsigned int integration_points_circumference_;

    //! Type of rotational coupling.
    INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling rotation_coupling_;

    //! Shape function for the mortar Lagrange-multiplicators
    INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions mortar_shape_function_rotation_;

    //! Penalty parameter for rotational coupling.
    double rotational_coupling_penalty_parameter_;

    //! Pointer to the vtk output parameters for beam to solid volume meshtying.
    Teuchos::RCP<BeamToSolidVolumeMeshtyingVtkOutputParams> output_params_ptr_;

    //! If the coupling terms should be evaluated with the restart configuration.
    bool couple_restart_state_;
  };

}  // namespace BEAMINTERACTION

#endif
