/*----------------------------------------------------------------------*/
/*! \file

\brief Data container holding all beam to solid surface contact input parameters.

\level 3
*/
// End doxygen header.


#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_SURFACE_CONTACT_PARAMS_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_SURFACE_CONTACT_PARAMS_HPP


#include "4C_config.hpp"

#include "4C_beaminteraction_beam_to_solid_params_base.hpp"

FOUR_C_NAMESPACE_OPEN


// Forward declaration.
namespace BEAMINTERACTION
{
  class BeamToSolidSurfaceVisualizationOutputParams;
}

namespace BEAMINTERACTION
{
  /**
   * \brief Class for beam to solid contact parameters.
   */
  class BeamToSolidSurfaceContactParams : public BeamToSolidParamsBase
  {
   public:
    /**
     * \brief Constructor.
     */
    BeamToSolidSurfaceContactParams();


    /**
     * \brief Initialize with the stuff coming from input file.
     */
    void init() override;

    /**
     * \brief Returns true if the coupling should be evaluated with FAD.
     */
    inline bool get_is_fad() const override { return true; }

    /**
     * \brief Returns the order of the FAD type.
     */
    int get_fad_order() const override;

    /**
     * \brief Returns the contact type.
     */
    inline Inpar::BeamToSolid::BeamToSolidSurfaceContact get_contact_type() const
    {
      return contact_type_;
    }

    /**
     * \brief Returns the type of penalty law.
     */
    inline Inpar::BeamToSolid::BeamToSolidSurfaceContactPenaltyLaw get_penalty_law() const
    {
      return penalty_law_;
    }

    /**
     * \brief Returns the regularization parameter of the penalty law.
     */
    inline double get_penalty_parameter_g0() const { return penalty_parameter_g0_; }

    /**
     * \brief Returns the configuration where the mortar contact is defined in.
     */
    inline Inpar::BeamToSolid::BeamToSolidSurfaceContactMortarDefinedIn
    get_beam_to_solid_surface_contact_mortar_defined_in() const
    {
      return mortar_contact_configuration_;
    }

    /**
     * \brief Returns a pointer to the visualization output parameters.
     * @return Pointer to visualization output parameters.
     */
    Teuchos::RCP<BeamToSolidSurfaceVisualizationOutputParams> get_visualization_output_params_ptr()
    {
      return output_params_ptr_;
    }

   private:
    //! Type of contact constraints.
    Inpar::BeamToSolid::BeamToSolidSurfaceContact contact_type_;

    //! Type of penalty law.
    Inpar::BeamToSolid::BeamToSolidSurfaceContactPenaltyLaw penalty_law_;

    //! Regularization parameter for the penalty law.
    double penalty_parameter_g0_;

    //! Configuration where the mortar contact is defined
    Inpar::BeamToSolid::BeamToSolidSurfaceContactMortarDefinedIn mortar_contact_configuration_;

    //! Pointer to the visualization output parameters for beam to solid surface problems.
    Teuchos::RCP<BeamToSolidSurfaceVisualizationOutputParams> output_params_ptr_;
  };

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
