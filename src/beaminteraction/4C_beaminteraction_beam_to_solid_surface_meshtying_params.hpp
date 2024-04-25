/*----------------------------------------------------------------------*/
/*! \file

\brief Data container holding all beam to solid surface meshtying input parameters.

\level 3
*/
// End doxygen header.


#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_SURFACE_MESHTYING_PARAMS_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_SURFACE_MESHTYING_PARAMS_HPP


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
   * \brief Class for beam to solid meshtying parameters.
   */
  class BeamToSolidSurfaceMeshtyingParams : public BeamToSolidParamsBase
  {
   public:
    /**
     * \brief Constructor.
     */
    BeamToSolidSurfaceMeshtyingParams();


    /**
     * \brief Initialize with the stuff coming from input file.
     */
    void Init() override;

    /**
     * \brief Returns the coupling type for beam-to-surface coupling.
     */
    inline INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling GetCouplingType() const
    {
      return coupling_type_;
    }

    /**
     * \brief Returns true if the coupling should be evaluated with FAD.
     */
    inline bool GetIsFAD() const override
    {
      switch (coupling_type_)
      {
        case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::reference_configuration_forced_to_zero:
        case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::displacement:
          return false;
          break;
        case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::
            reference_configuration_forced_to_zero_fad:
        case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::displacement_fad:
        case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::consistent_fad:
          return true;
          break;
        default:
          FOUR_C_THROW("Wrong coupling type.");
          break;
      }
      return false;
    }

    /**
     * \brief Returns the order of the FAD type.
     */
    inline int GetFADOrder() const override
    {
      if (GetIsFAD())
        return 2;
      else
        return 0;
    }

    /**
     * \brief Returns if rotational coupling is activated.
     */
    inline bool GetIsRotationalCoupling() const { return rotational_coupling_; }

    /**
     * \brief Returns the penalty parameter for rotational coupling.
     */
    inline double GetRotationalCouplingPenaltyParameter() const
    {
      return rotational_coupling_penalty_parameter_;
    }

    /**
     * \brief Returns the type of surface triad construction for rotational coupling.
     */
    inline INPAR::BEAMTOSOLID::BeamToSolidSurfaceRotationCoupling GetSurfaceTriadConstruction()
        const
    {
      return rotational_coupling_triad_construction_;
    }

    /**
     * \brief Returns a pointer to the visualization output parameters.
     * @return Pointer to visualization output parameters.
     */
    Teuchos::RCP<BeamToSolidSurfaceVisualizationOutputParams> GetVisualizationOutputParamsPtr();

   private:
    //! How the coupling should be evaluated.
    INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling coupling_type_;

    //! Pointer to the visualization output parameters for beam to solid volume meshtying.
    Teuchos::RCP<BeamToSolidSurfaceVisualizationOutputParams> output_params_ptr_;

    //! If rotational coupling should be activated.
    bool rotational_coupling_;

    //! Penalty parameter for rotational coupling.
    double rotational_coupling_penalty_parameter_;

    //! Type of surface triad construction.
    INPAR::BEAMTOSOLID::BeamToSolidSurfaceRotationCoupling rotational_coupling_triad_construction_;
  };

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
