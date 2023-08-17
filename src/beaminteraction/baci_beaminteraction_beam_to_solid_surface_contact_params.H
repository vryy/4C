/*----------------------------------------------------------------------*/
/*! \file

\brief Data container holding all beam to solid surface contact input parameters.

\level 3
*/
// End doxygen header.


#ifndef BACI_BEAMINTERACTION_BEAM_TO_SOLID_SURFACE_CONTACT_PARAMS_H
#define BACI_BEAMINTERACTION_BEAM_TO_SOLID_SURFACE_CONTACT_PARAMS_H


#include "baci_beaminteraction_beam_to_solid_params_base.H"


// Forward declaration.
namespace BEAMINTERACTION
{
  class BeamToSolidSurfaceVtkOutputParams;
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
     * \brief Destructor.
     */
    virtual ~BeamToSolidSurfaceContactParams(){};

    /**
     * \brief Initialize with the stuff coming from input file.
     */
    void Init() override;

    /**
     * \brief Returns true if the coupling should be evaluated with FAD.
     */
    inline bool GetIsFAD() const override { return true; }

    /**
     * \brief Returns the order of the FAD type.
     */
    int GetFADOrder() const override;

    /**
     * \brief Returns the contact type.
     */
    inline INPAR::BEAMTOSOLID::BeamToSolidSurfaceContact GetContactType() const
    {
      return contact_type_;
    }

    /**
     * \brief Returns the type of penalty law.
     */
    inline INPAR::BEAMTOSOLID::BeamToSolidSurfaceContactPenaltyLaw GetPenaltyLaw() const
    {
      return penalty_law_;
    }

    /**
     * \brief Returns the regularization parameter of the penalty law.
     */
    inline double GetPenaltyParameterG0() const { return penalty_parameter_g0_; }

    /**
     * \brief Returns a pointer to the vtk output parameters.
     * @return Pointer to vtk output parameters.
     */
    Teuchos::RCP<BeamToSolidSurfaceVtkOutputParams> GetVtkOuputParamsPtr()
    {
      return output_params_ptr_;
    }

   private:
    //! Type of contact constraints.
    INPAR::BEAMTOSOLID::BeamToSolidSurfaceContact contact_type_;

    //! Type of penalty law.
    INPAR::BEAMTOSOLID::BeamToSolidSurfaceContactPenaltyLaw penalty_law_;

    //! Regularization parameter for the penalty law.
    double penalty_parameter_g0_;

    //! Pointer to the vtk output parameters for beam to solid surface problems.
    Teuchos::RCP<BeamToSolidSurfaceVtkOutputParams> output_params_ptr_;
  };

}  // namespace BEAMINTERACTION

#endif
