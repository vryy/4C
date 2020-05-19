/*----------------------------------------------------------------------*/
/*! \file

\brief Data container holding all beam to solid surface meshtying input parameters.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_SURFACE_MESHTYING_PARAMS_H
#define BEAM_TO_SOLID_SURFACE_MESHTYING_PARAMS_H


#include "beam_to_solid_params_base.H"


// Forward declaration.
namespace BEAMINTERACTION
{
  class BeamToSolidSurfaceVtkOutputParams;
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
     * \brief Destructor.
     */
    virtual ~BeamToSolidSurfaceMeshtyingParams(){};

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
        case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::configurations_forced_to_zero:
        case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::displacements:
          return false;
          break;
        case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::configurations_forced_to_zero_fad:
        case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::displacements_fad:
        case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::consistent_fad:
          return true;
          break;
        default:
          dserror("Wrong coupling type.");
          break;
      }
      return false;
    }

    /**
     * \brief Returns a pointer to the vtk output parameters.
     * @return Pointer to vtk output parameters.
     */
    Teuchos::RCP<BeamToSolidSurfaceVtkOutputParams> GetVtkOuputParamsPtr();

   private:
    //! How the coupling should be evaluated.
    INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling coupling_type_;

    //! Pointer to the vtk output parameters for beam to solid volume meshtying.
    Teuchos::RCP<BeamToSolidSurfaceVtkOutputParams> output_params_ptr_;
  };

}  // namespace BEAMINTERACTION

#endif
