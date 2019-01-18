/*!
\file beam_to_solid_volume_meshtying_pair_mortar.H

\brief Meshtying element for meshtying between a 3D beam and a 3D solid element using mortar shape
functions for the traction.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_MORTAR_H_
#define BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_MORTAR_H_


#include "beam_to_solid_volume_meshtying_pair_base.H"

#include "../drt_geometry_pair/geometry_pair_element_types.H"
#include "../drt_inpar/inpar_beaminteraction.H"
#include "../drt_lib/drt_element.H"


// Forward declarations.
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace LINALG
namespace GEOMETRYPAIR
{
  template <typename scalar_type, typename line, typename volume>
  class GeometryPairLineToVolume;
}  // namespace GEOMETRYPAIR


namespace BEAMINTERACTION
{
  /**
   * \brief This structure "converts" the type of a mortar shape function to a discretization. This
   * is done in this way, so the type conversion can be done at compile time.
   */
  template <INPAR::BEAMINTERACTION::BeamToSolidVolumeMortarShapefunctions mortar_discretization>
  struct MortarDiscretizationToElementDiscretization
  {
  };
  template <>
  struct MortarDiscretizationToElementDiscretization<
      INPAR::BEAMINTERACTION::BeamToSolidVolumeMortarShapefunctions::line2>
  {
    static const DRT::Element::DiscretizationType discretization_ = DRT::Element::line2;
    static const unsigned int n_val_ = 1;
  };
  template <>
  struct MortarDiscretizationToElementDiscretization<
      INPAR::BEAMINTERACTION::BeamToSolidVolumeMortarShapefunctions::line3>
  {
    static const DRT::Element::DiscretizationType discretization_ = DRT::Element::line3;
    static const unsigned int n_val_ = 1;
  };


  /**
   * This type can be used to template the beam to solid mortar pairs on. It gives all relevant data
   * needed at compile time, so fixed size matrices can be used for the element calcualtions.
   */
  template <INPAR::BEAMINTERACTION::BeamToSolidVolumeMortarShapefunctions mortar_discretization>
  class BeamToSolidVolumeMeshtyingMortarShape
      : public GEOMETRYPAIR::ElementDiscretizationBase<
            MortarDiscretizationToElementDiscretization<mortar_discretization>::discretization_,
            MortarDiscretizationToElementDiscretization<mortar_discretization>::n_val_>
  {
   public:
    static const unsigned int n_nodes_per_element_ =
        GEOMETRYPAIR::ElementDiscretizationBase<
            MortarDiscretizationToElementDiscretization<mortar_discretization>::discretization_,
            MortarDiscretizationToElementDiscretization<mortar_discretization>::n_val_>::n_nodes_ -
        2;
  };


  /**
   * Shortcuts to element types are created here, so the explicit template initialisations are
   * better readable.
   */
  typedef BeamToSolidVolumeMeshtyingMortarShape<
      INPAR::BEAMINTERACTION::BeamToSolidVolumeMortarShapefunctions::line2>
      t_mortar_line2;
  typedef BeamToSolidVolumeMeshtyingMortarShape<
      INPAR::BEAMINTERACTION::BeamToSolidVolumeMortarShapefunctions::line3>
      t_mortar_line3;


  /**
   * \brief Class for beam to solid meshtying using mortar shape functions for the contact
   * tractions.
   * @param beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @param solid Type from GEOMETRYPAIR::ElementDiscretization... representing the solid.
   * @param mortar Type from BEAMINTERACTION::BeamToSolidVolumeMeshtyingMortarShape...
   * mortar shape functions.
   */
  template <typename beam, typename solid, typename mortar>
  class BeamToSolidVolumeMeshtyingPairMortar
      : public BeamToSolidVolumeMeshtyingPairBase<beam, solid>
  {
   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidVolumeMeshtyingPairMortar();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidVolumeMeshtyingPairMortar(){};

   private:
    //! Type to be used for scalar AD variables. This can not be inherited from the base class.
    typedef Sacado::ELRFad::SLFad<double, beam::n_dof_ + solid::n_dof_> TYPE_BTS_VMT_AD;
  };
}  // namespace BEAMINTERACTION

#endif
