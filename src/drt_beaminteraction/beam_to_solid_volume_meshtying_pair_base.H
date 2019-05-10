/*!
\file beam_to_solid_volume_meshtying_pair_base.H

\brief Base meshtying element for meshtying between a 3D beam and a 3D solid element.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_BASE_H_
#define BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_BASE_H_


#include "beam_contact_pair.H"

#include "beam3contact_defines.H"
#include "../drt_geometry_pair/geometry_pair_utility_classes.H"

#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"


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
   * \brief Class for beam to solid meshtying.
   * @param beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @param solid Type from GEOMETRYPAIR::ElementDiscretization... representing the solid.
   */
  template <typename beam, typename solid>
  class BeamToSolidVolumeMeshtyingPairBase : public BeamContactPair
  {
   protected:
    //! Type to be used for scalar AD variables. This can not be inherited from the base class.
    typedef Sacado::ELRFad::SLFad<double, beam::n_dof_ + solid::n_dof_> TYPE_BTS_VMT_AD;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidVolumeMeshtyingPairBase();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidVolumeMeshtyingPairBase(){};
    /**
     * \brief Initialize the contact pair.
     */
    void Init(const Teuchos::RCP<BEAMINTERACTION::BeamContactParams> params_ptr,
        const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> geometry_evaluation_data_ptr,
        std::vector<DRT::Element const*> elements);

    /**
     * \brief Setup the contact pair.
     */
    void Setup() override;

    /**
     * \brief Things that need to be done in a separate loop before the actual evaluation loop over
     * all contact pairs.
     */
    void PreEvaluate() override;

    /**
     * \brief Evaluate this contact element pair.
     * @param forcevec1 (out) Force vector on element 1.
     * @param forcevec2 (out) Force vector on element 2.
     * @param stiffmat11 (out) Stiffness contributions on element 1 - element 1.
     * @param stiffmat12 (out) Stiffness contributions on element 1 - element 2.
     * @param stiffmat21 (out) Stiffness contributions on element 2 - element 1.
     * @param stiffmat22 (out) Stiffness contributions on element 2 - element 2.
     * @return True if pair is in contact.
     */
    bool Evaluate(LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
        LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
        LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22) override
    {
      return false;
    };

    /**
     * \brief Update state of translational nodal DoFs (absolute positions and tangents) of both
     * elements.
     * @param beam_centerline_dofvec
     * @param solid_nodal_dofvec
     */
    void ResetState(const std::vector<double>& beam_centerline_dofvec,
        const std::vector<double>& solid_nodal_dofvec) override;

    /**
     * \brief Print information about this beam contact element pair to screen.
     */
    void Print(std::ostream& out) const override;

    /**
     * \brief Print this beam contact element pair to screen.
     */
    void PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const override;

    /**
     * \brief Check if this pair is in contact. The correct value is only returned after PreEvaluate
     * and Evaluate are run on the geometry pair.
     * @return true if it is in contact.
     */
    inline bool GetContactFlag() const override
    {
      // The element pair is assumed to be active when we have at least one active contact point
      if (line_to_volume_segments_.size() > 0)
        return true;
      else
        return false;
    };

    /**
     * \brief Get number of active contact point pairs on this element pair. Not yet implemented.
     */
    unsigned int GetNumAllActiveContactPointPairs() const override
    {
      dserror("not yet implemented!");
      return 0;
    };

    /**
     * \brief Get coordinates of all active contact points on element1. Not yet implemented.
     */
    void GetAllActiveContactPointCoordsElement1(
        std::vector<LINALG::TMatrix<double, 3, 1>>& coords) const override
    {
      dserror("not yet implemented!");
    }

    /**
     * \brief Get coordinates of all active contact points on element2. Not yet implemented.
     */
    void GetAllActiveContactPointCoordsElement2(
        std::vector<LINALG::TMatrix<double, 3, 1>>& coords) const override
    {
      dserror("not yet implemented!");
    }

    /**
     * \brief Get all (scalar) contact forces of this contact pair. Not yet implemented.
     */
    void GetAllActiveContactForces(std::vector<double>& forces) const override
    {
      dserror("not yet implemented!");
    }

    /**
     * \brief Get all (scalar) gap values of this contact pair. Not yet implemented.
     */
    void GetAllActiveContactGaps(std::vector<double>& gaps) const override
    {
      dserror("not yet implemented!");
    }

    /**
     * \brief Get energy of penalty contact. Not yet implemented.
     */
    double GetEnergy() const override
    {
      dserror("not implemented yet!");
      return 0.0;
    }

    /**
     * \brief Add the visualization of this pair to the beam to solid vtu output writer.
     *
     * This base class creates output of (if selected in the input file) the segmentation, the
     * integration points - and if implemented in the derived classes - the forces at the
     * integration points.
     *
     * @param visualization_writer (out) Object that manages all visualization related data for beam
     * to solid pairs.
     * @param visualization_params (in) Parameter list (not used in this class).
     */
    void GetPairVisualization(Teuchos::RCP<BeamToSolidVtuOutputWriterBase> visualization_writer,
        const Teuchos::ParameterList& visualization_params) const override;

   protected:
    /**
     * \brief Create the geometry pair for this contact pair.
     * @param geometry_evaluation_data_ptr Evaluation data that will be linked to the pair.
     */
    void CreateGeometryPair(const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal>
            geometry_evaluation_data_ptr);

    /**
     * \brief Return a cast of the geometry pair to the type for this contact pair.
     * @return RPC with the type of geometry pair for this beam contact pair.
     */
    inline Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToVolume<double, beam, solid>>
    CastGeometryPair() const
    {
      return Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::GeometryPairLineToVolume<double, beam, solid>>(
          geometry_pair_, true);
    };

    /**
     * \brief This function evaluates the position at an integration point for the pairs.
     *
     * This is needed because the cross section pairs have 3 parameter coordinates on the beam and
     * the other pairs have 1. This method is mainly used for visualization.
     *
     * @param integration_point (in) Integration where the position should be evaluated.
     * @param r_beam (out) Position on the beam.
     * @param reference (in) True -> the reference position is calculated, False -> the current
     * position is calculated.
     */
    virtual void EvaluateBeamPosition(
        const GEOMETRYPAIR::ProjectionPointLineToVolume<double>& integration_point,
        LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1>& r_beam, bool reference) const;

    /**
     * \brief This function evaluates the penalty force from a given beam position and a given solid
     * position.
     *
     * This method is mainly used for visualization.
     *
     * @param r_beam (in) Position on the beam.
     * @param r_solid (in) Position on the solid.
     * @param force (out) Force acting on the beam (the negative force acts on the solid).
     */
    virtual void EvaluatePenaltyForce(const LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1>& r_beam,
        const LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1>& r_solid,
        LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1>& force) const;

   protected:
    //! Flag if the meshtying has been evaluated already.
    bool meshtying_is_evaluated_;

    //! Current nodal positions (and tangents) of the two elements.
    LINALG::TMatrix<TYPE_BTS_VMT_AD, beam::n_dof_, 1> ele1pos_;
    LINALG::TMatrix<TYPE_BTS_VMT_AD, solid::n_dof_, 1> ele2pos_;

    //! Reference nodal positions (and tangents) of the two elements.
    LINALG::TMatrix<TYPEBTS, beam::n_dof_, 1> ele1posref_;
    LINALG::TMatrix<TYPEBTS, solid::n_dof_, 1> ele2posref_;

    //! Vector with the segments of the line to volume pair.
    std::vector<GEOMETRYPAIR::LineSegment<double>> line_to_volume_segments_;
  };
}  // namespace BEAMINTERACTION

#endif
