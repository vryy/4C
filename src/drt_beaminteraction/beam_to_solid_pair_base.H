/*----------------------------------------------------------------------*/
/*! \file

\brief Base element for interactions between a beam and a solid.

\level 3
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_PAIR_BASE_H_
#define BEAM_TO_SOLID_PAIR_BASE_H_


#include "beam_contact_pair.H"

#include "../drt_geometry_pair/geometry_pair_utility_classes.H"


namespace BEAMINTERACTION
{
  /**
   * \brief Base class for beam to solid interactions.
   * @tparam scalar_type_fad Scalar FAD type to be used in this pair.
   * @tparam beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @tparam solid Type from GEOMETRYPAIR::ElementDiscretization... representing the solid.
   */
  template <typename scalar_type_fad, typename beam, typename solid>
  class BeamToSolidPairBase : public BeamContactPair
  {
   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidPairBase();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidPairBase(){};

    /**
     * \brief Setup the contact pair (derived).
     *
     * This method sets the beam reference positions for this pair.
     */
    void Setup() override;

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
     * \brief Update state of translational nodal DoFs (absolute positions and tangents) of the beam
     * element.
     * @param beam_centerline_dofvec
     * @param solid_nodal_dofvec
     */
    void ResetState(const std::vector<double>& beam_centerline_dofvec,
        const std::vector<double>& solid_nodal_dofvec) override;

    /**
     * \brief Set the restart displacement in this pair.
     *
     * If coupling interactions should be evaluated w.r.t the restart state, this method will set
     * them in the pair accordingly.
     *
     * @param centerline_restart_vec_ (in) Vector with the centerline displacements at the restart
     * step, for all contained elements (Vector of vector).
     */
    virtual void SetRestartDisplacement(
        const std::vector<std::vector<double>>& centerline_restart_vec_) override;

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
      if (line_to_3D_segments_.size() > 0)
        return true;
      else
        return false;
    };

    /**
     * \brief Get number of active contact point pairs on this element pair. Not yet implemented.
     */
    unsigned int GetNumAllActiveContactPointPairs() const override
    {
      dserror("GetNumAllActiveContactPointPairs not yet implemented!");
      return 0;
    };

    /**
     * \brief Get coordinates of all active contact points on element1. Not yet implemented.
     */
    void GetAllActiveContactPointCoordsElement1(
        std::vector<LINALG::Matrix<3, 1, double>>& coords) const override
    {
      dserror("GetAllActiveContactPointCoordsElement1 not yet implemented!");
    }

    /**
     * \brief Get coordinates of all active contact points on element2. Not yet implemented.
     */
    void GetAllActiveContactPointCoordsElement2(
        std::vector<LINALG::Matrix<3, 1, double>>& coords) const override
    {
      dserror("GetAllActiveContactPointCoordsElement2 not yet implemented!");
    }

    /**
     * \brief Get all (scalar) contact forces of this contact pair. Not yet implemented.
     */
    void GetAllActiveContactForces(std::vector<double>& forces) const override
    {
      dserror("GetAllActiveContactForces not yet implemented!");
    }

    /**
     * \brief Get all (scalar) gap values of this contact pair. Not yet implemented.
     */
    void GetAllActiveContactGaps(std::vector<double>& gaps) const override
    {
      dserror("GetAllActiveContactGaps not yet implemented!");
    }

    /**
     * \brief Get energy of penalty contact. Not yet implemented.
     */
    double GetEnergy() const override
    {
      dserror("GetEnergy not implemented yet!");
      return 0.0;
    }

   protected:
    /**
     * \brief This function evaluates the beam position at an integration point for the pairs.
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
        const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& integration_point,
        LINALG::Matrix<3, 1, scalar_type_fad>& r_beam, bool reference) const;

   protected:
    //! Vector with the segments of the line to 3D pair.
    std::vector<GEOMETRYPAIR::LineSegment<double>> line_to_3D_segments_;

    //! Current nodal positions (and tangents) of the beam.
    LINALG::Matrix<beam::n_dof_, 1, scalar_type_fad> ele1pos_;

    //! Reference nodal positions (and tangents) of the beam.
    LINALG::Matrix<beam::n_dof_, 1, double> ele1posref_;
  };
}  // namespace BEAMINTERACTION

#endif
