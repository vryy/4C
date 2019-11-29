/*----------------------------------------------------------------------*/
/*! \file

\brief Base meshtying element for meshtying between a 1D beam and a 3D fluid element.

\level 2
\maintainer Nora Hagmeyer
*/

#ifndef BEAM_TO_FLUID_VOLUME_MESHTYING_PAIR_BASE_H_
#define BEAM_TO_FLUID_VOLUME_MESHTYING_PAIR_BASE_H_

// todo declutter
#include "../drt_beaminteraction/beam_contact_pair.H"
#include "../drt_beaminteraction/beam3contact_defines.H"
#include "../drt_geometry_pair/geometry_pair_utility_classes.H"

#include "../linalg/linalg_utils.H"
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
  class LineTo3DEvaluationData;
  class GeometryEvaluationDataBase;
}  // namespace GEOMETRYPAIR


namespace BEAMINTERACTION
{
  /**
   * \brief Class representing a pair of elements for beam to fluid meshtying
   *
   * \params[in] beam Type from GEOMETRYPAIR::ElementDiscretization representing the beam.
   * \params[in] fluid Type from GEOMETRYPAIR::ElementDiscretization representing the fluid.
   */
  template <typename beam, typename fluid>
  class BeamToFluidMeshtyingPairBase
      : public BeamContactPair  // todo can i just derice from beamtosolidmeshtyingpairbase?
  {
   public:
    /**
     * \brief empty Destructor.
     */
    virtual ~BeamToFluidMeshtyingPairBase(){};

    /**
     * \brief Initialize the contact pair.
     */
    void Init(const Teuchos::RCP<BEAMINTERACTION::BeamContactParams> params_ptr,
        std::vector<DRT::Element const*> elements);

    /**
     * \brief Setup the contact pair and set information on the current position of the elements in
     * the pair
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
     * \brief Update state of translational nodal DoFs (absolute positions, tangents and velocities)
     * of both elements.
     * @param beam_centerline_dofvec current nodal beam positions extracted from the element and
     * nodal velocities computed by the time integrator
     * @param fluid_nodal_dofvec current nodal fluid positions (only for ALE different from the
     * reference nodal values) and nodal velocities
     */
    void ResetState(const std::vector<double>& beam_centerline_dofvec,
        const std::vector<double>& fluid_nodal_dofvec) override;

    /**
     *\brief Unset meshtying_is evaluated in order to reevaluate the geometry pair todo maybe we
     *find something better than a public method?
     */
    void UnsetEvaluationFlag() { meshtying_is_evaluated_ = false; }
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
        std::vector<LINALG::Matrix<3, 1, double>>& coords) const override
    {
      dserror("not yet implemented!");
    }

    /**
     * \brief Get coordinates of all active contact points on element2. Not yet implemented.
     */
    void GetAllActiveContactPointCoordsElement2(
        std::vector<LINALG::Matrix<3, 1, double>>& coords) const override
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

    /**
     * \brief Create the geometry pair for this contact pair.
     * @param geometry_evaluation_data_ptr Evaluation data that will be linked to the pair.
     */
    void CreateGeometryPair(
        const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
        override;

   private:  //! Type to be used for scalar AD variables. This can not be inherited from the base
             //! class.
    typedef Sacado::ELRFad::SLFad<double, beam::n_dof_ + fluid::n_dof_> TYPE_BTS_VMT_AD;

   protected:
    /** \brief You will have to use the FBI::PairFactory
     *
     */
    BeamToFluidMeshtyingPairBase();

    virtual void EvaluatePenaltyForce(LINALG::Matrix<3, 1, TYPE_BTS_VMT_AD>& force,
        const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& projected_gauss_point,
        LINALG::Matrix<3, 1, TYPE_BTS_VMT_AD> v_beam) const
    {
      dserror("Not implemented yet!\n");
    };

    void EvaluateBeamPosition(const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& integration_point,
        LINALG::Matrix<3, 1, TYPE_BTS_VMT_AD>& r_beam, bool reference) const;

    /**
     * \brief Return a cast of the geometry pair to the type for this contact pair.
     * \returns RPC with the type of geometry pair for this beam contact pair.
     */
    inline Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToVolume<double, beam, fluid>>
    CastGeometryPair() const
    {
      return Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::GeometryPairLineToVolume<double, beam, fluid>>(
          geometry_pair_, true);
    };

    //! Flag if the meshtying has been evaluated in this newton step already.
    bool meshtying_is_evaluated_;

    //! Current nodal positions (and tangents) of the two elements.
    LINALG::Matrix<beam::n_dof_, 1, TYPE_BTS_VMT_AD> ele1pos_;
    LINALG::Matrix<fluid::n_dof_, 1, TYPE_BTS_VMT_AD> ele2pos_;

    //! Current nodal velocities of the two elements.
    LINALG::Matrix<beam::n_dof_, 1, TYPE_BTS_VMT_AD> ele1vel_;
    LINALG::Matrix<fluid::n_dof_, 1, TYPE_BTS_VMT_AD> ele2vel_;

    //! Reference nodal positions (and tangents) of the two elements.
    LINALG::Matrix<beam::n_dof_, 1, TYPEBTS> ele1posref_;
    LINALG::Matrix<fluid::n_dof_, 1, TYPEBTS> ele2posref_;

    //! Current nodal positions (and tangents) of the two elements.
    LINALG::Matrix<beam::n_dof_, 1, TYPEBTS> ele1poscur_;
    LINALG::Matrix<fluid::n_dof_, 1, TYPEBTS> ele2poscur_;

    //! Vector with the segments of the line to volume pair.
    std::vector<GEOMETRYPAIR::LineSegment<double>> line_to_volume_segments_;
  };
}  // namespace BEAMINTERACTION

#endif
