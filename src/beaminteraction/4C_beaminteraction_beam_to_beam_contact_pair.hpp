/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief One beam contact pair (two beam elements) consisting of several contact segments

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_BEAM_CONTACT_PAIR_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_BEAM_CONTACT_PAIR_HPP

#include "4C_config.hpp"

#include "4C_beam3_base.hpp"
#include "4C_beaminteraction_beam_to_beam_contact_defines.hpp"
#include "4C_beaminteraction_beam_to_beam_contact_tangentsmoothing.hpp"
#include "4C_beaminteraction_beam_to_beam_contact_utils.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_inpar_beamcontact.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_utils_fad.hpp"

#include <Sacado.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declaration ...
namespace CORE::LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace CORE::LINALG

namespace BEAMINTERACTION
{
  template <unsigned int numnodes, unsigned int numnodalvalues>
  class BeamToBeamContactVariables;

  /*!
   \brief class for contact between two 3D beam elements
   */
  template <unsigned int numnodes, unsigned int numnodalvalues>
  class BeamToBeamContactPair : public BeamContactPair
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}

    //! @name Constructors and destructors and related methods
    /*!
    \brief Standard Constructor
    */
    BeamToBeamContactPair();



    //! Setup
    void Setup() override;
    //@}

    //! @name Derived methods from base class
    /*!
    \brief things that need to be done in a separate loop before the actual evaluation loop
           over all contact pairs
    */
    void PreEvaluate() override;

    /*!
    \brief Evaluate this contact element pair, return value indicates whether pair is active,
           i.e. non-zero values for force and stiffmat are returned
    */
    bool Evaluate(CORE::LINALG::SerialDenseVector* forcevec1,
        CORE::LINALG::SerialDenseVector* forcevec2, CORE::LINALG::SerialDenseMatrix* stiffmat11,
        CORE::LINALG::SerialDenseMatrix* stiffmat12, CORE::LINALG::SerialDenseMatrix* stiffmat21,
        CORE::LINALG::SerialDenseMatrix* stiffmat22) override;

    /*
    \brief Update state of translational nodal DoFs (absolute positions and tangents) of both
    elements
    */
    void ResetState(const std::vector<double>& centerline_dofvec_ele1,
        const std::vector<double>& centerline_dofvec_ele2) override;

    /** \brief print information about this beam contact element pair to screen
     *
     *  \author grill
     *  \date 05/16 */
    void Print(std::ostream& out) const override;


    /** \brief print this beam contact element pair to screen
     *
     *  \author grill
     *  \date 12/16 */
    void print_summary_one_line_per_active_segment_pair(std::ostream& out) const override;
    //@}

    //! @name Access methods

    /*!
    \brief Get flag indicating whether contact is active (true) or inactive (false)
    */
    inline bool GetContactFlag() const override
    {
      // The element pair is assumed to be active when we have at least one active contact point
      return (cpvariables_.size() + gpvariables_.size() + epvariables_.size());
    }

    /*!
    \brief Get number of active contact point pairs on this element pair
    */
    unsigned int get_num_all_active_contact_point_pairs() const override
    {
      return (unsigned int)(cpvariables_.size() + gpvariables_.size() + epvariables_.size());
    }

    /*!
    \brief Get coordinates of all active contact points on element1 and element2
    */
    void get_all_active_contact_point_coords_element1(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& coords) const override;

    void get_all_active_contact_point_coords_element2(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& coords) const override;

    /*!
    \brief Get all (scalar) contact forces of this contact pair
    */
    void get_all_active_contact_forces(std::vector<double>& forces) const override;

    /*!
    \brief Get all (scalar) gap values of this contact pair
    */
    void get_all_active_contact_gaps(std::vector<double>& gaps) const override;

    /*!
    \brief Get energy of penalty contact.
    */
    double GetEnergy() const override;

    //@}

   private:
    //! @name member variables

    //! current node coordinates of the two elements
    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE> ele1pos_;
    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE> ele2pos_;

    //! neighbor elements of element 1
    Teuchos::RCP<BEAMINTERACTION::B3CNeighbor> neighbors1_;

    //! neighbor elements of element 2
    Teuchos::RCP<BEAMINTERACTION::B3CNeighbor> neighbors2_;

    //! cross section radius of first beam
    double r1_;

    //! cross section radius of second beam
    double r2_;

    //! Maximal gap at which a contact can become active
    double maxactivegap_;

    //! Maximal distance between a real segment on beam element 1 and its straight approximation
    double maxsegdist1_;

    //! Maximal distance between a real segment on beam element 2 and its straight approximation
    double maxsegdist2_;

    //! Number of segments on element1
    int numseg1_;

    //! Number of segments on element2
    int numseg2_;

    //! bound for search of large angle contact segment pairs
    double deltalargeangle_;

    //! bound for search of small angle contact segment pairs
    double deltasmallangle_;

    //! Indicates if the left / right node of the slave element 1 coincides with the endpoint of the
    //! physical beam (true) or not (false)
    std::pair<bool, bool> boundarynode1_;

    //! Indicates if the left / right node of the master element 2 coincides with the endpoint of
    //! the physical beam (true) or not (false)
    std::pair<bool, bool> boundarynode2_;

    //! Variables stored at the closest points of the large-angle-contact algorithm
    std::vector<Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues>>> cpvariables_;

    //! Variables stored at the Gauss points of the small-angle-contact algorithm
    std::vector<Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues>>> gpvariables_;

    //! Variables stored at the end points of the endpoint-contact algorithm
    std::vector<Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues>>> epvariables_;

    //@}

    //! @name Private evaluation methods

    /*!
    \brief Get active large angle pairs
    */
    void get_active_large_angle_pairs(std::vector<CORE::LINALG::Matrix<3, 1, double>>& endpoints1,
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& endpoints2,
        std::map<std::pair<int, int>, CORE::LINALG::Matrix<3, 1, double>>& closelargeanglesegments,
        const double& pp);

    /*!
    \brief Evaluate active large angle pairs
    */
    void evaluate_active_large_angle_pairs(CORE::LINALG::SerialDenseVector* forcevec1,
        CORE::LINALG::SerialDenseVector* forcevec2, CORE::LINALG::SerialDenseMatrix* stiffmat11,
        CORE::LINALG::SerialDenseMatrix* stiffmat12, CORE::LINALG::SerialDenseMatrix* stiffmat21,
        CORE::LINALG::SerialDenseMatrix* stiffmat22);

    /*!
    \brief Get active small angle pairs
    */
    void get_active_small_angle_pairs(
        std::map<std::pair<int, int>, CORE::LINALG::Matrix<3, 1, double>>& closesmallanglesegments,
        std::pair<int, int>* iminmax = nullptr,
        std::pair<bool, bool>* leftrightsolutionwithinsegment = nullptr,
        std::pair<double, double>* eta1_leftrightboundary = nullptr);

    /*!
    \brief Evaluate active small angle pairs
    */
    void evaluate_active_small_angle_pairs(CORE::LINALG::SerialDenseVector* forcevec1,
        CORE::LINALG::SerialDenseVector* forcevec2, CORE::LINALG::SerialDenseMatrix* stiffmat11,
        CORE::LINALG::SerialDenseMatrix* stiffmat12, CORE::LINALG::SerialDenseMatrix* stiffmat21,
        CORE::LINALG::SerialDenseMatrix* stiffmat22, std::pair<int, int>* iminmax = nullptr,
        std::pair<bool, bool>* leftrightsolutionwithinsegment = nullptr,
        std::pair<double, double>* eta1_leftrightboundary = nullptr);

    /*!
    \brief Get active endpoint pairs
    */
    void get_active_end_point_pairs(
        std::vector<std::pair<int, int>>& closeendpointsegments, const double pp);

    /*!
    \brief Evaluate active endpoint pairs
    */
    void evaluate_active_end_point_pairs(CORE::LINALG::SerialDenseVector* forcevec1,
        CORE::LINALG::SerialDenseVector* forcevec2, CORE::LINALG::SerialDenseMatrix* stiffmat11,
        CORE::LINALG::SerialDenseMatrix* stiffmat12, CORE::LINALG::SerialDenseMatrix* stiffmat21,
        CORE::LINALG::SerialDenseMatrix* stiffmat22);

    /*!
    \brief Find segments close to each other
    */
    void GetCloseSegments(const std::vector<CORE::LINALG::Matrix<3, 1, double>>& endpoints1,
        const std::vector<CORE::LINALG::Matrix<3, 1, double>>& endpoints2,
        std::map<std::pair<int, int>, CORE::LINALG::Matrix<3, 1, double>>& closesmallanglesegments,
        std::map<std::pair<int, int>, CORE::LINALG::Matrix<3, 1, double>>& closelargeanglesegments,
        std::vector<std::pair<int, int>>& closeendpointsegments, double maxactivedist);

    /*!
    \brief Find contact point via closest point projection
    */
    bool closest_point_projection(double& eta_left1, double& eta_left2, double& l1, double& l2,
        CORE::LINALG::Matrix<3, 1, double>& segmentdata, std::pair<TYPE, TYPE>& solutionpoints,
        int segid1, int segid2);

    /*!
    \brief Find closest point eta2_master on a line for a given slave point eta1_slave
    */
    bool point_to_line_projection(double& eta1_slave, double& eta_left2, double& l2,
        double& eta2_master, double& gap, double& alpha, bool& pairactive, bool smallanglepair,
        bool invertpairs = false, bool orthogonalprojection = false);

    /*!
    \brief Determine minimal distance and contact angle for unconverged segment pair
    */
    void check_unconverged_segment_pair(double& eta_left1, double& eta_left2, double& l1,
        double& l2, double& eta1_min, double& eta2_min, double& g_min, double& alpha_g_min,
        bool& pointtolinesolfound);

    /*!
    \brief Subdivide elements into segments for CPP
    */
    double CreateSegments(const DRT::Element* ele,
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& endpoints_final, int& numsegment, int i);

    /*!
    \brief Get maximal gap at which a contact can become active
    */
    double GetMaxActiveDist();

    /*!
    \brief Check, if segments are fine enough
    */
    bool CheckSegment(CORE::LINALG::Matrix<3, 1, double>& r1,
        CORE::LINALG::Matrix<3, 1, double>& t1, CORE::LINALG::Matrix<3, 1, double>& r2,
        CORE::LINALG::Matrix<3, 1, double>& t2, CORE::LINALG::Matrix<3, 1, double>& rm,
        double& segdist);

    /*!
    \brief Calculate scalar contact force
    */
    void CalcPenaltyLaw(
        Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues>> variables);

    /*!
    \brief Calculate angle-dependent penalty scale factor for large-angle-contact
    */
    void calc_perp_penalty_scale_fac(
        Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues>> cpvariables,
        CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi, CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const double shiftangle1, const double shiftangle2);

    /*!
     * Todo redundant ?!
    \brief Calculate angle-dependent penalty scale factor for small-angle-contact
    */
    void calc_par_penalty_scale_fac(
        Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues>> gpvariables,
        CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi, CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const double shiftangle1, const double shiftangle2);

    /*!
     \brief Compute contact forces
     */
    void EvaluateFcContact(CORE::LINALG::SerialDenseVector& forcevec1,
        CORE::LINALG::SerialDenseVector& forcevec2, const CORE::LINALG::Matrix<3, 1, TYPE>& r1,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2, const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
        Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues>> variables,
        const double& intfac, bool cpp, bool gp, bool fixedendpointxi, bool fixedendpointeta,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE>* fc1_FAD = nullptr,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE>* fc2_FAD = nullptr);

    /*!
    \brief Evaluate contact stiffness
    */
    void evaluate_stiffc_contact(CORE::LINALG::SerialDenseMatrix& stiffmat11,
        CORE::LINALG::SerialDenseMatrix& stiffmat12, CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22, const CORE::LINALG::Matrix<3, 1, TYPE>& r1,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2, const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xixi,
        Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues>> variables,
        const double& intfac, bool cpp, bool gp, bool fixedendpointxi, bool fixedendpointeta);

#ifdef ENDPOINTSEGMENTATION
    /*!
    \brief FAD-based Evaluation of contact stiffness in case of ENDPOINTSEGMENTATION
    */
    void evaluate_stiffc_contact_int_seg(CORE::LINALG::SparseMatrix& stiffmatrix,
        const CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi_bound,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1, const CORE::LINALG::Matrix<3, 1, TYPE>& r2,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
        Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues>> cpvariables,
        const double& intfac, const double& d_xi_ele_d_xi_bound, TYPE signed_jacobi_interval);
#endif

    /*!
    \brief Linearizations of contact point
    */
    void compute_lin_xi_and_lin_eta(
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi);

    /*!
    \brief Lin. of contact point coordinate eta with fixed xi
    */
    void ComputeLinEtaFixXi(
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi);

    /*!
    \brief Lin. of contact point coordinate xi with fixed eta
    */
    void ComputeLinXiFixEta(
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi);

    /*!
    \brief Compute linearization of integration interval bounds (necessary in case of
    ENDPOINTSEGMENTATION)
    */
    void ComputeLinXiBound(
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi_bound,
        TYPE& eta1_bound, TYPE eta2);

    /*!
    \brief Compute linearization of gap
    */
    void ComputeLinGap(CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_gap,
        const CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
        const CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r, const TYPE& norm_delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2);

    /*!
    \brief Compute linearization of cosine of contact angle
    */
    void compute_lin_cos_contact_angle(
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_coscontactangle,
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi);

    /*!
    \brief Compute linearization of normal vector
    */
    void ComputeLinNormal(
        CORE::LINALG::Matrix<3, 2 * 3 * numnodes * numnodalvalues, TYPE>& delta_normal,
        const CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
        const CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2);

    /*!
    \brief Calculate shape function values for given parameter values
    */
    void GetShapeFunctions(CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xixi,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xixi, const TYPE& eta1,
        const TYPE& eta2);

    /*!
    \brief Calculate one specified shape function value / derivative for given parameter value and
    element
    */
    void GetShapeFunctions(CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N,
        const TYPE& eta, int deriv, const DRT::Element* ele) const;

    /*!
    \brief compute coordinate at given curve point
    */
    CORE::LINALG::Matrix<3, 1, TYPE> r(const TYPE& eta, const DRT::Element* ele) const;

    /*!
    \brief compute derivative at given curve point
    */
    CORE::LINALG::Matrix<3, 1, TYPE> r_xi(const TYPE& eta, const DRT::Element* ele);

    /*!
    \brief Compute coordinates and their derivatives from the discretization
    */
    void compute_coords_and_derivs(CORE::LINALG::Matrix<3, 1, TYPE>& r1,
        CORE::LINALG::Matrix<3, 1, TYPE>& r2, CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi, CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xixi);

    /*!
    \brief Compute coordinates of contact points of last time step from the discretization
    */
    void compute_old_coords_and_derivs(CORE::LINALG::Matrix<3, 1, TYPE>& r1_old,
        CORE::LINALG::Matrix<3, 1, TYPE>& r2_old, CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi_old,
        CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi_old,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi);

    /*!
    \brief Utility method for CPP (evaluate nonlinear function f)
    */
    void evaluate_orthogonality_condition(CORE::LINALG::Matrix<2, 1, TYPE>& f,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r, const double norm_delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi, const CORE::LINALG::Matrix<3, 1, TYPE>& t1,
        const CORE::LINALG::Matrix<3, 1, TYPE>& t2);

    /*!
    \brief Utility method for CPP (evaluate Jacobian of nonlinear function f)
    */
    void evaluate_lin_orthogonality_condition(CORE::LINALG::Matrix<2, 2, TYPE>& df,
        CORE::LINALG::Matrix<2, 2, TYPE>& dfinv, const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const double norm_delta_r, const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi, const CORE::LINALG::Matrix<3, 1, TYPE>& t1,
        const CORE::LINALG::Matrix<3, 1, TYPE>& t2, const CORE::LINALG::Matrix<3, 1, TYPE>& t1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& t2_xi, bool& elementscolinear);

    /*!
    \brief Evaluate orthogonality cond. of point to line projeciton
    */
    void evaluate_ptl_orthogonality_condition(TYPE& f,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r, const double norm_delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi, bool orthogonalprojection);

    /*!
    \brief Evaluate Jacobian df of PTLOrthogonalityCondition
    */
    bool evaluate_lin_ptl_orthogonality_condition(TYPE& df,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r, const double norm_delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi, bool orthogonalprojection);

    /*!
    \brief Compute normal vector and gap function at contact point
    */
    void ComputeNormal(CORE::LINALG::Matrix<3, 1, TYPE>& r1, CORE::LINALG::Matrix<3, 1, TYPE>& r2,
        CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi, CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues>> variables,
        int contacttype);

    /*!
    \brief Check, if we have contact or not (e.g. gap < gmax [e.g. gmax=0]?)
    */
    bool CheckContactStatus(const double& gap);

    /*!
      \brief Get jacobi factor of beam element
    */
    double GetJacobi(const DRT::Element* element1);

    /** \brief get Jacobi factor of beam element at xi \in [-1;1]
     *
     *  \author grill
     *  \date 06/16 */
    inline double GetJacobiAtXi(const DRT::Element* element1, const double& xi)
    {
      const DRT::ELEMENTS::Beam3Base* ele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(element1);

      if (ele == nullptr) FOUR_C_THROW("Dynamic cast to Beam3Base failed");

      return ele->GetJacobiFacAtXi(xi);
    }

    /*!
      \brief clear class variables at the beginning of a Newton step
    */
    void ClearClassVariables();

    /*!
      \brief Linearization-check of coordinates xi and eta via FAD
    */
    void fad_check_lin_xi_and_lin_eta(const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi);

    /*!
      \brief Linearization-check for local Newton in CPP via FAD
    */
    void fad_check_lin_orthogonality_condition(const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const double& norm_delta_r, const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi, const CORE::LINALG::Matrix<3, 1, TYPE>& t1,
        const CORE::LINALG::Matrix<3, 1, TYPE>& t2);

    //  /*!
    //    \brief FD-Check of stiffness matrix
    //  */
    //  void FDCheck( CORE::LINALG::SparseMatrix& stiffmatrix,
    //                Epetra_Vector& fint,
    //                const double& pp,
    //                std::map<std::pair<int,int>, Teuchos::RCP<BeamContactPair> >& contactpairmap,
    //                Teuchos::ParameterList& timeintparams,
    //                bool fdcheck);

    //@}
  };
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
