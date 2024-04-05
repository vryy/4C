/*----------------------------------------------------------------------*/
/*! \file

\brief contact element for contact between two 3D beam elements

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_BEAMCONTACT_BEAM3CONTACTNEW_HPP
#define FOUR_C_BEAMCONTACT_BEAM3CONTACTNEW_HPP

#include "baci_config.hpp"

#include "baci_beamcontact_beam3contactinterface.hpp"
#include "baci_beaminteraction_beam_to_beam_contact_defines.hpp"
#include "baci_beaminteraction_beam_to_beam_contact_tangentsmoothing.hpp"
#include "baci_lib_element.hpp"
#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_linalg_sparsematrix.hpp"
#include "baci_utils_fad.hpp"

#include <Sacado.hpp>

BACI_NAMESPACE_OPEN


namespace CONTACT
{
  /*!
   \brief contact element for contact between two 3D beam elements

   Refer also to the Semesterarbeit of Matthias Mayr, 2010

   */

  template <const int numnodes, const int numnodalvalues>
  class Beam3contactnew : public Beam3contactinterface
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}

    //! @name Constructors and destructors and related methods
    /*!
    \brief Standard Constructor
    \param pdiscret (in): the problem discretization
    \param cdiscret (in): the beam contact discretization
    \param dofoffset (in): offset of dof between pdiscret and cdiscret
    \param element1 (in): first element of contact pair
    \param element2 (in): second element of contact pair
    \param ele1pos (in): nodal coordinates of first element
    \param ele2pos (in): nodal coordinates of second element
    */
    Beam3contactnew(const DRT::Discretization& pdiscret, const DRT::Discretization& cdiscret,
        const std::map<int, int>& dofoffsetmap, DRT::Element* element1, DRT::Element* element2,
        Teuchos::ParameterList& beamcontactparams);


    //@}

    //! @name Access methods
    /*!
    \brief Get problem discretization
    */
    inline const DRT::Discretization& ProblemDiscret() const override { return pdiscret_; };

    /*!
    \brief Get beam contact discretization
    */
    inline const DRT::Discretization& ContactDiscret() const override { return cdiscret_; };

    /*!
    \brief Get offset of dofs between cdiscret and pdiscret
    */
    inline const std::map<int, int>& DofOffset() const override { return dofoffsetmap_; };

    /*!
    \brief Get first element
    */
    inline const DRT::Element* Element1() override { return element1_; };

    /*!
    \brief Get first element
    */
    inline const DRT::Element* Element2() override { return element2_; };

    /*!
    \brief Get number of contact points on this element pair
    */
    int GetNumCps() override { return 1; };
    int GetNumGps() override { return 0; };
    int GetNumEps() override { return 0; };

    /*!
    \brief Get vector of type declarations (0=closest point contact, 1=gauss point contact, 2= end
    point contact) of all contact pairs
    */
    std::vector<int> GetContactType() override
    {
      std::vector<int> types(1, 0);

      return types;
    };

    /*!
    \brief Get gap of this contact pair
    */
    std::vector<double> GetGap() override
    {
      std::vector<double> gaps(1, CORE::FADUTILS::CastToDouble(gap_));

      return gaps;
    };

    /*!
    \brief Get contact force of this contact pair
    */
    std::vector<double> GetContactForce() override
    {
      std::vector<double> forces(1, CORE::FADUTILS::CastToDouble(fp_));

      return forces;
    };

    /*!
    \brief Get contact angle of this contact pair
    */
    std::vector<double> GetContactAngle() override
    {
      double angle = 0.0;
      double cosangle = CORE::FADUTILS::CastToDouble(tangentproduct_);
      if (cosangle < 1.0)
        angle = acos(cosangle);  // returns an angle \in [0;pi/2] since scalarproduct \in [0;1.0]
      else
        angle = 0;
      std::vector<double> angles(1, angle);

      return angles;
    };

    /*!
    \brief Get closest point of this contact pair
    */
    std::vector<std::pair<double, double>> GetClosestPoint() override
    {
      std::pair<double, double> cp(0.0, 0.0);
      cp.first = CORE::FADUTILS::CastToDouble(xi1_);
      cp.second = CORE::FADUTILS::CastToDouble(xi2_);
      std::vector<std::pair<double, double>> cps(1, cp);

      return cps;
    };

    /*!
    \brief Return number of individual contact segments on element pair
    */
    std::pair<int, int> GetNumSegments() override { return std::make_pair(1, 1); };

    /*!
    \brief Return ids of active segments
    */
    std::vector<std::pair<int, int>> GetSegmentIds() override
    {
      std::vector<std::pair<int, int>> ids(1, std::make_pair(0, 0));

      return ids;
    };

    /*!
    \brief Get flag indicating whether contact is active (true) or inactive (false)
    */
    bool GetContactFlag() override { return contactflag_; };

    /*!
    \brief Get coordinates of contact point of element1
    */
    std::vector<CORE::LINALG::Matrix<3, 1>> GetX1() override
    {
      std::vector<CORE::LINALG::Matrix<3, 1>> r1(1, CORE::LINALG::Matrix<3, 1>(true));

      for (int j = 0; j < 3; j++) r1[0](j) = CORE::FADUTILS::CastToDouble(r1_(j));

      return r1;
    };

    /*!
    \brief Get coordinates of contact point of element2
    */
    std::vector<CORE::LINALG::Matrix<3, 1>> GetX2() override
    {
      std::vector<CORE::LINALG::Matrix<3, 1>> r2(1, CORE::LINALG::Matrix<3, 1>(true));

      for (int j = 0; j < 3; j++) r2[0](j) = CORE::FADUTILS::CastToDouble(r2_(j));

      return r2;
    };

    /*!
    \brief Get normal vector
    */
    CORE::LINALG::SerialDenseVector GetNormal() override
    {
      CORE::LINALG::SerialDenseVector normal(3);

      if (GetNewGapStatus() == true)
      {
        for (int i = 0; i < 3; i++) normal(i) = -CORE::FADUTILS::CastToDouble(normal_(i));
      }
      else
      {
        for (int i = 0; i < 3; i++) normal(i) = CORE::FADUTILS::CastToDouble(normal_(i));
      }

      return normal;
    }

    /*!
    \brief Get flag indicating whether the nodal values of one element had been shifted due to r1=r2
    */
    bool GetShiftStatus() override { return shiftnodalvalues_; };

    /*!
      \Check, if there is a difference between the result of the new and old gap definition, i.e. if
      the beams centerlines have already crossed or not.
    */
    bool GetNewGapStatus() override;
    //@}

    /*!
      \Get energy of penalty contact.
    */
    double GetEnergy() override
    {
      TYPE energy = 0.5 * pp_ * gap_ * gap_;
      return CORE::FADUTILS::CastToDouble(energy);
    };

    /*!
      Not needed for beam3contactnew variables
    */
    double GetUnscaledPerpEnergy() override { return 0.0; };

    /*!
      Not needed for beam3contactnew variables
    */
    double GetUnscaledParallelEnergy() override { return 0.0; };

    /*!
    \brief Get normal vector of last time step
    */
    CORE::LINALG::Matrix<3, 1, TYPE>* GetNormalOld() override
    {
      if (firsttimestep_)
      {
        dserror(
            "Vector normal_old_requested but not available in the first time step the pair has "
            "been found: Choose larger search radius!!!");
      }

      // Only deliver normal_old_ when a valid closest point projection for last time step existed!
      if (CORE::FADUTILS::Norm(xi1_old_) < 1.0 + XIETATOL and
          CORE::FADUTILS::Norm(xi2_old_) < 1.0 + XIETATOL and !oldcppunconverged_)
      {
        return &normal_old_;
      }
      else
        return nullptr;
    };

    /*!
    \brief Check, if it is the first time step the element is in contact
    */
    bool FirstTimeStep() override { return firsttimestep_; };
    //@}

    /** \brief print this beam contact element pair to screen
     *
     *  \author grill
     *  \date 05/16 */
    void Print() const override{};

    //! @name Public evaluation methods
    /*!
    \brief Evaluate this contact element pair
    */
    bool Evaluate(CORE::LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint, const double& pp,
        std::map<std::pair<int, int>, Teuchos::RCP<Beam3contactinterface>>& contactpairmap,
        Teuchos::ParameterList& timeintparams, bool fdcheck = false) override;

    /*!
    \brief Change the sign of the normal vector: This has to be done at the end of a time step when
    the remainig penetration is larger that the sum of the beam radii (R1+R2). Otherwise, the beams
    could cross in the next time step when the new gap function definition (ngf_=true) for slender
    beams is applied!
    */
    void InvertNormal() override;

    /*!
      \brief Update of class variables at the end of a time step
    */
    void UpdateClassVariablesStep() override;

    /*
    \brief Update nodal coordinates of both elements at the beginning of a new time step!
    */
    void UpdateElePos(CORE::LINALG::SerialDenseMatrix& newele1pos,
        CORE::LINALG::SerialDenseMatrix& newele2pos) override;

    /*
    \brief Update interpolated nodal tangents for tangent smoothing
    */
    void UpdateEleSmoothTangents(
        std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions) override;


    //@}
   private:
    //! @name member variables

    //! reference to problem discretization
    const DRT::Discretization& pdiscret_;

    //! reference to beam contact discretization
    const DRT::Discretization& cdiscret_;

    //! dof offset between pdiscret and cdiscret
    const std::map<int, int>& dofoffsetmap_;

    //! first element of contact pair
    DRT::Element* element1_;

    //! second element of contact pair
    DRT::Element* element2_;

    //! beam contact parameter list
    Teuchos::ParameterList& bcparams_;

    //! current node coordinates of the two elements
    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE> ele1pos_;
    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE> ele2pos_;

    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, double> ele1pos_old_;
    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, double> ele2pos_old_;

    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, double> ele1pos_lastiter_;
    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, double> ele2pos_lastiter_;

    //! variable to check which smoothing type should be applied
    int smoothing_;

    //! sgn(normal*normal_old)
    double sgn_;

    //! variable to identify first call of a pair within considered time step (for initializing)
    bool firstcallofstep_;

    //! variable to identify the first time step in which the pair was found
    bool firsttimestep_;

    //! gap function according to original (ngf_==false) or modified (ngf_==true) definition
    TYPE gap_;

    //! gap function according to original definition
    TYPE gap_original_;

    //! flag indicating contact (active/inactive)
    bool contactflag_;

    //! flag indicating (active/inactive) damping forces
    bool dampingcontactflag_;

    //! flag indicating if contact was active in the last time step
    bool oldcontactflag_;

    //! flag indicating whether an element has at least once been active during the current time
    //! step
    bool currentlyincontact_;

    //! flag indicating if elements are collinear or not
    bool elementscolinear_;

    //! flag indicating if elements share the same contact point, i.e. r1_=r2_ --> evaluation not
    //! possible
    bool elementscrossing_;

    //! flag indicating if the element nodal positions have been shifted in order to avoid r1_=r2_
    bool shiftnodalvalues_;

    //! coordinates of contact point on center lines of beams
    CORE::LINALG::Matrix<3, 1, TYPE> r1_;
    CORE::LINALG::Matrix<3, 1, TYPE> r2_;
    CORE::LINALG::Matrix<3, 1, TYPE> r1_old_;
    CORE::LINALG::Matrix<3, 1, TYPE> r2_old_;

    //! derivative of contact point coordinates
    CORE::LINALG::Matrix<3, 1, TYPE> r1_xi_;
    CORE::LINALG::Matrix<3, 1, TYPE> r2_xi_;
    CORE::LINALG::Matrix<3, 1, TYPE> r1_xi_old_;
    CORE::LINALG::Matrix<3, 1, TYPE> r2_xi_old_;

    //! parameter values of contact point
    TYPE xi1_;
    TYPE xi2_;

    //! parameter values of contact point of last time step
    double xi1_old_;
    double xi2_old_;

    //! normal vector of current time step
    CORE::LINALG::Matrix<3, 1, TYPE> normal_;

    //! normal vector of last time step
    CORE::LINALG::Matrix<3, 1, TYPE> normal_old_;

    //! neighbor elements of element 1
    Teuchos::RCP<BEAMINTERACTION::B3CNeighbor> neighbors1_;

    //! neighbor elements of element 2
    Teuchos::RCP<BEAMINTERACTION::B3CNeighbor> neighbors2_;

    //! averaged nodal tangents, necessary for smoothed tangent fields of C^0 Reissner beams
    CORE::LINALG::Matrix<3 * numnodes, 1> nodaltangentssmooth1_;
    CORE::LINALG::Matrix<3 * numnodes, 1> nodaltangentssmooth2_;

    //! current penalty parameter
    double pp_;

    //! scalar value of current contact force, standard case: fp_ = -pp_ * gap_
    TYPE fp_;

    //! derivative of scalar value of current contact force, standard case: fp_ = -pp_
    TYPE dfp_;

    //! scalar value of current damping force and its derivative, standard case: fdamp_ = -d_ * fd_
    //! e.g. = -d_ * gap_t =
    TYPE fd_;
    TYPE dfd_;
    TYPE d_;
    TYPE dd_;

    //! current Newton iteration
    int iter_;

    //! current time step
    int numstep_;

    //! time step size
    double dt_;

    //! flag indicating when one beam slides over the end point of the second beam
    bool beamendcontactopened_;

    //! flag indicating that a pair is sorted out due to almost parallelism
    bool beamsalmostparallel_;

    //! flag indicating that closest point projection did not converge
    bool cppunconverged_;

    //! flag indicating that closest point projection of last time step did not converge
    bool oldcppunconverged_;

    //! initial element lengths
    double ele1length_;
    double ele2length_;

    //! Has the pair asked for the vector normal_old_ of an neighbor element during the current
    //! iteration step?
    bool neighbornormalrequired_;

    //! Scalar product of normalized tangents at contact point (needed in order to exclude almost
    //! parallel elements)
    TYPE tangentproduct_;

    //! Cross-section radius of beam 1
    double radius1_;

    //! Cross-section radius of beam 2
    double radius2_;

    //! Additive extrusion value of bounding box applied in contact search
    double searchboxinc_;

    //@}

    //! @name Private evaluation methods

    /*!
    \brief Evaluate contact forces
    */
    void EvaluateFcContact(const double& pp, Epetra_Vector* fint,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE>* fc1_FAD = nullptr,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE>* fc2_FAD = nullptr);

    /*!
    \brief Evaluate contact stiffness
    */
    void EvaluateStiffcContact(const double& pp, const TYPE& norm_delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r, CORE::LINALG::SparseMatrix& stiffmatrix,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1, const CORE::LINALG::Matrix<3, 1, TYPE>& r2,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xixi);

    /*!
    \brief Evaluate algorithmic forces that make convergence more robust / faster
    */
    void EvaluateAlgorithmicForce(const double& pp, Epetra_Vector* fint,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE>* fc1_FAD = nullptr,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE>* fc2_FAD = nullptr);

    /*!
    \brief Evaluate algorithmic stiffnesses that make convergence more robust / faster
    */
    void EvaluateAlgorithmicStiff(const double& pp, const TYPE& norm_delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r, CORE::LINALG::SparseMatrix& stiffmatrix,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1, const CORE::LINALG::Matrix<3, 1, TYPE>& r2,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xixi);

    /*!
    \brief Compute linearizations of contact point parameter coordinates xi and eta
    */
    void ComputeLinXiAndLinEta(
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
    \brief Compute linearization of gap g
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
    \brief Compute linearization of time derivative g_t of gap g
    */
    void ComputeLinGapt(
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_gap_t,
        const CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
        const CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
        const CORE::LINALG::Matrix<3, 2 * 3 * numnodes * numnodalvalues, TYPE>& delta_n,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_old,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_old,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi_old,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi_old);

    /*!
    \brief Compute linearization of normal vector n
    */
    void ComputeLinNormal(
        CORE::LINALG::Matrix<3, 2 * 3 * numnodes * numnodalvalues, TYPE>& delta_normal,
        const CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
        const CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
        const TYPE& norm_delta_r, const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2);

    /*!
    \brief Find contact point via closest point projection
    */
    void ClosestPointProjection();

    /*!
    \brief Calculate scalar contact force and linearization according to chosen penalty law
    */
    void CalcPenaltyLaw();

    /*!
    \brief Calculate scalar damping contribution to contact force and linearization
    */
    void CalcDampingLaw();

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
    \brief Assemble the shape functions into corresponding matrices
    */
    void AssembleShapefunctions(const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N_i,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N_i_xi,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N_i_xixi,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N_xi,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N_xixi);

    /*!
    \brief Compute coordinates and their derivatives from the discretization
    */
    void ComputeCoordsAndDerivs(CORE::LINALG::Matrix<3, 1, TYPE>& r1,
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
    void ComputeOldCoordsAndDerivs(CORE::LINALG::Matrix<3, 1, TYPE>& r1_old,
        CORE::LINALG::Matrix<3, 1, TYPE>& r2_old, CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi_old,
        CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi_old,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi);

    /*!
    \brief Utility method for CPP (evaluate nonlinear function f)
    */
    void EvaluateOrthogonalityCondition(CORE::LINALG::Matrix<2, 1, TYPE>& f,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r, const double norm_delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi, const CORE::LINALG::Matrix<3, 1, TYPE>& t1,
        const CORE::LINALG::Matrix<3, 1, TYPE>& t2);

    /*!
    \brief Utility method for CPP (evaluate Jacobian of nonlinear function f)
    */
    void EvaluateLinOrthogonalityCondition(CORE::LINALG::Matrix<2, 2, TYPE>& df,
        CORE::LINALG::Matrix<2, 2, TYPE>& dfinv, const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const double norm_delta_r, const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi, const CORE::LINALG::Matrix<3, 1, TYPE>& t1,
        const CORE::LINALG::Matrix<3, 1, TYPE>& t2, const CORE::LINALG::Matrix<3, 1, TYPE>& t1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& t2_xi);

    /*!
    \brief Compute normal vector and gap function at contact point
    */
    void ComputeNormal(CORE::LINALG::Matrix<3, 1, TYPE>& delta_r, TYPE& norm_delta_r,
        std::map<std::pair<int, int>, Teuchos::RCP<Beam3contactinterface>>& contactpairmap);

    /*!
    \brief Check, if we have contact or not (gap < 0 ???)
    */
    void CheckContactStatus(const double& pp);

    /*!
      \brief Get global dofs of a node

      Internally this method first extracts the dofs of the given node
      in the beam contact discretization (which has its own dofs) and
      then transfers these dofs to their actual GIDs in the underlying
      problem discretization by applying the pre-computed dofoffset_.
      */
    std::vector<int> GetGlobalDofs(const DRT::Node* node);

    /*!
      \brief These method shifts the nodal positions applied within the beam contact framework by a
      small pre-defined amount in order to enable contact evaluation in the case of two identical
      contact points, i.e r1=r2
    */
    void ShiftNodalPositions();

    /*!
      \brief Get the vector normal_old_ from the neighbor element
    */
    void GetNeighborNormalOld(
        std::map<std::pair<int, int>, Teuchos::RCP<Beam3contactinterface>>& contactpairmap);

    void FADCheckLinXiAndLinEta(const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi);

    void FADCheckLinOrthogonalityCondition(const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi);

    /*!
      \brief Check, if the current contact point is at one of the two ends of a physical beam
    */
    void CheckBoundaryContact();

    /*!
      \brief Get jacobi factor of beam element
    */
    double GetJacobi(DRT::Element* element1);


    /*!
      \brief Update of class variables at the end of a Newton step
    */
    void UpdateClassVariablesIter();

    /*!
      \brief Set class variables at the beginning of a Newton step
    */
    void SetClassVariables(const double& pp, Teuchos::ParameterList timeintparams);

    //@}

  };  // class Beam3contactnew
}  // namespace CONTACT

BACI_NAMESPACE_CLOSE

#endif  // BEAMCONTACT_BEAM3CONTACTNEW_H
