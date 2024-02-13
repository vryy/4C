/*----------------------------------------------------------------------------*/
/*! \file

\brief One beam and solid contact pair (two elements)

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef BACI_BEAMCONTACT_BEAM3TOSOLIDCONTACT_HPP
#define BACI_BEAMCONTACT_BEAM3TOSOLIDCONTACT_HPP

#include "baci_config.hpp"

#include "baci_beaminteraction_beam_to_beam_contact_defines.hpp"
#include "baci_beaminteraction_beam_to_beam_contact_tangentsmoothing.hpp"
#include "baci_beaminteraction_beam_to_beam_contact_utils.hpp"
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_linalg_sparsematrix.hpp"
#include "baci_utils_fad.hpp"

#include <Sacado.hpp>

BACI_NAMESPACE_OPEN

namespace CONTACT
{
  /*!
   \brief contact element for contact between a 3D beam end a 2D surface (belonging to a 3D solid)
   element

   */

  class Beam3tosolidcontactinterface
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
    Beam3tosolidcontactinterface() {}

    /*!
    \brief Destructor
    */
    virtual ~Beam3tosolidcontactinterface() = default;
    //@}

    //! @name Access methods
    /*!
    \brief Get problem discretization
    */
    virtual const DRT::Discretization& ProblemDiscret() const = 0;

    /*!
    \brief Get beam contact discretization
    */
    virtual const DRT::Discretization& ContactDiscret() const = 0;

    /*!
    \brief Get offset of dofs between cdiscret and pdiscret
    */
    virtual const std::map<int, int>& DofOffset() const = 0;

    /*!
    \brief Get first element
    */
    virtual const DRT::Element* Element1() = 0;
    // inline const DRT::Element* Element1() { return element1_;};

    /*!
    \brief Get first element
    */
    virtual const DRT::Element* Element2() = 0;

    /*!
    \brief Get gap of this contact pair
    */
    virtual double GetGap() = 0;

    /*!
    \brief Get flag ndicating whether contact is active (true) or inactive (false)
    */
    virtual bool GetContactFlag() = 0;

    /*!
    \brief Get coordinates of contact point of element1 and element2
    */
    virtual CORE::LINALG::SerialDenseVector GetX1() = 0;

    virtual CORE::LINALG::SerialDenseVector GetX2() = 0;

    /*!
      \Check, if there is a difference between the result of the new and old gap definition, i.e. if
      the beams centerlines have already crossed or not.
    */
    virtual bool GetNewGapStatus() = 0;

    /*!
    \brief Get flag indicating whether the nodal values of one element had been shifted due to r1=r2
    */
    virtual bool GetShiftStatus() = 0;
    //@}


    //! @name Public evaluation methods
    /*!
    \brief Evaluate this contact element pair
    */
    virtual bool Evaluate(
        CORE::LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint, const double& pp) = 0;

    //! return appropriate internal implementation class (acts as a simple factory)
    static Teuchos::RCP<Beam3tosolidcontactinterface> Impl(const int numnodessol,
        const int numnodes, const int numnodalvalues, const DRT::Discretization& pdiscret,
        const DRT::Discretization& cdiscret, const std::map<int, int>& dofoffsetmap,
        DRT::Element* element1, DRT::Element* element2, Teuchos::ParameterList beamcontactparams);

    /*!
    \brief Change the sign of the normal vector: This has to be done at the end of a time step when
    the remainig penetration is larger that the sum of the beam radii (R1+R2). Otherwise, the beams
    could cross in the next time step when the new gap function definition (ngf_=true) for slender
    beams is applied!
    */
    virtual void InvertNormal() = 0;

    /*!
      \brief Update of class variables at the end of a time step
    */
    virtual void UpdateClassVariablesStep() = 0;

    /*!
      \brief Shift current normal vector to old normal vector at the end of a time step. This is
      necessary when the new gap function definition (ngf_=true) for slender beams is applied!
    */
    virtual void ShiftNormal() = 0;

    /*
    \brief Update nodal coordinates of both elements at the beginning of a new time step!
    */
    virtual void UpdateElePos(CORE::LINALG::SerialDenseMatrix& newele1pos,
        CORE::LINALG::SerialDenseMatrix& newele2pos) = 0;

    /*
    \brief Update interpolated nodal tangents for tangent smoothing
    */
    virtual void UpdateEleSmoothTangents(
        std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions) = 0;

    //! brief Struct for debug data in Gmsh
    struct gmshDebugPoint
    {
      CORE::LINALG::Matrix<3, 1, double> r1;
      CORE::LINALG::Matrix<3, 1, double> x2;
      CORE::LINALG::Matrix<3, 1, double> n2;
      double gap;
      double fp;
      int type;
    };

    /*
    \ brief Get debug data for Gmsh
     */
    virtual std::vector<gmshDebugPoint> GetGmshDebugPoints() = 0;


  };  // class Beam3tosolidcontactinterface



  template <const int numnodessol, const int numnodes, const int numnodalvalues>
  class Beam3tosolidcontact : public Beam3tosolidcontactinterface
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
    Beam3tosolidcontact(const DRT::Discretization& pdiscret, const DRT::Discretization& cdiscret,
        const std::map<int, int>& dofoffsetmap, DRT::Element* element1, DRT::Element* element2,
        Teuchos::ParameterList beamcontactparams);

    /*!
    \brief Copy Constructor
    Makes a deep copy of this contact element pair
    */
    Beam3tosolidcontact(const Beam3tosolidcontact& old);


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
    \brief Get gap of this contact pair
    */
    double GetGap() override { return CORE::FADUTILS::CastToDouble(gap_); };

    /*!
    \brief Get flag indicating whether contact is active (true) or inactive (false)
    */
    bool GetContactFlag() override { return contactflag_; };

    /*!
    \brief Get coordinates of contact point of element1 and element2
    */
    CORE::LINALG::SerialDenseVector GetX1() override
    {
      CORE::LINALG::SerialDenseVector r1;
      r1.resize(3);
      for (int i = 0; i < 3; i++) r1(i) = CORE::FADUTILS::CastToDouble(r1_(i));

      return r1;
    };

    CORE::LINALG::SerialDenseVector GetX2() override
    {
      CORE::LINALG::SerialDenseVector r2;
      r2.resize(3);
      for (int i = 0; i < 3; i++) r2(i) = CORE::FADUTILS::CastToDouble(r2_(i));

      return r2;
    };

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


    //! @name Public evaluation methods
    /*!
    \brief Evaluate this contact element pair
    */
    bool Evaluate(
        CORE::LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint, const double& pp) override;

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

    /*!
      \brief Shift current normal vector to old normal vector at the end of a time step. This is
      necessary when the new gap function definition (ngf_=true) for slender beams is applied!
    */
    void ShiftNormal() override;

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

    /*!
    \brief Get debugging data at Gauss points for Gmsh
    */
    std::vector<gmshDebugPoint> GetGmshDebugPoints() override { return gmshDebugPoints_; };


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

    //! current node coordinates of the two elements
    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPEBTS> ele1pos_;
    CORE::LINALG::Matrix<3 * numnodessol, 1, TYPEBTS> ele2pos_;

    //! variable to check if old or modified gap function
    bool ngf_;

    //! variable to check which smoothing type should be applied
    int smoothing_;

    //! sgn(normal*normal_old)
    double sgn_;

    //! variable to identify first call of a pair (for initializing)
    bool firstcall_;

    //! gap function according to original (ngf_==false) or modified (ngf_==true) definition
    TYPEBTS gap_;

    //! gap function according to original definition
    TYPEBTS gap_original_;

    //! flag indicating contact (active/inactive)
    bool contactflag_;

    //! flag indicating if elements are collinear or not
    bool elementscolinear_;

    //! flag indicating if elements share the same contact point, i.e. r1_=r2_ --> evaluation not
    //! possible
    bool elementscrossing_;

    //! flag indicating if the element nodal positions have been shifted in order to avoid r1_=r2_
    bool shiftnodalvalues_;

    //! coordinates of contact points
    CORE::LINALG::Matrix<3, 1, TYPEBTS> r1_;
    CORE::LINALG::Matrix<3, 1, TYPEBTS> r2_;

    //! parameter values of contact point
    TYPEBTS xi1_;
    TYPEBTS xi2_;

    //! Vector containing pairs of unit distance vector nD and beam parameter eta of current time
    //! step
    std::vector<std::pair<TYPEBTS, CORE::LINALG::Matrix<3, 1, TYPEBTS>>> normalsets_;

    //! Vector containing pairs of unit distance vector nD and beam parameter eta of last time step
    std::vector<std::pair<TYPEBTS, CORE::LINALG::Matrix<3, 1, TYPEBTS>>> normalsets_old_;

    //! neighbor elements of element 1
    Teuchos::RCP<BEAMINTERACTION::B3CNeighbor> neighbors1_;

    //! neighbor elements of element 2
    Teuchos::RCP<BEAMINTERACTION::B3CNeighbor> neighbors2_;

    //! averaged nodal tangents, necessary for smoothed tangent fields of C^0 Reissner beams
    CORE::LINALG::Matrix<3 * numnodes, 1> nodaltangentssmooth1_;
    CORE::LINALG::Matrix<3 * numnodes, 1> nodaltangentssmooth2_;

    //! Comparator for comparing the beam parameter of two parameter sets
    static bool CompareParsets(
        const std::pair<CORE::LINALG::Matrix<3, 1, TYPEBTS>, CORE::LINALG::Matrix<2, 1, int>>& lhs,
        const std::pair<CORE::LINALG::Matrix<3, 1, TYPEBTS>, CORE::LINALG::Matrix<2, 1, int>>& rhs)
    {
      // Compare eta
      return lhs.first(2) < rhs.first(2);
    }

    //! Comparator for comparing the beam parameter of two normal sets
    static bool CompareNormalsets(
        const std::pair<TYPEBTS, CORE::LINALG::Matrix<3, 1, TYPEBTS>>& lhs, const TYPEBTS& rhs)
    {
      // Compare eta
      return lhs.first < rhs;
    }

    //! Vector containing structs for Gmsh debug
    std::vector<gmshDebugPoint> gmshDebugPoints_;


    //@}

    //! @name Private evaluation methods

    /*!
    \brief Evaluate contact forces and stiffness for one contact interval
    */
    void EvaluateContactInterval(const double& pp,
        const std::pair<CORE::LINALG::Matrix<3, 1, TYPEBTS>, int>& parset_a,
        const std::pair<CORE::LINALG::Matrix<3, 1, TYPEBTS>, int>& parset_b,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPEBTS>& fc1,
        CORE::LINALG::Matrix<3 * numnodessol, 1, TYPEBTS>& fc2,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues,
            3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& stiffc1,
        CORE::LINALG::Matrix<3 * numnodessol, 3 * numnodes * numnodalvalues + 3 * numnodessol,
            TYPEBTS>& stiffc2,
        bool& doAssembleContactInterval,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues,
            3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& stiffc1_FAD,
        CORE::LINALG::Matrix<3 * numnodessol, 3 * numnodes * numnodalvalues + 3 * numnodessol,
            TYPEBTS>& stiffc2_FAD);

    /*!
    \brief Evaluate penalty force law for different regularizations
    */
    void EvaluatePenaltyForceLaw(const double& pp, const TYPEBTS& gap, TYPEBTS& fp, TYPEBTS& dfp);

    /*!
    \brief Evaluate contact forces
    */
    void EvaluateFcContact(const TYPEBTS& fp,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPEBTS>& fc1,
        CORE::LINALG::Matrix<3 * numnodessol, 1, TYPEBTS>& fc2, const TYPEBTS& eta_a,
        const TYPEBTS& eta_b, const double& w_gp, const double& sgn,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& nD,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& n2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2, const double& jacobi);

    /*!
    \brief Evaluate contact stiffness
    */
    void EvaluateStiffcContact(const TYPEBTS& fp, const TYPEBTS& dfp, const TYPEBTS& gap,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues,
            3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& stiffc1,
        CORE::LINALG::Matrix<3 * numnodessol, 3 * numnodes * numnodalvalues + 3 * numnodessol,
            TYPEBTS>& stiffc2,
        const double& sgn, const TYPEBTS& eta_a, const TYPEBTS& eta_b,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& eta_d,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& eta_a_d,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& eta_b_d,
        const double& w_gp, const CORE::LINALG::Matrix<3, 1, TYPEBTS>& rD, const TYPEBTS& norm_rD,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& nD,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& a2, const TYPEBTS& norm_a2,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& n2,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& r1_eta,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1xi1,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2xi2,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1xi2,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2xi1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N1_eta,
        const CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2_xi1,
        const CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2_xi2, const double& jacobi,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues,
            3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& stiffc1_FAD,
        CORE::LINALG::Matrix<3 * numnodessol, 3 * numnodes * numnodalvalues + 3 * numnodessol,
            TYPEBTS>& stiffc2_FAD);

    /*!
    \brief Compute linearizations of element parameters xi1, xi2 and eta
    */
    void ComputeLinParameter(const int& fixed_par,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& xi1_d,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& xi2_d,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& eta_d,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& rD,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& r1_eta,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1xi1,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2xi2,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1xi2,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2xi1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2_xi1,
        const CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2_xi2);

    /*!
    \brief Compute linearization of gap
    */
    void ComputeLinGap(
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& gap_d,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            xi1_d,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            xi2_d,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            eta_d,
        const double sgn, const CORE::LINALG::Matrix<3, 1, TYPEBTS>& rD,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& nD,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& n2, const TYPEBTS& norm_rD,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& r1_eta,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& rD_d);

    /*!
    \brief Compute linearization of unit distance vector nD and surface unit normal vector n2
    */
    void ComputeLinNormal(
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& nD_d,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& nD, const TYPEBTS& norm_rD,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& n2_d,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& n2, const TYPEBTS& norm_a2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>&
            rD_d,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            xi1_d,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            xi2_d,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1xi1,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2xi2,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1xi2,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2xi1,
        const CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2_xi1,
        const CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2_xi2);

    /*!
    \brief Assemble contact forces and stiffness
    */
    void AssembleFcAndStiffcContact(
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPEBTS> fc1,
        const CORE::LINALG::Matrix<3 * numnodessol, 1, TYPEBTS> fc2, Epetra_Vector* fint,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues,
            3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>
            stiffc1,
        const CORE::LINALG::Matrix<3 * numnodessol, 3 * numnodes * numnodalvalues + 3 * numnodessol,
            TYPEBTS>
            stiffc2,
        CORE::LINALG::SparseMatrix& stiffmatrix);

    /*!
    \brief Find projection of surface edges on beam and projection of beam center line on surface
    (CPP)
    */
    void Projection(
        const int& fixed_par, TYPEBTS& xi1, TYPEBTS& xi2, TYPEBTS& eta, bool& proj_allowed);

    /*!
    \brief Find contact interval borders
    */
    void GetContactIntervalBorders(
        std::vector<std::pair<CORE::LINALG::Matrix<3, 1, TYPEBTS>, int>>& parsets);

    /*!
    \brief Calculate beam shape function values for given parameter value eta
    */
    void GetBeamShapeFunctions(CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N_eta,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N_etaeta,
        const TYPEBTS& eta);

    /*!
    \brief Calculate solid surface shape function values for given parameter values xi1 and xi2
    */
    void GetSurfShapeFunctions(CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1xi1,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2xi2,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1xi2,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2xi1, const TYPEBTS& xi1,
        const TYPEBTS& xi2);

    /*!
    \brief Assemble beam shape functions into corresponding matrices
    */
    void AssembleBeamShapefunctions(
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPEBTS>& N_i,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPEBTS>& N_i_eta,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPEBTS>& N_i_etaeta,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N_eta,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N_etaeta);

    /*!
    \brief Assemble solid surface shape functions into corresponding matrices
    */
    void AssembleSurfShapefunctions(const CORE::LINALG::Matrix<1, numnodessol, TYPEBTS>& N_i,
        const CORE::LINALG::Matrix<2, numnodessol, TYPEBTS>& N_i_xi,
        const CORE::LINALG::Matrix<3, numnodessol, TYPEBTS>& N_i_xixi,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1xi1,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2xi2,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1xi2,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2xi1);

    /*!
    \brief Compute beam coordinates and their derivatives from the discretization
    */
    void ComputeBeamCoordsAndDerivs(CORE::LINALG::Matrix<3, 1, TYPEBTS>& r,
        CORE::LINALG::Matrix<3, 1, TYPEBTS>& r_eta, CORE::LINALG::Matrix<3, 1, TYPEBTS>& r_etaeta,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N_eta,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N_etaeta);

    /*!
    \brief Compute solid surface coordinates and their derivatives from the discretization
    */
    void ComputeSurfCoordsAndDerivs(CORE::LINALG::Matrix<3, 1, TYPEBTS>& r,
        CORE::LINALG::Matrix<3, 1, TYPEBTS>& r_xi1, CORE::LINALG::Matrix<3, 1, TYPEBTS>& r_xi2,
        CORE::LINALG::Matrix<3, 1, TYPEBTS>& r_xi1xi1,
        CORE::LINALG::Matrix<3, 1, TYPEBTS>& r_xi2xi2,
        CORE::LINALG::Matrix<3, 1, TYPEBTS>& r_xi1xi2,
        CORE::LINALG::Matrix<3, 1, TYPEBTS>& r_xi2xi1,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1xi1,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2xi2,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1xi2,
        CORE::LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2xi1);

    /*!
    \brief Compute distance vector rD, its norm norm_rD and unit distance vector nD
    */
    void ComputeDistanceNormal(const CORE::LINALG::Matrix<3, 1, TYPEBTS>& r1,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2, CORE::LINALG::Matrix<3, 1, TYPEBTS>& rD,
        TYPEBTS& norm_rD, CORE::LINALG::Matrix<3, 1, TYPEBTS>& nD);

    /*!
    \brief Compute tangent cross product a2, its norm norm_a2 and surface unit normal vector n2
    */
    void ComputeSurfaceNormal(const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2, CORE::LINALG::Matrix<3, 1, TYPEBTS>& a2,
        TYPEBTS& norm_a2, CORE::LINALG::Matrix<3, 1, TYPEBTS>& n2);

    /*!
    \brief Utility method for CPP (evaluate nonlinear function f)
    */
    void EvaluateOrthogonalityCondition(CORE::LINALG::Matrix<2, 1, TYPEBTS>& f,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& delta_r, const double norm_delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& t1,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& t2);

    /*!
    \brief Utility method for CPP (evaluate Jacobian of nonlinear function f)
    */
    void EvaluateLinOrthogonalityCondition(CORE::LINALG::Matrix<2, 2, TYPEBTS>& df,
        CORE::LINALG::Matrix<2, 2, TYPEBTS>& dfinv,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& delta_r, const double norm_delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& t1,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& t2,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& t1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& t2_xi);

    /*!
    \brief Check, if we have contact or not
    */
    void CheckContactStatus(const double& pp, const TYPEBTS& gap, bool& contactflag);

    /*!
    \brief These method shifts the nodal positions applied within the beam contact framework py a
    small pre-defined amount in order to enable contact evaluation in the case of two identical
    contact points, i.e r1=r2
    */
    void ShiftNodalPositions();

    void FADCheckLinParameter(const int& fixed_par, const CORE::LINALG::Matrix<3, 1, TYPEBTS>& rD,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            xi1_d_FAD,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            xi2_d_FAD,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            eta_d_FAD,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            xi1_d,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            xi2_d,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            eta_d);

    void FADCheckLinOrthogonalityCondition(const int& fixed_par,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& rD,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2,
        CORE::LINALG::Matrix<2, 2, TYPEBTS>& J_FAD, const CORE::LINALG::Matrix<2, 2, TYPEBTS>& J);

    void FADCheckLinGapAndDistanceVector(const TYPEBTS& gap,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& rD,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            xi1_d,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            xi2_d,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            eta_d,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            gap_d_FAD,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& rD_d_FAD,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            gap_d,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>&
            rD_d);

    void FADCheckLinNormal(const CORE::LINALG::Matrix<3, 1, TYPEBTS>& nD,
        const CORE::LINALG::Matrix<3, 1, TYPEBTS>& n2,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            xi1_d,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            xi2_d,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>&
            eta_d,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& nD_d_FAD,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& n2_d_FAD,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>&
            nD_d,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>&
            n2_d);

    void FDCheckStiffness(const double& pp,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPEBTS>& fc1,
        const CORE::LINALG::Matrix<3 * numnodessol, 1, TYPEBTS>& fc2,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues,
            3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& stiffc1,
        CORE::LINALG::Matrix<3 * numnodessol, 3 * numnodes * numnodalvalues + 3 * numnodessol,
            TYPEBTS>& stiffc2);

    void FADCheckStiffness(const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues,
                               3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& stiffc1,
        const CORE::LINALG::Matrix<3 * numnodessol, 3 * numnodes * numnodalvalues + 3 * numnodessol,
            TYPEBTS>& stiffc2,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues,
            3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& stiffc1_FAD,
        const CORE::LINALG::Matrix<3 * numnodessol, 3 * numnodes * numnodalvalues + 3 * numnodessol,
            TYPEBTS>& stiffc2_FAD);


    /*!
    \brief Get global dofs of a node

    Internally this method first extracts the dofs of the given node
    in the beam contact discretization (which has its own dofs) and
    then transfers these dofs to their actual GIDs in the underlying
    problem discretization by applying the pre-computed dofoffset_.
    */
    std::vector<int> GetGlobalDofs(const DRT::Node* node);

    //@}

  };  // class Beam3tosolidcontact
}  // namespace CONTACT

BACI_NAMESPACE_CLOSE

#endif  // BEAMCONTACT_BEAM3TOSOLIDCONTACT_H
