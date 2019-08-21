/*----------------------------------------------------------------------------*/
/*!

\brief One beam and solid contact pair (two elements)

\level 3

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------------*/

#ifndef BEAM3TSC_H
#define BEAM3TSC_H

#include "../drt_beaminteraction/beam3contact_utils.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include <Sacado.hpp>

#include "../drt_beaminteraction/beam3contact_defines.H"
#include "../drt_beaminteraction/beam3contact_tangentsmoothing.H"
#include "../headers/FAD_utils.H"

namespace CONTACT
{
  /*!
   \brief contact element for contact between a 3D beam end a 2D surface (belonging to a 3D solid)
   element

   \author meier (meier@lnm.mw.tum.de)
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
    virtual ~Beam3tosolidcontactinterface() {}

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
    virtual const Epetra_SerialDenseVector GetX1() = 0;

    virtual const Epetra_SerialDenseVector GetX2() = 0;

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
        LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint, const double& pp) = 0;

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
    virtual void UpdateElePos(
        Epetra_SerialDenseMatrix& newele1pos, Epetra_SerialDenseMatrix& newele2pos) = 0;

    /*
    \brief Update interpolated nodal tangents for tangent smoothing
    */
    virtual void UpdateEleSmoothTangents(std::map<int, LINALG::Matrix<3, 1>>& currentpositions) = 0;

    //! brief Struct for debug data in Gmsh
    struct gmshDebugPoint
    {
      LINALG::Matrix<3, 1, double> r1;
      LINALG::Matrix<3, 1, double> x2;
      LINALG::Matrix<3, 1, double> n2;
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

    /*!
    \brief Destructor
    */
    virtual ~Beam3tosolidcontact(){};
    //@}


    //! @name Access methods
    /*!
    \brief Get problem discretization
    */
    inline const DRT::Discretization& ProblemDiscret() const { return pdiscret_; };

    /*!
    \brief Get beam contact discretization
    */
    inline const DRT::Discretization& ContactDiscret() const { return cdiscret_; };

    /*!
    \brief Get offset of dofs between cdiscret and pdiscret
    */
    inline const std::map<int, int>& DofOffset() const { return dofoffsetmap_; };

    /*!
    \brief Get first element
    */
    inline const DRT::Element* Element1() { return element1_; };

    /*!
    \brief Get first element
    */
    inline const DRT::Element* Element2() { return element2_; };

    /*!
    \brief Get gap of this contact pair
    */
    double GetGap() { return FADUTILS::CastToDouble(gap_); };

    /*!
    \brief Get flag indicating whether contact is active (true) or inactive (false)
    */
    bool GetContactFlag() { return contactflag_; };

    /*!
    \brief Get coordinates of contact point of element1 and element2
    */
    const Epetra_SerialDenseVector GetX1()
    {
      Epetra_SerialDenseVector r1;
      r1.Resize(3);
      for (int i = 0; i < 3; i++) r1(i) = FADUTILS::CastToDouble(r1_(i));

      return r1;
    };

    const Epetra_SerialDenseVector GetX2()
    {
      Epetra_SerialDenseVector r2;
      r2.Resize(3);
      for (int i = 0; i < 3; i++) r2(i) = FADUTILS::CastToDouble(r2_(i));

      return r2;
    };

    /*!
    \brief Get flag indicating whether the nodal values of one element had been shifted due to r1=r2
    */
    bool GetShiftStatus() { return shiftnodalvalues_; };

    /*!
      \Check, if there is a difference between the result of the new and old gap definition, i.e. if
      the beams centerlines have already crossed or not.
    */
    bool GetNewGapStatus();
    //@}


    //! @name Public evaluation methods
    /*!
    \brief Evaluate this contact element pair
    */
    bool Evaluate(LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint, const double& pp);

    /*!
    \brief Change the sign of the normal vector: This has to be done at the end of a time step when
    the remainig penetration is larger that the sum of the beam radii (R1+R2). Otherwise, the beams
    could cross in the next time step when the new gap function definition (ngf_=true) for slender
    beams is applied!
    */
    void InvertNormal();

    /*!
      \brief Update of class variables at the end of a time step
    */
    void UpdateClassVariablesStep();

    /*!
      \brief Shift current normal vector to old normal vector at the end of a time step. This is
      necessary when the new gap function definition (ngf_=true) for slender beams is applied!
    */
    void ShiftNormal();

    /*
    \brief Update nodal coordinates of both elements at the beginning of a new time step!
    */
    void UpdateElePos(Epetra_SerialDenseMatrix& newele1pos, Epetra_SerialDenseMatrix& newele2pos);

    /*
    \brief Update interpolated nodal tangents for tangent smoothing
    */
    void UpdateEleSmoothTangents(std::map<int, LINALG::Matrix<3, 1>>& currentpositions);

    /*!
    \brief Get debugging data at Gauss points for Gmsh
    */
    std::vector<gmshDebugPoint> GetGmshDebugPoints() { return gmshDebugPoints_; };


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
    LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPEBTS> ele1pos_;
    LINALG::Matrix<3 * numnodessol, 1, TYPEBTS> ele2pos_;

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
    LINALG::Matrix<3, 1, TYPEBTS> r1_;
    LINALG::Matrix<3, 1, TYPEBTS> r2_;

    //! parameter values of contact point
    TYPEBTS xi1_;
    TYPEBTS xi2_;

    //! Vector containing pairs of unit distance vector nD and beam parameter eta of current time
    //! step
    std::vector<std::pair<TYPEBTS, LINALG::Matrix<3, 1, TYPEBTS>>> normalsets_;

    //! Vector containing pairs of unit distance vector nD and beam parameter eta of last time step
    std::vector<std::pair<TYPEBTS, LINALG::Matrix<3, 1, TYPEBTS>>> normalsets_old_;

    //! neighbor elements of element 1
    Teuchos::RCP<B3CNeighbor> neighbors1_;

    //! neighbor elements of element 2
    Teuchos::RCP<B3CNeighbor> neighbors2_;

    //! averaged nodal tangents, necessary for smoothed tangent fields of C^0 Reissner beams
    LINALG::Matrix<3 * numnodes, 1> nodaltangentssmooth1_;
    LINALG::Matrix<3 * numnodes, 1> nodaltangentssmooth2_;

    //! Comparator for comparing the beam parameter of two parameter sets
    static bool CompareParsets(
        const std::pair<LINALG::Matrix<3, 1, TYPEBTS>, LINALG::Matrix<2, 1, int>>& lhs,
        const std::pair<LINALG::Matrix<3, 1, TYPEBTS>, LINALG::Matrix<2, 1, int>>& rhs)
    {
      // Compare eta
      return lhs.first(2) < rhs.first(2);
    }

    //! Comparator for comparing the beam parameter of two normal sets
    static bool CompareNormalsets(
        const std::pair<TYPEBTS, LINALG::Matrix<3, 1, TYPEBTS>>& lhs, const TYPEBTS& rhs)
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
        const std::pair<LINALG::Matrix<3, 1, TYPEBTS>, int>& parset_a,
        const std::pair<LINALG::Matrix<3, 1, TYPEBTS>, int>& parset_b,
        LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPEBTS>& fc1,
        LINALG::Matrix<3 * numnodessol, 1, TYPEBTS>& fc2,
        LINALG::Matrix<3 * numnodes * numnodalvalues,
            3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& stiffc1,
        LINALG::Matrix<3 * numnodessol, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>&
            stiffc2,
        bool& doAssembleContactInterval,
        LINALG::Matrix<3 * numnodes * numnodalvalues,
            3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& stiffc1_FAD,
        LINALG::Matrix<3 * numnodessol, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>&
            stiffc2_FAD);

    /*!
    \brief Evaluate penalty force law for different regularizations
    */
    void EvaluatePenaltyForceLaw(const double& pp, const TYPEBTS& gap, TYPEBTS& fp, TYPEBTS& dfp);

    /*!
    \brief Evaluate contact forces
    */
    void EvaluateFcContact(const TYPEBTS& fp,
        LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPEBTS>& fc1,
        LINALG::Matrix<3 * numnodessol, 1, TYPEBTS>& fc2, const TYPEBTS& eta_a,
        const TYPEBTS& eta_b, const double& w_gp, const double& sgn,
        const LINALG::Matrix<3, 1, TYPEBTS>& nD, const LINALG::Matrix<3, 1, TYPEBTS>& n2,
        const LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N1,
        const LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2, const double& jacobi);

    /*!
    \brief Evaluate contact stiffness
    */
    void EvaluateStiffcContact(const TYPEBTS& fp, const TYPEBTS& dfp, const TYPEBTS& gap,
        LINALG::Matrix<3 * numnodes * numnodalvalues,
            3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& stiffc1,
        LINALG::Matrix<3 * numnodessol, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>&
            stiffc2,
        const double& sgn, const TYPEBTS& eta_a, const TYPEBTS& eta_b,
        LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& eta_d,
        LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& eta_a_d,
        LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& eta_b_d,
        const double& w_gp, const LINALG::Matrix<3, 1, TYPEBTS>& rD, const TYPEBTS& norm_rD,
        const LINALG::Matrix<3, 1, TYPEBTS>& nD, const LINALG::Matrix<3, 1, TYPEBTS>& a2,
        const TYPEBTS& norm_a2, const LINALG::Matrix<3, 1, TYPEBTS>& n2,
        const LINALG::Matrix<3, 1, TYPEBTS>& r1_eta, const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2, const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1xi1,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2xi2,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1xi2,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2xi1,
        const LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N1,
        const LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N1_eta,
        const LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2,
        const LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2_xi1,
        const LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2_xi2, const double& jacobi,
        LINALG::Matrix<3 * numnodes * numnodalvalues,
            3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& stiffc1_FAD,
        LINALG::Matrix<3 * numnodessol, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>&
            stiffc2_FAD);

    /*!
    \brief Compute linearizations of element parameters xi1, xi2 and eta
    */
    void ComputeLinParameter(const int& fixed_par,
        LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& xi1_d,
        LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& xi2_d,
        LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& eta_d,
        const LINALG::Matrix<3, 1, TYPEBTS>& rD, const LINALG::Matrix<3, 1, TYPEBTS>& r1_eta,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1, const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1xi1,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2xi2,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1xi2,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2xi1,
        const LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N1,
        const LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2,
        const LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2_xi1,
        const LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2_xi2);

    /*!
    \brief Compute linearization of gap
    */
    void ComputeLinGap(
        LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& gap_d,
        const LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& xi1_d,
        const LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& xi2_d,
        const LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& eta_d,
        const double sgn, const LINALG::Matrix<3, 1, TYPEBTS>& rD,
        const LINALG::Matrix<3, 1, TYPEBTS>& nD, const LINALG::Matrix<3, 1, TYPEBTS>& n2,
        const TYPEBTS& norm_rD, const LINALG::Matrix<3, 1, TYPEBTS>& r1_eta,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1, const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2,
        const LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N1,
        const LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2,
        LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& rD_d);

    /*!
    \brief Compute linearization of unit distance vector nD and surface unit normal vector n2
    */
    void ComputeLinNormal(
        LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& nD_d,
        const LINALG::Matrix<3, 1, TYPEBTS>& nD, const TYPEBTS& norm_rD,
        LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& n2_d,
        const LINALG::Matrix<3, 1, TYPEBTS>& n2, const TYPEBTS& norm_a2,
        const LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& rD_d,
        const LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& xi1_d,
        const LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& xi2_d,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1, const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1xi1,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2xi2,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1xi2,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2xi1,
        const LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2_xi1,
        const LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N2_xi2);

    /*!
    \brief Assemble contact forces and stiffness
    */
    void AssembleFcAndStiffcContact(
        const LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPEBTS> fc1,
        const LINALG::Matrix<3 * numnodessol, 1, TYPEBTS> fc2, Epetra_Vector* fint,
        const LINALG::Matrix<3 * numnodes * numnodalvalues,
            3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>
            stiffc1,
        const LINALG::Matrix<3 * numnodessol, 3 * numnodes * numnodalvalues + 3 * numnodessol,
            TYPEBTS>
            stiffc2,
        LINALG::SparseMatrix& stiffmatrix);

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
        std::vector<std::pair<LINALG::Matrix<3, 1, TYPEBTS>, int>>& parsets);

    /*!
    \brief Calculate beam shape function values for given parameter value eta
    */
    void GetBeamShapeFunctions(LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N,
        LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N_eta,
        LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N_etaeta, const TYPEBTS& eta);

    /*!
    \brief Calculate solid surface shape function values for given parameter values xi1 and xi2
    */
    void GetSurfShapeFunctions(LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1xi1,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2xi2,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1xi2,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2xi1, const TYPEBTS& xi1,
        const TYPEBTS& xi2);

    /*!
    \brief Assemble beam shape functions into corresponding matrices
    */
    void AssembleBeamShapefunctions(
        const LINALG::Matrix<1, numnodes * numnodalvalues, TYPEBTS>& N_i,
        const LINALG::Matrix<1, numnodes * numnodalvalues, TYPEBTS>& N_i_eta,
        const LINALG::Matrix<1, numnodes * numnodalvalues, TYPEBTS>& N_i_etaeta,
        LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N,
        LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N_eta,
        LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N_etaeta);

    /*!
    \brief Assemble solid surface shape functions into corresponding matrices
    */
    void AssembleSurfShapefunctions(const LINALG::Matrix<1, numnodessol, TYPEBTS>& N_i,
        const LINALG::Matrix<2, numnodessol, TYPEBTS>& N_i_xi,
        const LINALG::Matrix<3, numnodessol, TYPEBTS>& N_i_xixi,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1xi1,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2xi2,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1xi2,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2xi1);

    /*!
    \brief Compute beam coordinates and their derivatives from the discretization
    */
    void ComputeBeamCoordsAndDerivs(LINALG::Matrix<3, 1, TYPEBTS>& r,
        LINALG::Matrix<3, 1, TYPEBTS>& r_eta, LINALG::Matrix<3, 1, TYPEBTS>& r_etaeta,
        const LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N,
        const LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N_eta,
        const LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPEBTS>& N_etaeta);

    /*!
    \brief Compute solid surface coordinates and their derivatives from the discretization
    */
    void ComputeSurfCoordsAndDerivs(LINALG::Matrix<3, 1, TYPEBTS>& r,
        LINALG::Matrix<3, 1, TYPEBTS>& r_xi1, LINALG::Matrix<3, 1, TYPEBTS>& r_xi2,
        LINALG::Matrix<3, 1, TYPEBTS>& r_xi1xi1, LINALG::Matrix<3, 1, TYPEBTS>& r_xi2xi2,
        LINALG::Matrix<3, 1, TYPEBTS>& r_xi1xi2, LINALG::Matrix<3, 1, TYPEBTS>& r_xi2xi1,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1xi1,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2xi2,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi1xi2,
        LINALG::Matrix<3, 3 * numnodessol, TYPEBTS>& N_xi2xi1);

    /*!
    \brief Compute distance vector rD, its norm norm_rD and unit distance vector nD
    */
    void ComputeDistanceNormal(const LINALG::Matrix<3, 1, TYPEBTS>& r1,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2, LINALG::Matrix<3, 1, TYPEBTS>& rD,
        TYPEBTS& norm_rD, LINALG::Matrix<3, 1, TYPEBTS>& nD);

    /*!
    \brief Compute tangent cross product a2, its norm norm_a2 and surface unit normal vector n2
    */
    void ComputeSurfaceNormal(const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2, LINALG::Matrix<3, 1, TYPEBTS>& a2,
        TYPEBTS& norm_a2, LINALG::Matrix<3, 1, TYPEBTS>& n2);

    /*!
    \brief Utility method for CPP (evaluate nonlinear function f)
    */
    void EvaluateOrthogonalityCondition(LINALG::Matrix<2, 1, TYPEBTS>& f,
        const LINALG::Matrix<3, 1, TYPEBTS>& delta_r, const double norm_delta_r,
        const LINALG::Matrix<3, 1, TYPEBTS>& r1_xi, const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi,
        const LINALG::Matrix<3, 1, TYPEBTS>& t1, const LINALG::Matrix<3, 1, TYPEBTS>& t2);

    /*!
    \brief Utility method for CPP (evaluate Jacobian of nonlinear function f)
    */
    void EvaluateLinOrthogonalityCondition(LINALG::Matrix<2, 2, TYPEBTS>& df,
        LINALG::Matrix<2, 2, TYPEBTS>& dfinv, const LINALG::Matrix<3, 1, TYPEBTS>& delta_r,
        const double norm_delta_r, const LINALG::Matrix<3, 1, TYPEBTS>& r1_xi,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi, const LINALG::Matrix<3, 1, TYPEBTS>& r1_xixi,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xixi, const LINALG::Matrix<3, 1, TYPEBTS>& t1,
        const LINALG::Matrix<3, 1, TYPEBTS>& t2, const LINALG::Matrix<3, 1, TYPEBTS>& t1_xi,
        const LINALG::Matrix<3, 1, TYPEBTS>& t2_xi);

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

    void FADCheckLinParameter(const int& fixed_par, const LINALG::Matrix<3, 1, TYPEBTS>& rD,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1, const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2,
        LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& xi1_d_FAD,
        LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& xi2_d_FAD,
        LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& eta_d_FAD,
        const LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& xi1_d,
        const LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& xi2_d,
        const LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& eta_d);

    void FADCheckLinOrthogonalityCondition(const int& fixed_par,
        const LINALG::Matrix<3, 1, TYPEBTS>& rD, const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi1,
        const LINALG::Matrix<3, 1, TYPEBTS>& x2_xi2, LINALG::Matrix<2, 2, TYPEBTS>& J_FAD,
        const LINALG::Matrix<2, 2, TYPEBTS>& J);

    void FADCheckLinGapAndDistanceVector(const TYPEBTS& gap,
        const LINALG::Matrix<3, 1, TYPEBTS>& rD,
        const LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& xi1_d,
        const LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& xi2_d,
        const LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& eta_d,
        LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& gap_d_FAD,
        LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& rD_d_FAD,
        const LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& gap_d,
        const LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& rD_d);

    void FADCheckLinNormal(const LINALG::Matrix<3, 1, TYPEBTS>& nD,
        const LINALG::Matrix<3, 1, TYPEBTS>& n2,
        const LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& xi1_d,
        const LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& xi2_d,
        const LINALG::Matrix<3 * numnodes * numnodalvalues + 3 * numnodessol, 1, TYPEBTS>& eta_d,
        LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& nD_d_FAD,
        LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& n2_d_FAD,
        const LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& nD_d,
        const LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& n2_d);

    void FDCheckStiffness(const double& pp,
        const LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPEBTS>& fc1,
        const LINALG::Matrix<3 * numnodessol, 1, TYPEBTS>& fc2,
        LINALG::Matrix<3 * numnodes * numnodalvalues,
            3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& stiffc1,
        LINALG::Matrix<3 * numnodessol, 3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>&
            stiffc2);

    void FADCheckStiffness(const LINALG::Matrix<3 * numnodes * numnodalvalues,
                               3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& stiffc1,
        const LINALG::Matrix<3 * numnodessol, 3 * numnodes * numnodalvalues + 3 * numnodessol,
            TYPEBTS>& stiffc2,
        const LINALG::Matrix<3 * numnodes * numnodalvalues,
            3 * numnodes * numnodalvalues + 3 * numnodessol, TYPEBTS>& stiffc1_FAD,
        const LINALG::Matrix<3 * numnodessol, 3 * numnodes * numnodalvalues + 3 * numnodessol,
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

#endif  // #ifndef BEAM3TSC_H
