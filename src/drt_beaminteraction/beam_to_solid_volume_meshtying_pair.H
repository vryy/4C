/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_to_solid_volume_meshtying_pair.H

\brief meshtying element for meshtying between a 3D beam and a 3D solid element

\level 3

\maintainer Alexander Popp
*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef BEAMTOSOLVOLMT_H
#define BEAMTOSOLVOLMT_H

#include "beam_contact_pair.H"

#include "beam3contact_defines.H"

#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"


// forward declaration ...
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace LINALG
namespace BEAMINTERACTION
{
  template <unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
  class BeamToSolidVolumeMeshtyingPair : public BeamContactPair
  {
   public:
    //! @name Constructors and destructors and related methods
    /*!
    \brief Standard Constructor
    */
    BeamToSolidVolumeMeshtyingPair();

    /*!
    \brief Destructor
    */
    virtual ~BeamToSolidVolumeMeshtyingPair(){};

    //! Setup
    void Setup() override;
    //@}


    //! @name Derived methods from base class

    /*!
    \brief things that need to be done in a separate loop before the actual evaluation loop
           over all contact pairs
    */
    void PreEvaluate() override;

    /** \brief Evaluate this contact element pair
     *
     * \author popp
     * \date 05/16 */
    bool Evaluate(LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
        LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
        LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22) override;

    /*!
     \brief Update state of translational nodal DoFs (absolute positions and tangents) of both
     elements
     */
    void ResetState(const std::vector<double>& beam_centerline_dofvec,
        const std::vector<double>& solid_nodal_dofvec) override;

    /*!
     \brief print information about this beam contact element pair to screen
    */
    void Print(std::ostream& out) const override;

    /*!
     \brief print this beam contact element pair to screen
    */
    void PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const override;

    /*!
    \brief Get flag indicating whether contact is active (true) or inactive (false)
    */
    inline bool GetContactFlag() const override
    {
      // The element pair is assumed to be active when we have at least one active contact point
      return contactflag_;
    };

    /*!
    \brief Get number of active contact point pairs on this element pair
    */
    unsigned int GetNumAllActiveContactPointPairs() const override
    {
      dserror("not yet implemented!");
      return 0;
    };

    /*!
    \brief Get coordinates of all active contact points on element1 and element2
    */
    void GetAllActiveContactPointCoordsElement1(
        std::vector<LINALG::TMatrix<double, 3, 1>>& coords) const override
    {
      dserror("not yet implemented!");
    }

    void GetAllActiveContactPointCoordsElement2(
        std::vector<LINALG::TMatrix<double, 3, 1>>& coords) const override
    {
      dserror("not yet implemented!");
    }

    /*!
    \brief Get all (scalar) contact forces of this contact pair
    */
    void GetAllActiveContactForces(std::vector<double>& forces) const override
    {
      dserror("not yet implemented!");
    }

    /*!
    \brief Get all (scalar) gap values of this contact pair
    */
    void GetAllActiveContactGaps(std::vector<double>& gaps) const override
    {
      dserror("not yet implemented!");
    }

    /*!
    \brief Get energy of penalty contact.
    */
    double GetEnergy() const override
    {
      dserror("not implemented yet!");
      return 0.0;
    }
    //@}

   private:
    //! @name Private evaluation methods

    /*!
    \brief Compute reference beam coordinates and their derivatives from the discretization
    */
    void ComputeBeamCoordsAndDerivsRef(LINALG::TMatrix<TYPEBTS, 3, 1>& r,
        LINALG::TMatrix<TYPEBTS, 3, 1>& r_eta, LINALG::TMatrix<TYPEBTS, 3, 1>& r_etaeta,
        const LINALG::TMatrix<TYPEBTS, 3, 3 * numnodes * numnodalvalues>& N,
        const LINALG::TMatrix<TYPEBTS, 3, 3 * numnodes * numnodalvalues>& N_eta,
        const LINALG::TMatrix<TYPEBTS, 3, 3 * numnodes * numnodalvalues>& N_etaeta);

    /*!
    \brief Compute beam coordinates and their derivatives from the discretization
    */
    void ComputeBeamCoordsAndDerivsMom(LINALG::TMatrix<TYPEBTS, 3, 1>& r,
        LINALG::TMatrix<TYPEBTS, 3, 1>& r_eta, LINALG::TMatrix<TYPEBTS, 3, 1>& r_etaeta,
        const LINALG::TMatrix<TYPEBTS, 3, 3 * numnodes * numnodalvalues>& N,
        const LINALG::TMatrix<TYPEBTS, 3, 3 * numnodes * numnodalvalues>& N_eta,
        const LINALG::TMatrix<TYPEBTS, 3, 3 * numnodes * numnodalvalues>& N_etaeta);

    /*!
    \brief Get beam Jacobi for given parameter value eta
    */
    void GetBeamJacobi(const TYPEBTS& eta, TYPEBTS& jac);

    /*!
    \brief Assemble beam shape functions into corresponding matrices
    */
    void AssembleBeamShapefunctions(
        const LINALG::TMatrix<TYPEBTS, 1, numnodes * numnodalvalues>& N_i,
        const LINALG::TMatrix<TYPEBTS, 1, numnodes * numnodalvalues>& N_i_eta,
        const LINALG::TMatrix<TYPEBTS, 1, numnodes * numnodalvalues>& N_i_etaeta,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodes * numnodalvalues>& N,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodes * numnodalvalues>& N_eta,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodes * numnodalvalues>& N_etaeta);

    /*!
    \brief Calculate beam shape function values for given parameter value eta
    */
    void GetBeamShapeFunctions(LINALG::TMatrix<TYPEBTS, 3, 3 * numnodes * numnodalvalues>& N,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodes * numnodalvalues>& N_eta,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodes * numnodalvalues>& N_etaeta, const TYPEBTS& eta);

    /*!
    \brief "Project" GP on beam into solid volume element
    */
    void Projection(TYPEBTS& xi1, TYPEBTS& xi2, TYPEBTS& xi3, TYPEBTS& eta, bool& proj_allowed);

    /*!
    \brief Calculate solid shape function values for given parameter values xi1, xi2, xi3
    */
    void GetSolShapeFunctions(LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N_xi1,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N_xi2,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N_xi3, const TYPEBTS& xi1, const TYPEBTS& xi2,
        const TYPEBTS& xi3);

    /*!
    \brief Evaluate boundary segmentation for that pair
    */
    void EvaluateBoundarySegmentation();

    /*!
    \brief Find beam center line on surface
    */
    void SurfProjection(TYPEBTS& xi1, TYPEBTS& xi2, TYPEBTS& xi3, TYPEBTS& eta, const int fixedPar,
        double fixedAt, bool& proj_allowed);

    /*!
    \brief Assemble solid shape functions into corresponding matrices
    */
    void AssembleSolShapefunctions(const LINALG::TMatrix<TYPEBTS, 1, numnodessol>& N_i,
        const LINALG::TMatrix<TYPEBTS, 3, numnodessol>& N_i_xi,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N_xi1,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N_xi2,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N_xi3);

    /*!
    \brief Evaluate meshtying stiffness matrix
    */
    void EvaluateStiffmtMeshtying(LINALG::SerialDenseMatrix& stiffmat11,
        LINALG::SerialDenseMatrix& stiffmat12, LINALG::SerialDenseMatrix& stiffmat21,
        LINALG::SerialDenseMatrix& stiffmat22, const double& w_gp,
        const LINALG::TMatrix<TYPEBTS, 3, 3 * numnodes * numnodalvalues>& N1,
        const LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N2, TYPEBTS& jacobi, const double& pp);

    /*!
    \brief Evaluate meshtying forces
    */
    void EvaluateFmtMeshtying(LINALG::SerialDenseVector& forcevec1,
        LINALG::SerialDenseVector& forcevec2, LINALG::SerialDenseMatrix& stiffmat11,
        LINALG::SerialDenseMatrix& stiffmat12, LINALG::SerialDenseMatrix& stiffmat21,
        LINALG::SerialDenseMatrix& stiffmat22);

    /*!
    \brief Fill the values of the local stiffness matrix into the element stiffness matrix
    */
    void FillElementStiffmtMeshtying(LINALG::SerialDenseMatrix& local_stiffmat11,
        LINALG::SerialDenseMatrix& local_stiffmat12, LINALG::SerialDenseMatrix& local_stiffmat21,
        LINALG::SerialDenseMatrix& local_stiffmat22, LINALG::SerialDenseMatrix& stiffmat11,
        LINALG::SerialDenseMatrix& stiffmat12, LINALG::SerialDenseMatrix& stiffmat21,
        LINALG::SerialDenseMatrix& stiffmat22);

    /*!
    \brief Compute solid coordinates and their derivatives in ref-conf from the discretization
    */
    void ComputeSolCoordsAndDerivsRef(LINALG::TMatrix<TYPEBTS, 3, 1>& r,
        LINALG::TMatrix<TYPEBTS, 3, 1>& r_xi1, LINALG::TMatrix<TYPEBTS, 3, 1>& r_xi2,
        LINALG::TMatrix<TYPEBTS, 3, 1>& r_xi3, LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N_xi1,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N_xi2,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N_xi3);

    /*!
    \brief Compute solid surface coordinates and their derivatives from the discretization
    */
    void ComputeSolCoordsAndDerivsMom(LINALG::TMatrix<TYPEBTS, 3, 1>& r,
        LINALG::TMatrix<TYPEBTS, 3, 1>& r_xi1, LINALG::TMatrix<TYPEBTS, 3, 1>& r_xi2,
        LINALG::TMatrix<TYPEBTS, 3, 1>& r_xi3, LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N_xi1,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N_xi2,
        LINALG::TMatrix<TYPEBTS, 3, 3 * numnodessol>& N_xi3);
    //@}

   private:
    //! @name member variables

    //! current nodal positions (and tangents) of the two elements
    LINALG::TMatrix<TYPEBTS, 3 * numnodes * numnodalvalues, 1> ele1pos_;
    LINALG::TMatrix<TYPEBTS, 3 * numnodessol, 1> ele2pos_;

    //! reference nodal positions (and tangents) of the two elements
    LINALG::TMatrix<TYPEBTS, 3 * numnodes * numnodalvalues, 1> ele1posref_;
    LINALG::TMatrix<TYPEBTS, 3 * numnodessol, 1> ele2posref_;

    //! Beam Jacobi Determinant at Gausspoints
    std::vector<TYPEBTS> beamjacobi_;

    //! Gausspoints in solid
    std::vector<TYPEBTS> xi1_;
    std::vector<TYPEBTS> xi2_;
    std::vector<TYPEBTS> xi3_;

    //! Gausspoints in beam and weights
    std::vector<TYPEBTS> eta_;
    std::vector<TYPEBTS> wgp_;

    //! flag indicating if at least one Gauss point of the beam element projects onto solid
    bool atleastone_;

    //! vector that keeps track of the Gauss points with valid projection of this pair
    std::vector<bool> localGaussPointTracker_;

    //! stores the number of Gauss points for that element
    int n_gp_;

    //! stores the global element id of the beam element
    int ele1gid_;

    //! flag indicating contact (active/inactive)
    bool contactflag_;

    //@}
  };
}  // namespace BEAMINTERACTION

#endif
