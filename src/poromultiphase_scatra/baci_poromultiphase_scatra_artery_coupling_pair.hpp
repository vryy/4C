/*----------------------------------------------------------------------*/
/*! \file
 \brief one pair consisting of exactly one artery element and one poro-
        multiphase-scatra element which might be tied to each other

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_PAIR_HPP
#define FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_PAIR_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "baci_inpar_bio.hpp"
#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"
#include "baci_porofluidmultiphase_ele_phasemanager.hpp"
#include "baci_porofluidmultiphase_ele_variablemanager.hpp"

#include <Epetra_Vector.h>
#include <Sacado.hpp>
#include <Teuchos_RCP.hpp>

// define Fad object for evaluation
typedef Sacado::Fad::DFad<double> FAD;

// forward declaration
class Epetra_MultiVector;

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Element;
  class Discretization;
}  // namespace DRT
namespace MAT
{
  class MatList;
  class Cnst1dArt;
}  // namespace MAT
namespace CORE
{
  namespace LINALG
  {
    class SerialDenseVector;
    class SerialDenseMatrix;
  }  // namespace LINALG
  namespace UTILS
  {
    class FunctionOfAnything;
  }
}  // namespace CORE
namespace POROMULTIPHASESCATRA
{
  class PoroMultiPhaseScatraArteryCouplingPairBase
  {
   public:
    //! constructor
    PoroMultiPhaseScatraArteryCouplingPairBase() { return; };

    //! destructor
    virtual ~PoroMultiPhaseScatraArteryCouplingPairBase() = default;

    //! Init
    virtual void Init(std::vector<DRT::Element const*> elements,
        const Teuchos::ParameterList& couplingparams,
        const Teuchos::ParameterList& fluidcouplingparams, const std::vector<int>& coupleddofs_cont,
        const std::vector<int>& coupleddofs_art, const std::vector<std::vector<int>>& scale_vec,
        const std::vector<std::vector<int>>& funct_vec, const std::string condnames,
        const double penalty, const std::string couplingtype = "", const int eta_ntp = 0) = 0;

    //! query if pair active
    virtual bool IsActive() = 0;

    //! things that need to be done in a separate loop before the actual evaluation loop
    //! over all coupling pairs
    virtual void PreEvaluate(Teuchos::RCP<Epetra_MultiVector> gp_vector) = 0;

    //! things that need to be done in a separate loop before the actual evaluation loop
    //! over all coupling pairs
    virtual void DeleteUnnecessaryGPs(Teuchos::RCP<Epetra_MultiVector> gp_vector) = 0;

    /*!
     * @brief Evaluate this pair
     *
     * @returns integral of diameter of the segment
     */
    virtual double Evaluate(CORE::LINALG::SerialDenseVector* forcevec1,
        CORE::LINALG::SerialDenseVector* forcevec2, CORE::LINALG::SerialDenseMatrix* stiffmat11,
        CORE::LINALG::SerialDenseMatrix* stiffmat12, CORE::LINALG::SerialDenseMatrix* stiffmat21,
        CORE::LINALG::SerialDenseMatrix* stiffmat22, CORE::LINALG::SerialDenseMatrix* D_ele,
        CORE::LINALG::SerialDenseMatrix* M_ele, CORE::LINALG::SerialDenseVector* Kappa_ele,
        const std::vector<double>& segmentlengths) = 0;

    //! evaluate additional linearization of (integrated) element diameter dependent terms
    //! (Hagen-Poiseuille)
    virtual void EvaluateAdditionalLinearizationofIntegratedDiam(
        CORE::LINALG::SerialDenseMatrix* stiffmat11,
        CORE::LINALG::SerialDenseMatrix* stiffmat12) = 0;

    //! flag if diameter function is active, i.e., varying diameter linearization need to be
    //! calculated
    virtual bool DiamFunctionActive() = 0;

    //! reset state
    virtual void ResetState(
        Teuchos::RCP<DRT::Discretization> contdis, Teuchos::RCP<DRT::Discretization> artdis) = 0;

    /**
     * Setup the porofluid-managers and the materials for later evaluation
     * @param[in] disname: name of continuous discretization
     * @param[in] timefacrhs_art: right hand side factor for artery time integration
     * @param[in] timefacrhs_cont: right hand side factor for time integration of 2D/3D
     * discretization
     */
    virtual void SetupFluidManagersAndMaterials(
        const std::string disname, const double& timefacrhs_art, const double& timefacrhs_cont) = 0;

    //! beginning of integration segment
    virtual double EtaA() const = 0;
    //! end of integration segment
    virtual double EtaB() const = 0;

    //! element 1 (= artery) GID
    virtual int Ele1GID() const = 0;
    //! element 2 (= cont) GID
    virtual int Ele2GID() const = 0;

    //! apply mesh movement on artery element
    virtual double ApplyMeshMovement(
        const bool firstcall, Teuchos::RCP<DRT::Discretization> contdis) = 0;

    //! set segment id
    virtual void SetSegmentID(const int& segmentid) = 0;
    //! get segment id
    virtual int GetSegmentID() const = 0;

    //! get the volume of the 2D/3D element
    virtual double CalculateVol2D3D() const = 0;

    //! get number of Gauss points
    virtual int NumGP() const = 0;

    //! type of coupling pair
    enum CouplingType
    {
      type_undefined,
      type_porofluid,  //!< porofluid
      type_scatra      //!< scatra
    };
  };

  //! the coupling pair
  template <CORE::FE::CellType distypeArt, CORE::FE::CellType distypeCont, int dim>
  class PoroMultiPhaseScatraArteryCouplingPair : public PoroMultiPhaseScatraArteryCouplingPairBase
  {
   public:
    //! constructor
    PoroMultiPhaseScatraArteryCouplingPair();

    //! Init
    void Init(std::vector<DRT::Element const*> elements,
        const Teuchos::ParameterList& couplingparams,
        const Teuchos::ParameterList& fluidcouplingparams, const std::vector<int>& coupleddofs_cont,
        const std::vector<int>& coupleddofs_art, const std::vector<std::vector<int>>& scale_vec,
        const std::vector<std::vector<int>>& funct_vec, const std::string condname,
        const double penalty, const std::string couplingtype = "", const int eta_ntp = 0) override;

    //! query if pair active
    bool IsActive() override { return isactive_; }

    //! things that need to be done in a separate loop before the actual evaluation loop
    //! over all coupling pairs
    void PreEvaluate(Teuchos::RCP<Epetra_MultiVector> gp_vector) override;

    //! things that need to be done in a separate loop before the actual evaluation loop
    //! over all coupling pairs
    void DeleteUnnecessaryGPs(Teuchos::RCP<Epetra_MultiVector> gp_vector) override;

    //! flag if diameter function is active, i.e., varying diameter linearization need to be
    //! calculated
    bool DiamFunctionActive() override { return diam_funct_active_; }

    //! reset state
    void ResetState(Teuchos::RCP<DRT::Discretization> contdis,
        Teuchos::RCP<DRT::Discretization> artdis) override;

    /**
     * Setup the porofluid-managers and the materials for later evaluation
     * @param[in] disname: name of continuous discretization
     * @param[in] timefacrhs_art: right hand side factor for artery time integration
     * @param[in] timefacrhs_cont: right hand side factor for time integration of 2D/3D
     * discretization
     */
    void SetupFluidManagersAndMaterials(const std::string disname, const double& timefacrhs_art,
        const double& timefacrhs_cont) override;

    /*!
     * @brief Evaluate this pair
     *
     * @returns integral of diameter of the segment
     */
    double Evaluate(CORE::LINALG::SerialDenseVector* forcevec1,
        CORE::LINALG::SerialDenseVector* forcevec2, CORE::LINALG::SerialDenseMatrix* stiffmat11,
        CORE::LINALG::SerialDenseMatrix* stiffmat12, CORE::LINALG::SerialDenseMatrix* stiffmat21,
        CORE::LINALG::SerialDenseMatrix* stiffmat22, CORE::LINALG::SerialDenseMatrix* D_ele,
        CORE::LINALG::SerialDenseMatrix* M_ele, CORE::LINALG::SerialDenseVector* Kappa_ele,
        const std::vector<double>& segmentlengths) override;

    //! evaluate additional linearization of (integrated) element diameter dependent terms
    //! (Hagen-Poiseuille)
    void EvaluateAdditionalLinearizationofIntegratedDiam(
        CORE::LINALG::SerialDenseMatrix* stiffmat11,
        CORE::LINALG::SerialDenseMatrix* stiffmat12) override;

    //! beginning and end of integration segment
    double EtaA() const override { return eta_a_; }
    double EtaB() const override { return eta_b_; }

    //! element 1 (= artery) GID
    int Ele1GID() const override;
    //! element 2 (= cont) GID
    int Ele2GID() const override;

    //! number of GP
    int NumGP() const override { return n_gp_; };

    //! apply mesh movement on artery element
    double ApplyMeshMovement(
        const bool firstcall, Teuchos::RCP<DRT::Discretization> contdis) override;

    //! set segment id
    void SetSegmentID(const int& segmentid) override;
    //! get segment id
    int GetSegmentID() const override;

    //! get the volume of the 2D/3D element
    double CalculateVol2D3D() const override;

   private:
    // static variables
    //! number of nodes of 1D artery element
    static constexpr unsigned numnodesart_ = CORE::FE::num_nodes<distypeArt>;
    //! number of nodes of 2D/3D continuous element
    static constexpr unsigned numnodescont_ = CORE::FE::num_nodes<distypeCont>;
    //! number of nodes of spatial dimensions
    static constexpr unsigned numdim_ = CORE::FE::dim<distypeCont>;

    //! set time factor needed for evaluation of right hand side (function coupling) terms
    void SetTimeFacRhs(const double& arterydensity, Teuchos::RCP<MAT::MatList> contscatramat,
        const double& timefacrhs_art, const double& timefacrhs_cont);

    //! pre-evaluate for lateral surface coupling
    void PreEvaluateLateralSurfaceCoupling(Teuchos::RCP<Epetra_MultiVector> gp_vector);

    //! pre-evaluate for centerline coupling
    void PreEvaluateCenterlineCoupling();

    //! pre-evaluate for node-to-point coupling
    void PreEvaluateNodeToPointCoupling();

    //! extract velocity of solid phase
    void ExtractSolidVel(Teuchos::RCP<DRT::Discretization> contdis);

    //! recompute if deformable arteries are assumed
    void RecomputeEtaAndXiInDeformedConfiguration(const std::vector<double>& segmentlengths,
        std::vector<double>& myEta, std::vector<std::vector<double>>& myXi, double& etaA,
        double& etaB);

    /**
     * \brief create segment [eta_a, eta_b]
     *
     * \note: the following algorithm works only for linear 1D elements and linear 2D/3D
     * elements where always 0,1 or 2 intersections can be found. For higher order elements with
     * special cases, it has to be re-thought. For instance, if we find two
     * intersections, it is always assumed that the integration segment lies between these two
     * intersections --> for a higher order 1D element this may not be the case
     */
    void CreateIntegrationSegment();

    //! get all intersections of artery element with 2D/3D element
    std::vector<double> GetAllInterSections();

    //! project a Gauss point on 1D element into 2D/3D element
    template <typename T>
    void Projection(
        CORE::LINALG::Matrix<numdim_, 1, T>& r1, std::vector<T>& xi, bool& projection_valid);

    //! Check for duplicate projections
    bool ProjectionNotYetFound(const std::vector<double>& intersections, const double& eta);

    //! Intersect artery element with edges (2D) or surfaces (3D) of element
    void InterSectWith2D3D(std::vector<double>& xi, double& eta, const int& fixedPar,
        const double& fixedAt, bool& projection_valid);

    //! get 1D shapefunctions at eta
    template <typename T>
    void Get1DShapeFunctions(CORE::LINALG::Matrix<1, numnodesart_, T>& N1,
        CORE::LINALG::Matrix<1, numnodesart_, T>& N1_eta, const T& eta);

    //! get 2D/3D shapefunctions at xi1, xi2 (, xi3)
    template <typename T>
    void Get2D3DShapeFunctions(CORE::LINALG::Matrix<1, numnodescont_, T>& N2,
        CORE::LINALG::Matrix<numdim_, numnodescont_, T>& N2_xi, const std::vector<T>& xi);

    //! compute artery coordinates and derivatives in reference configuration
    template <typename T>
    void ComputeArteryCoordsAndDerivsRef(CORE::LINALG::Matrix<numdim_, 1, T>& r1,
        CORE::LINALG::Matrix<numdim_, 1, T>& r1_eta,
        const CORE::LINALG::Matrix<1, numnodesart_, T>& N1,
        const CORE::LINALG::Matrix<1, numnodesart_, T>& N1_eta);

    //! compute 2D/3D coordinates and derivatives in reference configuration
    template <typename T>
    void Compute2D3DCoordsAndDerivsRef(CORE::LINALG::Matrix<numdim_, 1, T>& x2,
        CORE::LINALG::Matrix<numdim_, numdim_, T>& x2_xi,
        const CORE::LINALG::Matrix<1, numnodescont_, T>& N2,
        const CORE::LINALG::Matrix<numdim_, numnodescont_, T>& N2_xi);

    //! evaluate the function coupling (return integral of diameter of the segment)
    void EvaluateFunctionCoupling(const std::vector<double>& eta,
        const std::vector<std::vector<double>>& xi, const std::vector<double>& segmentlengths,
        CORE::LINALG::SerialDenseVector* forcevec1, CORE::LINALG::SerialDenseVector* forcevec2,
        CORE::LINALG::SerialDenseMatrix* stiffmat11, CORE::LINALG::SerialDenseMatrix* stiffmat12,
        CORE::LINALG::SerialDenseMatrix* stiffmat21, CORE::LINALG::SerialDenseMatrix* stiffmat22,
        double& integrated_diam);

    /**
     * evaluate derivative of 1D shape function times solid velocity (only porofluid has this term)
     * @param[in] eta: GP coordinates in artery element parameter space
     * @param[in] xi: GP coordinates in porofluid element parameter space
     * @param[in] segmentlengths: length of all segments of this artery element
     * @param[in] forcevec1: rhs-vector to assemble into
     * @param[in] etaA: beginning of segment in artery element parameter space
     * @param[in] etaB: end of segment in artery element parameter space
     */
    void EvaluatedNdsSolidVel(const std::vector<double>& eta,
        const std::vector<std::vector<double>>& xi, const std::vector<double>& segmentlengths,
        CORE::LINALG::SerialDenseVector& forcevec1, const double& etaA, const double& etaB);

    //! evaluate stiffness for GPTS case
    void EvaluateGPTSStiff(const double& w_gp, const CORE::LINALG::Matrix<1, numnodesart_>& N1,
        const CORE::LINALG::Matrix<1, numnodescont_>& N2, const double& jacobi, const double& pp);

    //! evaluate stiffness for NTP case
    void EvaluateNTPStiff(const CORE::LINALG::Matrix<1, numnodesart_>& N1,
        const CORE::LINALG::Matrix<1, numnodescont_>& N2, const double& pp);

    //! evaluate mortar coupling matrices D and M
    void EvaluateDMKappa(const double& w_gp, const CORE::LINALG::Matrix<1, numnodesart_>& N1,
        const CORE::LINALG::Matrix<1, numnodescont_>& N2, const double& jacobi);

    //! evalute GPTS
    void EvaluateGPTS(const std::vector<double>& eta, const std::vector<std::vector<double>>& xi,
        const std::vector<double>& segmentlengths, CORE::LINALG::SerialDenseVector* forcevec1,
        CORE::LINALG::SerialDenseVector* forcevec2, CORE::LINALG::SerialDenseMatrix* stiffmat11,
        CORE::LINALG::SerialDenseMatrix* stiffmat12, CORE::LINALG::SerialDenseMatrix* stiffmat21,
        CORE::LINALG::SerialDenseMatrix* stiffmat22);

    //! evalute NTP
    void EvaluateNTP(const std::vector<double>& eta, const std::vector<std::vector<double>>& xi,
        CORE::LINALG::SerialDenseVector* forcevec1, CORE::LINALG::SerialDenseVector* forcevec2,
        CORE::LINALG::SerialDenseMatrix* stiffmat11, CORE::LINALG::SerialDenseMatrix* stiffmat12,
        CORE::LINALG::SerialDenseMatrix* stiffmat21, CORE::LINALG::SerialDenseMatrix* stiffmat22);

    //! evaluate mortar coupling matrices D and M
    void EvaluateDMKappa(const std::vector<double>& eta, const std::vector<std::vector<double>>& xi,
        const std::vector<double>& segmentlengths, CORE::LINALG::SerialDenseMatrix* D_ele,
        CORE::LINALG::SerialDenseMatrix* M_ele, CORE::LINALG::SerialDenseVector* Kappa_ele);

    //! evaluate the function coupling
    void EvaluateFunctionCoupling(const double& w_gp,
        const CORE::LINALG::Matrix<1, numnodesart_>& N1,
        const CORE::LINALG::Matrix<1, numnodescont_>& N2, const double& jacobi,
        CORE::LINALG::SerialDenseVector& forcevec1, CORE::LINALG::SerialDenseVector& forcevec2,
        CORE::LINALG::SerialDenseMatrix& stiffmat11, CORE::LINALG::SerialDenseMatrix& stiffmat12,
        CORE::LINALG::SerialDenseMatrix& stiffmat21, CORE::LINALG::SerialDenseMatrix& stiffmat22,
        double& integrated_diam);

    //! evaluate the diameter function and derivative (for couplingtype porofluid)
    void EvaluateDiamFunctionAndDeriv(const double artpressnpAtGP, const double& w_gp,
        const CORE::LINALG::Matrix<1, numnodesart_>& N1,
        const CORE::LINALG::Matrix<1, numnodescont_>& N2, const double& jacobi);

    //! integrate in deformed configuration from eta_A to eta_s
    FAD IntegrateLengthToEtaS(const FAD& eta_s);

    //! get values of artery at GP
    void GetArteryValuesAtGP(const CORE::LINALG::Matrix<1, numnodesart_>& N1, double& artpress,
        std::vector<double>& artscalar);

    //! get scalar values of continuous discretization at GP
    void GetContScalarValuesAtGP(
        const CORE::LINALG::Matrix<1, numnodescont_>& N2, std::vector<double>& contscalarnp);

    //! assemble the function coupling into stiffness matrix (artery-part)
    void AssembleFunctionCouplingIntoForceStiffArt(const int& i_art, const double& w_gp,
        const CORE::LINALG::Matrix<1, numnodesart_>& N1,
        const CORE::LINALG::Matrix<1, numnodescont_>& N2, const double& jacobi, const int& scale,
        const double& functval, const std::vector<double>& artderivs,
        const std::vector<double>& contderivs, CORE::LINALG::SerialDenseVector& forcevec1,
        CORE::LINALG::SerialDenseMatrix& stiffmat11, CORE::LINALG::SerialDenseMatrix& stiffmat12);

    //! assemble the function coupling into stiffness matrix (2D/3D-part)
    void AssembleFunctionCouplingIntoForceStiffCont(const std::vector<int>& assembleInto,
        const double& w_gp, const CORE::LINALG::Matrix<1, numnodesart_>& N1,
        const CORE::LINALG::Matrix<1, numnodescont_>& N2, const double& jacobi, const int& scale,
        const double& timefacrhs_cont, const double& functval, const std::vector<double>& artderivs,
        const std::vector<double>& contderivs, CORE::LINALG::SerialDenseVector& forcevec2,
        CORE::LINALG::SerialDenseMatrix& stiffmat21, CORE::LINALG::SerialDenseMatrix& stiffmat22);

    //! evaluate function and its derivative

    void EvaluateFunctionAndDeriv(const CORE::UTILS::FunctionOfAnything& funct,
        const double& artpressnpAtGP, const std::vector<double>& artscalarnpAtGP,
        const std::vector<double>& scalarnpAtGP, double& functval, std::vector<double>& artderivs,
        std::vector<double>& contderivs);

    //! set scalar as constants into function
    void SetScalarValuesAsConstants(std::vector<std::pair<std::string, double>>& constants,
        const std::vector<double>& artscalarnpAtGP, const std::vector<double>& scalarnpAtGP);

    //! set fluid as variables into function
    void SetFluidValuesAsVariables(
        std::vector<std::pair<std::string, double>>& variables, const double& artpressnpAtGP);

    //! set fluid as constants into function
    void SetFluidValuesAsConstants(
        std::vector<std::pair<std::string, double>>& constants, const double& artpressnpAtGP);

    //! set scalar as variables into function
    void SetScalarValuesAsVariables(std::vector<std::pair<std::string, double>>& variables,
        const std::vector<double>& artscalarnpAtGP, const std::vector<double>& scalarnpAtGP);

    //! evaluate derivatives w.r.t. fluid of function
    void EvaluateFluidDerivs(std::vector<double>& artderivs, std::vector<double>& contderivs,
        const std::vector<double>& functderivs);

    //! evaluate derivatives w.r.t. scalar of function
    void EvaluateScalarDerivs(std::vector<double>& artderivs, std::vector<double>& contderivs,
        const std::vector<double>& functderivs);

    //! evaluate force for GPTS  or NTP case
    void EvaluateGPTSNTPForce(CORE::LINALG::SerialDenseVector& forcevec1,
        CORE::LINALG::SerialDenseVector& forcevec2,
        const CORE::LINALG::SerialDenseMatrix& stiffmat11,
        const CORE::LINALG::SerialDenseMatrix& stiffmat12,
        const CORE::LINALG::SerialDenseMatrix& stiffmat21,
        const CORE::LINALG::SerialDenseMatrix& stiffmat22);

    //! update the stiffness for GPTS or NTP
    void UpdateGPTSNTPStiff(CORE::LINALG::SerialDenseMatrix& stiffmat11,
        CORE::LINALG::SerialDenseMatrix& stiffmat12, CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22);

    /**
     * \brief coupling to additional porous network is only possible if we also have an
     * element with a valid volume fraction pressure, i.e., if we also have a smeared representation
     * of the neovasculature at this point if not ---> corresponding matrices are set to zero
     */
    void CheckValidVolumeFractionPressureCoupling(CORE::LINALG::SerialDenseMatrix& stiffmat11,
        CORE::LINALG::SerialDenseMatrix& stiffmat12, CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22);

    //! update the D, M and Kappa for MP
    void UpdateDMKappa(CORE::LINALG::SerialDenseMatrix& D_ele,
        CORE::LINALG::SerialDenseMatrix& M_ele, CORE::LINALG::SerialDenseVector& Kappa_ele);

    //! fill the function vector
    void FillFunctionVector(std::vector<const CORE::UTILS::FunctionOfAnything*>& my_funct_vec,
        const std::vector<int>& funct_vec, const std::vector<int>& scale_vec);


    //! initialize a function
    void InitializeFunction(const CORE::UTILS::FunctionOfAnything& funct);

    //! initialize names used in functions
    void InitializeFunctionNames();

    //! initialize vector where to assemble continuous DOF functions into
    void InitializeAssembleIntoContDofVector();

    //! coupling type
    CouplingType coupltype_;

    //! coupling method (either GPTS or MP)
    INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod couplmethod_;

    //! name of the condition
    std::string condname_;

    //! indicates if the Init() function has been called
    bool isinit_;

    //! indicates if the PreEvaluate() function has been called
    bool ispreevaluated_;

    //! indicates if mesh tying is active, i.e., if projection possible
    bool isactive_;

    //! indicates if function coupling is active, i.e., if functions are defined
    bool funct_coupl_active_;

    //! indicates if diameter function is active, i.e, if diameter function is defined
    bool diam_funct_active_;

    //! so far, it is assumed that artery elements always follow the deformation of the
    //! underlying porous medium hence, we actually have to evalute them in current
    //! configuration if this flag is set to true, artery elements will not move and are
    //! evaluated in reference configuration
    bool evaluate_in_ref_config_;

    //! evaluate 1D-3D coupling on lateral surface?
    bool evaluate_on_lateral_surface_;

    //! first element of interacting pair (artery element)
    const DRT::Element* element1_;

    //! second element of interacting pair (2D/3D element)
    const DRT::Element* element2_;

    //! reference nodal positions element 1 (1D artery)
    CORE::LINALG::Matrix<numdim_ * numnodesart_, 1> ele1posref_;
    //! reference nodal positions element 2 (2D/3D continuous element)
    CORE::LINALG::Matrix<numdim_, numnodescont_> ele2posref_;

    //! current position of element 2
    CORE::LINALG::Matrix<numdim_, numnodescont_> ele2pos_;
    //! current velocity of element 2
    CORE::LINALG::Matrix<numdim_, numnodescont_> ele2vel_;

    //! reference diameter of the artery element (constant across element)
    double arterydiamref_;
    //! current diameter of the artery element at the GP
    double arterydiamAtGP_;
    //! derivatives of the diameter function
    std::vector<double> diamderivs_;

    //! numdofs of 2D/3D element
    int numdof_cont_;
    //! numdofs of artery element
    int numdof_art_;
    //! number of dofs * number of nodes of artery element
    int dim1_;
    //! number of dofs * number of nodes of 2D/3D element
    int dim2_;

    //! coupled dofs of continuous (2D/3D) element
    std::vector<int> coupleddofs_cont_;
    //! coupled dofs of artery (1D) element
    std::vector<int> coupleddofs_art_;
    //! number of coupled dofs
    int numcoupleddofs_;

    //! the id of the volume fraction pressure phase
    std::vector<int> volfracpressid_;

    //! number of fluid phases
    int numfluidphases_;
    //! number of volume fractions
    int numvolfrac_;

    //! number of scalars (cont)
    int numscalcont_;
    //! number of scalars (art)
    int numscalart_;

    //! dof-set number of porofluid (either 0 or 2)
    int nds_porofluid_;

    //! stores the number of Gauss points for that element
    int n_gp_;

    //! stores the number of Gauss points for that element per patch (only required for
    //! surface-based formulation)
    int n_gp_per_patch_;

    //! stores the artery element length in reference configuration
    double arteryelelengthref_;

    //! artery element length in current configuration
    double arteryelelength_;

    //! stores initial direction of artery element
    CORE::LINALG::Matrix<numdim_, 1> lambda0_;

    //! Jacobian determinant for integration segment = L/2.0*(eta_a - eta_b)/2.0
    double jacobi_;

    //! Gausspoints in solid
    std::vector<std::vector<double>> xi_;

    //! Gausspoints in arteries
    std::vector<double> eta_;
    //! Gausspoint weigths in arteries
    std::vector<double> wgp_;

    //! eta_s
    std::vector<double> eta_s_;

    //! primary variables of cont element
    std::vector<double> contelephinp_;
    //! primary variables of artery element
    std::vector<double> artelephinp_;

    //! nodal artery pressure values
    CORE::LINALG::Matrix<numnodesart_, 1> earterypressurenp_;

    //! nodal artery-scalar values
    std::vector<CORE::LINALG::Matrix<numnodesart_, 1>> eartscalarnp_;

    //! nodal continuous-scalar values for scatra coupling
    std::vector<CORE::LINALG::Matrix<numnodescont_, 1>> econtscalarnp_;

    //! penalty parameter
    double pp_;

    //! start of integration segment
    double eta_a_;
    //! end of integration segment
    double eta_b_;

    //! length of integration segment int current configuration
    double curr_segment_length_;

    //! check if constant part (i.e. GPTS and MP part) has already been evaluated if integration
    //! in reference configuration is performed
    bool constant_part_evaluated_;

    //! 1D Cooupling element type (can be ARTERY or AIRWAY)
    std::string coupling_element_type_;

    //! GPTS/NTP stiffness matrix (artery-artery contribution)
    CORE::LINALG::SerialDenseMatrix GPTS_NTP_stiffmat11_;
    //! GPTS/NTP stiffness matrix (artery-cont contribution)
    CORE::LINALG::SerialDenseMatrix GPTS_NTP_stiffmat12_;
    //! GPTS/NTP stiffness matrix (cont-artery contribution)
    CORE::LINALG::SerialDenseMatrix GPTS_NTP_stiffmat21_;
    //! GPTS/NTP stiffness matrix (cont-cont contribution)
    CORE::LINALG::SerialDenseMatrix GPTS_NTP_stiffmat22_;

    //! (varying) diameter stiffness matrix (artery-artery contribution)
    CORE::LINALG::SerialDenseMatrix diam_stiffmat11_;
    //! (varying) diameter stiffness matrix (artery-cont contribution)
    CORE::LINALG::SerialDenseMatrix diam_stiffmat12_;

    //! mortar coupling matrix D
    CORE::LINALG::SerialDenseMatrix D_;
    //! mortar coupling matrix M
    CORE::LINALG::SerialDenseMatrix M_;
    //! mortar coupling vector kappa
    CORE::LINALG::SerialDenseVector Kappa_;

    //! (dX/dxi)^-1
    std::vector<CORE::LINALG::Matrix<numdim_, numdim_>> invJ_;

    //! phase manager of the fluid
    Teuchos::RCP<DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface> phasemanager_;

    //! variable manager of the fluid
    Teuchos::RCP<DRT::ELEMENTS::POROFLUIDMANAGER::VariableManagerInterface<numdim_, numnodescont_>>
        variablemanager_;

    //! scale vector
    std::vector<std::vector<int>> scale_vec_;
    //! function vector
    std::vector<std::vector<const CORE::UTILS::FunctionOfAnything*>> funct_vec_;

    //! diameter function
    const CORE::UTILS::FunctionOfAnything* artdiam_funct_;

    //! string name used for scalars in function parser
    std::vector<std::string> scalarnames_;
    //! string name used for pressure in function parser
    std::vector<std::string> pressurenames_;
    //! string name used for saturation in function parser
    std::vector<std::string> saturationnames_;
    //! string name used for porosity in function parser
    const std::string porosityname_;
    //! string name used for artery-pressure in function parser
    const std::string artpressname_;
    //! string name used for artery-scalars in function parser
    std::vector<std::string> artscalarnames_;
    //! string name used for volume fractions in function parser
    std::vector<std::string> volfracnames_;
    //! string name used for volume fraction pressures in function parser
    std::vector<std::string> volfracpressurenames_;

    //! dofset of artery pressure in scatra-dis
    //! TODO: find a better way to do this
    const int ndsscatra_artery_ = 2;

    //! dofset of scatra primary variable in artery-dis
    //! TODO: find a better way to do this
    const int ndsartery_scatra_ = 2;

    //! segment id
    int segmentid_;

    //! number of integration patches in axial direction (for surface-based coupling)
    int numpatch_axi_;
    //! number of integration patches in radial direction (for surface-based coupling)
    int numpatch_rad_;

    //! right hand side factor for artery time integration scaled with inverse density
    double timefacrhs_art_dens_;
    //! right hand side factor for time integration of 2D/3D discretization scaled with inverse
    //! density of specific phase or species
    std::vector<double> timefacrhs_cont_dens_;

    //! right hand side factor for artery time integration
    double timefacrhs_art_;
    //! right hand side factor for time integration of 2D/3D discretization
    double timefacrhs_cont_;

    //! vector where to assemble rhs-(function) coupling into
    //! summed up phase requires special treatment
    std::vector<std::vector<int>> cont_dofs_to_assemble_functions_into_;

    //! the artery material
    Teuchos::RCP<MAT::Cnst1dArt> arterymat_;
  };

}  // namespace POROMULTIPHASESCATRA


FOUR_C_NAMESPACE_CLOSE

#endif
