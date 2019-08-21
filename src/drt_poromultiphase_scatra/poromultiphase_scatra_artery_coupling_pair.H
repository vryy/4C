/*----------------------------------------------------------------------*/
/*!
 \brief one pair consisting of exactly one artery element and one poro-
        multiphase-scatra element which might be tied to each other

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_PAIR_H_
#define SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_PAIR_H_

#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_utils.H"
#include "../drt_porofluidmultiphase_ele/porofluid_variablemanager.H"
#include "../drt_porofluidmultiphase_ele/porofluid_phasemanager.H"
#include "../drt_inpar/inpar_bio.H"
#include "Sacado.hpp"


#include "Teuchos_RCP.hpp"
#include <Epetra_Vector.h>


// define Fad object for evaluation
typedef Sacado::Fad::DFad<double> FAD;

// forward declaration
namespace DRT
{
  class Element;
  class Discretization;

  namespace UTILS
  {
    class VariableExprFunction;
  }
}  // namespace DRT
namespace MAT
{
  class MatList;
}  // namespace MAT
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace LINALG
namespace POROMULTIPHASESCATRA
{
  class PoroMultiPhaseScatraArteryCouplingPairBase
  {
   public:
    //! constructor
    PoroMultiPhaseScatraArteryCouplingPairBase() { return; };

    //! destructor
    virtual ~PoroMultiPhaseScatraArteryCouplingPairBase(){};

    //! Init
    virtual void Init(std::vector<DRT::Element const*> elements,
        const Teuchos::ParameterList& meshtyingparams, const std::vector<int>& coupleddofs_cont,
        const std::vector<int>& coupleddofs_art, const std::vector<std::vector<int>>& scale_vec,
        const std::vector<std::vector<int>>& funct_vec) = 0;

    //! query if pair active
    virtual bool IsActive() = 0;

    //! things that need to be done in a separate loop before the actual evaluation loop
    //! over all coupling pairs
    virtual void PreEvaluate() = 0;

    //! evaluate this pair
    virtual void Evaluate(LINALG::SerialDenseVector* forcevec1,
        LINALG::SerialDenseVector* forcevec2, LINALG::SerialDenseMatrix* stiffmat11,
        LINALG::SerialDenseMatrix* stiffmat12, LINALG::SerialDenseMatrix* stiffmat21,
        LINALG::SerialDenseMatrix* stiffmat22, LINALG::SerialDenseMatrix* D_ele,
        LINALG::SerialDenseMatrix* M_ele, LINALG::SerialDenseVector* Kappa_ele,
        const std::vector<double>& segmentlengths) = 0;

    //! reset state
    virtual void ResetState(
        Teuchos::RCP<DRT::Discretization> contdis, Teuchos::RCP<DRT::Discretization> artdis) = 0;

    /**
     * Setup the porofluid-managers and the materials for later evaluation
     * @param disname: name of continuous discretization
     * @param timefacrhs_art: right hand side factor for artery time integration
     * @param timefacrhs_cont: right hand side factor for time integration of 2D/3D discretization
     */
    virtual void SetupFluidManagersAndMaterials(
        const std::string disname, const double& timefacrhs_art, const double& timefacrhs_cont) = 0;

    //! beginning and end of integration segment
    virtual double EtaA() const = 0;
    virtual double EtaB() const = 0;

    //! the element GIDs (1 = artery, 2 = cont)
    virtual int Ele1GID() const = 0;
    virtual int Ele2GID() const = 0;

    //! apply mesh movement on artery element
    virtual double ApplyMeshMovement(
        const bool firstcall, Teuchos::RCP<DRT::Discretization> contdis) = 0;

    //! set and get segment id
    virtual void SetSegmentID(const int& segmentid) = 0;
    virtual int GetSegmentID() const = 0;

    //! type of coupling pair
    enum CouplingType
    {
      type_undefined,
      type_porofluid,  //!< porofluid
      type_scatra      //!< scatra
    };
  };

  //! the coupling pair
  template <DRT::Element::DiscretizationType distypeArt,
      DRT::Element::DiscretizationType distypeCont>
  class PoroMultiPhaseScatraArteryCouplingPair : public PoroMultiPhaseScatraArteryCouplingPairBase
  {
   public:
    //! constructor
    PoroMultiPhaseScatraArteryCouplingPair();

    //! destructor
    virtual ~PoroMultiPhaseScatraArteryCouplingPair(){};

    //! Init
    void Init(std::vector<DRT::Element const*> elements,
        const Teuchos::ParameterList& meshtyingparams, const std::vector<int>& coupleddofs_cont,
        const std::vector<int>& coupleddofs_art, const std::vector<std::vector<int>>& scale_vec,
        const std::vector<std::vector<int>>& funct_vec);

    //! query if pair active
    virtual bool IsActive() { return isactive_; }

    //! things that need to be done in a separate loop before the actual evaluation loop
    //! over all coupling pairs
    void PreEvaluate();

    //! reset state
    void ResetState(
        Teuchos::RCP<DRT::Discretization> contdis, Teuchos::RCP<DRT::Discretization> artdis);

    /**
     * Setup the porofluid-managers and the materials for later evaluation
     * @param disname: name of continuous discretization
     * @param timefacrhs_art: right hand side factor for artery time integration
     * @param timefacrhs_cont: right hand side factor for time integration of 2D/3D discretization
     */
    void SetupFluidManagersAndMaterials(
        const std::string disname, const double& timefacrhs_art, const double& timefacrhs_cont);

    //! evaluate this pair
    void Evaluate(LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
        LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
        LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22,
        LINALG::SerialDenseMatrix* D_ele, LINALG::SerialDenseMatrix* M_ele,
        LINALG::SerialDenseVector* Kappa_ele, const std::vector<double>& segmentlengths);

    //! beginning and end of integration segment
    double EtaA() const { return eta_a_; }
    double EtaB() const { return eta_b_; }

    //! the element GIDs (1 = artery, 2 = cont)
    int Ele1GID() const;
    int Ele2GID() const;

    //! apply mesh movement on artery element
    double ApplyMeshMovement(const bool firstcall, Teuchos::RCP<DRT::Discretization> contdis);

    //! set and get segment id
    void SetSegmentID(const int& segmentid);
    int GetSegmentID() const;

   private:
    //! static variables
    static const unsigned numnodesart_ =
        DRT::UTILS::DisTypeToNumNodePerEle<distypeArt>::numNodePerElement;
    static const unsigned numnodescont_ =
        DRT::UTILS::DisTypeToNumNodePerEle<distypeCont>::numNodePerElement;
    static const unsigned numdim_ = DRT::UTILS::DisTypeToDim<distypeCont>::dim;

    //! set time factor needed for evaluation of right hand side (function coupling) terms
    void SetTimeFacRhs(const double& arterydensity, Teuchos::RCP<MAT::MatList> contscatramat,
        const double& timefacrhs_art, const double& timefacrhs_cont);

    //! extract velocity of solid phase
    void ExtractSolidVel(Teuchos::RCP<DRT::Discretization> contdis);

    //! recompute if deformable arteries are assumed
    void RecomputeEtaAndXiInDeformedConfiguration(const std::vector<double>& segmentlengths,
        std::vector<double>& myEta, std::vector<std::vector<double>>& myXi, double& etaA,
        double& etaB);

    //! create segment [eta_a, eta_b]
    void CreateIntegrationSegment(const std::vector<bool>& validprojections);

    //! get all intersections of artery element with 2D/3D element
    std::vector<double> GetAllInterSections();

    //! project a Gauss point on 1D element into 2D/3D element
    template <typename T>
    void Projection(const T& eta, std::vector<T>& xi, bool& projection_valid);

    //! Check for duplicate projections
    bool ProjectionNotYetFound(const std::vector<double>& intersections, const double& eta);

    //! Intersect artery element with edges (2D) or surfaces (3D) of element
    void InterSectWith2D3D(std::vector<double>& xi, double& eta, const int& fixedPar,
        const double& fixedAt, bool& projection_valid);

    //! get 1D shapefunctions at eta
    template <typename T>
    void Get1DShapeFunctions(LINALG::Matrix<1, numnodesart_, T>& N1,
        LINALG::Matrix<1, numnodesart_, T>& N1_eta, const T& eta);

    //! get 2D/3D shapefunctions at xi1, xi2 (, xi3)
    template <typename T>
    void Get2D3DShapeFunctions(LINALG::Matrix<1, numnodescont_, T>& N2,
        LINALG::Matrix<numdim_, numnodescont_, T>& N2_xi, const std::vector<T>& xi);

    //! compute artery coordinates and derivatives in reference configuration
    template <typename T>
    void ComputeArteryCoordsAndDerivsRef(LINALG::Matrix<numdim_, 1, T>& r1,
        LINALG::Matrix<numdim_, 1, T>& r1_eta, const LINALG::Matrix<1, numnodesart_, T>& N1,
        const LINALG::Matrix<1, numnodesart_, T>& N1_eta);

    //! compute 2D/3D coordinates and derivatives in reference configuration
    template <typename T>
    void Compute2D3DCoordsAndDerivsRef(LINALG::Matrix<numdim_, 1, T>& x2,
        LINALG::Matrix<numdim_, numdim_, T>& x2_xi, const LINALG::Matrix<1, numnodescont_, T>& N2,
        const LINALG::Matrix<numdim_, numnodescont_, T>& N2_xi);

    //! evaluate the function coupling
    void EvaluateFunctionCoupling(const std::vector<double>& eta,
        const std::vector<std::vector<double>>& xi, const std::vector<double>& segmentlengths,
        LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
        LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
        LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22);

    /**
     * evaluate divergence of solid velocity term (only porofluid has this term)
     * @param eta [in]: GP coordinates in artery element parameter space
     * @param xi [in]: GP coordinates in porofluid element parameter space
     * @param segmentlengths [in]: length of all segments of this artery element
     * @param forcevec1 [in]: rhs-vector to assemble into
     * @param etaA [in]: beginning of segment in artery element parameter space
     * @param etaB [in]: end of segment in artery element parameter space
     */
    void EvaluateDivSolidVel(const std::vector<double>& eta,
        const std::vector<std::vector<double>>& xi, const std::vector<double>& segmentlengths,
        LINALG::SerialDenseVector& forcevec1, const double& etaA, const double& etaB);

    //! evaluate stiffness for GPTS case
    void EvaluateGPTSStiff(const double& w_gp, const LINALG::Matrix<1, numnodesart_>& N1,
        const LINALG::Matrix<1, numnodescont_>& N2, const double& jacobi, const double& pp);

    //! evaluate mortar coupling matrices D and M
    void EvaluateDMKappa(const double& w_gp, const LINALG::Matrix<1, numnodesart_>& N1,
        const LINALG::Matrix<1, numnodescont_>& N2, const double& jacobi);

    //! evalute GPTS
    void EvaluateGPTS(const std::vector<double>& eta, const std::vector<std::vector<double>>& xi,
        const std::vector<double>& segmentlengths, const int& dim1, const int& dim2,
        LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
        LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
        LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22);

    //! evaluate mortar coupling matrices D and M
    void EvaluateDMKappa(const std::vector<double>& eta, const std::vector<std::vector<double>>& xi,
        const std::vector<double>& segmentlengths, const int& dim1, const int& dim2,
        LINALG::SerialDenseMatrix* D_ele, LINALG::SerialDenseMatrix* M_ele,
        LINALG::SerialDenseVector* Kappa_ele);

    //! evaluate the function coupling
    void EvaluateFunctionCoupling(const double& w_gp, const LINALG::Matrix<1, numnodesart_>& N1,
        const LINALG::Matrix<1, numnodescont_>& N2, const double& jacobi,
        LINALG::SerialDenseVector& forcevec1, LINALG::SerialDenseVector& forcevec2,
        LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
        LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22);

    //! integrate in deformed configuration from eta_A to eta_s
    FAD IntegrateLengthToEtaS(const FAD& eta_s);

    //! get values of artery at GP
    void GetArteryValuesAtGP(const LINALG::Matrix<1, numnodesart_>& N1, double& artpress,
        std::vector<double>& artscalar);

    //! get scalar values of continuous discretization at GP
    void GetContScalarValuesAtGP(
        const LINALG::Matrix<1, numnodescont_>& N2, std::vector<double>& contscalarnp);

    //! assemble the function coupling into stiffness matrix (artery-part)
    void AssembleFunctionCouplingIntoForceStiffArt(const int& i_art, const double& w_gp,
        const LINALG::Matrix<1, numnodesart_>& N1, const LINALG::Matrix<1, numnodescont_>& N2,
        const double& jacobi, const int& scale, const double& functval,
        const std::vector<double>& artderivs, const std::vector<double>& contderivs,
        LINALG::SerialDenseVector& forcevec1, LINALG::SerialDenseMatrix& stiffmat11,
        LINALG::SerialDenseMatrix& stiffmat12);

    //! assemble the function coupling into stiffness matrix (2D/3D-part)
    void AssembleFunctionCouplingIntoForceStiffCont(const std::vector<int>& assembleInto,
        const double& w_gp, const LINALG::Matrix<1, numnodesart_>& N1,
        const LINALG::Matrix<1, numnodescont_>& N2, const double& jacobi, const int& scale,
        const double& timefacrhs_cont, const double& functval, const std::vector<double>& artderivs,
        const std::vector<double>& contderivs, LINALG::SerialDenseVector& forcevec2,
        LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22);

    //! evaluate function and its derivative
    void EvaluateFunctionAndDeriv(DRT::UTILS::VariableExprFunction* funct,
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

    //! evaluate force for GPTS case
    void EvaluateGPTSForce(LINALG::SerialDenseVector& forcevec1,
        LINALG::SerialDenseVector& forcevec2, const LINALG::SerialDenseMatrix& stiffmat11,
        const LINALG::SerialDenseMatrix& stiffmat12, const LINALG::SerialDenseMatrix& stiffmat21,
        const LINALG::SerialDenseMatrix& stiffmat22);

    //! update the stiffness for GPTS
    void UpdateGPTSStiff(LINALG::SerialDenseMatrix& stiffmat11,
        LINALG::SerialDenseMatrix& stiffmat12, LINALG::SerialDenseMatrix& stiffmat21,
        LINALG::SerialDenseMatrix& stiffmat22);

    //! update the D, M and Kappa for MP
    void UpdateDMKappa(LINALG::SerialDenseMatrix& D_ele, LINALG::SerialDenseMatrix& M_ele,
        LINALG::SerialDenseVector& Kappa_ele);

    //! fill the function vector
    void FillFunctionVector(std::vector<DRT::UTILS::VariableExprFunction*>* my_funct_vec,
        const std::vector<int>& funct_vec, const std::vector<int>& scale_vec);

    //! initialize a function
    void InitializeFunction(DRT::UTILS::VariableExprFunction* funct);

    //! initialize names used in functions
    void InitializeFunctionNames();

    //! initialize vector where to assemble continuous DOF functions into
    void InitializeAssembleIntoContDofVector();

    //! get function
    DRT::UTILS::VariableExprFunction* Function(int functnum) const;

    //! coupling type
    CouplingType coupltype_;

    //! coupling method (either GPTS or MP)
    INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod couplmethod_;

    //! indicates if the Init() function has been called
    bool isinit_;

    //! indicates if the PreEvaluate() function has been called
    bool ispreevaluated_;

    //! indicates if mesh tying is active, i.e., if projection possible
    bool isactive_;

    //! indicates if function coupling is active, i.e., if functions are defined
    bool funct_coupl_active_;

    //! so far, it is assumed that artery elements always follow the deformation of the underlying
    //! porous medium
    //! hence, we actually have to evalute them in current configuration
    //! if this flag is set to true, artery elements will not move and are evaluated in reference
    //! configuration
    bool evaluate_in_ref_config_;

    //! first element of interacting pair (artery element)
    const DRT::Element* element1_;

    //! second element of interacting pair (2D/3D element)
    const DRT::Element* element2_;

    //! reference nodal positions of the two elements
    LINALG::Matrix<numdim_ * numnodesart_, 1> ele1posref_;
    LINALG::Matrix<numdim_, numnodescont_> ele2posref_;

    //! current position and velocity of element 2
    LINALG::Matrix<numdim_, numnodescont_> ele2pos_;
    LINALG::Matrix<numdim_, numnodescont_> ele2vel_;

    //! diameter of the artery element
    double arterydiam_;

    //! numdofs of artery and 2D/3D element
    int numdof_cont_;
    int numdof_art_;

    //! coupled dofs
    std::vector<int> coupleddofs_cont_;
    std::vector<int> coupleddofs_art_;
    int numcoupleddofs_;

    //! the id of the volume fraction pressure phase
    std::vector<int> volfracpressid_;

    //! number of fluid phases and number of volume fractions
    int numfluidphases_;
    int numvolfrac_;

    //! number of scalars (cont and art)
    int numscalcont_;
    int numscalart_;

    //! dof-set number of porofluid (either 0 or 2)
    int nds_porofluid_;

    //! stores the number of Gauss points for that element
    int n_gp_;

    //! stores the artery element length in reference configuration
    double arteryelelengthref_;

    //! stores initial direction of artery element
    LINALG::Matrix<numdim_, 1> lambda0_;

    //! Jacobian determinant for integration segment = L/2.0*(eta_a - eta_b)/2.0
    double jacobi_;

    //! Gausspoints in solid
    std::vector<std::vector<double>> xi_;

    //! Gausspoints in arteries and weights
    std::vector<double> eta_;
    std::vector<double> wgp_;

    //! eta_s
    std::vector<double> eta_s_;

    //! primary variables of cont and artery elements
    std::vector<double> contelephinp_;
    std::vector<double> artelephinp_;

    //! nodal artery pressure values
    LINALG::Matrix<numnodesart_, 1> earterypressurenp_;

    //! nodal artery-scalar values
    std::vector<LINALG::Matrix<numnodesart_, 1>> eartscalarnp_;

    //! nodal continuous-scalar values for scatra coupling
    std::vector<LINALG::Matrix<numnodescont_, 1>> econtscalarnp_;

    //! penalty parameter
    double pp_;

    //! start and end of integration segment
    double eta_a_;
    double eta_b_;

    //! length of integration segment int current configuration
    double curr_segment_length_;

    //! check if constant part (i.e. GPTS and MP part) has already been evaluated if integration
    //! in reference configuration is performed
    bool constant_part_evaluated_;

    //! GPTS stiffness matrices
    LINALG::SerialDenseMatrix GPTS_stiffmat11_;
    LINALG::SerialDenseMatrix GPTS_stiffmat12_;
    LINALG::SerialDenseMatrix GPTS_stiffmat21_;
    LINALG::SerialDenseMatrix GPTS_stiffmat22_;

    //! mortar coupling matrices
    LINALG::SerialDenseMatrix D_;
    LINALG::SerialDenseMatrix M_;
    LINALG::SerialDenseVector Kappa_;

    //! (dX/dxi)^-1
    std::vector<LINALG::Matrix<numdim_, numdim_>> invJ_;

    //! phase manager of the fluid
    Teuchos::RCP<DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface> phasemanager_;

    //! variable manager of the fluid
    Teuchos::RCP<DRT::ELEMENTS::POROFLUIDMANAGER::VariableManagerInterface<numdim_, numnodescont_>>
        variablemanager_;

    //! scale and function vectors
    std::vector<std::vector<int>> scale_vec_;
    std::vector<std::vector<DRT::UTILS::VariableExprFunction*>> funct_vec_;

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
  };

}  // namespace POROMULTIPHASESCATRA


#endif /* SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_PAIR_H_ */
