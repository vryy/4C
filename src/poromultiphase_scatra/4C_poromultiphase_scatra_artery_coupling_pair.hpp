/*----------------------------------------------------------------------*/
/*! \file
 \brief one pair consisting of exactly one artery element and one poro-
        multiphase-scatra element which might be tied to each other

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_PAIR_HPP
#define FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_PAIR_HPP

#include "4C_config.hpp"

#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_inpar_bio.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_porofluidmultiphase_ele_phasemanager.hpp"
#include "4C_porofluidmultiphase_ele_variablemanager.hpp"

#include <Epetra_Vector.h>
#include <Sacado.hpp>
#include <Teuchos_RCP.hpp>

// define Fad object for evaluation
typedef Sacado::Fad::DFad<double> FAD;

// forward declaration
class Epetra_MultiVector;

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Elements
{
  class Element;
}

namespace Mat
{
  class MatList;
  class Cnst1dArt;
}  // namespace Mat
namespace Core
{
  namespace LinAlg
  {
    class SerialDenseVector;
    class SerialDenseMatrix;
  }  // namespace LinAlg
  namespace UTILS
  {
    class FunctionOfAnything;
  }
}  // namespace Core
namespace PoroMultiPhaseScaTra
{
  class PoroMultiPhaseScatraArteryCouplingPairBase
  {
   public:
    //! constructor
    PoroMultiPhaseScatraArteryCouplingPairBase() { return; };

    //! destructor
    virtual ~PoroMultiPhaseScatraArteryCouplingPairBase() = default;

    //! Init
    virtual void init(std::vector<Core::Elements::Element const*> elements,
        const Teuchos::ParameterList& couplingparams,
        const Teuchos::ParameterList& fluidcouplingparams, const std::vector<int>& coupleddofs_cont,
        const std::vector<int>& coupleddofs_art, const std::vector<std::vector<int>>& scale_vec,
        const std::vector<std::vector<int>>& funct_vec, const std::string condnames,
        const double penalty, const std::string couplingtype = "", const int eta_ntp = 0) = 0;

    //! query if pair active
    virtual bool is_active() = 0;

    //! things that need to be done in a separate loop before the actual evaluation loop
    //! over all coupling pairs
    virtual void pre_evaluate(Teuchos::RCP<Epetra_MultiVector> gp_vector) = 0;

    //! things that need to be done in a separate loop before the actual evaluation loop
    //! over all coupling pairs
    virtual void delete_unnecessary_g_ps(Teuchos::RCP<Epetra_MultiVector> gp_vector) = 0;

    /*!
     * @brief Evaluate this pair
     *
     * @returns integral of diameter of the segment
     */
    virtual double evaluate(Core::LinAlg::SerialDenseVector* forcevec1,
        Core::LinAlg::SerialDenseVector* forcevec2, Core::LinAlg::SerialDenseMatrix* stiffmat11,
        Core::LinAlg::SerialDenseMatrix* stiffmat12, Core::LinAlg::SerialDenseMatrix* stiffmat21,
        Core::LinAlg::SerialDenseMatrix* stiffmat22, Core::LinAlg::SerialDenseMatrix* D_ele,
        Core::LinAlg::SerialDenseMatrix* M_ele, Core::LinAlg::SerialDenseVector* Kappa_ele,
        const std::vector<double>& segmentlengths) = 0;

    //! evaluate additional linearization of (integrated) element diameter dependent terms
    //! (Hagen-Poiseuille)
    virtual void evaluate_additional_linearizationof_integrated_diam(
        Core::LinAlg::SerialDenseMatrix* stiffmat11,
        Core::LinAlg::SerialDenseMatrix* stiffmat12) = 0;

    //! flag if diameter function is active, i.e., varying diameter linearization need to be
    //! calculated
    virtual bool DiamFunctionActive() = 0;

    //! reset state
    virtual void ResetState(Teuchos::RCP<Core::FE::Discretization> contdis,
        Teuchos::RCP<Core::FE::Discretization> artdis) = 0;

    /**
     * Setup the porofluid-managers and the materials for later evaluation
     * @param[in] disname: name of continuous discretization
     * @param[in] timefacrhs_art: right hand side factor for artery time integration
     * @param[in] timefacrhs_cont: right hand side factor for time integration of 2D/3D
     * discretization
     */
    virtual void setup_fluid_managers_and_materials(
        const std::string disname, const double& timefacrhs_art, const double& timefacrhs_cont) = 0;

    //! beginning of integration segment
    virtual double Etadata() const = 0;
    //! end of integration segment
    virtual double EtaB() const = 0;

    //! element 1 (= artery) GID
    virtual int Ele1GID() const = 0;
    //! element 2 (= cont) GID
    virtual int Ele2GID() const = 0;

    //! apply mesh movement on artery element
    virtual double ApplyMeshMovement(
        const bool firstcall, Teuchos::RCP<Core::FE::Discretization> contdis) = 0;

    //! set segment id
    virtual void SetSegmentID(const int& segmentid) = 0;
    //! get segment id
    virtual int GetSegmentID() const = 0;

    //! get the volume of the 2D/3D element
    virtual double CalculateVol2D3D() const = 0;

    //! get number of Gauss points
    virtual int num_gp() const = 0;

    //! type of coupling pair
    enum CouplingType
    {
      type_undefined,
      type_porofluid,  //!< porofluid
      type_scatra      //!< scatra
    };
  };

  //! the coupling pair
  template <Core::FE::CellType distype_art, Core::FE::CellType distype_cont, int dim>
  class PoroMultiPhaseScatraArteryCouplingPair : public PoroMultiPhaseScatraArteryCouplingPairBase
  {
   public:
    //! constructor
    PoroMultiPhaseScatraArteryCouplingPair();

    //! Init
    void init(std::vector<Core::Elements::Element const*> elements,
        const Teuchos::ParameterList& couplingparams,
        const Teuchos::ParameterList& fluidcouplingparams, const std::vector<int>& coupleddofs_cont,
        const std::vector<int>& coupleddofs_art, const std::vector<std::vector<int>>& scale_vec,
        const std::vector<std::vector<int>>& funct_vec, const std::string condname,
        const double penalty, const std::string couplingtype = "", const int eta_ntp = 0) override;

    //! query if pair active
    bool is_active() override { return isactive_; }

    //! things that need to be done in a separate loop before the actual evaluation loop
    //! over all coupling pairs
    void pre_evaluate(Teuchos::RCP<Epetra_MultiVector> gp_vector) override;

    //! things that need to be done in a separate loop before the actual evaluation loop
    //! over all coupling pairs
    void delete_unnecessary_g_ps(Teuchos::RCP<Epetra_MultiVector> gp_vector) override;

    //! flag if diameter function is active, i.e., varying diameter linearization need to be
    //! calculated
    bool DiamFunctionActive() override { return diam_funct_active_; }

    //! reset state
    void ResetState(Teuchos::RCP<Core::FE::Discretization> contdis,
        Teuchos::RCP<Core::FE::Discretization> artdis) override;

    /**
     * Setup the porofluid-managers and the materials for later evaluation
     * @param[in] disname: name of continuous discretization
     * @param[in] timefacrhs_art: right hand side factor for artery time integration
     * @param[in] timefacrhs_cont: right hand side factor for time integration of 2D/3D
     * discretization
     */
    void setup_fluid_managers_and_materials(const std::string disname, const double& timefacrhs_art,
        const double& timefacrhs_cont) override;

    /*!
     * @brief Evaluate this pair
     *
     * @returns integral of diameter of the segment
     */
    double evaluate(Core::LinAlg::SerialDenseVector* forcevec1,
        Core::LinAlg::SerialDenseVector* forcevec2, Core::LinAlg::SerialDenseMatrix* stiffmat11,
        Core::LinAlg::SerialDenseMatrix* stiffmat12, Core::LinAlg::SerialDenseMatrix* stiffmat21,
        Core::LinAlg::SerialDenseMatrix* stiffmat22, Core::LinAlg::SerialDenseMatrix* D_ele,
        Core::LinAlg::SerialDenseMatrix* M_ele, Core::LinAlg::SerialDenseVector* Kappa_ele,
        const std::vector<double>& segmentlengths) override;

    //! evaluate additional linearization of (integrated) element diameter dependent terms
    //! (Hagen-Poiseuille)
    void evaluate_additional_linearizationof_integrated_diam(
        Core::LinAlg::SerialDenseMatrix* stiffmat11,
        Core::LinAlg::SerialDenseMatrix* stiffmat12) override;

    //! beginning and end of integration segment
    double Etadata() const override { return eta_a_; }
    double EtaB() const override { return eta_b_; }

    //! element 1 (= artery) GID
    int Ele1GID() const override;
    //! element 2 (= cont) GID
    int Ele2GID() const override;

    //! number of GP
    int num_gp() const override { return n_gp_; };

    //! apply mesh movement on artery element
    double ApplyMeshMovement(
        const bool firstcall, Teuchos::RCP<Core::FE::Discretization> contdis) override;

    //! set segment id
    void SetSegmentID(const int& segmentid) override;
    //! get segment id
    int GetSegmentID() const override;

    //! get the volume of the 2D/3D element
    double CalculateVol2D3D() const override;

   private:
    // static variables
    //! number of nodes of 1D artery element
    static constexpr unsigned numnodesart_ = Core::FE::num_nodes<distype_art>;
    //! number of nodes of 2D/3D continuous element
    static constexpr unsigned numnodescont_ = Core::FE::num_nodes<distype_cont>;
    //! number of nodes of spatial dimensions
    static constexpr unsigned numdim_ = Core::FE::dim<distype_cont>;

    //! set time factor needed for evaluation of right hand side (function coupling) terms
    void set_time_fac_rhs(const double& arterydensity, Teuchos::RCP<Mat::MatList> contscatramat,
        const double& timefacrhs_art, const double& timefacrhs_cont);

    //! pre-evaluate for lateral surface coupling
    void pre_evaluate_lateral_surface_coupling(Teuchos::RCP<Epetra_MultiVector> gp_vector);

    //! pre-evaluate for centerline coupling
    void pre_evaluate_centerline_coupling();

    //! pre-evaluate for node-to-point coupling
    void pre_evaluate_node_to_point_coupling();

    //! extract velocity of solid phase
    void extract_solid_vel(Teuchos::RCP<Core::FE::Discretization> contdis);

    //! recompute if deformable arteries are assumed
    void recompute_eta_and_xi_in_deformed_configuration(const std::vector<double>& segmentlengths,
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
    void create_integration_segment();

    //! get all intersections of artery element with 2D/3D element
    std::vector<double> get_all_inter_sections();

    //! project a Gauss point on 1D element into 2D/3D element
    template <typename T>
    void projection(
        Core::LinAlg::Matrix<numdim_, 1, T>& r1, std::vector<T>& xi, bool& projection_valid);

    //! Check for duplicate projections
    bool projection_not_yet_found(const std::vector<double>& intersections, const double& eta);

    //! Intersect artery element with edges (2D) or surfaces (3D) of element
    void inter_sect_with2_d3_d(std::vector<double>& xi, double& eta, const int& fixedPar,
        const double& fixedAt, bool& projection_valid);

    //! get 1D shapefunctions at eta
    template <typename T>
    void get1_d_shape_functions(Core::LinAlg::Matrix<1, numnodesart_, T>& N1,
        Core::LinAlg::Matrix<1, numnodesart_, T>& N1_eta, const T& eta);

    //! get 2D/3D shapefunctions at xi1, xi2 (, xi3)
    template <typename T>
    void get2_d3_d_shape_functions(Core::LinAlg::Matrix<1, numnodescont_, T>& N2,
        Core::LinAlg::Matrix<numdim_, numnodescont_, T>& N2_xi, const std::vector<T>& xi);

    //! compute artery coordinates and derivatives in reference configuration
    template <typename T>
    void compute_artery_coords_and_derivs_ref(Core::LinAlg::Matrix<numdim_, 1, T>& r1,
        Core::LinAlg::Matrix<numdim_, 1, T>& r1_eta,
        const Core::LinAlg::Matrix<1, numnodesart_, T>& N1,
        const Core::LinAlg::Matrix<1, numnodesart_, T>& N1_eta);

    //! compute 2D/3D coordinates and derivatives in reference configuration
    template <typename T>
    void compute2_d3_d_coords_and_derivs_ref(Core::LinAlg::Matrix<numdim_, 1, T>& x2,
        Core::LinAlg::Matrix<numdim_, numdim_, T>& x2_xi,
        const Core::LinAlg::Matrix<1, numnodescont_, T>& N2,
        const Core::LinAlg::Matrix<numdim_, numnodescont_, T>& N2_xi);

    //! evaluate the function coupling (return integral of diameter of the segment)
    void evaluate_function_coupling(const std::vector<double>& eta,
        const std::vector<std::vector<double>>& xi, const std::vector<double>& segmentlengths,
        Core::LinAlg::SerialDenseVector* forcevec1, Core::LinAlg::SerialDenseVector* forcevec2,
        Core::LinAlg::SerialDenseMatrix* stiffmat11, Core::LinAlg::SerialDenseMatrix* stiffmat12,
        Core::LinAlg::SerialDenseMatrix* stiffmat21, Core::LinAlg::SerialDenseMatrix* stiffmat22,
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
    void evaluated_nds_solid_vel(const std::vector<double>& eta,
        const std::vector<std::vector<double>>& xi, const std::vector<double>& segmentlengths,
        Core::LinAlg::SerialDenseVector& forcevec1, const double& etaA, const double& etaB);

    //! evaluate stiffness for GPTS case
    void evaluate_gpts_stiff(const double& w_gp, const Core::LinAlg::Matrix<1, numnodesart_>& N1,
        const Core::LinAlg::Matrix<1, numnodescont_>& N2, const double& jacobi, const double& pp);

    //! evaluate stiffness for NTP case
    void evaluate_ntp_stiff(const Core::LinAlg::Matrix<1, numnodesart_>& N1,
        const Core::LinAlg::Matrix<1, numnodescont_>& N2, const double& pp);

    //! evaluate mortar coupling matrices D and M
    void evaluate_dm_kappa(const double& w_gp, const Core::LinAlg::Matrix<1, numnodesart_>& N1,
        const Core::LinAlg::Matrix<1, numnodescont_>& N2, const double& jacobi);

    //! evalute GPTS
    void evaluate_gpts(const std::vector<double>& eta, const std::vector<std::vector<double>>& xi,
        const std::vector<double>& segmentlengths, Core::LinAlg::SerialDenseVector* forcevec1,
        Core::LinAlg::SerialDenseVector* forcevec2, Core::LinAlg::SerialDenseMatrix* stiffmat11,
        Core::LinAlg::SerialDenseMatrix* stiffmat12, Core::LinAlg::SerialDenseMatrix* stiffmat21,
        Core::LinAlg::SerialDenseMatrix* stiffmat22);

    //! evalute NTP
    void evaluate_ntp(const std::vector<double>& eta, const std::vector<std::vector<double>>& xi,
        Core::LinAlg::SerialDenseVector* forcevec1, Core::LinAlg::SerialDenseVector* forcevec2,
        Core::LinAlg::SerialDenseMatrix* stiffmat11, Core::LinAlg::SerialDenseMatrix* stiffmat12,
        Core::LinAlg::SerialDenseMatrix* stiffmat21, Core::LinAlg::SerialDenseMatrix* stiffmat22);

    //! evaluate mortar coupling matrices D and M
    void evaluate_dm_kappa(const std::vector<double>& eta,
        const std::vector<std::vector<double>>& xi, const std::vector<double>& segmentlengths,
        Core::LinAlg::SerialDenseMatrix* D_ele, Core::LinAlg::SerialDenseMatrix* M_ele,
        Core::LinAlg::SerialDenseVector* Kappa_ele);

    //! evaluate the function coupling
    void evaluate_function_coupling(const double& w_gp,
        const Core::LinAlg::Matrix<1, numnodesart_>& N1,
        const Core::LinAlg::Matrix<1, numnodescont_>& N2, const double& jacobi,
        Core::LinAlg::SerialDenseVector& forcevec1, Core::LinAlg::SerialDenseVector& forcevec2,
        Core::LinAlg::SerialDenseMatrix& stiffmat11, Core::LinAlg::SerialDenseMatrix& stiffmat12,
        Core::LinAlg::SerialDenseMatrix& stiffmat21, Core::LinAlg::SerialDenseMatrix& stiffmat22,
        double& integrated_diam);

    //! evaluate the diameter function and derivative (for couplingtype porofluid)
    void evaluate_diam_function_and_deriv(const double artpressnpAtGP, const double& w_gp,
        const Core::LinAlg::Matrix<1, numnodesart_>& N1,
        const Core::LinAlg::Matrix<1, numnodescont_>& N2, const double& jacobi);

    //! integrate in deformed configuration from eta_A to eta_s
    FAD integrate_length_to_eta_s(const FAD& eta_s);

    //! get values of artery at GP
    void get_artery_values_at_gp(const Core::LinAlg::Matrix<1, numnodesart_>& N1, double& artpress,
        std::vector<double>& artscalar);

    //! get scalar values of continuous discretization at GP
    void get_cont_scalar_values_at_gp(
        const Core::LinAlg::Matrix<1, numnodescont_>& N2, std::vector<double>& contscalarnp);

    //! assemble the function coupling into stiffness matrix (artery-part)
    void assemble_function_coupling_into_force_stiff_art(const int& i_art, const double& w_gp,
        const Core::LinAlg::Matrix<1, numnodesart_>& N1,
        const Core::LinAlg::Matrix<1, numnodescont_>& N2, const double& jacobi, const int& scale,
        const double& functval, const std::vector<double>& artderivs,
        const std::vector<double>& contderivs, Core::LinAlg::SerialDenseVector& forcevec1,
        Core::LinAlg::SerialDenseMatrix& stiffmat11, Core::LinAlg::SerialDenseMatrix& stiffmat12);

    //! assemble the function coupling into stiffness matrix (2D/3D-part)
    void assemble_function_coupling_into_force_stiff_cont(const std::vector<int>& assembleInto,
        const double& w_gp, const Core::LinAlg::Matrix<1, numnodesart_>& N1,
        const Core::LinAlg::Matrix<1, numnodescont_>& N2, const double& jacobi, const int& scale,
        const double& timefacrhs_cont, const double& functval, const std::vector<double>& artderivs,
        const std::vector<double>& contderivs, Core::LinAlg::SerialDenseVector& forcevec2,
        Core::LinAlg::SerialDenseMatrix& stiffmat21, Core::LinAlg::SerialDenseMatrix& stiffmat22);

    //! evaluate function and its derivative

    void evaluate_function_and_deriv(const Core::UTILS::FunctionOfAnything& funct,
        const double& artpressnpAtGP, const std::vector<double>& artscalarnpAtGP,
        const std::vector<double>& scalarnpAtGP, double& functval, std::vector<double>& artderivs,
        std::vector<double>& contderivs);

    //! set scalar as constants into function
    void set_scalar_values_as_constants(std::vector<std::pair<std::string, double>>& constants,
        const std::vector<double>& artscalarnpAtGP, const std::vector<double>& scalarnpAtGP);

    //! set fluid as variables into function
    void set_fluid_values_as_variables(
        std::vector<std::pair<std::string, double>>& variables, const double& artpressnpAtGP);

    //! set fluid as constants into function
    void set_fluid_values_as_constants(
        std::vector<std::pair<std::string, double>>& constants, const double& artpressnpAtGP);

    //! set scalar as variables into function
    void set_scalar_values_as_variables(std::vector<std::pair<std::string, double>>& variables,
        const std::vector<double>& artscalarnpAtGP, const std::vector<double>& scalarnpAtGP);

    //! evaluate derivatives w.r.t. fluid of function
    void evaluate_fluid_derivs(std::vector<double>& artderivs, std::vector<double>& contderivs,
        const std::vector<double>& functderivs);

    //! evaluate derivatives w.r.t. scalar of function
    void evaluate_scalar_derivs(std::vector<double>& artderivs, std::vector<double>& contderivs,
        const std::vector<double>& functderivs);

    //! evaluate force for GPTS  or NTP case
    void evaluate_gptsntp_force(Core::LinAlg::SerialDenseVector& forcevec1,
        Core::LinAlg::SerialDenseVector& forcevec2,
        const Core::LinAlg::SerialDenseMatrix& stiffmat11,
        const Core::LinAlg::SerialDenseMatrix& stiffmat12,
        const Core::LinAlg::SerialDenseMatrix& stiffmat21,
        const Core::LinAlg::SerialDenseMatrix& stiffmat22);

    //! update the stiffness for GPTS or NTP
    void update_gptsntp_stiff(Core::LinAlg::SerialDenseMatrix& stiffmat11,
        Core::LinAlg::SerialDenseMatrix& stiffmat12, Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22);

    /**
     * \brief coupling to additional porous network is only possible if we also have an
     * element with a valid volume fraction pressure, i.e., if we also have a smeared representation
     * of the neovasculature at this point if not ---> corresponding matrices are set to zero
     */
    void check_valid_volume_fraction_pressure_coupling(Core::LinAlg::SerialDenseMatrix& stiffmat11,
        Core::LinAlg::SerialDenseMatrix& stiffmat12, Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22);

    //! update the D, M and Kappa for MP
    void update_dm_kappa(Core::LinAlg::SerialDenseMatrix& D_ele,
        Core::LinAlg::SerialDenseMatrix& M_ele, Core::LinAlg::SerialDenseVector& Kappa_ele);

    //! fill the function vector
    void fill_function_vector(std::vector<const Core::UTILS::FunctionOfAnything*>& my_funct_vec,
        const std::vector<int>& funct_vec, const std::vector<int>& scale_vec);


    //! initialize a function
    void initialize_function(const Core::UTILS::FunctionOfAnything& funct);

    //! initialize names used in functions
    void initialize_function_names();

    //! initialize vector where to assemble continuous DOF functions into
    void initialize_assemble_into_cont_dof_vector();

    //! coupling type
    CouplingType coupltype_;

    //! coupling method (either GPTS or MP)
    Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod couplmethod_;

    //! name of the condition
    std::string condname_;

    //! indicates if the init() function has been called
    bool isinit_;

    //! indicates if the pre_evaluate() function has been called
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
    const Core::Elements::Element* element1_;

    //! second element of interacting pair (2D/3D element)
    const Core::Elements::Element* element2_;

    //! reference nodal positions element 1 (1D artery)
    Core::LinAlg::Matrix<numdim_ * numnodesart_, 1> ele1posref_;
    //! reference nodal positions element 2 (2D/3D continuous element)
    Core::LinAlg::Matrix<numdim_, numnodescont_> ele2posref_;

    //! current position of element 2
    Core::LinAlg::Matrix<numdim_, numnodescont_> ele2pos_;
    //! current velocity of element 2
    Core::LinAlg::Matrix<numdim_, numnodescont_> ele2vel_;

    //! reference diameter of the artery element (constant across element)
    double arterydiamref_;
    //! current diameter of the artery element at the GP
    double arterydiam_at_gp_;
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
    Core::LinAlg::Matrix<numdim_, 1> lambda0_;

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
    Core::LinAlg::Matrix<numnodesart_, 1> earterypressurenp_;

    //! nodal artery-scalar values
    std::vector<Core::LinAlg::Matrix<numnodesart_, 1>> eartscalarnp_;

    //! nodal continuous-scalar values for scatra coupling
    std::vector<Core::LinAlg::Matrix<numnodescont_, 1>> econtscalarnp_;

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
    Core::LinAlg::SerialDenseMatrix gpts_ntp_stiffmat11_;
    //! GPTS/NTP stiffness matrix (artery-cont contribution)
    Core::LinAlg::SerialDenseMatrix gpts_ntp_stiffmat12_;
    //! GPTS/NTP stiffness matrix (cont-artery contribution)
    Core::LinAlg::SerialDenseMatrix gpts_ntp_stiffmat21_;
    //! GPTS/NTP stiffness matrix (cont-cont contribution)
    Core::LinAlg::SerialDenseMatrix gpts_ntp_stiffmat22_;

    //! (varying) diameter stiffness matrix (artery-artery contribution)
    Core::LinAlg::SerialDenseMatrix diam_stiffmat11_;
    //! (varying) diameter stiffness matrix (artery-cont contribution)
    Core::LinAlg::SerialDenseMatrix diam_stiffmat12_;

    //! mortar coupling matrix D
    Core::LinAlg::SerialDenseMatrix d_;
    //! mortar coupling matrix M
    Core::LinAlg::SerialDenseMatrix m_;
    //! mortar coupling vector kappa
    Core::LinAlg::SerialDenseVector kappa_;

    //! (dX/dxi)^-1
    std::vector<Core::LinAlg::Matrix<numdim_, numdim_>> inv_j_;

    //! phase manager of the fluid
    Teuchos::RCP<Discret::ELEMENTS::PoroFluidManager::PhaseManagerInterface> phasemanager_;

    //! variable manager of the fluid
    Teuchos::RCP<
        Discret::ELEMENTS::PoroFluidManager::VariableManagerInterface<numdim_, numnodescont_>>
        variablemanager_;

    //! scale vector
    std::vector<std::vector<int>> scale_vec_;
    //! function vector
    std::vector<std::vector<const Core::UTILS::FunctionOfAnything*>> funct_vec_;

    //! diameter function
    const Core::UTILS::FunctionOfAnything* artdiam_funct_;

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
    Teuchos::RCP<Mat::Cnst1dArt> arterymat_;
  };

}  // namespace PoroMultiPhaseScaTra


FOUR_C_NAMESPACE_CLOSE

#endif
