/*----------------------------------------------------------------------*/
/*! \file
\brief This file is used to manage the homogenized constraint mixture during growth and remodeling

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_GROWTHREMODEL_ELASTHYPER_HPP
#define FOUR_C_MAT_GROWTHREMODEL_ELASTHYPER_HPP


#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_linalg_FADmatrix_utils.hpp"
#include "4C_mat_anisotropy.hpp"
#include "4C_mat_membrane_material_interfaces.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration due to avoid header definition
namespace MAT
{
  namespace ELASTIC
  {
    class Summand;
    class RemodelFiber;
  }  // namespace ELASTIC

  // forward declaration
  class GrowthRemodelElastHyper;

  namespace PAR
  {
    class GrowthRemodelElastHyper : public CORE::MAT::PAR::Parameter
    {
      friend class MAT::GrowthRemodelElastHyper;

     public:
      /// standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      GrowthRemodelElastHyper(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{

      /// length of remodelfiber material list
      const int nummat_remodelfiber_;

      /// length of 3d elastin matrix material list
      const int nummat_elastiniso_;

      /// length of membrane elastin matrix material list
      const int nummat_elastinmem_;

      /// the list of remodelfiber material IDs
      const std::vector<int> matids_remodelfiber_;

      /// the list of 3d elastin matrix material IDs
      const std::vector<int> matids_elastiniso_;

      /// the list of membrane elastin matrix material IDs
      const std::vector<int> matids_elastinmem_;

      /// material ID of growth penalty material
      const int matid_penalty_;

      /// initial mass fraction of elastin in constraint mixture
      const double init_w_el_;

      /// material mass density
      const double density_;

      /// circumferential prestretch of elastin matrix
      const double lamb_prestretch_cir_;

      /// axial prestretch of elastin matrix
      const double lamb_prestretch_ax_;

      /// reference wall thickness of the idealized cylindrical aneurysm
      const double t_ref_;

      /// mean blood pressure
      const double p_mean_;

      /// inner radius of the idealized cylindrical aneurysm
      const double ri_;

      /// flag for turning elastin damage on or off: 1: damage on; 0: damage off
      const int damage_;

      /// flag to decide what type of collagen growth is used: 1: anisotropic growth; 0: isotropic
      /// growth
      const int growthtype_;

      /// flag to decide what type of local time integration scheme is used: 1: Backward Euler
      /// Method; 0: Forward Euler Method
      const int loctimeint_;

      /// Indicator whether Hex or Membrane elements are used ( Membrane: 1, Hex: Everything else )
      const int membrane_;

      /// Indicator whether geometry is a cylinder and that we want to calculate the AXI, CIR and
      /// RAD-directions with the help of the location of the center of each element in the
      /// reference configuration (cylinder_ == 1: cylinder aligned in x-direction; cylinder_ == 2:
      /// cylinder aligned in y-direction cylinder_ == 3: cylinder aligned in z-direction)
      const int cylinder_;

      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<CORE::MAT::Material> CreateMaterial() override;

    };  // class GrowthRemodelElastHyper

  }  // namespace PAR

  class GrowthRemodelElastHyperType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "GrowthRemodel_ElastHyperType"; }

    static GrowthRemodelElastHyperType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static GrowthRemodelElastHyperType instance_;
  };


  /*----------------------------------------------------------------------*/
  class Material;

  class GrowthRemodelElastHyper : public So3Material,
                                  public MembraneMaterialGlobalCoordinates,
                                  public MembraneMaterialInelasticThickness
  {
   public:
    /// construct empty material object
    GrowthRemodelElastHyper();

    /// construct the material object given material parameters
    explicit GrowthRemodelElastHyper(MAT::PAR::GrowthRemodelElastHyper* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of parobject.H (this file) and should return it in this method.
    int UniqueParObjectId() const override
    {
      return GrowthRemodelElastHyperType::Instance().UniqueParObjectId();
    }

    /// \brief Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
    ///
    /// \param data (in/out): char vector to store class information
    void Pack(CORE::COMM::PackBuffer& data) const override;

    /// \brief Unpack data from a char vector into this class
    ///
    /// The vector data contains all information to rebuild the
    /// exact copy of an instance of a class on a different processor.
    /// The first entry in data has to be an integer which is the unique
    /// parobject id defined at the top of this file and delivered by
    /// UniqueParObjectId().
    ///
    /// \param data (in) : vector storing all data to be unpacked into this
    ///                    instance.
    void Unpack(const std::vector<char>& data) override;

    //@}

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::KinemType::nonlinearTotLag))
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    /// material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_growthremodel_elasthyper;
    }

    /// return copy of this material object
    Teuchos::RCP<CORE::MAT::Material> Clone() const override
    {
      return Teuchos::rcp(new GrowthRemodelElastHyper(*this));
    }

    /// material mass density
    double Density() const override { return params_->density_; }

    /// hyperelastic stress response plus elasticity tensor
    void Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const CORE::LINALG::Matrix<6, 1>* glstrain,          ///< Green-Lagrange strain
        Teuchos::ParameterList& params,      ///< Container for additional information
        CORE::LINALG::Matrix<6, 1>* stress,  ///< 2nd Piola-Kirchhoff stresses
        CORE::LINALG::Matrix<6, 6>* cmat,    ///< Constitutive matrix
        int gp,                              ///< Gauss point
        int eleGID) override;                ///< Element ID


    /// setup
    void Setup(int numgp, INPUT::LineDefinition* linedef) override;

    /*!
     * \brief Post setup routine called before the first Evaluate call
     *
     * \param params Container for additional information
     * \param eleGID Global element id
     */
    void post_setup(Teuchos::ParameterList& params, int eleGID) override;

    /// This material uses the extended update call
    bool UsesExtendedUpdate() override { return true; }

    /// update
    void Update(CORE::LINALG::Matrix<3, 3> const& defgrd,  ///< Deformation gradient
        int const gp,                                      ///< Current Gauss-Point
        Teuchos::ParameterList& params,                    ///< Container for additional information
        int const eleGID) override;                        ///< Element ID

    /// Return quick accessible material parameter data
    CORE::MAT::PAR::Parameter* Parameter() const override { return params_; }

    /// hyperelastic stress response plus elasticity tensor for membrane element (membrane
    /// formulation)
    void EvaluateMembrane(CORE::LINALG::Matrix<3, 3> const&
                              defgrd_glob,      ///< Deformation gradient in global coordinates
        Teuchos::ParameterList& params,         ///< Container for additional information
        CORE::LINALG::Matrix<3, 3>& pk2M_glob,  ///< 2nd Piola-Kirchhoff stress global coordinates
        CORE::LINALG::Matrix<6, 6>& cmat_glob,  ///< Elasticity tensor in global coordinates
        int gp,                                 ///< Gauss point
        int eleGID) override;

    /// hyperelastic stress response plus elasticity tensor for membrane element (membrane
    /// formulation)
    double evaluate_membrane_thickness_stretch(
        CORE::LINALG::Matrix<3, 3> const&
            defgrd_glob,                 ///< Deformation gradient in global coordinates
        Teuchos::ParameterList& params,  ///< Container for additional information
        int gp,                          ///< Gauss point
        int eleGID) override;            ///< Element ID

    /// Return names of visualization data
    void VisNames(std::map<std::string, int>& names) override;

    /// Return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID) override;

   private:
    /// Setup circumferential, radial and axial structural tensor
    void setup_axi_cir_rad_structural_tensor(INPUT::LineDefinition* linedef);

    /// Setup prestretch (optional: setup element axi-, circ-, and rad-directions) for 3D elements
    void setup_g_r3_d(CORE::LINALG::Matrix<3, 3> const* const defgrd,  ///< Deformation gradient
        Teuchos::ParameterList& params,  ///< Container for additional information
        const double dt,                 ///< Time step size
        const int gp,                    ///< Current Gauss-Point
        const int eleGID);               ///< Element ID

    /// Setup prestretch (optional: setup element axi-, circ-, and rad-directions) for 2D elements
    /// -> membrane
    void setup_g_r2_d(Teuchos::ParameterList& params,  ///< Container for additional information
        const double dt,                               ///< Time step size
        const int gp);                                 ///< Current Gauss-Point

    /// Calculates AXI, CIR and RAD structural tensors and sets new fiber directions in the case of
    /// a cylinder
    void setup_axi_cir_rad_cylinder(CORE::LINALG::Matrix<3, 1> elecenter, double const dt);

    /// Setup anisotropic growth tensors. Here we assume that the growth direction corresponds with
    /// the radial/ thickness direction
    void setup_aniso_growth_tensors();

    /// Read AXI CIR RAD direction
    void read_dir(INPUT::LineDefinition* linedef, const std::string& specifier,
        CORE::LINALG::Matrix<3, 1>& dir);

    /// Evaluate Prestretches
    void evaluate_prestretch(
        CORE::LINALG::Matrix<3, 3> const* const defgrd, int const gp, int const eleGID);

    /// Internal Newton to implicitly solve for current mass density and inelastic remodeling
    /// stretch of each fiber family
    void solve_for_rho_lambr(
        CORE::LINALG::SerialDenseMatrix& K_T,  ///< Tangent stiffness matrix of internal Newton
        CORE::LINALG::Matrix<3, 3>& FgM,       ///< Growth deformation gradient
        CORE::LINALG::Matrix<3, 3>& iFgM,      ///< Inverse growth deformation gradient
        CORE::LINALG::Matrix<3, 3>&
            dFgdrhoM,  ///< Derivative of growth deformation gradient w.r.t. the mass density
        CORE::LINALG::Matrix<3, 3>& diFgdrhoM,  ///< Derivative of inverse growth deformation
                                                ///< gradient w.r.t. the mass density
        CORE::LINALG::Matrix<3, 3> const* const defgrd,  ///< Deformation gradient
        double const& dt,                                ///< Time step size
        int const gp,                                    ///< Current Gauss-Point
        int const eleGID);                               ///< Element ID

    /// Internal Newton to implicitly solve for current mass density and inelastic remodeling
    /// stretch of each fiber family
    void solve_fordrhod_cdlambrd_c(
        std::vector<CORE::LINALG::Matrix<1, 6>>&
            drhodC,  ///< Derivative of collagen mass density w.r.t. right Cauchy Green tensor
        std::vector<CORE::LINALG::Matrix<1, 6>>&
            dlambrdC,  ///< Derivative of inelastic remodel stretch w.r.t. right Cauchy Green tensor
        CORE::LINALG::Matrix<1, 6>& sum_drhodC,  ///< Sum of derivatives of all collagen mass
                                                 ///< density w.r.t. right Cauchy Green tensor
        CORE::LINALG::SerialDenseMatrix& K_T,    ///< Tangent stiffness matrix of internal Newton
        CORE::LINALG::Matrix<3, 3> const& iFgM,  ///< Inverse growth deformation gradient
        CORE::LINALG::Matrix<3, 3> const* const defgrd,  ///< Deformation gradient
        double const& dt,                                ///< Time step size
        int const gp,                                    ///< Current Gauss-Point
        int const eleGID) const;                         ///< Element ID

    /// Evaluates some kinematic quantities which are used in stress and elasticity tensor
    /// calculation
    static void evaluate_kin_quant_elast(
        CORE::LINALG::Matrix<3, 3> const* const defgrd,  ///< Deformation gradient
        CORE::LINALG::Matrix<3, 3> const& iFinM,         ///< Inverse inelastic deformation gradient
        CORE::LINALG::Matrix<6, 1>& iCinv,        ///< Inverse inelastic right Cauchy-Green tensor
        CORE::LINALG::Matrix<6, 1>& iCinCiCinv,   ///< C_{in}^{-1} * C * C_{in}^{-1}
        CORE::LINALG::Matrix<6, 1>& iCv,          ///< Inverse right Cauchy-Green tensor
        CORE::LINALG::Matrix<3, 3>& iCinCM,       ///< C_{in}^{-1} * C
        CORE::LINALG::Matrix<3, 3>& iFinCeM,      ///< F_{in}^{-1} * C_e
        CORE::LINALG::Matrix<9, 1>& CiFin9x1,     ///< C * F_{in}^{-1}
        CORE::LINALG::Matrix<9, 1>& CiFinCe9x1,   ///< C * F_{in}^{-1} * C_e
        CORE::LINALG::Matrix<9, 1>& CiFiniCe9x1,  ///< C * F_{in}^{-1} * C_e^{-1}
        CORE::LINALG::Matrix<3, 1>&
            prinv,      ///< Principal invariants of elastic right Cauchy-Green tensor
        int const gp);  ///< Current Gauss-Point

    /// calculate modified invariants
    void invariants_modified(CORE::LINALG::Matrix<3, 1>& modinv,  ///< modified invariants
        CORE::LINALG::Matrix<3, 1> const& prinv) const;           ///< principal invariants

    /// calculates the derivatives of the hyperelastic laws with respect to the (modified)
    /// invariants
    void evaluate_invariant_derivatives(
        CORE::LINALG::Matrix<3, 1> const&
            prinv,                         ///< Principal invariants of right Cauchy-Green tensor
        CORE::LINALG::Matrix<3, 1>& dPIw,  ///< First derivative with respect to invariants weighted
                                           ///< with the corresponding volume fraction
        CORE::LINALG::Matrix<6, 1>& ddPIIw,  ///< Second derivative with respect to invariants
                                             ///< weighted with the corresponding volume fraction
        int const gp,                        ///< Current Gauss-Point
        int const eleGID) const;             ///< Element ID

    /// converts the derivatives with respect to the modified invariants in derivatives with respect
    /// to principal invariants. Uses the following conversions:
    ///\f[
    ///  \overline{I}_{\boldsymbol{C}} = J^{-2/3} I_{\boldsymbol{C}},
    ///\f]
    ///\f[
    ///  \overline{II}_{\boldsymbol{C}} = J^{-4/3} II_{\boldsymbol{C}},
    ///\f]
    ///\f[
    ///  J = \sqrt{III_{\boldsymbol{C}}}
    ///\f]
    void convert_mod_to_princ(CORE::LINALG::Matrix<3, 1> const&
                                  prinv,  ///< principal invariants of right Cauchy-Green tensor
        CORE::LINALG::Matrix<3, 1> const&
            dPmodI,  ///< first derivative with respect to modified invariants
        CORE::LINALG::Matrix<6, 1> const&
            ddPmodII,                     ///< second derivative with respect to modified invariants
        CORE::LINALG::Matrix<3, 1>& dPI,  ///< first derivative with respect to invariants
        CORE::LINALG::Matrix<6, 1>& ddPII) const;  ///< second derivative with respect to invariants

    /// calculates the isotropic stress and elasticity tensor for coupled configuration
    void evaluate_isotropic_princ_elast(
        CORE::LINALG::Matrix<6, 1>& stressisoprinc,  ///< 2nd Piola Kirchhoff stress
        CORE::LINALG::Matrix<6, 6>& cmatisoprinc,    ///< Elasticity tensor
        CORE::LINALG::Matrix<6, 1> const& iCinv,  ///< Inverse inelastic right Cauchy-Green tensor
        CORE::LINALG::Matrix<6, 1> const& iCinCiCinv,  ///< C_{in}^{-1} * C * C_{in}^{-1}
        CORE::LINALG::Matrix<6, 1> const& iCv,         ///< Inverse right Cauchy-Green tensor
        CORE::LINALG::Matrix<3, 1> const& gamma,       ///< Factors for stress calculation
        CORE::LINALG::Matrix<8, 1> const& delta)
        const;  ///< Factors for elasticity tensor calculation

    /// Evaluate derivative of 2nd Piola Kirchhoff stress w.r.t. the inelastic deformation gradient
    void evaluated_sdi_fg(CORE::LINALG::Matrix<6, 9>& dSdiFg,  ///< Output
        CORE::LINALG::Matrix<3, 1> const& gamma,               ///< Factors for stress calculation
        CORE::LINALG::Matrix<8, 1> const& delta,   ///< Factors for elasticity tensor calculation
        CORE::LINALG::Matrix<3, 3> const& iFinM,   ///< Inverse inelastic deformation gradient
        CORE::LINALG::Matrix<3, 3> const& iCinCM,  ///< C_{in}^{-1} * C
        CORE::LINALG::Matrix<6, 1> const& iCinv,   ///< Inverse inelastic right Cauchy-Green tensor
        CORE::LINALG::Matrix<9, 1> const& CiFin9x1,     ///< C * F_{in}^{-1}
        CORE::LINALG::Matrix<9, 1> const& CiFinCe9x1,   ///< C * F_{in}^{-1} * C_e
        CORE::LINALG::Matrix<6, 1> const& iCinCiCinv,   ///< C_{in}^{-1} * C * C_{in}^{-1}
        CORE::LINALG::Matrix<9, 1> const& CiFiniCe9x1,  ///< C * F_{in}^{-1} * C_e^{-1}
        CORE::LINALG::Matrix<6, 1> const& iCv,          ///< Inverse right Cauchy-Green tensor
        CORE::LINALG::Matrix<3, 3> const& iFinCeM,      ///< F_{in}^{-1} * C_e
        int const gp) const;                            ///< Current Gauss-Point

    /// Evaluate additional terms for the elasticity tensor
    /// Additional terms are caused by implicit implementation of the evolution equation
    void evaluate_additional_cmat(
        CORE::LINALG::Matrix<6, 6>& cmatadd,          ///< Additional elasticity tensor
        CORE::LINALG::Matrix<3, 3> const& diFgdrhoM,  ///< Derivative of inverse growth deformation
                                                      ///< gradient w.r.t. the mass density
        CORE::LINALG::Matrix<1, 6> const&
            sum_drhodC,  ///< Sum over all fiber families of the derivatives of the individual mass
                         ///< density w.r.t. right Cauchy Green tensor
        CORE::LINALG::Matrix<6, 9> const&
            dSdiFg,           ///< Derivative of 2nd Piola Kirchhoff stress w.r.t.
                              ///< the inverse inelastic growth deformation gradient
        int const gp) const;  ///< Current Gauss-Point

    /// Evaluate current individual volume of elastin
    void evaluate_elastin_damage(
        double const& dt_pre);  ///< time step size of first time step in simulation (initialization
                                ///< and prestretching)

    /// build growth deformation gradient
    void evaluate_growth_def_grad(CORE::LINALG::Matrix<3, 3>& FgM,  ///< Growth deformation gradient
        CORE::LINALG::Matrix<3, 3>& iFgM,  ///< Inverse growth deformation gradient
        CORE::LINALG::Matrix<3, 3>&
            dFgdrhoM,  ///< Derivative of growth deformation gradient w.r.t. the mass density
        CORE::LINALG::Matrix<3, 3>& diFgdrhoM,  ///< Derivative of inverse growth deformation
                                                ///< gradient w.r.t. the mass density
        const int gp);                          ///< Current Gauss-Point

    /// Evaluates stress and cmat during prestressing
    void evaluate_stress_cmat_iso(
        CORE::LINALG::Matrix<3, 3> const* const defgrd,  ///< Deformation gradient
        CORE::LINALG::Matrix<3, 3> const& iFinM,         ///< Inverse inelastic deformation gradient
        CORE::LINALG::Matrix<6, 1>& stressiso,           ///< Isotropic stress tensor
        CORE::LINALG::Matrix<6, 6>& cmatiso,             ///< Isotropic stiffness matrix
        CORE::LINALG::Matrix<6, 9>& dSdiFg,  ///< Derivative of 2nd Piola Kirchhoff stress w.r.t.
                                             ///< the inverse inelastic growth deformation gradient
        int const gp,                        ///< Current Gauss-Point
        int const eleGID) const;             ///< Element ID

    /// Contribution of membrane material to 2nd Piola-Kirchhoff stress and elasticity tensor
    void evaluate_stress_cmat_membrane(
        CORE::LINALG::Matrix<3, 3> const& CM,    ///< Right Cauchy Green tensor
        CORE::LINALG::Matrix<3, 3> const& iFgM,  ///< Inelastic growth deformation gradient
        CORE::LINALG::Matrix<6, 1>& stress,      ///< 2nd Piola-Kirchhoff stress
        CORE::LINALG::Matrix<6, 6>& cmat,        ///< Elasticity tensor
        CORE::LINALG::Matrix<6, 9>& dSdiFg,  ///< Derivative of 2nd Piola Kirchhoff stress w.r.t.
                                             ///< the inelastic growth deformation gradient
        const int gp,                        ///< Current Gauss-Point
        const int eleGID) const;             ///< Element ID

    /// My material parameters
    MAT::PAR::GrowthRemodelElastHyper* params_;

    /// Map to remodelfiber material summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::RemodelFiber>> potsumrf_;

    /// Map to elastin 3d matrix material summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>> potsumeliso_;

    /// Map to membrane elastin material summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>> potsumelmem_;

    /// Map to penalty elastin matrix material summands
    Teuchos::RCP<MAT::ELASTIC::Summand> potsumelpenalty_;

    /// Current individual mass density of elastin
    std::vector<double> cur_rho_el_;

    /// Initial individual mass density of elastin
    std::vector<double> init_rho_el_;

    /// Current volume change induced by growth
    std::vector<double> v_;

    /// Circumferential structural tensor
    CORE::LINALG::Matrix<3, 3> acir_m_;

    /// Axial structural tensor
    CORE::LINALG::Matrix<3, 3> aax_m_;

    /// Radial structural tensor
    CORE::LINALG::Matrix<3, 3> arad_m_;

    /// Radial structural tensor in "stress-like" Voigt notation
    CORE::LINALG::Matrix<6, 1> aradv_;

    /// Prestretch of elastin matrix in axial, circumferential and radial direction (used for
    /// prestressing)
    std::vector<CORE::LINALG::Matrix<3, 3>> gm_;

    /// Total simulation time
    double t_tot_;

    /// Total number of remodel fibers
    unsigned nr_rf_tot_;

    ///                       (rad_x, axi_x, circ_x)
    /// Cylinder coordinates  (rad_y, axi_y, circ_y)
    ///                       (rad_z, axi_z, circ_z)
    CORE::LINALG::Matrix<3, 3> radaxicirc_;

    /// Axial coordinate of each Gauss-Point
    std::vector<double> gp_ax_;

    /// Radial coordinate of each Gauss-Point
    std::vector<double> gp_rad_;

    /// Structural tensor of growth direction
    CORE::LINALG::Matrix<3, 3> ag_m_;

    /// Structural tensor of the plane in which all fibers are located
    CORE::LINALG::Matrix<3, 3> apl_m_;

    /// Fraction of 2D material parameter of elastin
    std::vector<double> mue_frac_;

    /// First call of Evaluate()
    std::vector<int> setup_;

    /// Holder for anisotropic behavior
    Anisotropy anisotropy_;
  };

}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
