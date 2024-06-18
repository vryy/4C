/*----------------------------------------------------------------------*/
/*! \file
\brief evaluation of a generic material whose deformation gradient is modeled to be split
multiplicatively into elastic and inelastic parts

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_MULTIPLICATIVE_SPLIT_DEFGRAD_ELASTHYPER_HPP
#define FOUR_C_MAT_MULTIPLICATIVE_SPLIT_DEFGRAD_ELASTHYPER_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_anisotropy.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN


// forward declaration
namespace Mat
{
  class Anisotropy;
  class InelasticDefgradFactors;

  namespace Elastic
  {
    class Summand;
  }

  namespace PAR
  {
    enum class InelasticSource
    {
      none,
      concentration,
      temperature
    };

    class MultiplicativeSplitDefgradElastHyper : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      explicit MultiplicativeSplitDefgradElastHyper(const Core::Mat::PAR::Parameter::Data& matdata);

      Teuchos::RCP<Core::Mat::Material> create_material() override;

      /// length of elastic material list
      const int nummat_elast_;

      /// the list of elastic material IDs
      const std::vector<int> matids_elast_;

      /// number of factors of inelastic deformation gradient F_{in} = F_{in,1} . F_{in,2}. ... .
      /// F_{in,n} (n factors)
      const int numfac_inel_;

      /// IDs of inelastic deformation gradient factors (i-th ID specifies calculation of F_{in,i})
      const std::vector<int> inel_defgradfacids_;

      /// material mass density
      const double density_;

    };  // class MultiplicativeSplitDefgrad_ElastHyper

  }  // namespace PAR

  class MultiplicativeSplitDefgradElastHyperType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "MultiplicativeSplitDefgrad_ElastHyperType"; }

    static MultiplicativeSplitDefgradElastHyperType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static MultiplicativeSplitDefgradElastHyperType instance_;
  };

  /*----------------------------------------------------------------------*/
  /*! \class InelasticFactorsHandler

     All factors (class InelasticDefgradFactors) contributing to inelastic deformation are stored
     within this class. They are initially classified according to the source of inelastic
     deformation (e.g. concentration). All contributions from one source can be returned.
  */
  class InelasticFactorsHandler
  {
   public:
    /*!
     * @brief Evaluate inelastic deformation gradient from all factors and assign to different
     * sources
     *
     * @param[in]  defgrad  Deformation gradient
     * @param[out] iFinM    inverse inelastic deformation gradient
     */
    void evaluate_inverse_inelastic_def_grad(
        const Core::LinAlg::Matrix<3, 3>* defgrad, Core::LinAlg::Matrix<3, 3>& iFinM);

    /// Returns all inelastic factors as a vector
    const std::vector<std::pair<PAR::InelasticSource, Teuchos::RCP<Mat::InelasticDefgradFactors>>>&
    FacDefGradIn() const
    {
      return facdefgradin_;
    }

    /// Return vector of all inelastic deformation gradients
    const std::vector<std::pair<PAR::InelasticSource, Core::LinAlg::Matrix<3, 3>>>& GetiFinj() const
    {
      return i_finj_;
    }

    /// total number of inelastic contributions
    int NumInelasticDefGrad() const { return static_cast<int>(facdefgradin_.size()); }

    /// Assigns the different inelastic factors to different sources
    void Setup(Mat::PAR::MultiplicativeSplitDefgradElastHyper* params);

   private:
    /// vector that holds pairs of inelastic contribution and respective source
    std::vector<std::pair<PAR::InelasticSource, Teuchos::RCP<Mat::InelasticDefgradFactors>>>
        facdefgradin_;

    /// vector that holds pairs of inelastic deformation gradients and respective source
    std::vector<std::pair<PAR::InelasticSource, Core::LinAlg::Matrix<3, 3>>> i_finj_;
  };

  /*----------------------------------------------------------------------*/
  /*! @class MultiplicativeSplitDefgrad_ElastHyper

    In this class the deformation gradient is modeled to be split multiplicatively
    in elastic and inelastic deformation gradients (F = F_{el} * F_{in}).
    The elastic contribution can be any elastic material law provided in elast_summand.
    Only elastic strains cause stresses!
    The inelastic deformation gradient itself can be a product of different inelastic
    deformation gradients, i.e. F_{in} = F_{in,1} * F_{in,2} * ... * F_{in,n}.
    The inverse inelastic deformation gradients and their derivative w.r.t. the primary variables
    that are needed to set up the system to be solved are evaluated in the derived classes
    of the interface class 'InelasticDefgradFactors'.
*/
  class MultiplicativeSplitDefgradElastHyper : public So3Material
  {
   public:
    /// construct empty material object
    MultiplicativeSplitDefgradElastHyper();

    /// construct the material object given material parameters
    explicit MultiplicativeSplitDefgradElastHyper(
        Mat::PAR::MultiplicativeSplitDefgradElastHyper* params);

    int UniqueParObjectId() const override
    {
      return MultiplicativeSplitDefgradElastHyperType::Instance().UniqueParObjectId();
    }

    void Pack(Core::Communication::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    void ValidKinematics(Inpar::STR::KinemType kinem) override
    {
      if (kinem != Inpar::STR::KinemType::nonlinearTotLag)
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_multiplicative_split_defgrad_elasthyper;
    }

    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new MultiplicativeSplitDefgradElastHyper(*this));
    }

    double Density() const override { return params_->density_; }

    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

    void evaluate(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmat, int gp,
        int eleGID) override;

    void evaluate_cauchy_n_dir_and_derivatives(const Core::LinAlg::Matrix<3, 3>& defgrd,
        const Core::LinAlg::Matrix<3, 1>& n, const Core::LinAlg::Matrix<3, 1>& dir,
        double& cauchy_n_dir, Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dn,
        Core::LinAlg::Matrix<3, 1>* d_cauchyndir_ddir, Core::LinAlg::Matrix<9, 1>* d_cauchyndir_dF,
        Core::LinAlg::Matrix<9, 9>* d2_cauchyndir_dF2,
        Core::LinAlg::Matrix<9, 3>* d2_cauchyndir_dF_dn,
        Core::LinAlg::Matrix<9, 3>* d2_cauchyndir_dF_ddir, int gp, int eleGID,
        const double* concentration, const double* temp, double* d_cauchyndir_dT,
        Core::LinAlg::Matrix<9, 1>* d2_cauchyndir_dF_dT) override;

    void evaluate_linearization_od(const Core::LinAlg::Matrix<3, 3>& defgrd, double concentration,
        Core::LinAlg::Matrix<9, 1>* d_F_dx) override;

    void Setup(int numgp, Input::LineDefinition* linedef) override;

    void Update() override;

    /*!
     * @brief Evaluate off-diagonal stiffness matrix (required for monolithic algorithms)
     *
     * @param[in] source   Source of inelastic deformation
     * @param[in] defgrad  Deformation gradient
     * @param[in] dSdiFin  Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
     *                     deformation gradient
     * @param[out] dstressdx  Derivative of 2nd Piola Kirchhoff stresses w.r.t. primary variable of
     *                        different field
     */
    void EvaluateODStiffMat(PAR::InelasticSource source, const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<6, 9>& dSdiFin, Core::LinAlg::Matrix<6, 1>& dstressdx);

    /*!
     * @brief Evaluate additional terms of the elasticity tensor
     *
     * \f[
     *   cmat_{add} = \frac{\partial\mathbf{S}}{\partial \mathbf{F}_\text{in}^{-1}}
     *             : \frac{\partial\mathbf{F}_\text{in}^{-1}}{\partial \mathbf{C}}
     * \f]
     *
     * @param[in] defgrad  Deformation gradient
     * @param[in] iCV      Inverse right Cauchy-Green tensor stored as 6x1 vectors
     * @param[in] dSdiFin  Derivative of 2nd Piola Kirchhoff stress w.r.t. the inverse inelastic
     *                     deformation gradient
     * @param[out] cmatadd Additional elasticity tensor \f$ cmat_{add} \f$
     */
    void evaluate_additional_cmat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<6, 1>& iCV, const Core::LinAlg::Matrix<6, 9>& dSdiFin,
        Core::LinAlg::Matrix<6, 6>& cmatadd);

    /*!
     * @brief  Evaluate derivative of 2nd Piola Kirchhoff stresses w.r.t. the inelastic deformation
     * gradient
     *
     * @param[in] gamma  Factors for stress calculation (see e.g. Holzapfel - Nonlinear Solid
     *                   Mechanics)
     * @param[in] delta  Factors for elasticity tensor calculation (see e.g. Holzapfel  - Nonlinear
     *                   Solid Mechanics)
     * @param[in] iFinM   Inverse inelastic deformation gradient
     * @param[in] iCinCM  \f$ \mathbf{C}_\text{in}^{-1} \cdot \mathbf{C} \f$
     * @param[in] iCinV   Inverse inelastic right Cauchy-Green tensor
     * \f$\mathbf{C}_\text{in}^{-1}\f$ stored as 6x1 vector
     * @param[in] CiFin9x1    \f$ \mathbf{C} \cdot \mathbf{F}_\text{in}^{-1} \f$ stored as 9x1
     * vector
     * @param[in] CiFinCe9x1  \f$ \mathbf{C} \cdot \mathbf{F}_\text{in}^{-1} \cdot
     * \mathbf{C}_\text{el} \f$ stored as 9x1 vector
     * @param[in] iCinCiCinV    \f$ \mathbf{C}_\text{in}^{-1} \cdot \mathbf{C} \cdot
     *                          \mathbf{C}_\text{in}^{-1} \f$ stored as 6x1 vector
     * @param[in] CiFiniCe9x1   \f$ \mathbf{C} \cdot \mathbf{F}_\text{in}^{-1} \cdot
     *                          \mathbf{C}_\text{el}^{-1} \f$ stored as 9x1 vector
     * @param[in] iCV           Inverse right Cauchy-Green tensor \f$ \mathbf{C}^{-1} \f$ stored as
     *                          6x1 vector
     * @param[in] iFinCeM       \f$ \mathbf{F}_\text{in}^{-1} \cdot \mathbf{C}_\text{el} \f$
     * @param[in] detFin        determinant of inelastic deformation gradient
     * @param[out] dSdiFin      derivative of 2nd Piola Kirchhoff stresses w.r.t. inverse inelastic
     *                          deformation gradient
     */
    void EvaluatedSdiFin(const Core::LinAlg::Matrix<3, 1>& gamma,
        const Core::LinAlg::Matrix<8, 1>& delta, const Core::LinAlg::Matrix<3, 3>& iFinM,
        const Core::LinAlg::Matrix<3, 3>& iCinCM, const Core::LinAlg::Matrix<6, 1>& iCinV,
        const Core::LinAlg::Matrix<9, 1>& CiFin9x1, const Core::LinAlg::Matrix<9, 1>& CiFinCe9x1,
        const Core::LinAlg::Matrix<6, 1>& iCinCiCinV, const Core::LinAlg::Matrix<9, 1>& CiFiniCe9x1,
        const Core::LinAlg::Matrix<6, 1>& iCV, const Core::LinAlg::Matrix<3, 3>& iFinCeM,
        double detFin, Core::LinAlg::Matrix<6, 9>& dSdiFin) const;

    /*!
     * @brief Evaluates stress and cmat
     *
     * \f[
     *   cmat = 2 \frac{\partial \mathbf{S}}{\partial \mathbf{C}}
     * \f]
     *
     *  cf. Holzapfel [1], p. 216 (6.32) and p. 248
     *  The formulation in Holzapfel is slightly changed as only the elastic strains cause stresses
     *  with \f$\mathbf{C}_\text{in} = \mathbf{F}_\text{in}^T \cdot \mathbf{F}_\text{in}\f$
     *  and \f$\mathbf{F}_\text{in}\f$ denoting the inelastic deformation gradient
     * \f[
     *   \mathbf{S} = 2 \det(\mathbf{F}_\text{in}) \mathbf{F}_\text{in}^{-1} \frac{\partial
     *   \Psi_\text{e}}{\partial \mathbf{C}_\text{e}} \mathbf{F}_\text{in}^{-T}
     * \f]
     * \f[
     * \mathbf{S} = \det(\mathbf{F}_\text{in}) \left(\gamma_1 \ \mathbf{C}_\text{in}^{-1} + \gamma_2
     * \
     * \mathbf{C}_\text{in}^{-1} \cdot \mathbf{C} \cdot \mathbf{C}_\text{in}^{-1} + \gamma_3 \
     * \mathbf{C}^{-1} \right)
     * \f] with \f$\gamma_i\f$ as defined in Mat::CalculateGammaDelta()
     *
     * \f[
     *   \mathbb{C} = \det(\mathbf{F}_\text{in}) \left( \delta_1 \left( \mathbf{C}_\text{in}^{-1}
     *                \otimes \mathbf{C}_\text{in}^{-1} \right)
     *              + \delta_2 \left( \mathbf{C}_\text{in}^{-1} \otimes \mathbf{C}_\text{in}^{-1}
     *                \mathbf{C} \mathbf{C}_\text{in}^{-1}  + \mathbf{C}_\text{in}^{-1} \mathbf{C}
     *                \mathbf{C}_\text{in}^{-1} \otimes \mathbf{C}_\text{in}^{-1} \right)
     *              + \delta_3 \left( \mathbf{C}_\text{in}^{-1} \otimes \mathbf{C}^{-1} +
     * \mathbf{C}^{-1} \otimes \mathbf{C}_\text{in}^{-1} \right)
     *              + \delta_4 \left( \mathbf{C}_\text{in}^{-1} \mathbf{C} \mathbf{C}_\text{in}^{-1}
     *                \otimes \mathbf{C}_\text{in}^{-1} \mathbf{C} \mathbf{C}_\text{in}^{-1} \right)
     *              + \delta_5 \left( \mathbf{C}_\text{in}^{-1} \mathbf{C} \mathbf{C}_\text{in}^{-1}
     *                \otimes \mathbf{C}^{-1} + \mathbf{C}^{-1} \otimes \mathbf{C}_\text{in}^{-1}
     * \mathbf{C} \mathbf{C}_\text{in}^{-1} \right)
     *              + \delta_6 \left( \mathbf{C}^{-1} \otimes \mathbf{C}^{-1} \right)
     *              + \delta_7 \left( \mathbf{C}^{-1} \odot \mathbf{C}^{-1} \right)
     *              + \delta_8 \left( \mathbf{C}_\text{in}^{-1} \odot \mathbf{C}_\text{in}^{-1}
     * \right) \right) \f]
     *
     * \f$ \odot \f$ is defined as follows:
     * \f[
     *   - (A^{-1} \odot A^{-1})_{abcd} = 1/2 (A^{-1}_{ac} A^{-1}_{bd} +
     *   A^{-1}_{ad} A^{-1}_{bc}); \text{ for } A_{ab} = A_{ba}
     * \f] meaning a symmetric 2-tensor \f$A\f$
     *
     * @param[in] iCV     Inverse right Cauchy-Green tensor \f$ \mathbf{C}^{-1} \f$ stored as 6x1
     *                    vector
     * @param[in] iCinV   Inverse inelastic right Cauchy-Green tensor
     * \f$\mathbf{C}_\text{in}^{-1}\f$ stored as 6x1 vector
     * @param[in] iCinCiCinV  \f$ \mathbf{C}_\text{in}^{-1} \cdot \mathbf{C} \cdot
     *                        \mathbf{C}_\text{in}^{-1} \f$ stored as 6x1 vector
     * @param[in] gamma  Factors for stress calculation (see e.g. Holzapfel - Nonlinear Solid
     *                   Mechanics)
     * @param[in] delta  Factors for elasticity tensor calculation (see e.g. Holzapfel  - Nonlinear
     *                   Solid Mechanics)
     * @param[in] detFin    determinant of inelastic deformation gradient
     * @param[out] stress   2nd Piola--Kirchhoff stress tensor
     * @param[out] cmatiso  part of the elasticity tensor as shown above
     */
    void evaluate_stress_cmat_iso(const Core::LinAlg::Matrix<6, 1>& iCV,
        const Core::LinAlg::Matrix<6, 1>& iCinV, const Core::LinAlg::Matrix<6, 1>& iCinCiCinV,
        const Core::LinAlg::Matrix<3, 1>& gamma, const Core::LinAlg::Matrix<8, 1>& delta,
        double detFin, Core::LinAlg::Matrix<6, 1>& stress,
        Core::LinAlg::Matrix<6, 6>& cmatiso) const;

    /*!
     * @brief Evaluates some kinematic quantities that are used in stress and elasticity tensor
     * calculation
     *
     * @param[in] defgrad  Deformation gradient
     * @param[in] iFinM    Inverse inelastic deformation gradient
     * @param[out] iCinV   Inverse inelastic right Cauchy-Green tensor
     * \f$\mathbf{C}_\text{in}^{-1}\f$ stored as 6x1 vector
     * @param[out] iCinCiCinV  \f$ \mathbf{C}_\text{in}^{-1} \cdot \mathbf{C} \cdot
     *                         \mathbf{C}_\text{in}^{-1} \f$ stored as 6x1 vector
     * @param[out] iCV         Inverse right Cauchy-Green tensor \f$ \mathbf{C}^{-1} \f$ stored as
     * 6x1 vector
     * @param[out] iCinCM      \f$ \mathbf{C}_\text{in}^{-1} \cdot \mathbf{C} \f$
     * @param[out] iFinCeM     \f$ \mathbf{F}_\text{in}^{-1} \cdot \mathbf{C}_\text{el} \f$
     * @param[out] CiFin9x1    \f$ \mathbf{C} \cdot \mathbf{F}_\text{in}^{-1} \f$ stored as 9x1
     * vector
     * @param[out] CiFinCe9x1  \f$ \mathbf{C} \cdot \mathbf{F}_\text{in}^{-1} \cdot
     *                         \mathbf{C}_\text{el} \f$ stored as 9x1 vector
     * @param[out] CiFiniCe9x1   \f$ \mathbf{C} \cdot \mathbf{F}_\text{in}^{-1} \cdot
     *                           \mathbf{C}_\text{el}^{-1} \f$ stored as 9x1 vector
     * @param[out] prinv         Principal invariants of the elastic right Cauchy-Green tensor
     */
    void evaluate_kin_quant_elast(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinM, Core::LinAlg::Matrix<6, 1>& iCinV,
        Core::LinAlg::Matrix<6, 1>& iCinCiCinV, Core::LinAlg::Matrix<6, 1>& iCV,
        Core::LinAlg::Matrix<3, 3>& iCinCM, Core::LinAlg::Matrix<3, 3>& iFinCeM,
        Core::LinAlg::Matrix<9, 1>& CiFin9x1, Core::LinAlg::Matrix<9, 1>& CiFinCe9x1,
        Core::LinAlg::Matrix<9, 1>& CiFiniCe9x1, Core::LinAlg::Matrix<3, 1>& prinv) const;

    /*!
     * @brief calculates the derivatives of the hyper-elastic laws with respect to the invariants
     *
     * @param[in] prinv   Principal invariants of the elastic right Cauchy-Green tensor
     * @param[in] gp      current gauss point
     * @param[in] eleGID  Element ID
     * @param[out] dPI    First derivative w.r.t. principle invariants
     * @param[out] ddPII  Second derivative w.r.t. principle invariants
     */
    void evaluate_invariant_derivatives(const Core::LinAlg::Matrix<3, 1>& prinv, int gp, int eleGID,
        Core::LinAlg::Matrix<3, 1>& dPI, Core::LinAlg::Matrix<6, 1>& ddPII) const;

    /*!
     * @brief pre-evaluation, intended to be used for stuff that have to be done only once
     *
     * @param[in] params  parameter list as handed in from the element
     * @param[in] gp      current gauss point
     */
    void pre_evaluate(Teuchos::ParameterList& params, int gp) const;

    /*!
     * @brief set the gauss point concentration to the respective parameter class of the inelastic
     * material
     *
     * @param[in]  concentration concentration at gauss point
     */
    void SetConcentrationGP(double concentration);

   private:
    /// Holder for anisotropy
    Teuchos::RCP<Anisotropy> anisotropy_;

    /// Holds and classifies all inelastic factors
    Teuchos::RCP<InelasticFactorsHandler> inelastic_;

    /// My material parameters
    Mat::PAR::MultiplicativeSplitDefgradElastHyper* params_;

    /// map to elastic materials/potential summands
    std::vector<Teuchos::RCP<Mat::Elastic::Summand>> potsumel_;
  };

}  // namespace Mat
FOUR_C_NAMESPACE_CLOSE

#endif
