/*----------------------------------------------------------------------*/
/*! \file
\brief Contains the functions to establish local material law
stress-strain law for isotropic material for a 3D element
following Drucker Prager plasticity model
Reference:
EA de Souza Neto, D Peric, DRJ Owen. Computational Methods of Plasticity: Theory and Applications,
John Wiley & Sons, Ltd, 2008
\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_PLASTICDRUCKERPRAGER_HPP
#define FOUR_C_MAT_PLASTICDRUCKERPRAGER_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_inpar_material.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_utils_local_newton.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace PAR
  {
    using FAD = Sacado::Fad::DFad<double>;
    /**
     * \brief Elasto-plasitc Drucker-Prager Material Model
     *
     * This material model simulates the elasto-plastic behaviour of materials including soil
     * and concrete based on the yield surface criteria and the plastic hardening of the material.
     *
     * Following the approach provided by EA de Souza Neto, D Peric, DRJ Owen.
     * Computational Methods of Plasticity: Theory and Applications, Page 338-339
     */
    class PlasticDruckerPrager : public Core::Mat::PAR::Parameter
    {
     public:
      //! standard constructor
      PlasticDruckerPrager(const Core::Mat::PAR::Parameter::Data& matdata);
      //! @name material parameters
      //@{
      //! Young's modulus
      const double youngs_;
      //! Possion's ratio
      const double poissonratio_;
      //! Density
      const double density_;
      //! linear isotropic hardening modulus
      const double isohard_;
      //! tolerance for local Newton iteration
      const double abstol_;
      //! initial cohesion
      const double cohesion_;
      //! friction angle variable
      const double eta_;
      //! hardening factor
      const double xi_;
      //! dilatancy angle variable
      const double etabar_;
      //! method to compute the material tangent
      const std::string tang_;
      //! maximum iterations for local system solve
      const int itermax_;
      //@}
      //! create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;
    };
  }  // namespace PAR
  class PlasticDruckerPragerType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "PlasticDruckerPragerType"; }
    static PlasticDruckerPragerType& Instance() { return instance_; };
    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static PlasticDruckerPragerType instance_;
  };

  class PlasticDruckerPrager : public So3Material
  {
   public:
    PlasticDruckerPrager();

    explicit PlasticDruckerPrager(Mat::PAR::PlasticDruckerPrager* params);
    int UniqueParObjectId() const override
    {
      return PlasticDruckerPragerType::Instance().UniqueParObjectId();
    }
    void pack(Core::Communication::PackBuffer& data) const override;
    void unpack(const std::vector<char>& data) override;
    Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_pldruckprag;
    }
    void ValidKinematics(Inpar::STR::KinemType kinem) override
    {
      if (kinem != Inpar::STR::KinemType::linear)
        FOUR_C_THROW(
            "The plastic Drucker Prager material model is only compatible with linear kinematics.");
    }
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new PlasticDruckerPrager(*this));
    }
    void setup(int numgp, Input::LineDefinition* linedef) override;
    void update() override;
    /**
     * \brief Evaulate the stresses from the strains in the material
     *
     * \param defgrad :deformation gradient
     * \param linstrain :linear total strains
     * \param params :parameter list for communication
     * \param stress :2nd PK-stress
     * \param cmat :material stiffness matrix
     * \param gp :Gauss point
     * \param eleGID :element global identifier
     */
    void evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>* linstrain, Teuchos::ParameterList& params,
        Core::LinAlg::Matrix<NUM_STRESS_3D, 1>* stress,
        Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat, int gp, int eleGID) override
    {
      this->EvaluateFAD(defgrd, linstrain, params, stress, cmat, gp, eleGID);
    };

    template <typename ScalarT>
    void evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT>* linstrain,
        Teuchos::ParameterList& params, Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT>* stress,
        Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat, int gp, int eleGID)
    {
      this->EvaluateFAD(defgrd, linstrain, params, stress, cmat, gp, eleGID);
    };
    template <typename ScalarT>
    void EvaluateFAD(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT>* linstrain,
        Teuchos::ParameterList& params, Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT>* stress,
        Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat, int gp, int eleGID);
    template <typename T>
    void Stress(const T p, const Core::LinAlg::Matrix<NUM_STRESS_3D, 1, T>& devstress,
        Core::LinAlg::Matrix<NUM_STRESS_3D, 1, T>& stress);

    void setup_cmat(Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat);
    /**
     * \brief setup the elastoplasticity tensor in matrix notation for 3d return to cone
     *
     * \param cmat :elasto-plastic tangent modulus (out)
     * \param Dgamma :plastic multiplier
     * \param G :shear modulus
     * \param Kappa :Bulk modulus
     * \param devstrain :deviatoric strain
     * \param xi :Mohr-Coulomb parameter
     * \param Hiso :isotropic hardening modulus
     * \param eta :Mohr-Coulomb parameter
     * \param etabar :Mohr-Coulomb parameter
     */
    void setup_cmat_elasto_plastic_cone(Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
        double Dgamma, double G, double Kappa, Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& devstrain,
        double xi, double Hiso, double eta, double etabar);
    /**
     * \brief setup the elastoplasticity tensor in matrix notation for 3d return to apex
     *
     * \param cmat :elasto-plastic tangent modulus (out)
     * \param Kappa :Bulk modulus
     * \param devstrain :deviatoric strain
     * \param xi :Mohr-Coulomb parameter
     * \param Hiso :isotropic hardening modulus
     * \param eta :Mohr-Coulomb parameter
     * \param etabar :Mohr-Coulomb parameter
     */
    void setup_cmat_elasto_plastic_apex(Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>&
                                            cmat,           // elasto-plastic tangent modulus (out)
        double Kappa,                                       // Bulk modulus
        Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& devstrain,  // deviatoric strain
        double xi,                                          // Mohr-Coulomb parameter
        double Hiso,                                        // isotropic hardening modulus
        double eta,                                         // Mohr-Coulomb parameter
        double etabar                                       // Mohr-Coulomb parameter
    );
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }
    double Density() const override { return params_->density_; }
    template <typename T>
    std::pair<T, T> return_to_cone_funct_and_deriv(T Dgamma, T G, T kappa, T Phi_trial);
    template <typename T>
    std::pair<T, T> return_to_apex_funct_and_deriv(T dstrainv, T p, T kappa, T strainbar_p);
    bool Initialized() const { return (isinit_ and !strainplcurr_.empty()); }
    void register_output_data_names(
        std::unordered_map<std::string, int>& names_and_size) const override;

    bool EvaluateOutputData(
        const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const override;

   private:
    Mat::PAR::PlasticDruckerPrager* params_;
    std::vector<Core::LinAlg::Matrix<NUM_STRESS_3D, 1>> strainpllast_;
    std::vector<Core::LinAlg::Matrix<NUM_STRESS_3D, 1>> strainplcurr_;
    std::vector<double> strainbarpllast_;
    std::vector<double> strainbarplcurr_;
    bool isinit_;
  };
}  // namespace Mat
FOUR_C_NAMESPACE_CLOSE

#endif
