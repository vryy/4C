/*----------------------------------------------------------------------*/
/*! \file
\brief
St. Venant-Kirchhoff material

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_STVENANTKIRCHHOFF_HPP
#define FOUR_C_MAT_STVENANTKIRCHHOFF_HPP


#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for St. Venant--Kirchhoff
    class StVenantKirchhoff : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      StVenantKirchhoff(const Core::Mat::PAR::Parameter::Data& matdata);

      /// @name material parameters
      //@{

      /// Young's modulus
      const double youngs_;
      /// Possion's ratio
      const double poissonratio_;
      /// mass density
      const double density_;

      //@}

      Teuchos::RCP<Core::Mat::Material> create_material() override;

    };  // class StVenantKirchhoff
  }     // namespace PAR

  class StVenantKirchhoffType : public Core::Communication::ParObjectType
  {
   public:
    [[nodiscard]] std::string name() const override { return "StVenantKirchhoffType"; }

    static StVenantKirchhoffType& instance() { return instance_; };

    Core::Communication::ParObject* create(const std::vector<char>& data) override;

   private:
    static StVenantKirchhoffType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for St.-Venant-Kirchhoff material
  class StVenantKirchhoff : public So3Material
  {
   public:
    /// construct empty material object
    StVenantKirchhoff();

    /// construct the material object given material parameters
    explicit StVenantKirchhoff(Mat::PAR::StVenantKirchhoff* params);

    [[nodiscard]] int unique_par_object_id() const override
    {
      return StVenantKirchhoffType::instance().unique_par_object_id();
    }

    void pack(Core::Communication::PackBuffer& data) const override;

    void unpack(const std::vector<char>& data) override;

    //@}

    //! @name Access methods

    [[nodiscard]] Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::m_stvenant;
    }

    void valid_kinematics(Inpar::Solid::KinemType kinem) override
    {
      if (kinem != Inpar::Solid::KinemType::linear &&
          kinem != Inpar::Solid::KinemType::nonlinearTotLag)
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    [[nodiscard]] Teuchos::RCP<Core::Mat::Material> clone() const override
    {
      return Teuchos::rcp(new StVenantKirchhoff(*this));
    }

    /// Young's modulus
    [[nodiscard]] double youngs() const { return params_->youngs_; }

    /// Poisson's ratio
    [[nodiscard]] double poisson_ratio() const { return params_->poissonratio_; }

    [[nodiscard]] double density() const override { return params_->density_; }

    /// shear modulus
    [[nodiscard]] double shear_mod() const
    {
      return 0.5 * params_->youngs_ / (1.0 + params_->poissonratio_);
    }

    [[nodiscard]] Core::Mat::PAR::Parameter* parameter() const override { return params_; }

    //@}

    //! @name Evaluation methods

    /// evaluates material law
    void evaluate(const Core::LinAlg::SerialDenseVector* glstrain_e,
        Core::LinAlg::SerialDenseMatrix* cmat_e, Core::LinAlg::SerialDenseVector* stress_e);

    void evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmat, int gp,
        int eleGID) override;

    void strain_energy(
        const Core::LinAlg::Matrix<6, 1>& glstrain, double& psi, int gp, int eleGID) override;

    // computes isotropic elasticity tensor in matrix notion for 3d
    void setup_cmat(Core::LinAlg::Matrix<6, 6>& cmat);
    //@}

    //! general setup of constitutive tensor based on Young's and poisson's ratio
    static void fill_cmat(Core::LinAlg::Matrix<6, 6>& cmat, double Emod, double nu);

   private:
    /// my material parameters
    Mat::PAR::StVenantKirchhoff* params_;
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
