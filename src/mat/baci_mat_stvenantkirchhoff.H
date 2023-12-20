/*----------------------------------------------------------------------*/
/*! \file
\brief
St. Venant-Kirchhoff material

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_MAT_STVENANTKIRCHHOFF_H
#define BACI_MAT_STVENANTKIRCHHOFF_H


#include "baci_config.H"

#include "baci_comm_parobjectfactory.H"
#include "baci_mat_par_parameter.H"
#include "baci_mat_so3_material.H"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for St. Venant--Kirchhoff
    class StVenantKirchhoff : public Parameter
    {
     public:
      /// standard constructor
      StVenantKirchhoff(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{

      /// Young's modulus
      const double youngs_;
      /// Possion's ratio
      const double poissonratio_;
      /// mass density
      const double density_;

      //@}

      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class StVenantKirchhoff
  }     // namespace PAR

  class StVenantKirchhoffType : public CORE::COMM::ParObjectType
  {
   public:
    [[nodiscard]] std::string Name() const override { return "StVenantKirchhoffType"; }

    static StVenantKirchhoffType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

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
    explicit StVenantKirchhoff(MAT::PAR::StVenantKirchhoff* params);

    [[nodiscard]] int UniqueParObjectId() const override
    {
      return StVenantKirchhoffType::Instance().UniqueParObjectId();
    }

    void Pack(CORE::COMM::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    //@}

    //! @name Access methods

    [[nodiscard]] INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_stvenant;
    }

    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (kinem != INPAR::STR::kinem_linear && kinem != INPAR::STR::kinem_nonlinearTotLag)
        dserror("element and material kinematics are not compatible");
    }

    [[nodiscard]] Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new StVenantKirchhoff(*this));
    }

    /// Young's modulus
    [[nodiscard]] double Youngs() const { return params_->youngs_; }

    /// Poisson's ratio
    [[nodiscard]] double PoissonRatio() const { return params_->poissonratio_; }

    [[nodiscard]] double Density() const override { return params_->density_; }

    /// shear modulus
    [[nodiscard]] double ShearMod() const
    {
      return 0.5 * params_->youngs_ / (1.0 + params_->poissonratio_);
    }

    [[nodiscard]] MAT::PAR::Parameter* Parameter() const override { return params_; }

    //@}

    //! @name Evaluation methods

    /// evaluates material law
    void Evaluate(const CORE::LINALG::SerialDenseVector* glstrain_e,
        CORE::LINALG::SerialDenseMatrix* cmat_e, CORE::LINALG::SerialDenseVector* stress_e);

    void Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat, int gp,
        int eleGID) override;

    void StrainEnergy(
        const CORE::LINALG::Matrix<6, 1>& glstrain, double& psi, int gp, int eleGID) override;

    void EvaluateGEMM(CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* stress,
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>* cmat,
        const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_m,
        const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_new,
        const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_old,
        const CORE::LINALG::Matrix<3, 3>* rcg_new, const CORE::LINALG::Matrix<3, 3>* rcg_old,
        int gp, int eleGID) override;

    // computes isotropic elasticity tensor in matrix notion for 3d
    void SetupCmat(CORE::LINALG::Matrix<6, 6>& cmat);
    //@}

    //! general setup of constitutive tensor based on Young's and poisson's ratio
    static void FillCmat(CORE::LINALG::Matrix<6, 6>& cmat, double Emod, double nu);

   private:
    /// my material parameters
    MAT::PAR::StVenantKirchhoff* params_;
  };
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif  // MAT_STVENANTKIRCHHOFF_H
