/*----------------------------------------------------------------------*/
/*! \file
\brief
St. Venant-Kirchhoff material

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef MAT_STVENANTKIRCHHOFF_H
#define MAT_STVENANTKIRCHHOFF_H


#include "mat_par_parameter.H"
#include "mat_so3_material.H"
#include "lib_parobjectfactory.H"

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
      /// linear coefficient of thermal expansion
      const double thermexpans_;

      //@}

      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class StVenantKirchhoff
  }     // namespace PAR

  class StVenantKirchhoffType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "StVenantKirchhoffType"; }

    static StVenantKirchhoffType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

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

    int UniqueParObjectId() const override
    {
      return StVenantKirchhoffType::Instance().UniqueParObjectId();
    }

    void Pack(DRT::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    //@}

    //! @name Access methods

    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_stvenant; }

    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::kinem_linear || kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new StVenantKirchhoff(*this));
    }

    /// Young's modulus
    double Youngs() const { return params_->youngs_; }

    /// Poisson's ratio
    double PoissonRatio() const { return params_->poissonratio_; }

    double Density() const override { return params_->density_; }

    /// shear modulus
    double ShearMod() const { return 0.5 * params_->youngs_ / (1.0 + params_->poissonratio_); }

    MAT::PAR::Parameter* Parameter() const override { return params_; }

    //@}

    //! @name Evaluation methods

    /// evaluates material law
    void Evaluate(const Epetra_SerialDenseVector* glstrain_e, Epetra_SerialDenseMatrix* cmat_e,
        Epetra_SerialDenseVector* stress_e);

    void Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat, int gp,
        int eleGID) override;

    void StrainEnergy(
        const CORE::LINALG::Matrix<6, 1>& glstrain, double& psi, int gp, int eleGID) override;

    void EvaluateGEMM(CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* stress,
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>* cmat, double* density,
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_m,
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_new,
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_old,
        CORE::LINALG::Matrix<3, 3>* rcg_new, CORE::LINALG::Matrix<3, 3>* rcg_old, int gp,
        int eleGID) override;

    // computes isotropic elasticity tensor in matrix notion for 3d
    void SetupCmat(CORE::LINALG::Matrix<6, 6>& cmat);

    // computes isotropic elasticity tensor in matrix notion for 2d
    // plane strain, rotational symmetry
    void SetupCmat2d(Epetra_SerialDenseMatrix* cmat);

    //@}

    //! general setup of constitutive tensor based on Young's and poisson's ratio
    static void FillCmat(CORE::LINALG::Matrix<6, 6>& cmat, double Emod, double nu);

   private:
    /// my material parameters
    MAT::PAR::StVenantKirchhoff* params_;
  };
}  // namespace MAT

#endif  // MAT_STVENANTKIRCHHOFF_H
