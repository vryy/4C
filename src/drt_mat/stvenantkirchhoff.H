/*----------------------------------------------------------------------*/
/*! \file
\brief
St. Venant-Kirchhoff material

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef STVENANTKIRCHHOFF_H_
#define STVENANTKIRCHHOFF_H_


#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Teuchos_ParameterList.hpp>


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

    /*----------------------------------------------------------------------*/
    /// material parameters for St. Venant--Kirchhoff growth
    class StVKGrowth : public StVenantKirchhoff
    {
     public:
      /// standard constructor
      StVKGrowth(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{

      /// reference concentration without inelastic deformation
      const double c0_;

      /// order of polynomial for inelastic growth
      const int poly_num_;

      /// parameters of polynomial for inelastic growth
      const std::vector<double> poly_params_;

      /// growth proportional to amount of substance (true) or porportional to concentration (false)
      const bool amount_prop_growth_;

      //@}

      Teuchos::RCP<MAT::Material> CreateMaterial() override;
    };
  }  // namespace PAR

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

    void Evaluate(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
        int gp, int eleGID) override;

    void StrainEnergy(
        const LINALG::Matrix<6, 1>& glstrain, double& psi, int gp, int eleGID) override;

    void EvaluateGEMM(LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* stress,
        LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>* cmat, double* density,
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_m,
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_new,
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_old, LINALG::Matrix<3, 3>* rcg_new,
        LINALG::Matrix<3, 3>* rcg_old, int gp, int eleGID) override;

    // computes isotropic elasticity tensor in matrix notion for 3d
    void SetupCmat(LINALG::Matrix<6, 6>& cmat);

    // computes isotropic elasticity tensor in matrix notion for 2d
    // plane strain, rotational symmetry
    void SetupCmat2d(Epetra_SerialDenseMatrix* cmat);

    //@}

    //! general setup of constitutive tensor based on Young's and poisson's ratio
    static void FillCmat(LINALG::Matrix<6, 6>& cmat, double Emod, double nu);

   private:
    /// my material parameters
    MAT::PAR::StVenantKirchhoff* params_;
  };


  class StVKGrowthType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "StVKGrowthType"; }

    static StVKGrowthType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static StVKGrowthType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for St.-Venant-Kirchhoff material with growth properties
  class StVKGrowth : public StVenantKirchhoff
  {
   public:
    /// construct empty material object
    StVKGrowth();

    /// construct the material object given material parameters
    explicit StVKGrowth(MAT::PAR::StVKGrowth* params);

    int UniqueParObjectId() const override
    {
      return StVKGrowthType::Instance().UniqueParObjectId();
    }

    void Unpack(const std::vector<char>& data) override;

    void Pack(DRT::PackBuffer& data) const override;

    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_stvenant_growth; }

    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new StVKGrowth(*this)); }

    /// reference concentration without inelastic deformation
    double C0() const { return growth_params_->c0_; }

    /// parameters of polynomial for inelastic growth
    int PolyNum() const { return growth_params_->poly_num_; }

    /// parameters of polynomial for inelastic growth
    std::vector<double> PolyParams() const { return growth_params_->poly_params_; }

    /// growth proportional to amount of substance or to concentration
    bool AmountPropGrowth() const { return growth_params_->amount_prop_growth_; }

    MAT::PAR::Parameter* Parameter() const override { return growth_params_; }

   private:
    /// my material parameters
    MAT::PAR::StVKGrowth* growth_params_;
  };
}  // namespace MAT

#endif /*STVENANTKIRCHHOFF_H_*/
