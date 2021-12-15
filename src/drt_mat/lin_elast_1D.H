/*----------------------------------------------------------------------*/
/*! \file
\brief
Linear elastic material in one dimension and material that supports growth due to an external
quantity (e.g. concentration)

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef LINELAST1D_H_
#define LINELAST1D_H_


#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    class LinElast1D : public Parameter
    {
     public:
      LinElast1D(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{
      /// Young's modulus
      const double youngs_;

      /// mass density
      const double density_;
      //@}

      Teuchos::RCP<MAT::Material> CreateMaterial() override;
    };

    class LinElast1DGrowth : public LinElast1D
    {
     public:
      LinElast1DGrowth(Teuchos::RCP<MAT::PAR::Material> matdata);

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

  class LinElast1DType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "LinElast1DType"; }

    static LinElast1DType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static LinElast1DType instance_;
  };

  class LinElast1D : public Material
  {
   public:
    explicit LinElast1D(MAT::PAR::LinElast1D* params);

    int UniqueParObjectId() const override
    {
      return LinElast1DType::Instance().UniqueParObjectId();
    }

    void Pack(DRT::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_linelast1D; }

    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new LinElast1D(*this)); }

    /// Young's modulus
    double Youngs() const { return params_->youngs_; }

    /// mass density
    double Density() const override { return params_->density_; }

    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    MAT::PAR::LinElast1D* params_;
  };


  class LinElast1DGrowthType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "LinElast1DGrowthType"; }

    static LinElast1DGrowthType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static LinElast1DGrowthType instance_;
  };

  class LinElast1DGrowth : public LinElast1D
  {
   public:
    explicit LinElast1DGrowth(MAT::PAR::LinElast1DGrowth* params);

    int UniqueParObjectId() const override
    {
      return LinElast1DGrowthType::Instance().UniqueParObjectId();
    }

    void Unpack(const std::vector<char>& data) override;

    void Pack(DRT::PackBuffer& data) const override;

    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_linelast1D_growth;
    }

    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new LinElast1DGrowth(*this));
    }

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
    MAT::PAR::LinElast1DGrowth* growth_params_;
  };
}  // namespace MAT

#endif /*LINELAST1D_H_*/
