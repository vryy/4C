/*----------------------------------------------------------------------*/
/*! \file
\brief
Linear elastic material in one dimension and material that supports growth due to an external
quantity (e.g. concentration)

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef MAT_LIN_ELAST_1D_H
#define MAT_LIN_ELAST_1D_H


#include "mat_par_parameter.H"
#include "mat_so3_material.H"
#include "lib_parobjectfactory.H"

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

    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new LinElast1D(*this)); }

    /// mass density
    double Density() const override { return params_->density_; }

    /// elastic energy based on @p epsilon
    double EvaluateElasticEnergy(const double epsilon) const
    {
      return 0.5 * EvaluatePK2(epsilon) * epsilon;
    }

    /// evaluate 2nd Piola-Kirchhoff stress based on @param epsilon (Green-Lagrange strain)
    double EvaluatePK2(const double epsilon) const { return params_->youngs_ * epsilon; }

    /// evaluate stiffness of material i.e. derivative of 2nd Piola Kirchhoff stress w.r.t.
    /// Green-Lagrange strain
    double EvaluateStiffness() const { return params_->youngs_; }

    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_linelast1D; }

    void Pack(DRT::PackBuffer& data) const override;

    MAT::PAR::Parameter* Parameter() const override { return params_; }

    int UniqueParObjectId() const override
    {
      return LinElast1DType::Instance().UniqueParObjectId();
    }

    void Unpack(const std::vector<char>& data) override;

   private:
    /// my material parameters
    MAT::PAR::LinElast1D* params_;
  };


  class LinElast1DGrowthType : public DRT::ParObjectType
  {
   public:
    DRT::ParObject* Create(const std::vector<char>& data) override;

    static LinElast1DGrowthType& Instance() { return instance_; }

    std::string Name() const override { return "LinElast1DGrowthType"; }

   private:
    static LinElast1DGrowthType instance_;
  };

  class LinElast1DGrowth : public LinElast1D
  {
   public:
    explicit LinElast1DGrowth(MAT::PAR::LinElast1DGrowth* params);

    /// growth proportional to amount of substance or to concentration
    bool AmountPropGrowth() const { return growth_params_->amount_prop_growth_; }

    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new LinElast1DGrowth(*this));
    }
    /// elastic energy based on @p def_grad and @p conc
    double EvaluateElasticEnergy(double def_grad, double conc) const;

    /// 2nd Piola-Kirchhoff stress based on @p def_grad and @p conc
    double EvaluatePK2(double def_grad, double conc) const;

    /// stiffness, i.e. derivative of 2nd Piola-Kirchhoff stress w.r.t. @p def_grad based on @p
    /// def_grad and @p conc
    double EvaluateStiffness(double def_grad, double conc) const;

    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_linelast1D_growth;
    }

    void Pack(DRT::PackBuffer& data) const override;

    MAT::PAR::Parameter* Parameter() const override { return growth_params_; }

    int UniqueParObjectId() const override
    {
      return LinElast1DGrowthType::Instance().UniqueParObjectId();
    }

    void Unpack(const std::vector<char>& data) override;

   private:
    /// polynomial growth factor based on amount of substance (@p conc * @p def_grad)
    double GetGrowthFactorAoSProp(double conc, double def_grad) const;

    /// derivative of polynomial growth factor based on amount of substance w.r.t @p def_grad
    double GetGrowthFactorAoSPropDeriv(double conc, double def_grad) const;

    /// polynomial growth factor based on concentration (@p conc)
    double GetGrowthFactorConcProp(double conc) const;

    /// my material parameters
    MAT::PAR::LinElast1DGrowth* growth_params_;
  };
}  // namespace MAT

#endif  // MAT_LIN_ELAST_1D_H
