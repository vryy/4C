/*----------------------------------------------------------------------*/
/*! \file
\brief
Linear elastic material in one dimension and material that supports growth due to an external
quantity (e.g. concentration)

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_LIN_ELAST_1D_HPP
#define FOUR_C_MAT_LIN_ELAST_1D_HPP


#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace PAR
  {
    class LinElast1D : public Core::Mat::PAR::Parameter
    {
     public:
      LinElast1D(Teuchos::RCP<Core::Mat::PAR::Material> matdata);

      /// @name material parameters
      //@{
      /// Young's modulus
      const double youngs_;

      /// mass density
      const double density_;
      //@}

      Teuchos::RCP<Core::Mat::Material> create_material() override;
    };

    class LinElast1DGrowth : public LinElast1D
    {
     public:
      LinElast1DGrowth(Teuchos::RCP<Core::Mat::PAR::Material> matdata);

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

      Teuchos::RCP<Core::Mat::Material> create_material() override;
    };
  }  // namespace PAR

  class LinElast1DType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "LinElast1DType"; }

    static LinElast1DType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static LinElast1DType instance_;
  };

  class LinElast1D : public Core::Mat::Material
  {
   public:
    explicit LinElast1D(Mat::PAR::LinElast1D* params);

    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new LinElast1D(*this));
    }

    /// mass density
    double Density() const override { return params_->density_; }

    /// elastic energy based on @p epsilon
    double evaluate_elastic_energy(const double epsilon) const
    {
      return 0.5 * EvaluatePK2(epsilon) * epsilon;
    }

    /// evaluate 2nd Piola-Kirchhoff stress based on @param epsilon (Green-Lagrange strain)
    double EvaluatePK2(const double epsilon) const { return params_->youngs_ * epsilon; }

    /// evaluate stiffness of material i.e. derivative of 2nd Piola Kirchhoff stress w.r.t.
    /// Green-Lagrange strain
    double EvaluateStiffness() const { return params_->youngs_; }

    Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_linelast1D;
    }

    void Pack(Core::Communication::PackBuffer& data) const override;

    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

    int UniqueParObjectId() const override
    {
      return LinElast1DType::Instance().UniqueParObjectId();
    }

    void Unpack(const std::vector<char>& data) override;

   private:
    /// my material parameters
    Mat::PAR::LinElast1D* params_;
  };


  class LinElast1DGrowthType : public Core::Communication::ParObjectType
  {
   public:
    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

    static LinElast1DGrowthType& Instance() { return instance_; }

    std::string Name() const override { return "LinElast1DGrowthType"; }

   private:
    static LinElast1DGrowthType instance_;
  };

  class LinElast1DGrowth : public LinElast1D
  {
   public:
    explicit LinElast1DGrowth(Mat::PAR::LinElast1DGrowth* params);

    /// growth proportional to amount of substance or to concentration
    bool AmountPropGrowth() const { return growth_params_->amount_prop_growth_; }

    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new LinElast1DGrowth(*this));
    }
    /// elastic energy based on @p def_grad and @p conc
    double evaluate_elastic_energy(double def_grad, double conc) const;

    /// 2nd Piola-Kirchhoff stress based on @p def_grad and @p conc
    double EvaluatePK2(double def_grad, double conc) const;

    /// stiffness, i.e. derivative of 2nd Piola-Kirchhoff stress w.r.t. @p def_grad based on @p
    /// def_grad and @p conc
    double EvaluateStiffness(double def_grad, double conc) const;

    Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_linelast1D_growth;
    }

    void Pack(Core::Communication::PackBuffer& data) const override;

    Core::Mat::PAR::Parameter* Parameter() const override { return growth_params_; }

    int UniqueParObjectId() const override
    {
      return LinElast1DGrowthType::Instance().UniqueParObjectId();
    }

    void Unpack(const std::vector<char>& data) override;

   private:
    /// polynomial growth factor based on amount of substance (@p conc * @p def_grad)
    double get_growth_factor_ao_s_prop(double conc, double def_grad) const;

    /// derivative of polynomial growth factor based on amount of substance w.r.t @p def_grad
    double get_growth_factor_ao_s_prop_deriv(double conc, double def_grad) const;

    /// polynomial growth factor based on concentration (@p conc)
    double get_growth_factor_conc_prop(double conc) const;

    /// my material parameters
    Mat::PAR::LinElast1DGrowth* growth_params_;
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
