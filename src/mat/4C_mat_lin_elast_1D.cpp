/*----------------------------------------------------------------------*/
/*! \file
\brief
Linear elastic material in one dimension and material that supports growth due to an external
quantity (e.g. concentration)

\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_mat_lin_elast_1D.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_utils_function_library.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mat::PAR::LinElast1D::LinElast1D(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      youngs_(matdata.parameters.get<double>("YOUNG")),
      density_(matdata.parameters.get<double>("DENS"))
{
  if (youngs_ <= 0.) FOUR_C_THROW("Young's modulus must be greater zero");
  if (density_ <= 0.) FOUR_C_THROW("Density must be greater zero");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::LinElast1D::create_material()
{
  return Teuchos::rcp(new Mat::LinElast1D(this));
}

Mat::LinElast1DType Mat::LinElast1DType::instance_;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Mat::LinElast1DType::Create(const std::vector<char>& data)
{
  auto* stvenantk = new Mat::LinElast1D(nullptr);
  stvenantk->unpack(data);
  return stvenantk;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mat::LinElast1D::LinElast1D(Mat::PAR::LinElast1D* params) : params_(params) {}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::LinElast1D::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::LinElast1D::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::LinElast1D*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mat::PAR::LinElast1DGrowth::LinElast1DGrowth(const Core::Mat::PAR::Parameter::Data& matdata)
    : LinElast1D(matdata),
      c0_(matdata.parameters.get<double>("C0")),
      poly_num_(matdata.parameters.get<int>("POLY_PARA_NUM")),
      poly_params_(matdata.parameters.get<std::vector<double>>("POLY_PARAMS")),
      amount_prop_growth_(matdata.parameters.get<bool>("AOS_PROP_GROWTH"))
{
  if (c0_ <= 0.0) FOUR_C_THROW("Reference concentration must be greater than zero");
  if (poly_num_ <= 0) FOUR_C_THROW("Polynomial order must be greater than zero");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::LinElast1DGrowth::create_material()
{
  return Teuchos::rcp(new Mat::LinElast1DGrowth(this));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mat::LinElast1DGrowthType Mat::LinElast1DGrowthType::instance_;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Mat::LinElast1DGrowthType::Create(const std::vector<char>& data)
{
  auto* stvk_growth = new Mat::LinElast1DGrowth(nullptr);
  stvk_growth->unpack(data);
  return stvk_growth;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mat::LinElast1DGrowth::LinElast1DGrowth(Mat::PAR::LinElast1DGrowth* params)
    : LinElast1D(static_cast<Mat::PAR::LinElast1D*>(params)), growth_params_(params)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::LinElast1DGrowth::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  Mat::LinElast1D::pack(data);

  // matid
  int matid = -1;
  if (growth_params_ != nullptr)
    matid = growth_params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::LinElast1DGrowth::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Mat::LinElast1D::unpack(basedata);

  // matid and recover params_
  int matid;
  extract_from_pack(position, data, matid);
  growth_params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        growth_params_ = static_cast<Mat::PAR::LinElast1DGrowth*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::LinElast1DGrowth::EvaluatePK2(const double def_grad, const double conc) const
{
  const double def_grad_inel = AmountPropGrowth() ? get_growth_factor_ao_s_prop(conc, def_grad)
                                                  : get_growth_factor_conc_prop(conc);

  const double def_grad_el = def_grad / def_grad_inel;
  const double epsilon_el = 0.5 * (def_grad_el * def_grad_el - 1.0);

  return growth_params_->youngs_ * epsilon_el / def_grad_inel;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::LinElast1DGrowth::EvaluateStiffness(const double def_grad, const double conc) const
{
  // F_in
  const double def_grad_inel = AmountPropGrowth() ? get_growth_factor_ao_s_prop(conc, def_grad)
                                                  : get_growth_factor_conc_prop(conc);

  // F_el
  const double def_grad_el = def_grad / def_grad_inel;

  // E_el
  const double epsilon_el = 0.5 * (def_grad_el * def_grad_el - 1.0);

  // dF_in/dF
  const double d_def_grad_inel_d_def_grad =
      AmountPropGrowth() ? get_growth_factor_ao_s_prop_deriv(conc, def_grad) : 0.0;

  // dF_el_dF
  const double d_def_grad_el_d_def_grad =
      (def_grad_inel - def_grad * d_def_grad_inel_d_def_grad) / (def_grad_inel * def_grad_inel);

  // dE_el_dFel
  const double d_epsilon_el_d_def_grad_el = def_grad_el;

  // dE_el_dF
  const double d_epsilon_el_d_def_grad = d_epsilon_el_d_def_grad_el * d_def_grad_el_d_def_grad;

  return growth_params_->youngs_ *
         (d_epsilon_el_d_def_grad * def_grad_inel - epsilon_el * d_def_grad_inel_d_def_grad) /
         (def_grad_inel * def_grad_inel);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::LinElast1DGrowth::evaluate_elastic_energy(double def_grad, double conc) const
{
  const double def_grad_inel = AmountPropGrowth() ? get_growth_factor_ao_s_prop(conc, def_grad)
                                                  : get_growth_factor_conc_prop(conc);

  const double def_grad_el = def_grad / def_grad_inel;
  const double epsilon_el = 0.5 * (def_grad_el * def_grad_el - 1.0);

  return 0.5 * (2.0 * growth_params_->youngs_ * epsilon_el / def_grad_inel) * epsilon_el;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::LinElast1DGrowth::get_growth_factor_conc_prop(const double conc) const
{
  return Core::FE::Polynomial(growth_params_->poly_params_).evaluate(conc - growth_params_->c0_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::LinElast1DGrowth::get_growth_factor_ao_s_prop(
    const double conc, const double def_grad) const
{
  return Core::FE::Polynomial(growth_params_->poly_params_)
      .evaluate(conc * def_grad - growth_params_->c0_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::LinElast1DGrowth::get_growth_factor_ao_s_prop_deriv(
    const double conc, const double def_grad) const
{
  const double first_deriv = Core::FE::Polynomial(growth_params_->poly_params_)
                                 .EvaluateDerivative(conc * def_grad - growth_params_->c0_, 1);

  return first_deriv * conc;
}

FOUR_C_NAMESPACE_CLOSE
