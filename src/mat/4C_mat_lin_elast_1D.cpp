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
MAT::PAR::LinElast1D::LinElast1D(Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : Parameter(matdata),
      youngs_(matdata->Get<double>("YOUNG")),
      density_(matdata->Get<double>("DENS"))
{
  if (youngs_ <= 0.) FOUR_C_THROW("Young's modulus must be greater zero");
  if (density_ <= 0.) FOUR_C_THROW("Density must be greater zero");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::MAT::Material> MAT::PAR::LinElast1D::create_material()
{
  return Teuchos::rcp(new MAT::LinElast1D(this));
}

MAT::LinElast1DType MAT::LinElast1DType::instance_;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::COMM::ParObject* MAT::LinElast1DType::Create(const std::vector<char>& data)
{
  auto* stvenantk = new MAT::LinElast1D(nullptr);
  stvenantk->Unpack(data);
  return stvenantk;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::LinElast1D::LinElast1D(MAT::PAR::LinElast1D* params) : params_(params) {}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::LinElast1D::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::LinElast1D::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      CORE::MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::LinElast1D*>(mat);
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
MAT::PAR::LinElast1DGrowth::LinElast1DGrowth(Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : LinElast1D(matdata),
      c0_(matdata->Get<double>("C0")),
      poly_num_(matdata->Get<int>("POLY_PARA_NUM")),
      poly_params_(matdata->Get<std::vector<double>>("POLY_PARAMS")),
      amount_prop_growth_(matdata->Get<bool>("AOS_PROP_GROWTH"))
{
  if (c0_ <= 0.0) FOUR_C_THROW("Reference concentration must be greater than zero");
  if (poly_num_ <= 0) FOUR_C_THROW("Polynomial order must be greater than zero");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::MAT::Material> MAT::PAR::LinElast1DGrowth::create_material()
{
  return Teuchos::rcp(new MAT::LinElast1DGrowth(this));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::LinElast1DGrowthType MAT::LinElast1DGrowthType::instance_;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::COMM::ParObject* MAT::LinElast1DGrowthType::Create(const std::vector<char>& data)
{
  auto* stvk_growth = new MAT::LinElast1DGrowth(nullptr);
  stvk_growth->Unpack(data);
  return stvk_growth;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::LinElast1DGrowth::LinElast1DGrowth(MAT::PAR::LinElast1DGrowth* params)
    : LinElast1D(static_cast<MAT::PAR::LinElast1D*>(params)), growth_params_(params)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::LinElast1DGrowth::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  MAT::LinElast1D::Pack(data);

  // matid
  int matid = -1;
  if (growth_params_ != nullptr)
    matid = growth_params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::LinElast1DGrowth::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  MAT::LinElast1D::Unpack(basedata);

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  growth_params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      CORE::MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        growth_params_ = static_cast<MAT::PAR::LinElast1DGrowth*>(mat);
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
double MAT::LinElast1DGrowth::EvaluatePK2(const double def_grad, const double conc) const
{
  const double def_grad_inel = AmountPropGrowth() ? get_growth_factor_ao_s_prop(conc, def_grad)
                                                  : get_growth_factor_conc_prop(conc);

  const double def_grad_el = def_grad / def_grad_inel;
  const double epsilon_el = 0.5 * (def_grad_el * def_grad_el - 1.0);

  return growth_params_->youngs_ * epsilon_el / def_grad_inel;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::LinElast1DGrowth::EvaluateStiffness(const double def_grad, const double conc) const
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
double MAT::LinElast1DGrowth::evaluate_elastic_energy(double def_grad, double conc) const
{
  const double def_grad_inel = AmountPropGrowth() ? get_growth_factor_ao_s_prop(conc, def_grad)
                                                  : get_growth_factor_conc_prop(conc);

  const double def_grad_el = def_grad / def_grad_inel;
  const double epsilon_el = 0.5 * (def_grad_el * def_grad_el - 1.0);

  return 0.5 * (2.0 * growth_params_->youngs_ * epsilon_el / def_grad_inel) * epsilon_el;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::LinElast1DGrowth::get_growth_factor_conc_prop(const double conc) const
{
  return CORE::FE::Polynomial(growth_params_->poly_params_).Evaluate(conc - growth_params_->c0_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::LinElast1DGrowth::get_growth_factor_ao_s_prop(
    const double conc, const double def_grad) const
{
  return CORE::FE::Polynomial(growth_params_->poly_params_)
      .Evaluate(conc * def_grad - growth_params_->c0_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::LinElast1DGrowth::get_growth_factor_ao_s_prop_deriv(
    const double conc, const double def_grad) const
{
  const double first_deriv = CORE::FE::Polynomial(growth_params_->poly_params_)
                                 .EvaluateDerivative(conc * def_grad - growth_params_->c0_, 1);

  return first_deriv * conc;
}

FOUR_C_NAMESPACE_CLOSE
