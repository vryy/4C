/*----------------------------------------------------------------------*/
/*! \file
\brief
Linear elastic material in one dimension and material that supports growth due to an external
quantity (e.g. concentration)

\level 2

*/
/*----------------------------------------------------------------------*/

#include "lin_elast_1D.H"

#include "../drt_lib/drt_function_library.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::LinElast1D::LinElast1D(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata), youngs_(matdata->GetDouble("YOUNG")), density_(matdata->GetDouble("DENS"))
{
  if (youngs_ <= 0.) dserror("Young's modulus must be greater zero");
  if (density_ <= 0.) dserror("Density must be greater zero");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::LinElast1D::CreateMaterial()
{
  return Teuchos::rcp(new MAT::LinElast1D(this));
}

MAT::LinElast1DType MAT::LinElast1DType::instance_;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::LinElast1DType::Create(const std::vector<char>& data)
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
void MAT::LinElast1D::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
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
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::LinElast1D*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::LinElast1DGrowth::LinElast1DGrowth(Teuchos::RCP<MAT::PAR::Material> matdata)
    : LinElast1D(matdata),
      c0_(matdata->GetDouble("C0")),
      poly_num_(matdata->GetInt("POLY_PARA_NUM")),
      poly_params_(*matdata->Get<std::vector<double>>("POLY_PARAMS")),
      amount_prop_growth_(static_cast<bool>(matdata->GetInt("AOS_PROP_GROWTH")))
{
  if (c0_ <= 0.0) dserror("Reference concentration must be greater than zero");
  if (poly_num_ <= 0) dserror("Polynomial order must be greater than zero");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::LinElast1DGrowth::CreateMaterial()
{
  return Teuchos::rcp(new MAT::LinElast1DGrowth(this));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::LinElast1DGrowthType MAT::LinElast1DGrowthType::instance_;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::LinElast1DGrowthType::Create(const std::vector<char>& data)
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
void MAT::LinElast1DGrowth::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
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
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  MAT::LinElast1D::Unpack(basedata);

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  growth_params_ = nullptr;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        growth_params_ = static_cast<MAT::PAR::LinElast1DGrowth*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::LinElast1DGrowth::EvaluatePK2(const double def_grad, const double conc) const
{
  const double def_grad_inel =
      AmountPropGrowth() ? GetGrowthFactorAoSProp(conc, def_grad) : GetGrowthFactorConcProp(conc);

  const double def_grad_el = def_grad / def_grad_inel;
  const double epsilon_el = 0.5 * (def_grad_el * def_grad_el - 1.0);

  return growth_params_->youngs_ * epsilon_el / def_grad_inel;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::LinElast1DGrowth::EvaluateStiffness(const double def_grad, const double conc) const
{
  // F_in
  const double def_grad_inel =
      AmountPropGrowth() ? GetGrowthFactorAoSProp(conc, def_grad) : GetGrowthFactorConcProp(conc);

  // F_el
  const double def_grad_el = def_grad / def_grad_inel;

  // E_el
  const double epsilon_el = 0.5 * (def_grad_el * def_grad_el - 1.0);

  // dF_in/dF
  const double d_def_grad_inel_d_def_grad =
      AmountPropGrowth() ? GetGrowthFactorAoSPropDeriv(conc, def_grad) : 0.0;

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
double MAT::LinElast1DGrowth::EvaluateElasticEnergy(double def_grad, double conc) const
{
  const double def_grad_inel =
      AmountPropGrowth() ? GetGrowthFactorAoSProp(conc, def_grad) : GetGrowthFactorConcProp(conc);

  const double def_grad_el = def_grad / def_grad_inel;
  const double epsilon_el = 0.5 * (def_grad_el * def_grad_el - 1.0);

  return 0.5 * (2.0 * growth_params_->youngs_ * epsilon_el / def_grad_inel) * epsilon_el;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::LinElast1DGrowth::GetGrowthFactorConcProp(const double conc) const
{
  return DRT::UTILS::Polynomial(growth_params_->poly_params_).Evaluate(conc - growth_params_->c0_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::LinElast1DGrowth::GetGrowthFactorAoSProp(const double conc, const double def_grad) const
{
  return DRT::UTILS::Polynomial(growth_params_->poly_params_)
      .Evaluate(conc * def_grad - growth_params_->c0_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::LinElast1DGrowth::GetGrowthFactorAoSPropDeriv(
    const double conc, const double def_grad) const
{
  LINALG::Matrix<2, 1> derivs(false);
  DRT::UTILS::Polynomial(growth_params_->poly_params_)
      .Evaluate(conc * def_grad - growth_params_->c0_, derivs);

  return derivs(1) * conc;
}
