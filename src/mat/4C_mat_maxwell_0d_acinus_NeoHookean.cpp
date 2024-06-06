/*----------------------------------------------------------------------*/
/*! \file

\brief Four-element Maxwell material model for reduced dimensional acinus elements simplified to
linear spring (Stiffness1) only (Neo Hookean), inherits from Maxwell_0d_acinus

Necessary input lines:
(material section)
MAT 3 MAT_0D_MAXWELL_ACINUS_OGDEN Stiffness1 1.0 Stiffness2 5249.1 Viscosity1 3221.86 Viscosity2
1000.0 // acinus properties;
(element section)
1 RED_ACINUS LINE2 2 3 MAT 3 TYPE NeoHookean AcinusVolume 300 AlveolarDuctVolume 0.03711


\level 3
*/
/*----------------------------------------------------------------------*/


#include "4C_mat_maxwell_0d_acinus_NeoHookean.hpp"

#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_red_airways_elem_params.hpp"
#include "4C_red_airways_elementbase.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::Maxwell0dAcinusNeoHookean::Maxwell0dAcinusNeoHookean(
    Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : Maxwell0dAcinus(matdata)
{
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::Maxwell0dAcinusNeoHookean::create_material()
{
  return Teuchos::rcp(new Mat::Maxwell0dAcinusNeoHookean(this));
}


Mat::Maxwell0dAcinusNeoHookeanType Mat::Maxwell0dAcinusNeoHookeanType::instance_;


Core::Communication::ParObject* Mat::Maxwell0dAcinusNeoHookeanType::Create(
    const std::vector<char>& data)
{
  Mat::Maxwell0dAcinusNeoHookean* mxwll_0d_acin = new Mat::Maxwell0dAcinusNeoHookean();
  mxwll_0d_acin->Unpack(data);
  return mxwll_0d_acin;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::Maxwell0dAcinusNeoHookean::Maxwell0dAcinusNeoHookean() : Maxwell0dAcinus() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::Maxwell0dAcinusNeoHookean::Maxwell0dAcinusNeoHookean(Mat::PAR::Maxwell0dAcinus* params)
    : Maxwell0dAcinus(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::Maxwell0dAcinusNeoHookean::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // Pack type of this instance of ParObject
  int type = UniqueParObjectId();

  AddtoPack(data, type);

  // Pack matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::Maxwell0dAcinusNeoHookean::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // Extract matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::Maxwell0dAcinusNeoHookean*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*
 | Setup routine for NeoHookean material                                |
 |                                                          roth 10/2014|
 *----------------------------------------------------------------------*/
void Mat::Maxwell0dAcinusNeoHookean::Setup(Input::LineDefinition* linedef)
{
  // do nothing, all parameters are read by base class already
}


/*----------------------------------------------------------------------*
 | Evaluate NeoHookean material and build system matrix and rhs         |
 |                                                          roth 10/2014|
 *----------------------------------------------------------------------*/
void Mat::Maxwell0dAcinusNeoHookean::Evaluate(Core::LinAlg::SerialDenseVector& epnp,
    Core::LinAlg::SerialDenseVector& epn, Core::LinAlg::SerialDenseVector& epnm,
    Core::LinAlg::SerialDenseMatrix& sysmat, Core::LinAlg::SerialDenseVector& rhs,
    const Discret::ReducedLung::ElemParams& params, const double NumOfAcini, const double Vo,
    double time, double dt)
{
  // Set sysmat and rhs to zero
  sysmat.putScalar(0.0);
  rhs.putScalar(0.0);

  // Get acini pressure and beginning and end of acinus element
  double p1n = epn(0);
  double p2n = epn(1);

  // Safety check for NumOfAcini
  if (NumOfAcini < 1.0)
  {
    FOUR_C_THROW("Acinus condition at node (%d) has zero acini");
  }

  // Linear branches of the Maxwell model (Stiffness2(), B=R_t, B_a=R_a), notation according to
  // interacinar dependency paper
  const double Kp_np = 1.0 / (Stiffness1() * dt);
  const double Kp_n = 1.0 / (Stiffness1() * dt);

  // Build the system matrix for \boldsymbol{K} * \boldsymbol{P} = \boldsymbol{Q}
  sysmat(0, 0) = -1.0 * (Kp_np)*NumOfAcini;
  sysmat(0, 1) = 1.0 * (Kp_np)*NumOfAcini;
  sysmat(1, 0) = 1.0 * (Kp_np)*NumOfAcini;
  sysmat(1, 1) = -1.0 * (Kp_np)*NumOfAcini;

  // Build the corresponding right hand side
  rhs(0) = -1.0 * (Kp_n * (p1n - p2n)) * NumOfAcini;
  rhs(1) = 1.0 * (Kp_n * (p1n - p2n)) * NumOfAcini;
}

FOUR_C_NAMESPACE_CLOSE
