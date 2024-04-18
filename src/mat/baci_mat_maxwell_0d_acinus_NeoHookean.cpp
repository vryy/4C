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


#include "baci_mat_maxwell_0d_acinus_NeoHookean.hpp"

#include "baci_global_data.hpp"
#include "baci_io_linedefinition.hpp"
#include "baci_mat_par_bundle.hpp"
#include "baci_red_airways_elem_params.hpp"
#include "baci_red_airways_elementbase.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Maxwell0dAcinusNeoHookean::Maxwell0dAcinusNeoHookean(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Maxwell0dAcinus(matdata)
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::Maxwell0dAcinusNeoHookean::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Maxwell0dAcinusNeoHookean(this));
}


MAT::Maxwell0dAcinusNeoHookeanType MAT::Maxwell0dAcinusNeoHookeanType::instance_;


CORE::COMM::ParObject* MAT::Maxwell0dAcinusNeoHookeanType::Create(const std::vector<char>& data)
{
  MAT::Maxwell0dAcinusNeoHookean* mxwll_0d_acin = new MAT::Maxwell0dAcinusNeoHookean();
  mxwll_0d_acin->Unpack(data);
  return mxwll_0d_acin;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Maxwell0dAcinusNeoHookean::Maxwell0dAcinusNeoHookean() : Maxwell0dAcinus() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Maxwell0dAcinusNeoHookean::Maxwell0dAcinusNeoHookean(MAT::PAR::Maxwell0dAcinus* params)
    : Maxwell0dAcinus(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Maxwell0dAcinusNeoHookean::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
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
void MAT::Maxwell0dAcinusNeoHookean::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // Extract matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Maxwell0dAcinusNeoHookean*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*
 | Setup routine for NeoHookean material                                |
 |                                                          roth 10/2014|
 *----------------------------------------------------------------------*/
void MAT::Maxwell0dAcinusNeoHookean::Setup(INPUT::LineDefinition* linedef)
{
  // do nothing, all parameters are read by base class already
}


/*----------------------------------------------------------------------*
 | Evaluate NeoHookean material and build system matrix and rhs         |
 |                                                          roth 10/2014|
 *----------------------------------------------------------------------*/
void MAT::Maxwell0dAcinusNeoHookean::Evaluate(CORE::LINALG::SerialDenseVector& epnp,
    CORE::LINALG::SerialDenseVector& epn, CORE::LINALG::SerialDenseVector& epnm,
    CORE::LINALG::SerialDenseMatrix& sysmat, CORE::LINALG::SerialDenseVector& rhs,
    const DRT::REDAIRWAYS::ElemParams& params, const double NumOfAcini, const double Vo,
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
    dserror("Acinus condition at node (%d) has zero acini");
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
