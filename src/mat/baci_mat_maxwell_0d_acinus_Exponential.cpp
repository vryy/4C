/*----------------------------------------------------------------------*/
/*! \file

\brief Four-element Maxwell material model for reduced dimensional acinus elements with non-linear
spring with exponential behaviour, inherits from Maxwell_0d_acinus

The originally linear spring (Stiffness1) of the 4-element Maxwell model is substituted by an
exponential pressure-volume relation (derivation: see Ismail Mahmoud's dissertation, chapter 3.4)

Necessary input lines:
(material section)
MAT 3 MAT_0D_MAXWELL_ACINUS_OGDEN Stiffness1 1.0 Stiffness2 5249.1 Viscosity1 3221.86 Viscosity2
1000.0 // acinus properties;
(element section)
1 RED_ACINUS  MAT 3 TYPE Exponential AcinusVolume 300 AlveolarDuctVolume 0.03711 E1_0 0 E1_EXP 0
E1_LIN 0 TAU 0


\level 3
*/
/*----------------------------------------------------------------------*/


#include "baci_mat_maxwell_0d_acinus_Exponential.hpp"

#include "baci_global_data.hpp"
#include "baci_io_linedefinition.hpp"
#include "baci_mat_par_bundle.hpp"
#include "baci_red_airways_elem_params.hpp"
#include "baci_red_airways_elementbase.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Maxwell0dAcinusExponential::Maxwell0dAcinusExponential(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Maxwell0dAcinus(matdata)
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::Maxwell0dAcinusExponential::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Maxwell0dAcinusExponential(this));
}


MAT::Maxwell0dAcinusExponentialType MAT::Maxwell0dAcinusExponentialType::instance_;


CORE::COMM::ParObject* MAT::Maxwell0dAcinusExponentialType::Create(const std::vector<char>& data)
{
  MAT::Maxwell0dAcinusExponential* mxwll_0d_acin = new MAT::Maxwell0dAcinusExponential();
  mxwll_0d_acin->Unpack(data);
  return mxwll_0d_acin;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Maxwell0dAcinusExponential::Maxwell0dAcinusExponential() : Maxwell0dAcinus() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Maxwell0dAcinusExponential::Maxwell0dAcinusExponential(MAT::PAR::Maxwell0dAcinus* params)
    : Maxwell0dAcinus(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Maxwell0dAcinusExponential::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // Pack type of this instance of ParObject
  int type = UniqueParObjectId();

  AddtoPack(data, type);
  AddtoPack(data, e1_0_);
  AddtoPack(data, e1_lin_);
  AddtoPack(data, e1_exp_);
  AddtoPack(data, tau_);

  // Pack matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Maxwell0dAcinusExponential::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // Extract e1_0_, e1_lin_, e1_exp_, tau_
  ExtractfromPack(position, data, e1_0_);
  ExtractfromPack(position, data, e1_lin_);
  ExtractfromPack(position, data, e1_exp_);
  ExtractfromPack(position, data, tau_);

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
        params_ = static_cast<MAT::PAR::Maxwell0dAcinusExponential*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*
 | Setup routine to add Exponential material specific parameters E1_0   |
 | E1_LIN, E1_EXP, TAU to material                          roth 10/2014|
 *----------------------------------------------------------------------*/
void MAT::Maxwell0dAcinusExponential::Setup(INPUT::LineDefinition* linedef)
{
  linedef->ExtractDouble("E1_0", e1_0_);
  linedef->ExtractDouble("E1_LIN", e1_lin_);
  linedef->ExtractDouble("E1_EXP", e1_exp_);
  linedef->ExtractDouble("TAU", tau_);
  // TODO bool -variable init, in Evaluate abfragen ob init=true
}


/*----------------------------------------------------------------------*
 | Evaluate Exponential material and build system matrix and rhs.       |
 |                                                          roth 10/2014|
 *----------------------------------------------------------------------*/
void MAT::Maxwell0dAcinusExponential::Evaluate(CORE::LINALG::SerialDenseVector& epnp,
    CORE::LINALG::SerialDenseVector& epn, CORE::LINALG::SerialDenseVector& epnm,
    CORE::LINALG::SerialDenseMatrix& sysmat, CORE::LINALG::SerialDenseVector& rhs,
    const DRT::REDAIRWAYS::ElemParams& params, const double NumOfAcini, const double Vo,
    double time, double dt)
{
  // Set sysmat and rhs to zero
  sysmat.putScalar(0.0);
  rhs.putScalar(0.0);

  // Get acinar volume in next and current timestep
  double acin_vnp = params.acin_vnp;
  double acin_vn = params.acin_vn;

  // Get flow in next and current timestep
  double qnp = params.qin_np;
  double qn = params.qin_n;

  // Get acini pressure and beginning and end of acinus element
  double p1n = epn(0);
  double p2n = epn(1);

  // Safety check for NumOfAcini
  if (NumOfAcini < 1.0)
  {
    dserror("Acinus condition at node (%d) has zero acini");
  }

  // Calculate volume and flow difference per acinuar duct
  double dvnp = (acin_vnp / NumOfAcini) - Vo;
  double dvn = (acin_vn / NumOfAcini) - Vo;

  //------------------------------------------------------------
  // V  = A + B*exp(-K*P)
  //
  // The P-V curve is fitted to create the following
  // P1 = E1.(V-Vo)
  //
  // E1 = e1_0_ + e1_lin_.(V-Vo) + e1_exp_.exp(tau_.(V-Vo))
  //------------------------------------------------------------
  double kp_np = Viscosity1() / (Stiffness2() * dt) + 1;
  double kp_n = -Viscosity1() / (Stiffness2() * dt);
  double kq_np = Viscosity1() * Viscosity2() / (Stiffness2() * dt) + (Viscosity2() + Viscosity1());
  double kq_n = -Viscosity1() * Viscosity2() / (Stiffness2() * dt);

  // Get the terms assosciated with the nonlinear behavior of E1
  double term_nonlin = 0.0;
  double pnpi = 0.0;
  double pnpi2 = 0.0;
  double dpnpi_dt = 0.0;
  double dpnpi2_dt = 0.0;

  // Components of linearized E1
  pnpi = (e1_0_ + e1_lin_ * dvnp + e1_exp_ * exp(tau_ * dvnp)) * dvnp;
  pnpi2 = (e1_0_ + 2 * e1_lin_ * dvnp + e1_exp_ * exp(tau_ * dvnp) * (tau_ * dvnp + 1));

  // Components of linearized d(E1)/dt
  dpnpi_dt = (e1_0_ + 2 * e1_lin_ * dvnp + e1_exp_ * exp(tau_ * dvnp) * (1 + tau_ * dvnp)) *
             (dvnp - dvn) / dt;
  dpnpi2_dt = (2 * e1_lin_ + tau_ * e1_exp_ * exp(tau_ * dvnp) * (1 + tau_ * dvnp) +
                  e1_exp_ * tau_ * exp(tau_ * dvnp)) *
                  (dvnp - dvn) / dt +
              (e1_0_ + 2 * e1_lin_ * dvnp + e1_exp_ * exp(tau_ * dvnp) * (1 + tau_ * dvnp)) / dt;

  term_nonlin = pnpi + pnpi2 * (-(dvnp) + (qn / NumOfAcini) * dt / 2 + dvn);
  kq_np = kq_np + pnpi2 / 2 * dt;
  term_nonlin =
      term_nonlin + dpnpi_dt * Viscosity1() / Stiffness2() +
      dpnpi2_dt * Viscosity1() / Stiffness2() * (-(dvnp) + (qnp / NumOfAcini) * dt / 2 + dvn);
  kq_np = kq_np + dpnpi2_dt * Viscosity1() / Stiffness2() / 2 * dt;

  // Build the system matrix for \boldsymbol{K} * \boldsymbol{P} = \boldsymbol{Q}
  sysmat(0, 0) = -1.0 * (kp_np / kq_np) * NumOfAcini;
  sysmat(0, 1) = 1.0 * (kp_np / kq_np) * NumOfAcini;
  sysmat(1, 0) = 1.0 * (kp_np / kq_np) * NumOfAcini;
  sysmat(1, 1) = -1.0 * (kp_np / kq_np) * NumOfAcini;

  // Build the corresponding right hand side
  rhs(0) = -1.0 * (-(kp_n * (p1n - p2n) - term_nonlin) * NumOfAcini / kq_np + (kq_n * qn) / kq_np);
  rhs(1) = 1.0 * (-(kp_n * (p1n - p2n) - term_nonlin) * NumOfAcini / kq_np + (kq_n * qn) / kq_np);
}

FOUR_C_NAMESPACE_CLOSE
