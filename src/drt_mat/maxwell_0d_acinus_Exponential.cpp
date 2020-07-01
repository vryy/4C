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


#include <vector>
#include "maxwell_0d_acinus_Exponential.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_red_airways/red_airway.H"
#include "../drt_lib/drt_linedefinition.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Maxwell_0d_acinus_Exponential::Maxwell_0d_acinus_Exponential(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Maxwell_0d_acinus(matdata)
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::Maxwell_0d_acinus_Exponential::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Maxwell_0d_acinus_Exponential(this));
}


MAT::Maxwell_0d_acinusExponentialType MAT::Maxwell_0d_acinusExponentialType::instance_;


DRT::ParObject* MAT::Maxwell_0d_acinusExponentialType::Create(const std::vector<char>& data)
{
  MAT::Maxwell_0d_acinus_Exponential* mxwll_0d_acin = new MAT::Maxwell_0d_acinus_Exponential();
  mxwll_0d_acin->Unpack(data);
  return mxwll_0d_acin;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Maxwell_0d_acinus_Exponential::Maxwell_0d_acinus_Exponential() : Maxwell_0d_acinus() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Maxwell_0d_acinus_Exponential::Maxwell_0d_acinus_Exponential(
    MAT::PAR::Maxwell_0d_acinus* params)
    : Maxwell_0d_acinus(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Maxwell_0d_acinus_Exponential::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
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
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Maxwell_0d_acinus_Exponential::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // Extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  ;
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // Extract e1_0_, e1_lin_, e1_exp_, tau_
  ExtractfromPack(position, data, e1_0_);
  ExtractfromPack(position, data, e1_lin_);
  ExtractfromPack(position, data, e1_exp_);
  ExtractfromPack(position, data, tau_);

  // Extract matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Maxwell_0d_acinus_Exponential*>(mat);
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
void MAT::Maxwell_0d_acinus_Exponential::Setup(DRT::INPUT::LineDefinition* linedef)
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
void MAT::Maxwell_0d_acinus_Exponential::Evaluate(Epetra_SerialDenseVector& epnp,
    Epetra_SerialDenseVector& epn, Epetra_SerialDenseVector& epnm, Epetra_SerialDenseMatrix& sysmat,
    Epetra_SerialDenseVector& rhs, Teuchos::ParameterList& params, const double NumOfAcini,
    const double Vo, double time, double dt)
{
  // Set sysmat and rhs to zero
  sysmat.Scale(0.0);
  rhs.Scale(0.0);

  // Get acinar volume in next and current timestep
  double acin_vnp = params.get<double>("acin_vnp");
  double acin_vn = params.get<double>("acin_vn");

  // Get flow in next and current timestep
  double qnp = params.get<double>("qin_np");
  double qn = params.get<double>("qin_n");

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
