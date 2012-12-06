/*!----------------------------------------------------------------------
\file myocard.cpp

<pre>
Maintainer: Lasse Jagschies
            lasse.jagschies@tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                              dano 09/09 |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  headers                                                  dano 09/09 |
 *----------------------------------------------------------------------*/

#include <vector>
#include "myocard.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::Myocard::Myocard( Teuchos::RCP<MAT::PAR::Material> matdata )
: Parameter(matdata),
  maindirdiffusivity(matdata->GetDouble("MAIN_DIFFUSIVITY")),
  offdirdiffusivity(matdata->GetDouble("OFF_DIFFUSIVITY")),
  dt_deriv(matdata->GetDouble("PERTUBATION_DERIV")),
  u_o(matdata->GetDouble("U_O")), //
  u_u(matdata->GetDouble("U_U")),
  Theta_v(matdata->GetDouble("THETA_V")),
  Theta_w(matdata->GetDouble("THETA_W")),
  Theta_vm(matdata->GetDouble("THETA_VM")),
  Theta_o(matdata->GetDouble("THETA_O")),
  Tau_v1m(matdata->GetDouble("TAU_V1M")),
  Tau_v2m(matdata->GetDouble("TAU_V2M")),
  Tau_vp(matdata->GetDouble("TAU_VP")),
  Tau_w1m(matdata->GetDouble("TAU_W1M")),
  Tau_w2m(matdata->GetDouble("TAU_W2M")),
  k_wm(matdata->GetDouble("K_WM")),
  u_wm(matdata->GetDouble("U_WM")),
  Tau_wp(matdata->GetDouble("TAU_WP")),
  Tau_fi(matdata->GetDouble("TAU_FI")),
  Tau_o1(matdata->GetDouble("TAU_O1")),
  Tau_o2(matdata->GetDouble("TAU_O2")),
  Tau_so1(matdata->GetDouble("TAU_SO1")),
  Tau_so2(matdata->GetDouble("TAU_SO2")),
  k_so(matdata->GetDouble("K_SO")),
  u_so(matdata->GetDouble("U_SO")),
  Tau_s1(matdata->GetDouble("TAU_S1")),
  Tau_s2(matdata->GetDouble("TAU_S2")),
  k_s(matdata->GetDouble("K_S")),
  u_s(matdata->GetDouble("U_S")),
  Tau_si(matdata->GetDouble("TAU_SI")),
  Tau_winf(matdata->GetDouble("TAU_WINF")),
  w_infs(matdata->GetDouble("W_INFS"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::Myocard::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Myocard(this));
}


MAT::MyocardType MAT::MyocardType::instance_;


DRT::ParObject* MAT::MyocardType::Create( const std::vector<char> & data )
{
  MAT::Myocard* myocard = new MAT::Myocard();
  myocard->Unpack(data);
  return myocard;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::Myocard::Myocard()
  : params_(NULL),
    difftensor_(true)
{
    v0_ = 1.0;
    w0_ = 1.0;
    s0_ = 0.0;
}

/*----------------------------------------------------------------------*
 |  Constructor                                  (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::Myocard::Myocard(MAT::PAR::Myocard* params)
  : params_(params),
    difftensor_(true)
{
    v0_ = 1.0;
    w0_ = 1.0;
    s0_ = 0.0;
}

/*----------------------------------------------------------------------*
 |  Pack                                          (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
void MAT::Myocard::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);

  // pack history data
  AddtoPack(data,v0_);
  AddtoPack(data,w0_);
  AddtoPack(data,s0_);
  AddtoPack(data,difftensor_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack                                        (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
void MAT::Myocard::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Myocard*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  ExtractfromPack(position,data,v0_);
  ExtractfromPack(position,data,w0_);
  ExtractfromPack(position,data,s0_);
  ExtractfromPack(position,data,difftensor_);


  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

/*----------------------------------------------------------------------*
 |  Setup conductivity tensor                                ljag 06/12 |
 *----------------------------------------------------------------------*/

void MAT::Myocard::Setup(DRT::INPUT::LineDefinition* linedef)
    {
    // get conductivity values of main fibre direction and
    const double maindirdiffusivity = params_->maindirdiffusivity;
    const double offdirdiffusivity  = params_->maindirdiffusivity;

    // read local eigenvectors of diffusion tensor at current element
    std::vector<double> rad;
    std::vector<double> axi;
    std::vector<double> cir;
    linedef->ExtractDoubleVector("RAD",rad);
    linedef->ExtractDoubleVector("AXI",axi);
    linedef->ExtractDoubleVector("CIR",cir);

    // eigenvector matrix
    double radnorm=0.; double axinorm=0.; double cirnorm=0.;
    LINALG::Matrix<3,3> evmat(true);
    LINALG::Matrix<3,3> evmatinv(true);
    for (int i = 0; i < 3; ++i) {
	radnorm += rad[i]*rad[i]; axinorm += axi[i]*axi[i]; cirnorm += cir[i]*cir[i];
    }
    radnorm = sqrt(radnorm); axinorm = sqrt(axinorm); cirnorm = sqrt(cirnorm);
    for (int i=0; i<3; ++i){
	evmat(i,0) = rad[i]/radnorm;
	evmat(i,1) = axi[i]/axinorm;
	evmat(i,2) = cir[i]/cirnorm;
    }

    // determinant of eigenvector matrix
    double const evmatdet = evmat(0,0)*evmat(1,1)*evmat(2,2)
				    + evmat(0,1)*evmat(1,2)*evmat(2,0)
				    + evmat(0,2)*evmat(1,0)*evmat(2,1)
				    - evmat(0,2)*evmat(1,1)*evmat(2,0)
				    - evmat(0,1)*evmat(1,0)*evmat(2,2)
				    - evmat(0,0)*evmat(1,2)*evmat(2,1);
    // double const evmatdet = evmat.Determinant();

    // inverse of eigenvector matrix
    evmatinv(0,0) = evmat(1,1)*evmat(2,2)-evmat(1,2)*evmat(2,1);
    evmatinv(0,1) = evmat(0,2)*evmat(2,1)-evmat(0,1)*evmat(2,2);
    evmatinv(0,2) = evmat(0,1)*evmat(1,2)-evmat(0,2)*evmat(1,1);
    evmatinv(1,0) = evmat(1,2)*evmat(2,0)-evmat(1,0)*evmat(2,2);
    evmatinv(1,1) = evmat(0,0)*evmat(2,2)-evmat(0,2)*evmat(2,0);
    evmatinv(1,2) = evmat(0,2)*evmat(1,0)-evmat(0,0)*evmat(1,2);
    evmatinv(2,0) = evmat(1,0)*evmat(2,1)-evmat(1,1)*evmat(2,0);
    evmatinv(2,1) = evmat(0,1)*evmat(2,0)-evmat(0,0)*evmat(2,1);
    evmatinv(2,2) = evmat(0,0)*evmat(1,1)-evmat(0,1)*evmat(1,0);
    evmatinv.Scale(1/evmatdet);
    //const LINALG::Matrix<3,3> ematinv;
    //evmat.Invert(ematinv);
    //cout << "EVMAT: " << evmat << "   EVMATinv: " << evmatinv << endl;
    // Conductivity matrix D = EVmat*DiagonalConductivityMatrix*EVmatinv
    for (int i = 0; i<3; i++){
	evmatinv(0,i) *= maindirdiffusivity;
	evmatinv(1,i) *= offdirdiffusivity;
	evmatinv(2,i) *= offdirdiffusivity;
    }

    difftensor_.Multiply(evmat, evmatinv);

    // done
  return;
}

void MAT::Myocard::ComputeDiffusivity(LINALG::Matrix<1,1>& diffus3) const
{
    diffus3(0,0) = difftensor_(0,0); return;
}
void MAT::Myocard::ComputeDiffusivity(LINALG::Matrix<2,2>& diffus3) const
{
    for (int i=0; i<2; i++){for (int j=0; j<2; j++){diffus3(i,j) = difftensor_(i,j);}}
    return;
}
void MAT::Myocard::ComputeDiffusivity(LINALG::Matrix<3,3>& diffus3) const
{
    for (int i=0; i<3; i++){for (int j=0; j<3; j++){diffus3(i,j) = difftensor_(i,j);}}
    return;
}


double MAT::Myocard::ComputeReactionCoeff(const double phi, const double dt) const
{
    const double p = 1000;
    // Phenomenological model [5]-[8]
    // calculate voltage dependent time constants ([7] page 545)
     const double Tau_vm = GatingFunction(params_->Tau_v1m, params_->Tau_v2m, p, phi, params_->Theta_vm);
     const double Tau_wm = GatingFunction(params_->Tau_w1m, params_->Tau_w2m, params_->k_wm, phi, params_->u_wm);
     const double Tau_so = GatingFunction(params_->Tau_so1, params_->Tau_so2, params_->k_so, phi, params_->u_so);
     const double Tau_s = GatingFunction(params_->Tau_s1, params_->Tau_s2, p, phi, params_->Theta_w);
     const double Tau_o = GatingFunction(params_->Tau_o1, params_->Tau_o2, p, phi, params_->Theta_o);

     // calculate infinity values ([7] page 545)
     const double v_inf = GatingFunction(1.0, 0.0, p, phi, params_->Theta_vm);
     const double w_inf = GatingFunction(1.0 - phi/params_->Tau_winf, params_->w_infs, p, phi, params_->Theta_o);

     // calculate gating variables according to [8]
     double rhs = GatingFunction(v_inf/Tau_vm, 0.0, p, phi, params_->Theta_v);
     double A = -GatingFunction(1.0/Tau_vm, 1/params_->Tau_vp, p, phi, params_->Theta_v);
     const double v = 1/(1/dt-A)*(1/dt*v0_ + rhs);

     rhs = GatingFunction(w_inf/Tau_wm, 0.0, p, phi, params_->Theta_w);
     A = -GatingFunction(1.0/Tau_wm, 1/params_->Tau_wp, p, phi, params_->Theta_w);
     const double w = 1/(1/dt-A)*(1/dt*w0_ + rhs);

     rhs = GatingFunction(0.0, 1.0/Tau_s, params_->k_s, phi, params_->u_s);
     A = -1.0/Tau_s;
     const double s = 1/(1/dt-A)*(1/dt*s0_ + rhs);

     // calculate currents J_fi, J_so and J_si ([7] page 545)
     const double J_fi = -GatingFunction(0.0, v*(phi - params_->Theta_v)*(params_->u_u - phi)/params_->Tau_fi, p, phi, params_->Theta_v); // fast inward current
     const double J_so = GatingFunction((phi - params_->u_o)/Tau_o, 1.0/Tau_so, p, phi, params_->Theta_w);// slow outward current
     const double J_si = -GatingFunction(0.0, w*s/params_->Tau_si, p, phi, params_->Theta_w); // slow inward current

     const double reacoeff = (J_fi + J_so + J_si);

     return reacoeff;
}

/*----------------------------------------------------------------------*
 |                                                           ljag 06/12 |
 *----------------------------------------------------------------------*/
double MAT::Myocard::ComputeReactionCoeffDeriv(const double phi, const double dt) const
{
    const double ReaCoeff_t1 = ComputeReactionCoeff(phi, dt);
    const double ReaCoeff_t2 = ComputeReactionCoeff((phi+params_->dt_deriv), dt);

    const double ReaCoeffDeriv = (ReaCoeff_t2 - ReaCoeff_t1)/(params_->dt_deriv);

    return ReaCoeffDeriv;
}

/*----------------------------------------------------------------------*
 |                                                           ljag 06/12 |
 *----------------------------------------------------------------------*/
double MAT::Myocard::GatingFunction(const double Gate1, const double Gate2, const double p, const double var, const double thresh) const
{
    double Erg = Gate1+(Gate2-Gate1)*(1.0+tanh(p*(var-thresh)))/2;
    return Erg;
}


/*----------------------------------------------------------------------*
 |  update of material at the end of a time step             ljag 07/12 |
 *----------------------------------------------------------------------*/
void MAT::Myocard::Update(const double phi, const double dt)
{

    const double p = 1000;
    // calculate voltage dependent time constants ([7] page 545)
    const double Tau_vm = GatingFunction(params_->Tau_v1m, params_->Tau_v2m, p, phi, params_->Theta_vm);
    const double Tau_wm = GatingFunction(params_->Tau_w1m, params_->Tau_w2m, params_->k_wm, phi, params_->u_wm);
    //const double Tau_so = params_->Tau_so1 + (params_->Tau_so2 - params_->Tau_so1)*(1.0 + tanh(params_->k_so*(phi - params_->u_so)))/2.0;
    const double Tau_s = GatingFunction(params_->Tau_s1, params_->Tau_s2, p, phi, params_->Theta_w);
    //const double Tau_o = GatingFunction(params_->Tau_o1, params_->Tau_o2, phi, params_->Theta_o);

    // calculate infinity values ([7] page 545)
    const double v_inf = GatingFunction(1.0, 0.0, p, phi, params_->Theta_vm);
    const double w_inf = GatingFunction(1.0 - phi/params_->Tau_winf, params_->w_infs, p, phi, params_->Theta_o);

    // calculate gating variables according to [8]
    double rhs = v_inf/Tau_vm*GatingFunction(1.0, 0.0, p, phi, params_->Theta_v);
    double A = -GatingFunction(1.0/Tau_vm, 1/params_->Tau_vp, p, phi, params_->Theta_v);
    const double v = 1/(1/dt-A)*(1/dt*v0_ + rhs);

    rhs = w_inf/Tau_wm*GatingFunction(1.0, 0.0, p, phi, params_->Theta_w);
    A = -GatingFunction(1.0/Tau_wm, 1/params_->Tau_wp, p, phi, params_->Theta_w);
    const double w = 1/(1/dt-A)*(1/dt*w0_ + rhs);

    rhs = GatingFunction(0.0, 1.0/Tau_s, params_->k_s, phi, params_->u_s);
    A = -1.0/Tau_s;
    const double s = 1/(1/dt-A)*(1/dt*s0_ + rhs);

    // update initial values according to [8]
    v0_ = v;
    w0_ = w;
    s0_ = s;

    return;
}

