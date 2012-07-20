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

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::Myocard::Myocard( Teuchos::RCP<MAT::PAR::Material> matdata )
: Parameter(matdata),
  diffusivity(matdata->GetDouble("DIFFUSIVITY")),
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
  : params_(NULL)
{
    v0_ = 1.0;
    w0_ = 1.0;
    s0_ = 0.0;
}

/*----------------------------------------------------------------------*
 |  Constructor                                  (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::Myocard::Myocard(MAT::PAR::Myocard* params)
  : params_(params)
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

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack                                        (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
void MAT::Myocard::Unpack(const std::vector<char>& data)
{
  vector<char>::size_type position = 0;
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


  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

/*----------------------------------------------------------------------*
 |  calculate cmat*gradphi                                   ljag 06/12 |
 *----------------------------------------------------------------------*/

/*void MAT::Myocard::Evaluate(
  const LINALG::Matrix<3,1>& gradphi,
  LINALG::Matrix<3,3>& cmat,
  LINALG::Matrix<3,1>& electricflux
  ) const
{
  // conductivity tensor
  cmat.Clear();
  for (int i=0; i<3; ++i)
      {
      for (int j=0; j<3; ++j)
	  {
	  cmat(i,j) = params_->conduct_; // TODO How to get the parameter from the .dat-file?
	  }
      }

  // electric flux
  electricflux.MultiplyNN(cmat,gradphi);

  // done
  return;
}*/

/*----------------------------------------------------------------------*
 |  calculate reaction coefficient                           ljag 06/12 |
 *----------------------------------------------------------------------*/
// TODO find and implement an appropriate model for the cell ionic currents during depolarisation.
// [1] AL Hodgkin and AF Huxley - A quantatice description of membrane current and its application to conduction and excitation in nerve
// [2] D Noble - A modification of the Hodgkin-Huxley equations applicable to purkinje fibre action and pace-maker potentials
// [3] CH Luo and Y Rudy - A model of the ventricular cardiac action potential. Depolarization, repolarization and their interaction
// [4] AG Kleber and Y Rudy - Basic mechanisms of cardiac impulse propagation and associated arrhythmias
// [5] Fenton and Karma - Vortex dynamics in three-dimensional continuous myocardium with fiber rotation: Filament instability and fibrillation
// [6] Cherry and Fenton - Suppression of alternans and conduction blocks despite steep APD restitution:  electrotonic, memory and conduction velocity restitution effects
// [7] Bueno-Orovio et. al. - Minimal model for human ventricular action potentials in tissue
// [8] Rush and Larsen - A practical algorithm for solving dynamic membrane equations
// [9] ten Tusscher et. al. - A model for human ventricular tissue

double MAT::Myocard::ComputeReactionCoeff(const double phi, const double dt)
{
    // Phenomenological model [5]-[8]

    // calculate voltage dependent time constants ([7] page 545)
     const double Tau_vm = GatingFunction(params_->Tau_v1m, params_->Tau_v2m, phi, params_->Theta_vm);
     const double Tau_wm = params_->Tau_w1m + (params_->Tau_w2m - params_->Tau_w1m)*(1.0 + tanh(params_->k_wm*(phi - params_->u_wm)))/2.0;
     const double Tau_so = params_->Tau_so1 + (params_->Tau_so2 - params_->Tau_so1)*(1.0 + tanh(params_->k_so*(phi - params_->u_so)))/2.0;
     const double Tau_s = GatingFunction(params_->Tau_s1, params_->Tau_s2, phi, params_->Theta_w);
     const double Tau_o = GatingFunction(params_->Tau_o1, params_->Tau_o2, phi, params_->Theta_o);

     // calculate infinity values ([7] page 545)
     const double v_inf = GatingFunction(1.0, 0.0, phi, params_->Theta_vm);
     const double w_inf = GatingFunction(1.0 - phi/params_->Tau_winf, params_->w_infs, phi, params_->Theta_o);

     // calculate gating variables according to [8]
     const double exp_v = -GatingFunction(dt/Tau_vm, dt/params_->Tau_vp, phi, params_->Theta_v);
     const double exp_w = -GatingFunction(dt/Tau_wm, dt/params_->Tau_wp, phi, params_->Theta_w);

     const double v = GatingFunction(v_inf, 0.0, phi, params_->Theta_v) + GatingFunction(v0_ - v_inf, v0_, phi, params_->Theta_v)*exp(exp_v);
     const double w = GatingFunction(w_inf, 0.0, phi, params_->Theta_w) + GatingFunction(w0_ - w_inf, w0_, phi, params_->Theta_w)*exp(exp_w);
     const double s = (1.0 + tanh(params_->k_s*(phi - params_->u_s)))/2.0 + (s0_ - (1.0 + tanh(params_->k_s*(phi - params_->u_s)))/2.0)*exp(-dt/Tau_s);

     // update initial values according to [8]
     v0_ = v;
     w0_ = w;
     s0_ = s;

     // calculate currents J_fi, J_so and J_si ([7] page 545)
     const double J_fi = -GatingFunction(0.0, v*(phi - params_->Theta_v)*(params_->u_u - phi)/params_->Tau_fi, phi, params_->Theta_v); // fast inward current
     const double J_so = GatingFunction((phi - params_->u_o)/Tau_o, 1.0/Tau_so, phi, params_->Theta_w);// slow outward current
     const double J_si = -GatingFunction(0.0, w*s/params_->Tau_si, phi, params_->Theta_w); // slow inward current

     const double reacoeff = (J_fi + J_so + J_si);

     if (reacoeff<1E4) return reacoeff;
     else return 1E4;
}

/*----------------------------------------------------------------------*
 |                                                           ljag 06/12 |
 *----------------------------------------------------------------------*/
double MAT::Myocard::ComputeReactionCoeffDeriv(const double phi, const double dt)
{
    // Phenomenological model [5]-[8]

    // calculate voltage dependent time constants ([7] page 545)
     const double Tau_vm = GatingFunction(params_->Tau_v1m, params_->Tau_v2m, phi, params_->Theta_vm);
     const double Tau_wm = params_->Tau_w1m + (params_->Tau_w2m - params_->Tau_w1m)*(1.0 + tanh(params_->k_wm*(phi - params_->u_wm)))/2.0;
//     const double Tau_so = params_->Tau_so1 + (params_->Tau_so2 - params_->Tau_so1)*(1.0 + tanh(params_->k_so*(phi - params_->u_so)))/2.0;
     const double Tau_s = GatingFunction(params_->Tau_s1, params_->Tau_s2, phi, params_->Theta_w);
     const double Tau_o = GatingFunction(params_->Tau_o1, params_->Tau_o2, phi, params_->Theta_o);

     // calculate infinity values ([7] page 545)
     const double v_inf = GatingFunction(1.0, 0.0, phi, params_->Theta_vm);
     const double w_inf = GatingFunction((1.0 - phi/params_->Tau_winf), params_->w_infs, phi, params_->Theta_o);

     // calculate gating variables according to [8]
     const double exp_v = -dt*GatingFunction(1.0/Tau_vm, 1.0/params_->Tau_vp, phi, params_->Theta_v);
     const double exp_w = -dt*GatingFunction(1.0/Tau_wm, 1.0/params_->Tau_wp, phi, params_->Theta_w);

     const double v = GatingFunction(v_inf, 0.0, phi, params_->Theta_v) + GatingFunction(v0_ - v_inf, v0_, phi, params_->Theta_v)*exp(exp_v);
     const double w = GatingFunction(w_inf, 0.0, phi, params_->Theta_w) + GatingFunction(w0_ - w_inf, w0_, phi, params_->Theta_w)*exp(exp_w);
     const double s = (1.0 + tanh(params_->k_s*(phi - params_->u_s)))/2.0 - (s0_ + tanh(params_->k_s*(phi - params_->u_s)))/2.0*exp(-dt/Tau_s);

    // calculate derivatives of variables
     const double dw_inf = -GatingFunction(1.0/params_->Tau_winf, 0.0, phi, params_->Theta_o);

     const double dw = GatingFunction(dw_inf, 0.0, phi, params_->Theta_w) - GatingFunction(1.0 - dw_inf, dw_inf, phi, params_->Theta_w)*exp(-dt/GatingFunction(Tau_wm, params_->Tau_wp, phi, params_->Theta_w));
     const double ds = (1.0 - pow(tanh(params_->k_s*(phi - params_->u_s)), 2))*(1.0 - exp(dt/Tau_s));

     const double dJ_fi = -v*GatingFunction(0.0, 1.0/params_->Tau_fi, phi, params_->Theta_v)*(params_->Theta_v + params_->u_u - 2*phi);
     const double dJ_so = GatingFunction(1.0/Tau_o, 0.0, phi, params_->Theta_w);
     const double dJ_si = -GatingFunction(0.0, 1.0/params_->Tau_si, phi, params_->Theta_w)*(w*ds + dw*s);

     const double reacoeffderiv = (dJ_fi + dJ_so + dJ_si);

    if (reacoeffderiv<1E4) return reacoeffderiv;
    else return 1E4;
}

/*----------------------------------------------------------------------*
 |                                                           ljag 06/12 |
 *----------------------------------------------------------------------*/
double MAT::Myocard::GatingFunction(const double Gate1, const double Gate2, const double var, const double thresh) const
{
    if (var<thresh) return Gate1;
    else return Gate2;
}


/*----------------------------------------------------------------------*
 |  update of material at the end of a time step             ljag 07/12 |
 *----------------------------------------------------------------------*/
void MAT::Myocard::Update()
{
  // update initial values according to [8]
  //v0_ = v;
  //w0_ = w;
  //s0_ = s;
  return;
}

