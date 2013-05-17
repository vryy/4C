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
 |  definitions                                              ljag 07/12 |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  headers                                                  ljag 07/12 |
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
  model(matdata->Get<string>("MODEL")),
  tissue(matdata->Get<string>("TISSUE"))
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
 |  Constructor                                    (public)  ljag 07/12 |
 *----------------------------------------------------------------------*/
MAT::Myocard::Myocard()
  : params_(NULL),
  difftensor_(true)
  {
        v0_ = 1.0;
        w0_ = 1.0;
        s0_ = 0.0;

        Na_i_ = 11.6;
        Ca_i_ = 0.000045;
        K_i_ = 138.3;
        m_ = 0;
        h_ = 0.75;
        j_ = 0.75;
        d_ = 0;
        f_= 1;
        f_Ca_ = 1;
        s_ = 1;
        r_ = 0;
        x_s_ = 0;
        x_r1_ = 0;
        x_r2_ = 1;
        Ca_sr_ = 0.2;
        g_ = 1;
  }

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)   ljag 07/12 |
 *----------------------------------------------------------------------*/
MAT::Myocard::Myocard(MAT::PAR::Myocard* params)
  : params_(params),
  difftensor_(true)
  {

        v0_ = 1.0;
        w0_ = 1.0;
        s0_ = 0.0;

        Na_i_ = 11.6;
        Ca_i_ = 0.000045;
        K_i_ = 138.3;
        m_ = 0;
        h_ = 0.75;
        j_ = 0.75;
        d_ = 0;
        f_= 1;
        f_Ca_ = 1;
        s_ = 1;
        r_ = 0;
        x_s_ = 0;
        x_r1_ = 0;
        x_r2_ = 1;
        Ca_sr_ = 0.2;
        g_ = 1;
  }

/*----------------------------------------------------------------------*
 |  Pack                                           (public)  ljag 07/12 |
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
  AddtoPack(data, v0_);
  AddtoPack(data, w0_);
  AddtoPack(data, s0_);
  AddtoPack(data, m_);
  AddtoPack(data, h_);
  AddtoPack(data, j_);
  AddtoPack(data, r_);
  AddtoPack(data, s_);
  AddtoPack(data, x_r1_);
  AddtoPack(data, x_r2_);
  AddtoPack(data, x_s_);
  AddtoPack(data, d_);
  AddtoPack(data, f_);
  AddtoPack(data, f_Ca_);
  AddtoPack(data, g_);
  AddtoPack(data, Na_i_);
  AddtoPack(data, K_i_);
  AddtoPack(data, Ca_i_);
  AddtoPack(data, Ca_sr_);
  AddtoPack(data, difftensor_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack                                         (public)  ljag 07/12 |
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

    ExtractfromPack(position, data, v0_);
    ExtractfromPack(position, data, w0_);
    ExtractfromPack(position, data, s0_);
    ExtractfromPack(position, data, m_);
    ExtractfromPack(position, data, h_);
    ExtractfromPack(position, data, j_);
    ExtractfromPack(position, data, r_);
    ExtractfromPack(position, data, s_);
    ExtractfromPack(position, data, x_r1_);
    ExtractfromPack(position, data, x_r2_);
    ExtractfromPack(position, data, x_s_);
    ExtractfromPack(position, data, d_);
    ExtractfromPack(position, data, f_);
    ExtractfromPack(position, data, f_Ca_);
    ExtractfromPack(position, data, g_);
    ExtractfromPack(position, data, Na_i_);
    ExtractfromPack(position, data, K_i_);
    ExtractfromPack(position, data, Ca_i_);
    ExtractfromPack(position, data, Ca_sr_);
    ExtractfromPack(position, data, difftensor_);

    if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
  }


/*----------------------------------------------------------------------*
 |  Setup conductivity tensor                                cbert 02/13 |
 *----------------------------------------------------------------------*/
void MAT::Myocard::Setup(const std::vector<double> &fiber1)
{
  SetupDiffusionTensor(fiber1);
}

void MAT::Myocard::Setup(DRT::INPUT::LineDefinition* linedef)
{
  std::vector<double> fiber1(3);
  linedef->ExtractDoubleVector("FIBER1",fiber1);
  SetupDiffusionTensor(fiber1);
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
  double reacoeff = 0.0;
  const double p = 1000.0;


  if (*(params_->model) == "MV")
  {
    // Model parameter
    const double u_o = 0.0;
    const double u_u = 1.55;//1.58;
    const double Theta_v = 0.3;
    const double Theta_w = 0.13;//0.015;
    const double Theta_vm = 0.006;//0.015;
    const double Theta_o = 0.006;
    const double Tau_v1m = 60.0;
    const double Tau_v2m = 1150.0;
    const double Tau_vp = 1.4506;
    const double Tau_w1m = 60.0;//70.0;
    const double Tau_w2m = 15.0;//20.0;
    const double k_wm = 65.0;
    const double u_wm = 0.03;
    const double Tau_wp = 200.0;//280.0;
    const double Tau_fi = 0.11;
    const double Tau_o1 = 400.0;//6.0;
    const double Tau_o2 = 6.0;
    const double Tau_so1 = 30.0181;//43.0;
    const double Tau_so2 = 0.9957;//0.2;
    const double k_so = 2.0458;//2.0;
    const double u_so = 0.65;
    const double Tau_s1 = 2.7342;
    const double Tau_s2 = 16.0;//3.0;
    const double k_s = 2.0994;
    const double u_s = 0.9087;
    const double Tau_si = 1.8875;//2.8723;
    const double Tau_winf = 0.07;
    const double w_infs = 0.94;

    // calculate voltage dependent time constants ([7] page 545)
    const double Tau_vm = GatingFunction(Tau_v1m, Tau_v2m, p   , phi, Theta_vm);
    const double Tau_wm = GatingFunction(Tau_w1m, Tau_w2m, k_wm, phi, u_wm    );
    const double Tau_so = GatingFunction(Tau_so1, Tau_so2, k_so, phi, u_so    );
    const double Tau_s  = GatingFunction(Tau_s1 , Tau_s2 , p   , phi, Theta_w );
    const double Tau_o  = GatingFunction(Tau_o1 , Tau_o2 , p   , phi, Theta_o );

    // calculate infinity values ([7] page 545)
    const double v_inf = GatingFunction(1.0, 0.0, p, phi, Theta_vm);
    const double w_inf = GatingFunction(1.0 - phi/Tau_winf, w_infs, p, phi, Theta_o);

    // calculate gating variables according to [8]
    const double Tau_v    = GatingFunction(Tau_vm, Tau_vp, p, phi, Theta_v);
    const double v_inf_GF = GatingFunction(v_inf , 0.0   , p, phi, Theta_v);
    const double v = GatingVarCalc(dt, v0_, v_inf_GF, Tau_v);

    const double Tau_w    = GatingFunction(Tau_wm, Tau_wp, p, phi, Theta_w);
    const double w_inf_GF = GatingFunction(w_inf , 0.0   , p, phi, Theta_w);
    const double w = GatingVarCalc(dt, w0_, w_inf_GF, Tau_w);

    const double s_inf = GatingFunction(0.0, 1.0, k_s, phi, u_s);
    const double s = GatingVarCalc(dt, s0_, s_inf, Tau_s);



    // calculate currents J_fi, J_so and J_si ([7] page 545)
    const double J_fi = -GatingFunction(0.0, v*(phi - Theta_v)*(u_u - phi)/Tau_fi, p, phi, Theta_v); // fast inward current
    const double J_so =  GatingFunction((phi - u_o)/Tau_o, 1.0/Tau_so, p, phi, Theta_w);// slow outward current
    const double J_si = -GatingFunction(0.0, w*s/Tau_si, p, phi, Theta_w); // slow inward current

    reacoeff = (J_fi + J_so + J_si);
  }
  else if (*(params_->model) == "TNNP")
  {
    // Model parameter
    const double R      = 8.3143; //
    const double T      = 310.0; //
    const double F      = 96.4867; //
    const double K_O    = 5.4; //
    const double Na_O   = 140.0; //
    const double Ca_O   = 2.0; //
    const double G_Na   = 14.838; //
    const double G_K1   = 5.405; //
    double G_to = 0.294;
    if (*(params_->tissue) == "Endo") { G_to = 0.073; }
    const double G_Kr = 0.096; //
    double G_Ks = 0.245;
    if (*(params_->tissue) == "M") { G_Ks = 0.062; } //
    const double p_KNa  = 0.03; //
    const double G_CaL  = 0.000175; //
    const double k_NaCa = 1000; //
    const double Gamma  = 0.35; //
    const double K_mCa  = 1.38; //
    const double K_mNai = 87.5; //
    const double k_sat  = 0.1; //
    const double Alpha  = 2.5; //
    const double P_NaK  = 1.362; //
    const double K_mK   = 1.0; //
    const double K_mNa  = 40.0; //
    const double G_pK   = 0.0146; //
    const double G_pCa  = 0.025; //
    const double K_pCa  = 0.0005; //
    const double G_bNa  = 0.00029; //
    const double G_bCa  = 0.000592; //

    // Calculate reverse potentials
    const double E_Ks = R*T/(F)*log(K_O + p_KNa*Na_O/(K_i_ + p_KNa*Na_i_));
    const double E_K  = R*T/(F)*log(K_O/K_i_);
    const double E_Na = R*T/(F)*log(Na_O/Na_i_) ;
    const double E_Ca = R*T/(2*F)*log(Ca_O/Ca_i_) ;

    // Calculate ion channel currents
    // ------------------------------

    // fast Natrium channel
    const double a_m   = 1.0/(1.0 + exp( (-60.0 - phi)/5.0 ));
    const double b_m   = 0.1/(1.0 + exp( ( 35.0 + phi)/5.0 )) + 0.1/(1.0 + exp( (-50.0 + phi)/200.0 ));
    const double a_h   = GatingFunction( 0.057*exp( (-80.0 - phi)/6.8 ), 0.0, p, phi, -40.0);
    const double b_h   = GatingFunction( 2.7*exp( 0.079*phi ) + 3.1E5*exp( 0.3485*phi ), 0.77/(0.13*(1.0+exp(-(10.66+phi)/11.1))), p, phi, -40.0);
    const double a_j   = GatingFunction( ((-2.5428E4)*exp( 0.2444*phi )-(6.948E-6)*exp( -0.0439*phi ))*(phi + 37.78)/(1.0 + exp( 0.311*(79.23 + phi) )), 0.0, p, phi, -40.0);
    const double b_j   = GatingFunction( 0.02424*exp( -0.01052*phi )/(1.0 + exp( -0.1378*(40.14+phi) )), 0.6*exp( 0.057*phi )/(1.0 + exp( -0.1*(32.0 + phi) )), p, phi, -40.0);
    const double m_inf = 1.0/pow(1.0 + exp((-56.86 - phi)/9.03),2);
    const double h_inf = 1.0/pow(1.0 + exp((71.55 + phi)/7.43),2);
    const double j_inf = 1.0/pow(1.0 + exp((71.55 + phi)/7.43),2);
    const double tau_m = a_m*b_m;
    const double tau_h = 1/(a_h+b_h);
    const double tau_j = 1/(a_j+b_j);
    const double m     = GatingVarCalc( dt, m_, m_inf, tau_m);
    const double h     = GatingVarCalc( dt, h_, h_inf, tau_h);
    const double j     = GatingVarCalc( dt, j_, j_inf, tau_j);
    const double I_Na  = G_Na*pow( m, 3 )*h*j*(phi - E_Na);

    // Inward rectifier K+ current
    const double a_K1     = 0.1/(1.0+exp( 0.06*(-E_K - 200.0 + phi)/20.0 ));
    const double b_K1     = (3.0*exp( 0.0002*(-E_K + 100.0 + phi)) + exp( 0.1*(-E_K - 10.0 + phi) ))/(1.0 + exp( -0.5*(-E_K + phi) ));
    const double x_K1_inf = a_K1/(a_K1 + b_K1);
    const double I_K1     = G_K1*sqrt(K_O/5.4)*x_K1_inf*(phi - E_K);

    // Transient outward current
    const double r_inf = 1.0/(1.0 + exp((20.0 - phi)/6.0));

    double s_inf = 0;
    double tau_s = 0;
    if (*(params_->tissue) == "Endo")
      {
      s_inf = 1.0/(1.0+exp( (phi + 28.0)/5.0) );
      tau_s = 1000.0*exp( -pow( phi + 67.0, 2 )/1000.0 ) + 8.0;
      }
    else
      {
      s_inf = 1.0/(1.0+exp((phi+20.0)/5.0));
      tau_s = 85.0*exp( -pow( phi + 45.0, 2 )/320.0 ) + 5.0/(1.0 + exp( (phi - 20.0)/5.0) ) + 3.0;
      }
    const double tau_r = 9.5*exp( -pow( 40.0 + phi, 2 )/1800.0 ) + 0.8;
    const double r     = GatingVarCalc( dt, r_, r_inf, tau_r);
    const double s     = GatingVarCalc( dt, s_, s_inf, tau_s);
    const double I_to  = G_to*r*s*(phi - E_K);

    // Rapid delayed rectifier current
    const double a_xr1    = 450.0/(1.0 + exp( (-45.0 - phi)/10.0) );
    const double b_xr1    =   6.0/(1.0 + exp( ( 30.0 + phi)/11.5) );
    const double a_xr2    =   3.0/(1.0 + exp( (-60.0 - phi)/20.0) );
    const double b_xr2    =  1.12/(1.0 + exp( (-60.0 + phi)/20.0) );
    const double x_r1_inf = 450.0/(1.0 + exp( (-26.0 - phi)/ 7.0) );
    const double x_r2_inf =   1.0/(1.0 + exp( ( 88.0 + phi)/24.0) );
    const double tau_xr1  = a_xr1*b_xr1;
    const double tau_xr2  = a_xr2*b_xr2;
    const double x_r1     = GatingVarCalc(dt, x_r1_, x_r1_inf, tau_xr1);
    const double x_r2     = GatingVarCalc(dt, x_r2_, x_r2_inf, tau_xr2);
    const double I_Kr     = G_Kr*sqrt( K_O/5.4 )*x_r1*x_r2*(phi - E_K);

    // Slow delayed rectifier current
    const double a_xs    = 1100.0/(1.0 + exp( (-10.0 - phi)/ 6.0) );
    const double b_xs    =    1.0/(1.0 + exp( (-60.0 + phi)/20.0) );
    const double x_s_inf =    1.0/(1.0 + exp( ( -5.0 - phi)/14.0) );
    const double tau_xs  = a_xs*b_xs;
    const double x_s     = GatingVarCalc(dt, x_s_, x_s_inf, tau_xs);
    const double I_Ks    = G_Ks*pow(x_s, 2)*(phi - E_Ks);

    // L-type Ca2+ current
    const double a_d      = 1.4/(1.0 + exp( (-35.0 - phi)/13.0) ) + 0.25;
    const double b_d      = 1.4/(1.0 + exp( (  5.0 + phi)/ 5.0) );
    const double g_d      = 1.0/(1.0 + exp( ( 50.0 - phi)/20.0) );
    const double a_fca    = 1.0/(1.0 + pow(  Ca_i_/0.000325, 8) );
    const double b_fca    = 0.1/(1.0 + exp( (Ca_i_ - 0.0005)/0.0001) );
    const double g_fca    = 0.2/(1.0 + exp( (Ca_i_ - 0.00075)/0.0008) );
    const double d_inf    = 1.0/(1.0 + exp( ( -5.0 - phi)/7.5) );
    const double f_inf    = 1.0/(1.0 + exp( ( 20.0 + phi)/7.0) );
    const double f_Ca_inf = (a_fca + b_fca + g_fca + 0.23)/1.46;
    const double tau_d    = a_d*b_d + g_d;
    const double tau_f    = 1125.0*exp( -pow( 27.0 + phi, 2 )/240 ) + 165.0/(1.0+exp( (25.0 - phi)/10.0) ) + 80.0;
    const double tau_f_Ca = 2.0; // [ms]
    const double d        = GatingVarCalc(dt, d_, d_inf, tau_d);
    const double f        = GatingVarCalc(dt, f_, f_inf, tau_f);
    double f_Ca = f_Ca_;
    if (f_Ca_inf < f_Ca_ || phi < -60.0)
    { f_Ca = GatingVarCalc(dt, f_Ca_, f_Ca_inf, tau_f_Ca); }
    const double I_CaL  = G_CaL*d*f*f_Ca*4*phi*pow(F, 2)/(R*T)*(Ca_i_*exp(2*phi*F/(R*T)) - 0.341*Ca_O)/(exp( 2.0*phi*F/(R*T) ) - 1.0);

    // Na+/Ca2+ exchanger current
    const double I_NaCa = k_NaCa*(exp( Gamma*phi*F/(R*T) )*pow( Na_i_, 3 )*Ca_O - exp( (Gamma - 1)*phi*F/(R*T) )*pow( Na_O, 3 )*Ca_i_*Alpha)/((pow( K_mNai, 3 ) + pow( Na_O, 3) )*(K_mCa + Ca_O)*(1 + k_sat*exp( (Gamma - 1)*phi*F/(R*T) )));

    // Na+/K+ pump current
    const double I_NaK  = P_NaK*K_O*Na_i_/((K_O + K_mK)*(Na_i_ + K_mNa)*(1+0.1245*exp(-0.1*phi*F/(R*T)) + 0.0353*exp(-phi*F/(R*T))));

    // Plateau currents
    const double I_pCa  = G_pCa*Ca_i_/(K_pCa + Ca_i_);
    const double I_pK   = G_pK*(phi-E_K)/(1 + exp((25.0 - phi)/5.98));

    // Background currents
    const double I_bCa  = G_bCa*(phi - E_Ca);
    const double I_bNa  = G_bNa*(phi - E_Na);

    // Compute reaction coefficient as the sum of all ion currents
    reacoeff = (I_Na + I_K1 + I_to + I_Kr + I_Ks + I_CaL + I_NaCa + I_NaK + I_pCa + I_pK + I_bCa + I_bNa);// /C_m;
  }
  else dserror("Myocard cell model type not found!");

  return reacoeff;
}


/*----------------------------------------------------------------------*
 |                                                           ljag 07/12 |
 *----------------------------------------------------------------------*/
double MAT::Myocard::ComputeReactionCoeffDeriv(const double phi, const double dt) const
{
  const double ReaCoeff_t1 = ComputeReactionCoeff(phi, dt);
  const double ReaCoeff_t2 = ComputeReactionCoeff((phi+params_->dt_deriv), dt);
  const double ReaCoeffDeriv = (ReaCoeff_t2 - ReaCoeff_t1)/(params_->dt_deriv);

  return ReaCoeffDeriv;

}


/*----------------------------------------------------------------------*
 |                                                           ljag 07/12 |
 *----------------------------------------------------------------------*/
double MAT::Myocard::GatingFunction(const double Gate1, const double Gate2, const double p, const double var, const double thresh) const
{
  double Erg = Gate1+(Gate2-Gate1)*(1.0+tanh(p*(var-thresh)))/2;
  return Erg;
}


/*----------------------------------------------------------------------*
 |                                                           ljag 07/12 |
 *----------------------------------------------------------------------*/
double MAT::Myocard::GatingVarCalc(const double dt, double y_0, const double y_inf, const double y_tau) const
{
  double Erg =  1.0/(1.0/dt + 1.0/y_tau)*(y_0/dt + y_inf/y_tau);
  return Erg;



}

/*----------------------------------------------------------------------*
 |  update of material at the end of a time step             ljag 07/12 |
 *----------------------------------------------------------------------*/
void MAT::Myocard::Update(const double phi, const double dt)
{
  const double p = 1000;

  if (*(params_->model) == "MV")
  {
    // Model parameter
    const double Theta_v = 0.3;
    const double Theta_w = 0.13;//0.015;
    const double Theta_vm = 0.006;//0.015;
    const double Theta_o = 0.006;
    const double Tau_v1m = 60.0;
    const double Tau_v2m = 1150.0;
    const double Tau_vp = 1.4506;
    const double Tau_w1m = 60.0;//70.0;
    const double Tau_w2m = 15.0;//20.0;
    const double k_wm = 65.0;
    const double u_wm = 0.03;
    const double Tau_wp = 200.0;//280.0;
    const double Tau_s1 = 2.7342;
    const double Tau_s2 = 16.0;//3.0;
    const double k_s = 2.0994;
    const double u_s = 0.9087;
    const double Tau_winf = 0.07;
    const double w_infs = 0.94;

    // calculate voltage dependent time constants ([7] page 545)
    const double Tau_vm = GatingFunction(Tau_v1m, Tau_v2m, p   , phi, Theta_vm);
    const double Tau_wm = GatingFunction(Tau_w1m, Tau_w2m, k_wm, phi, u_wm    );
    const double Tau_s  = GatingFunction(Tau_s1 , Tau_s2 , p   , phi, Theta_w);

    // calculate infinity values ([7] page 545)
    const double v_inf = GatingFunction(1.0, 0.0, p, phi, Theta_vm);
    const double w_inf = GatingFunction(1.0 - phi/Tau_winf, w_infs, p, phi, Theta_o);
    const double s_inf = GatingFunction(0.0, 1.0, k_s, phi, u_s);

    // calculate gating variables according to [8]
    const double Tau_v    = GatingFunction(Tau_vm, Tau_vp, p, phi, Theta_v);
    const double v_inf_GF = GatingFunction(v_inf , 0.0   , p, phi, Theta_v);
    const double v = GatingVarCalc(dt, v0_, v_inf_GF, Tau_v);

    const double Tau_w    = GatingFunction(Tau_wm, Tau_wp, p, phi, Theta_w);
    const double w_inf_GF = GatingFunction(w_inf , 0.0   , p, phi, Theta_w);
    const double w = GatingVarCalc(dt, w0_, w_inf_GF, Tau_w);

    const double s = GatingVarCalc(dt, s0_, s_inf, Tau_s);

    // update initial values for next time step
    v0_ = v;
    w0_ = w;
    s0_ = s;

  }

  else if (*(params_->model) == "TNNP")
  {
    // Model parameter
    const double R       = 8.3143; //
    const double T       = 310.0; //
    const double F       = 96.4867; //
    const double V_C     = 16404.0; //
    const double V_SR    = 1094.0; //
    const double K_O     = 5.4; //
    const double Na_O    = 140.0; //
    const double Ca_O    = 2.0; //
    const double G_Na    = 14.838; //
    const double G_K1    = 5.405; //
    double G_to = 0.294;
    if (*(params_->tissue) == "Endo") { G_to = 0.073; }
    const double G_Kr = 0.096; //
    double G_Ks = 0.245;
    if (*(params_->tissue) == "M") { G_Ks = 0.062; } //
    const double p_KNa   = 0.03; //
    const double G_CaL   = 1.75; //
    const double k_NaCa  = 1000; //
    const double Gamma   = 0.35; //
    const double K_mCa   = 1.38; //
    const double K_mNai  = 87.5; //
    const double k_sat   = 0.1; //
    const double Alpha   = 2.5; //
    const double P_NaK   = 1.362; //
    const double K_mK    = 1.0; //
    const double K_mNa   = 40.0; //
    const double G_pK    = 0.0146; //
    const double G_pCa   = 0.025; //
    const double K_pCa   = 0.0005; //
    const double G_bNa   = 0.00029; //
    const double G_bCa   = 0.000592; //
    const double V_maxup = 0.000425; //
    const double K_up    = 0.00025; //
    const double a_rel   = 16.464; //
    const double b_rel   = 0.25; //
    const double c_rel   = 8.232; //
    const double V_leak  = 0.00008; //
    const double Buf_c   = 0.15; //
    const double K_bufc  = 0.001; //
    const double Buf_sr  = 10.0; //
    const double K_bufsr = 0.3; //

    // Calculate reverse potentials
    const double E_Ks = R*T/(F)*log(K_O + p_KNa*Na_O/(K_i_ + p_KNa*Na_i_));
    const double E_K  = R*T/(F)*log(K_O/K_i_);
    const double E_Na = R*T/(F)*log(Na_O/Na_i_) ;
    const double E_Ca = R*T/(2*F)*log(Ca_O/Ca_i_) ;

    // Calculate ion channel currents

    // fast Natrium channel
    const double a_m   = 1.0/(1.0 + exp( (-60.0-phi)/5.0) );
    const double b_m   = 0.1/(1.0 + exp( ( 35.0+phi)/5.0)) + 0.1/(1.0 + exp( (-50.0 + phi)/200.0) );
    const double a_h   = GatingFunction( 0.057*exp( (-80.0 - phi)/6.8 ), 0.0, p, phi, -40.0);
    const double b_h   = GatingFunction( 2.7*exp(0.079*phi) + 3.1E5*exp( 0.3485*phi ), 0.77/(0.13*(1.0 + exp( -(10.66 + phi)/11.1 ))), p, phi, -40.0);
    const double a_j   = GatingFunction( ((-2.5428E4)*exp( 0.2444*phi ) - (6.948E-6)*exp( -0.0439*phi ))*(phi + 37.78)/(1 + exp( 0.311*(79.23 + phi) )), 0.0, p, phi, -40.0);
    const double b_j   = GatingFunction( 0.02424*exp( -0.01052*phi )/(1.0 + exp( -0.1378*(40.14+phi) )), 0.6*exp( 0.057*phi )/(1.0 + exp( -0.1*(32.0 + phi) )), p, phi, -40.0);
    const double m_inf = 1.0/pow( 1.0 + exp((-56.86 - phi)/9.03), 2 );
    const double h_inf = 1.0/pow( 1.0 + exp(( 71.55 + phi)/7.43), 2 );
    const double j_inf = 1.0/pow( 1.0 + exp(( 71.55 + phi)/7.43), 2 );
    const double tau_m = a_m*b_m;
    const double tau_h = 1/(a_h + b_h);
    const double tau_j = 1/(a_j + b_j);
    m_ = GatingVarCalc(dt, m_, m_inf, tau_m);
    h_ = GatingVarCalc(dt, h_, h_inf, tau_h);
    j_ = GatingVarCalc(dt, j_, j_inf, tau_j);
    const double I_Na  = G_Na*pow( m_, 3 )*h_*j_*(phi - E_Na);

    // Inward rectifier K+ current
    const double a_K1     = 0.1/(1.0 + exp( 0.06*(-E_K - 200.0 + phi)/20.0 ));
    const double b_K1     = (3.0*exp( 0.0002*(-E_K + 100.0 + phi) ) + exp( 0.1*(-E_K - 10.0 + phi) ))/(1.0 + exp( -0.5*(-E_K + phi) ));
    const double x_K1_inf = a_K1/(a_K1 + b_K1);
    const double I_K1     = G_K1*sqrt( K_O/5.4 )*x_K1_inf*(phi - E_K);


    // Transient outward current
    const double I_to   = G_to*r_*s_*(phi - E_K);

    // Rapid delayed rectifier current
    const double I_Kr     = G_Kr*sqrt( K_O/5.4 )*x_r1_*x_r2_*(phi - E_K);

    // Slow delayed rectifier current
    const double I_Ks   = G_Ks*pow( x_s_, 2 )*(phi - E_Ks);

    // L-type Ca2+ current
    const double a_d      = 1.4/(1.0 + exp( (-35.0 - phi)/13.0 ))+0.25;
    const double b_d      = 1.4/(1.0 + exp( (  5.0 + phi)/ 5.0 ));
    const double g_d      = 1.0/(1.0 + exp( ( 50.0 - phi)/20.0 ));
    const double a_fca    = 1.0/(1.0 + pow( Ca_i_/0.000325,8 ));
    const double b_fca    = 0.1/(1.0 + exp( (Ca_i_ - 0.0005 )/0.0001 ));
    const double g_fca    = 0.2/(1.0 + exp( (Ca_i_ - 0.00075)/0.0008 ));
    const double d_inf    = 1.0/(1.0 + exp( (-5.0 - phi)/7.5 ));
    const double f_inf    = 1.0/(1.0 + exp( (20.0 + phi)/7.0 ));
    const double f_Ca_inf = (a_fca + b_fca + g_fca + 0.23)/1.46;
    const double tau_d    = a_d*b_d + g_d;
    const double tau_f    = 1125.0*exp(-pow( 27.0+phi, 2 )/240) + 165.0/(1.0 + exp( (25.0 - phi)/10.0 )) + 80.0;
    const double tau_f_Ca = 2.0; // [ms]

    d_ = GatingVarCalc(dt, d_, d_inf, tau_d);
    f_ = GatingVarCalc(dt, f_, f_inf, tau_f);
    if (f_Ca_inf < f_Ca_ || phi < -60)
    { f_Ca_ = GatingVarCalc(dt, f_Ca_, f_Ca_inf, tau_f_Ca); }

    const double I_CaL  = G_CaL*d_*f_*f_Ca_*4*phi*pow( F, 2 )/(R*T)*(Ca_i_*exp( 2*phi*F/(R*T) ) - 0.341*Ca_O)/(exp( 2*phi*F/(R*T) ) - 1);

    // Na+/Ca2+ exchanger current
    const double I_NaCa = k_NaCa*(exp(Gamma*phi*F/(R*T))*pow(Na_i_, 3)*Ca_O - exp((Gamma - 1)*phi*F/(R*T))*pow(Na_O, 3)*Ca_i_*Alpha)/((pow(K_mNai, 3) + pow(Na_O, 3))*(K_mCa + Ca_O)*(1 + k_sat*exp((Gamma - 1)*phi*F/(R*T))));

    // Na+/K+ pump current
    const double I_NaK  = P_NaK*K_O*Na_i_/((K_O + K_mK)*(Na_i_ + K_mNa)*(1+0.1245*exp(-0.1*phi*F/(R*T)) + 0.0353*exp(-phi*F/(R*T))));

    // I_pCa
    const double I_pCa  = G_pCa*Ca_i_/(K_pCa + Ca_i_);

    // I_pK
    const double I_pK   = G_pK*(phi-E_K)/(1 + exp((25.0 - phi)/5.98));

    // Background currents
    const double I_bCa  = G_bCa*(phi - E_Ca);
    const double I_bNa  = G_bNa*(phi - E_Na);

    // Ionic concentrations
    Na_i_ += -dt*(I_Na + I_bNa + 3*I_NaK + 3*I_NaCa)/(V_C*F);
    K_i_  += -dt*(I_K1 + I_to + I_Kr + I_Ks - 2*I_NaK + I_pK /*+ I_stim - I_ax*/)/(V_C*F);

    // Calcium dynamics
    const double tau_g    = 2.0; // [ms]
    const double g_inf    = GatingFunction(1/(1+pow(Ca_i_/0.00035,6)), 1/(1+pow(Ca_i_/0.00035,16)), p, Ca_i_, 0.00035);
    g_ = GatingVarCalc(dt, g_, g_inf, tau_g);
    const double I_leak   = V_leak*(Ca_sr_ - Ca_i_); // Ca2+ leakage current from SR into cytoplasm
    const double I_up     = V_maxup/(1.0 + pow(K_up, 2)/pow(Ca_i_, 2)); // pump current taking up calcium in the SR
    const double I_rel    = (a_rel*pow(Ca_sr_,2)/(pow(b_rel,2) + pow(Ca_sr_,2)) + c_rel)*d_*g_; // Calcium induced calcium current

    const double Ca_ibufc   = Ca_i_*Buf_c/(Ca_i_ + K_bufc); // Buffered calcium in cytoplasm
    const double dCa_itotal = -dt*(I_CaL + I_bCa + I_pCa - 2*I_NaCa)/(2*V_C*F) + I_leak + I_up + I_rel; // Total calcium in cytoplasm
    const double bc         = Buf_c - Ca_ibufc - dCa_itotal - Ca_i_ + K_bufc;
    const double cc         = K_bufc*(Ca_ibufc + dCa_itotal + Ca_i_);
    Ca_i_ = (sqrt(bc*bc+4*cc)-bc)/2; // Free Ca2+ in cytoplasm

    const double Ca_srbufsr  = Ca_sr_*Buf_sr/(Ca_sr_ + K_bufsr); // Buffered calcium in SR
    const double dCa_srtotal = -dt*V_C*(I_leak - I_up + I_rel)/(V_SR); // Total calcium in SR
    const double bjsr        = Buf_sr - Ca_srbufsr - dCa_srtotal - Ca_sr_ + K_bufsr;
    const double cjsr        = K_bufsr*(Ca_srbufsr + dCa_srtotal + Ca_sr_);
    Ca_sr_ = (sqrt(bjsr*bjsr + 4*cjsr) - bjsr)/2; // Free Ca2+ in SR
  }
  else dserror("Myocard model type not found!");
  return;
}

void MAT::Myocard::SetupDiffusionTensor(const std::vector<double> &fiber1)
  {
  // get conductivity values of main fiber direction and perpendicular to fiber direction (rot symmetry)
   const double maindirdiffusivity = params_->maindirdiffusivity;
   const double offdirdiffusivity  = params_->offdirdiffusivity;

   // read local eigenvectors of diffusion tensor at current element
   std::vector<double> fiber2(3);
   std::vector<double> fiber3(3);

   // eigenvector matrix
   double fiber1norm=0.;
   double fiber2norm=0.;
   double fiber3norm=0.;
   LINALG::Matrix<3,3> evmat(true);
   LINALG::Matrix<3,3> evmatinv(true);

   // construction of ONB for fibers
   fiber2[0] = -fiber1[1];
   fiber2[1] = fiber1[0];
   fiber2[2] = 0;
   fiber3[0] = -fiber1[2]*fiber2[1];
   fiber3[1] = fiber1[2]*fiber2[0];
   fiber3[2] = fiber1[0]*fiber2[1]-fiber1[1]*fiber2[0];

   for (int i = 0; i < 3; ++i)
     {
       fiber1norm += fiber1[i]*fiber1[i];
       fiber2norm += fiber2[i]*fiber2[i];
       fiber3norm += fiber3[i]*fiber3[i];
     }
   fiber1norm = sqrt(fiber1norm);
   fiber2norm = sqrt(fiber2norm);
   fiber3norm = sqrt(fiber3norm);
   for (int i=0; i<3; ++i)
     {
       evmat(i,0) = fiber1[i]/fiber1norm;
       evmat(i,1) = fiber2[i]/fiber2norm;
       evmat(i,2) = fiber3[i]/fiber3norm;
     }

   //paranoia check if we really computed an ONB
   double test1 = 0.0;
   double test2 = 0.0;
   double test3 = 0.0;
   for (int i = 0; i < 3; ++i)
     {
       test1 += fiber1[i]*fiber2[i];
       test2 += fiber1[i]*fiber3[i];
       test3 += fiber2[i]*fiber3[i];
     }
   if (test1>1e-10 or test2>1e-10 or test3>1e-10)
     {
       dserror("ONB Calculation in Myocard Material failed");
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
   for (int i = 0; i<3; i++)
     {
       evmatinv(0,i) *= maindirdiffusivity;
       evmatinv(1,i) *= offdirdiffusivity;
       evmatinv(2,i) *= offdirdiffusivity;
     }

   difftensor_.Multiply(evmat, evmatinv);
  // cout << "Diffusion tensor" << endl;
   // cout << difftensor_ << endl;


   // done
   return;

  }
