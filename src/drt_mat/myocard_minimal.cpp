/*!----------------------------------------------------------------------
\file myocard.cpp

<pre>
Maintainer: Cristobal Bertogloi
            bertoglio@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*/

/*----------------------------------------------------------------------*
 |  headers                                                  ljag 07/12 |
 *----------------------------------------------------------------------*/

#include <vector>
#include "myocard_minimal.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |  Constructor                                    (public)  cbert 08/13 |
 *----------------------------------------------------------------------*/
Myocard_Minimal::Myocard_Minimal()
{

}


/*----------------------------------------------------------------------*
 |  Constructor                                    (public)  cbert 08/13 |
 *----------------------------------------------------------------------*/
Myocard_Minimal::Myocard_Minimal(const double eps_deriv_myocard,const std::string tissue)
:
  tools_()
{
  // Initial condition
  v0_ = 1.0;
  w0_ = 1.0;
  s0_ = 0.0;
  v_ = v0_;
  w_ = w0_;
  s_ = s0_;
  Jfi_ = 0.0; /// fast inward current
  Jso_ = 0.0; /// slow outward current
  Jsi_ = 0.0; /// slow inward current


  eps_deriv_ = eps_deriv_myocard;
  // Model parameter (To pack later in the parameter list of inputfile!)
  if (tissue == "M"){
    u_o_ = 0.0;
    u_u_ = 1.61;
    Theta_v_ = 0.3;
    Theta_w_ = 0.13;
    Theta_vm_ = 0.1;
    Theta_o_ = 0.005;
    Tau_v1m_ = 80.0;
    Tau_v2m_ = 1.4506;
    Tau_vp_ = 1.4506;
    Tau_w1m_ = 70.0;
    Tau_w2m_ = 8.0;
    k_wm_ = 200.0;
    u_wm_ = 0.016;
    Tau_wp_ = 280.0;
    Tau_fi_ = 0.078;
    Tau_o1_ = 410.0;
    Tau_o2_ = 7.0;
    Tau_so1_ = 91;
    Tau_so2_ = 0.8;
    k_so_ = 2.1;
    u_so_ = 0.6;
    Tau_s1_ = 2.7342;
    Tau_s2_ = 4.0;
    k_s_ = 2.0994;
    u_s_ = 0.9087;
    Tau_si_ = 3.3849;
    Tau_winf_ = 0.01;
    w_infs_ = 0.5;
  }else if (tissue == "Atria")
  {
    u_o_ = 0.0;
    u_u_ = 1.02;//1.58;
    Theta_v_ = 0.302;
    Theta_w_ = 0.33;//0.015;
    Theta_vm_ = 0.172;//0.015;
    Theta_o_ = 0.06;
    Tau_v1m_ = 65.6;
    Tau_v2m_ = 1150.0;
    Tau_vp_ = 0.95;
    Tau_w1m_ = 170.8;//70.0;
    Tau_w2m_ = 112.4;//20.0;
    k_wm_ = 135.0;
    u_wm_ = 0.0744;
    Tau_wp_ = 217.0;//280.0;
    Tau_fi_ = 0.0678;
    Tau_o1_ = 100.0;//6.0;
    Tau_o2_ = 64.87;
    Tau_so1_ = 53.54;//43.0;
    Tau_so2_ = 8.03;//0.2;
    k_so_ = 1.748;//2.0;
    u_so_ = 0.644;
    Tau_s1_ = 5.406;
    Tau_s2_ = 52.91;//3.0;
    k_s_ = 1.008;
    u_s_ = 0.814;
    Tau_si_ = 6.978;//2.8723;
    Tau_winf_ = 4.97;
    w_infs_ = 1.0;
  }
  else
  {
    dserror("Parameters for tissue type not supported for minimal model (only M, Atria)");
  }

  // Variables for electromechanical coupling
   mechanical_activation_ = 0.0;; // to store the variable for activation (phi in this case=)
   act_thres_ = 0.5; // activation threshold (so that activation = 1.0 if mechanical_activation_ >= act_thres_)


}

double Myocard_Minimal::ReaCoeff(const double phi, const double dt)
{

  double reacoeff = 0.0;
  const double p = 1000.0;

    // calculate voltage dependent time constants ([7] page 545)
    double Tau_vm = tools_.GatingFunction(Tau_v1m_, Tau_v2m_, p   , phi, Theta_vm_);
    double Tau_wm = tools_.GatingFunction(Tau_w1m_, Tau_w2m_, k_wm_, phi, u_wm_    );
    double Tau_so = tools_.GatingFunction(Tau_so1_, Tau_so2_, k_so_, phi, u_so_    );
    double Tau_s  = tools_.GatingFunction(Tau_s1_, Tau_s2_, p   , phi, Theta_w_);
    double Tau_o  = tools_.GatingFunction(Tau_o1_, Tau_o2_, p   , phi, Theta_o_);

    // calculate infinity values ([7] page 545)
    double v_inf = tools_.GatingFunction(1.0, 0.0, p, phi, Theta_vm_);
    double w_inf = tools_.GatingFunction(1.0 - phi/Tau_winf_, w_infs_, p, phi, Theta_o_);

    // calculate gating variables according to [8]
    double Tau_v    = tools_.GatingFunction(Tau_vm, Tau_vp_, p, phi, Theta_v_);
    double v_inf_GF = tools_.GatingFunction(v_inf , 0.0   , p, phi, Theta_v_);
    v_ = tools_.GatingVarCalc(dt, v0_, v_inf_GF, Tau_v);

    double Tau_w    = tools_.GatingFunction(Tau_wm, Tau_wp_, p, phi, Theta_w_);
    double w_inf_GF = tools_.GatingFunction(w_inf , 0.0   , p, phi, Theta_w_);
    w_ = tools_.GatingVarCalc(dt, w0_, w_inf_GF, Tau_w);

    const double s_inf = tools_.GatingFunction(0.0, 1.0, k_s_, phi, u_s_);
    s_ = tools_.GatingVarCalc(dt, s0_, s_inf, Tau_s);

    // calculate currents J_fi, J_so and J_si ([7] page 545)
    Jfi_ = -tools_.GatingFunction(0.0, v_*(phi - Theta_v_)*(u_u_ - phi)/Tau_fi_, p, phi, Theta_v_); // fast inward current
    Jso_ =  tools_.GatingFunction((phi - u_o_)/Tau_o, 1.0/Tau_so, p, phi, Theta_w_);// slow outward current
    Jsi_ = -tools_.GatingFunction(0.0, w_*s_/Tau_si_, p, phi, Theta_w_); // slow inward current

    reacoeff = (Jfi_ + Jso_ + Jsi_);

    // Store necessary variables for mechanical activation and electromechanical coupling
    mechanical_activation_ = phi;

  return reacoeff;
}

/*----------------------------------------------------------------------*
 |  returns number of internal state variables of the material  cbert 08/13 |
 *----------------------------------------------------------------------*/
int Myocard_Minimal::GetNumberOfInternalStateVariables() const
{
  return 3;
}

/*----------------------------------------------------------------------*
 |  returns current internal state of the material          cbert 08/13 |
 *----------------------------------------------------------------------*/
double Myocard_Minimal::GetInternalState(const int k) const
{
  double val=0.0;
  switch(k){
    case -1:{ // Compute activation function for electromechanical coupling
      if(mechanical_activation_>=act_thres_) val = 1.0;
      break;
    }
    case 0: {val=v_; break;}
    case 1: {val=w_; break;}
    case 2: {val=s_; break;}
  }
  return val;
}

/*----------------------------------------------------------------------*
 |  set  internal state of the material                     cbert 08/13 |
 *----------------------------------------------------------------------*/
void Myocard_Minimal::SetInternalState(const int k, const double val)
{
  switch(k){
    case 0: {v0_ = val; v_ = val; break;}
    case 1: {w0_ = val; w_ = val; break;}
    case 2: {s0_ = val; s_ = val; break;}
    default: {dserror("There are only 3 internal variables in this material!"); break;}
  }
}


/*----------------------------------------------------------------------*
 |  returns number of internal state variables of the material  cbert 08/13 |
 *----------------------------------------------------------------------*/
int Myocard_Minimal::GetNumberOfIonicCurrents() const
{
  return 3;
}

/*----------------------------------------------------------------------*
 |  returns current internal currents          cbert 08/13 |
 *----------------------------------------------------------------------*/
double Myocard_Minimal::GetIonicCurrents(const int k) const
{
  double val=0.0;
  switch(k){
    case 0: {val=Jfi_; break;}
    case 1: {val=Jso_; break;}
    case 2: {val=Jsi_; break;}
  }
  return val;
}

/*----------------------------------------------------------------------*
 |  update of material at the end of a time step             ljag 07/12 |
 *----------------------------------------------------------------------*/
void Myocard_Minimal::Update(const double phi, const double dt)
{
    // update initial values for next time step
    v0_ = v_;
    w0_ = w_;
    s0_ = s_;

    return;
}
