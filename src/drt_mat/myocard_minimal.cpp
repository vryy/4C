/*!----------------------------------------------------------------------
\file myocard_minimal.cpp

\brief minimal model for myocard material

<pre>
\level 3

\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
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
Myocard_Minimal::Myocard_Minimal(const double eps_deriv_myocard,const std::string tissue, int num_gp)
:
  tools_(),
  v0_(num_gp),
  w0_(num_gp),
  s0_(num_gp),
  v_(num_gp),
  w_(num_gp),
  s_(num_gp),
  Jfi_(num_gp),
  Jso_(num_gp),
  Jsi_(num_gp)

{

  for (int i=0; i<num_gp; ++i)
  {
    // Initial condition
    v0_[i] = 1.0;
    w0_[i] = 1.0;
    s0_[i] = 0.0;
    v_[i] = v0_[i];
    w_[i] = w0_[i];
    s_[i] = s0_[i];
    Jfi_[i] = 0.0; /// fast inward current
    Jso_[i] = 0.0; /// slow outward current
    Jsi_[i] = 0.0; /// slow inward current
  }


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
    }
    else if(tissue == "ENDO"){
      u_o_ = 0.0;
      u_u_ = 1.56;
      Theta_v_ = 0.3;
      Theta_w_ = 0.13;
      Theta_vm_ = 0.2;
      Theta_o_ = 0.006;
      Tau_v1m_ = 75.0;
      Tau_v2m_ = 10.0;
      Tau_vp_ = 1.4506;
      Tau_w1m_ = 6.0;
      Tau_w2m_ = 140.0;
      k_wm_ = 200.0;
      u_wm_ = 0.016;
      Tau_wp_ = 280.0;
      Tau_fi_ = 0.1;
      Tau_o1_ = 470.0;
      Tau_o2_ = 6.0;
      Tau_so1_ = 40.0;
      Tau_so2_ = 1.2;
      k_so_ = 2.0;
      u_so_ = 0.65;
      Tau_s1_ = 2.7342;
      Tau_s2_ = 2.0;
      k_s_ = 2.0994;
      u_s_ = 0.9087;
      Tau_si_ = 2.9013;
      Tau_winf_ = 0.0273;
      w_infs_ = 0.78;
    }
    else if(tissue == "EPI"){
      u_o_ = 0.0;
      u_u_ = 1.55;
      Theta_v_ = 0.3;
      Theta_w_ = 0.13;
      Theta_vm_ = 0.006;
      Theta_o_ = 0.006;
      Tau_v1m_ = 60.0;
      Tau_v2m_ = 1150;
      Tau_vp_ = 1.4506;
      Tau_w1m_ = 60.0;
      Tau_w2m_ = 15.0;
      k_wm_ = 65.0;
      u_wm_ = 0.03;
      Tau_wp_ = 200.0;
      Tau_fi_ = 0.11;
      Tau_o1_ = 400.0;
      Tau_o2_ = 6.0;
      Tau_so1_ = 30.0181;
      Tau_so2_ = 0.9957;
      k_so_ = 2.0458;
      u_so_ = 0.65;
      Tau_s1_ = 2.7342;
      Tau_s2_ = 16.0;
      k_s_ = 2.0994;
      u_s_ = 0.9087;
      Tau_si_ = 1.8875;
      Tau_winf_ = 0.07;
      w_infs_ = 0.94;
    }
    else if (tissue == "Atria")
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
    else if (tissue == "pAtria") // Lenk et al. 2015, "Initation of atrial fibrillation..." (JTB)
    {
      u_o_ = 0.0;
      u_u_ = 0.9205;
      Theta_v_ = 0.35;
      Theta_w_ = 0.328;
      Theta_vm_ = 0.126;
      Theta_o_ =  0.00005;
      Tau_v1m_ = 41.857;
      Tau_v2m_ = 1150;
      Tau_vp_ =  1.7;
      Tau_w1m_ = 138.69;
      Tau_w2m_ = 62.341;
      k_wm_ = 202.66;
      u_wm_ = 0.055;
      Tau_wp_ = 177.41;
      Tau_fi_ = 0.045;
      Tau_o1_ = 410;
      Tau_o2_ = 64.914;
      Tau_so1_ = 115;
      Tau_so2_ = 6.5;
      k_so_ = 1.386;
      u_so_ = 0.332;
      Tau_s1_ = 11.457;
      Tau_s2_ = 53.902;
      k_s_ = 1.226;
      u_s_ = 0.792;
      Tau_si_ = 7.802;
      Tau_winf_ = 0.05;
      w_infs_ = 1.0;
    }
    else if (tissue == "pAtriaRe") // Lenk et al. 2015, "Initation of atrial fibrillation..." (JTB)
    {
      u_o_ = 0.0;
      u_u_ = 0.9205;
      Theta_v_ = 0.35;
      Theta_w_ = 0.328;
      Theta_vm_ = 0.126;
      Theta_o_ =  0.00005;
      Tau_v1m_ = 41.857;
      Tau_v2m_ = 1150;
      Tau_vp_ =  1.7;
      Tau_w1m_ = 138.69;
      Tau_w2m_ = 62.341;
      k_wm_ = 202.66;
      u_wm_ = 0.055;
      Tau_wp_ = 177.41;
      Tau_fi_ = 0.11;
      Tau_o1_ = 410;
      Tau_o2_ = 64.914;
      Tau_so1_ = 115;
      Tau_so2_ = 6.5;
      k_so_ = 1.386;
      u_so_ = 0.332;
      Tau_s1_ = 11.457;
      Tau_s2_ = 53.902;
      k_s_ = 1.226;
      u_s_ = 0.792;
      Tau_si_ = 7.802;
      Tau_winf_ = 0.05;
      w_infs_ = 1.0;
    }
    else if (tissue == "rAtria")
      // Lenk et al. 2015, "Initation of atrial fibrillation..." (JTB);
      //Richter et al. 2016, "Anatomical and spiral wave..." (JTB)
    {
      u_o_ = 0.00;
      u_u_ = 1.0089;
      Theta_v_ = 0.3;
      Theta_w_ = 0.18171;
      Theta_vm_ = 0.1007;
      Theta_o_ = 0.015473;
      Tau_v1m_ = 16.3;
      Tau_v2m_ = 1150;
      Tau_vp_ = 1.7026;
      Tau_w1m_ = 79.963;
      Tau_w2m_ = 28.136;
      k_wm_ = 60.219;
      u_wm_ = 0.00991;
      Tau_wp_ = 213.55;
      Tau_fi_ = 0.083536;
      Tau_o1_ = 250.03;
      Tau_o2_ = 16.632;
      Tau_so1_ = 73.675;
      Tau_so2_ = 6.5537;
      k_so_ = 2.9748;
      u_so_ = 0.592093;
      Tau_s1_ = 9.876;
      Tau_s2_ = 4.2036;
      k_s_ = 2.2268;
      u_s_ = 0.81568;
      Tau_si_ = 10.699;
      Tau_winf_ = 0.2233;
      w_infs_ = 0.902;
    }
    else
    {
      dserror("Parameters for tissue type not supported for minimal model (only M, EPI, ENDO, Atria, pAtria, pAtriaRe, rAtria)");
    }

  // Variables for electromechanical coupling
   mechanical_activation_ = 0.0; // to store the variable for activation (phi in this case=)

}

double Myocard_Minimal::ReaCoeff(const double phi, const double dt)
{
  return Myocard_Minimal::ReaCoeff(phi, dt, 0);
}

double Myocard_Minimal::ReaCoeff(const double phi, const double dt, int gp)
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
    v_[gp] = tools_.GatingVarCalc(dt, v0_[gp], v_inf_GF, Tau_v);

    double Tau_w    = tools_.GatingFunction(Tau_wm, Tau_wp_, p, phi, Theta_w_);
    double w_inf_GF = tools_.GatingFunction(w_inf , 0.0   , p, phi, Theta_w_);
    w_[gp] = tools_.GatingVarCalc(dt, w0_[gp], w_inf_GF, Tau_w);

    const double s_inf = tools_.GatingFunction(0.0, 1.0, k_s_, phi, u_s_);
    s_[gp] = tools_.GatingVarCalc(dt, s0_[gp], s_inf, Tau_s);

    // calculate currents J_fi, J_so and J_si ([7] page 545)
    Jfi_[gp] = -tools_.GatingFunction(0.0, v_[gp]*(phi - Theta_v_)*(u_u_ - phi)/Tau_fi_, p, phi, Theta_v_); // fast inward current
    Jso_[gp] =  tools_.GatingFunction((phi - u_o_)/Tau_o, 1.0/Tau_so, p, phi, Theta_w_);// slow outward current
    Jsi_[gp] = -tools_.GatingFunction(0.0, w_[gp]*s_[gp]/Tau_si_, p, phi, Theta_w_); // slow inward current

    reacoeff = (Jfi_[gp] + Jso_[gp] + Jsi_[gp]);


    // Store necessary variables for mechanical activation and electromechanical coupling
    mechanical_activation_ = phi;

  return reacoeff;
}


double Myocard_Minimal::ReaCoeffN(const double phi, const double dt, int gp)
{
  double reacoeff = 0.0;
  const double p = 1000.0;

    // calculate voltage dependent time constants ([7] page 545)

    double Tau_so = tools_.GatingFunction(Tau_so1_, Tau_so2_, k_so_, phi, u_so_    );
    double Tau_o  = tools_.GatingFunction(Tau_o1_, Tau_o2_, p   , phi, Theta_o_);

    // calculate currents J_fi, J_so and J_si ([7] page 545)
    Jfi_[gp] = -tools_.GatingFunction(0.0, v0_[gp]*(phi - Theta_v_)*(u_u_ - phi)/Tau_fi_, p, phi, Theta_v_); // fast inward current
    Jso_[gp] =  tools_.GatingFunction((phi - u_o_)/Tau_o, 1.0/Tau_so, p, phi, Theta_w_);// slow outward current
    Jsi_[gp] = -tools_.GatingFunction(0.0, w0_[gp]*s0_[gp]/Tau_si_, p, phi, Theta_w_); // slow inward current

    reacoeff = (Jfi_[gp] + Jso_[gp] + Jsi_[gp]);


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
  return GetInternalState(k, 0);
}

/*----------------------------------------------------------------------*
 |  returns current internal state of the material       hoermann 09/15 |
 |  for multiple points per element                                     |
 *----------------------------------------------------------------------*/
double Myocard_Minimal::GetInternalState(const int k, int gp) const
{
  double val=0.0;
  switch(k){
    case -1: {val =mechanical_activation_; break;}
    case 0: {val=v_[gp]; break;}
    case 1: {val=w_[gp]; break;}
    case 2: {val=s_[gp]; break;}
  }
  return val;
}

/*----------------------------------------------------------------------*
 |  set  internal state of the material                     cbert 08/13 |
 *----------------------------------------------------------------------*/
void Myocard_Minimal::SetInternalState(const int k, const double val)
{
  SetInternalState(k, val, 0);
  return;
}

/*----------------------------------------------------------------------*
 |  set  internal state of the material                  hoermann 09/15 |
 |  for multiple points per element                                     |
 *----------------------------------------------------------------------*/
void Myocard_Minimal::SetInternalState(const int k, const double val, int gp)
{
  if (v0_.size() < (unsigned) gp )
    dserror("Number of gp does not match");
  switch(k){
    case -1: {mechanical_activation_ = val; break;}
    case 0: {v0_[gp] = val; v_[gp] = val; break;}
    case 1: {w0_[gp] = val; w_[gp] = val; break;}
    case 2: {s0_[gp] = val; s_[gp] = val; break;}
    default: {dserror("There are only 3 internal variables in this myocard material!"); break;}
  }
  return;
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
  return GetIonicCurrents(k, 0);
}

/*----------------------------------------------------------------------*
 |  returns current internal currents                    hoermann 09/15 |
 |  for multiple points per element                                     |
 *----------------------------------------------------------------------*/
double Myocard_Minimal::GetIonicCurrents(const int k, int gp) const
{
  double val=0.0;
  switch(k){
    case 0: {val=Jfi_[gp]; break;}
    case 1: {val=Jso_[gp]; break;}
    case 2: {val=Jsi_[gp]; break;}
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

/*----------------------------------------------------------------------*
 |  resize internal state variables                      hoermann 12/16 |
 *----------------------------------------------------------------------*/
void Myocard_Minimal::ResizeInternalStateVariables(int gp)
{
  v0_.resize(gp);
  w0_.resize(gp);
  s0_.resize(gp);
  v_.resize(gp);
  w_.resize(gp);
  s_.resize(gp);
  Jfi_.resize(gp);
  Jso_.resize(gp);
  Jsi_.resize(gp);
  return;
}

/*----------------------------------------------------------------------*
 |  get number of Gauss points                           hoermann 12/16 |
 *----------------------------------------------------------------------*/
int Myocard_Minimal::GetNumberOfGP() const
{
  return v_.size();
};
