/*!----------------------------------------------------------------------
\file art_terminal_bc.cpp
\brief evaluation of 1d-artery inlet bc

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <stdio.h>

#include "art_terminal_bc.H"

#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_condition_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_ana.H"

#include "../drt_lib/drt_globalproblem.H"


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Solve (public)                                          ismail 08/09|
 |                                                                      |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::UTILS::SolvePrescribedTerminalBC(RefCountPtr<DRT::Discretization> actdis,
                                           const DRT::Condition *condition,
                                           ParameterList & params)
{

  // define BC name string (e.g: BC   = "inflow")
  string BC;
  // define BC type string (e.g: Type = "forced")
  string Type;

  // -------------------------------------------------------------------
  // Read in Condition type and name
  // -------------------------------------------------------------------
  Type = *(condition->Get<string>("type"));
  BC   = *(condition->Get<string>("boundarycond"));
  double time = params.get<double>("total time");

  // Define the reflection cooficient
  double Rf;

  // -------------------------------------------------------------------
  // Read in the bc curve information
  // -------------------------------------------------------------------
  const  vector<int>*    curve  = condition->Get<vector<int>    >("curve");
  double curvefac = 1.0;
  const  vector<double>* vals   = condition->Get<vector<double> >("val");

  // -------------------------------------------------------------------
  // Check whether the BC is absorbing or forced
  // -------------------------------------------------------------------
  if(Type == "absorbing")  // => without Reflection
  {
    Rf = 0.0;
  }
  else if(Type == "forced")// => with Reflection
  {
    // If forced curve exists => Rf = curve
    if ((*curve)[1]>=0)
    {
      curvefac = DRT::Problem::Instance()->Curve((*curve)[1]).f(time);
      Rf = (*vals)[1]*curvefac;
    }
    // else the BC is totally forced
    else
    {
      Rf = 1.0;
    }
    if (Rf <0.0 || Rf > 1.0)
    dserror("forced reflection (Rf = %f) should always belong to the range :[0  1.0]",Rf);
  }
  else
    dserror("%s is not defined as a 1D artery's inlet BC type",Type.c_str());

  // -------------------------------------------------------------------
  // Read in the value of the applied BC
  // -------------------------------------------------------------------
  double BCin;
  if((*curve)[0]>=0)
  {
    curvefac = DRT::Problem::Instance()->Curve((*curve)[0]).f(time);
    BCin = (*vals)[0]*curvefac;
  }
  else
    dserror("no inlet boundary condition defined!");

  // -------------------------------------------------------------------
  // Read in the parameters asosciated with the artery terminal to whom
  // the BC will be applied
  // -------------------------------------------------------------------

  double Wfnp, Wbnp;
  // IO = -1 if terminal is an inlet
  // IO =  1 if terminal is an outlet
  const int    IO    =  params.get<int>("in out flag");
  const double beta  =  params.get<double>("artery beta");
  const double Ao    =  params.get<double>("artery area");
  const double dens  =  params.get<double>("blood density");
  const double pext  =  params.get<double>("external pressure");

  if(IO == -1) // If BC is prescribed at the inlet
  {
    // read in backward characteristic wave speed
    Wbnp  =  params.get<double>("backward characteristic wave speed");

    // Initial backward characteristic speed at terminal 1
    const double Wbo   = -4.0*sqrt(beta/(2.0*dens*sqrt(Ao)));
    // backward characteristic wave, 
    const double Wb    = (Rf*Wbnp + (1.0-Rf)*Wbo);

    if(BC=="inflow")
    {

      /*
       Prescribed Volumetric flow rate:
                                                                                 
              /2.rho.Ao\2  /Wf - Wb\4  /Wf - Wb\                                
       Q    = |--------| . |-------| . |-------|                            
              \ beta   /   \   8   /   \   2   /                            
                                                                       
              /2.rho.Ao\2  /Wf - Wb\4  /Wf - Wb\                       
       f    = |--------| . |-------| . |-------|  - Q = 0               
              \ beta   /   \   8   /   \   2   /                        
                                                                          
        df    /2.rho.Ao\2  /Wf - Wb\3  /5*Wf - 3*Wb\                 
       ---- = |--------| . |-------| . |-----------|                 
       dWf    \ beta   /   \   8   /   \     16    /              
             
       The nonlinear equation: f could be solve using Newton-Raphson
       method as following:
                                                                                 
         1- U(first guess) = Q*(Ao) => W1(first guess) = 2Q/Ao - W2
         2- Calculate df/dWf
         3- Find Wf,i+1 = Wf,i - f,i/(df/dWf),i
         4- Calculate the Error f
         5- if Tolarance is Lager than Tolarance go to step (2)
      */
      double f;
      double dfdw;
      int itrs;
      itrs = 0;
      //step 1
      Wfnp   = 2.0*BCin/Ao - Wb;
      f      = pow(2.0*dens*Ao/beta,2)
              *pow((Wfnp - Wb)/8.0,4)*(  Wfnp +   Wb)/2.0;
      while(fabs(f)>0.00000001)
      {
        //step 2
        dfdw   =  pow(2.0*dens*Ao/beta,2)
                * pow((Wfnp - Wb)/8.0,3)*(5.0*Wfnp + 3.0*Wb)/16.0;
        //step 3
        Wfnp  = Wfnp - f/dfdw;
        //step 4
        f      =  pow(2.0*dens*Ao/beta,2)
                * pow((Wfnp - Wb)/8.0,4)*(  Wfnp +   Wb)/2.0
                - BCin;
        // an escape routine to prevent infinit loop
        itrs++;
        if(itrs>=30)
          dserror("Inflow boundary condition for Newton-Raphson exceeded the maximum allowed iterations");
      }
    }
    else if(BC == "velocity")
    {
      /*
       Prescribed Inlet Velocity
       Wf1 = 2*U - Wb1
      */
      Wfnp = 2.0*BCin - Wb;
    }
    else if(BC == "pressure")
    {
      /*
       Prescribed Inlet Pressure
       Wf1 = Wb1 + 8*sqrt((p-pext) + beta/sqrt(Ao))/(2.rho))
      */
      Wfnp = Wb + 8.0*sqrt((BCin-pext + beta/sqrt(Ao))/(2.0*dens));
    }
    else if(BC == "area")
    {
      /*
       Prescribed Inlet Area
       Wf1 = Wb1 + 8*sqrt(beta.A/(2.rho.Ao))
      */
      Wfnp = Wb + 8.0*sqrt(beta*sqrt(BCin)/(2.0*dens*Ao));
    }
    else if(BC == "characteristicWave")
    {
      /*
       Charachteristic wave
      */
      Wfnp = BCin;
    }
    else
      dserror("%s is not defined!",BC.c_str());
  }//If BC is prescribed at the inlet
  else if(IO == 1) // If BC is prescribed at the outlet
  {
    Wfnp  =  params.get<double>("forward characteristic wave speed");
    // Initial forward characteristic speed at terminal 2
    const double Wfo   =  4.0*sqrt(beta/(2.0*dens*sqrt(Ao)));
    // forward characteristic wave, 
    const double Wf    = (Rf*Wfnp + (1.0-Rf)*Wfo);

    if(BC=="inflow")
    {
      /*
       Prescribed Volumetric flow rate:
                                                                              
               /2.rho.Ao\2  /Wf - Wb\4  /Wf + Wb\                                
       Q    =  |--------| . |-------| . |-------|                            
               \ beta   /   \   8   /   \   2   /                            
                                                                       
               /2.rho.Ao\2  /Wf - Wb\4  /Wf + Wb\                       
       f    =  |--------| . |-------| . |-------|  - Q = 0               
               \ beta   /   \   8   /   \   2   /                        
                                                                          
        df      /2.rho.Ao\2  /Wf - Wb\3  /3*Wf + 5*Wb\                 
       ---- = - |--------| . |-------| . |-----------|                 
       dWb      \ beta   /   \   8   /   \     16    /              
             
       The nonlinear equation: f could be solve using Newton-Raphson
       method as following:
                                                                                    
         1- U(first guess) = Q*(Ao) => W2(first guess) = 2Q/Ao - W1
         2- Calculate df/dWb
         3- Find Wb,i+1 = Wb,i - f,i/(df/dWb),i
         4- Calculate the Error f
         5- if Tolarance is Lager than Tolarance go to step (2)

       */
      double f;
      double dfdw;
      int itrs;
      itrs = 0;
      //step 1
      Wbnp  = 2.0*BCin/Ao - Wf;
      f     = pow(2.0*dens*Ao/beta,2)
             *pow((Wf - Wbnp)/8.0,4)*(  Wf + Wbnp)/2.0;
      while(fabs(f)>0.000001)
      {
        //step 2
        dfdw   =-pow(2.0*dens*Ao/beta,2)
                *pow((Wf - Wbnp)/8.0,3)*(3.0*Wf + 5.0*Wbnp)/16.0;
        //step 3
        Wbnp   = Wbnp - f/dfdw;
        //step 4
        f      =  pow(2.0*dens*Ao/beta,2)
                * pow((Wf - Wbnp)/8.0,4)*( Wf +   Wbnp)/2.0
                - BCin;
  
        // a small routine to prevent infinit loop
        itrs++;
        if(itrs>=30)
          dserror("Inflow boundary condition for Newton-Raphson exceeded the maximum allowed iterations");
      }
    }
    else if(BC == "velocity")
    {
      /*
       Prescribed Inlet Velocity
       Wb2 = 2*U - Wf2
      */
      Wbnp = 2.0*BCin - Wf;
    }
    else if(BC == "pressure")
    {
      /*
       Prescribed Inlet Pressure
       Wb2 = Wf2 - 8*sqrt((p-pext) + beta/sqrt(Ao))/(2.rho))
      */
      Wbnp = Wf - 8.0*sqrt((BCin-pext + beta/sqrt(Ao))/(2.0*dens));
    }
    else if(BC == "area")
    {
      /*
       Prescribed Inlet Area
       Wb2 = Wf2 - 8*sqrt(beta.A/(2.rho.Ao))
      */
      Wbnp = Wf - 8.0*sqrt(beta*sqrt(BCin)/(2.0*dens*Ao));
    }
    else if(BC == "characteristicWave")
    {
      /*
       Charachteristic wave
      */
     Wbnp = BCin;
    }
  }//If BC is prescribed at the outlet
  else
    dserror("IO flag must be either 1 (for outlets) or 0 (for inlets)\n");

  // -------------------------------------------------------------------
  // finally return the updated value if the characteristic speeds
  // -------------------------------------------------------------------
  params.set<double>("forward characteristic wave speed",Wfnp);
  params.set<double>("backward characteristic wave speed",Wbnp);

}// void ART::UTILS::SolvePrescribedTerminalBC


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  SolveReflectiveTerminal (public)                        ismail 08/09|
 |  Solves the BC value for a certain Reflective value                  |
 |                                                                      |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::UTILS::SolveReflectiveTerminal(RefCountPtr<DRT::Discretization> actdis,
                                         const DRT::Condition *condition,
                                         ParameterList & params)
{

  // Define the reflection cooficient
  double Rf;

  // -------------------------------------------------------------------
  // Read in the bc curve information
  // -------------------------------------------------------------------
  const  vector<int>*    curve  = condition->Get<vector<int>    >("curve");
  double curvefac = 1.0;
  const  vector<double>* vals   = condition->Get<vector<double> >("val");

  // if the curve exist => Rf = val*curve(time)
  if ((*curve)[0]>=0)
  {
    double time = params.get<double>("total time");
    curvefac = DRT::Problem::Instance()->Curve((*curve)[0]).f(time);
    Rf = (*vals)[0]*curvefac;
  }
  // else Rf = val
  else
  {
    Rf = (*vals)[0];
  }

  // -------------------------------------------------------------------
  // Read in the parameters asosciated with the artery terminal to whom
  // the BC will be applied
  // -------------------------------------------------------------------
  const double beta  =  params.get<double>("artery beta");
  const double Ao    =  params.get<double>("artery area");
  const double dens  =  params.get<double>("blood density");
  const double Wfnp  =  params.get<double>("forward characteristic wave speed");


  // -------------------------------------------------------------------
  // Calculate the BC results
  // -------------------------------------------------------------------

  // Initial backward characteristic speed at terminal 1
  const double Wbo  = -4.0*sqrt(beta/(2.0*dens*sqrt(Ao)));
  const double Wfo  =  4.0*sqrt(beta/(2.0*dens*sqrt(Ao)));
  // backward characteristic wave, 
  const double Wbnp = -Rf*(Wfnp-Wfo) + Wbo;

  // -------------------------------------------------------------------
  // finally return the updated value if the characteristic speeds
  // -------------------------------------------------------------------
  params.set<double>("forward characteristic wave speed",Wfnp);
  params.set<double>("backward characteristic wave speed",Wbnp);

} //void ART::UTILS::SolveReflectiveTerminal



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  SolveExplWindkesselBC (public)                          ismail 08/09|
 |  Solves the windkessel boundary condition                            |
 |                                                                      |
 |  This code                                                           |
 |                                                                      |
 |                                                                      |
 |                                                                      |
 |                                                                      |
 |                                                                      |
 |                                                                      |
 |                                                                      |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::UTILS::SolveExplWindkesselBC(RefCountPtr<DRT::Discretization> actdis,
                                       const DRT::Condition *condition,
                                       ParameterList & params)
{
  // define BC windkessel inigration type string (e.g: BC   = "inflow")
  string int_type = *(condition->Get<string>("intigrationType"));
  // define windkessel BC type string (e.g: Type = "forced")
  string wk_type  = *(condition->Get<string>("windkesselType"));


  // -------------------------------------------------------------------
  // Read in the bc curve information
  // -------------------------------------------------------------------
  const  vector<int>*    curve  = condition->Get<vector<int>    >("curve");
  double curvefac = 1.0;
  const  vector<double>* vals   = condition->Get<vector<double> >("val");

  double Wb;
  if(int_type == "ExplicitWindkessel")
  {
    const double time = params.get<double>("total time");
    const double Wf   = params.get<double>("forward characteristic wave speed");
    const double Ao   = params.get<double>("artery area");
    const double beta = params.get<double>("artery beta");
    const double dens = params.get<double>("blood density");
    const double dt   = params.get<double>("time step size");
    const double Pext = params.get<double>("external pressure");
    double       Q    = params.get<double>("terminal volumetric flow rate");
    double       A    = params.get<double>("terminal cross-sectional area");
    double       P    = beta/Ao*(sqrt(A)-sqrt(Ao)) + Pext;


    // calculate the initial wave speed
    const double co = sqrt(beta/(2.0*dens*Ao));

    // -----------------------------------------------------------------
    // Find the type of windkessel model
    // -----------------------------------------------------------------
    if(wk_type == "R")  // a resister with a peripheral Pressure (Pout)
    {
      dserror("So far, only the 3 element windkessel model is implimented\n");
      // ---------------------------------------------------------------
      // Read in the wind kessel parameters
      // ---------------------------------------------------------------
      // define a terminal resistance
      double R, Rf;
      // read in the reflection value
      if ((*curve)[1]>=0)
      {
        curvefac = DRT::Problem::Instance()->Curve((*curve)[0]).f(time);
        R = (*vals)[1]*curvefac;
      }
      else
      {
        R = (*vals)[1];
      }

      if (R<0.0)
      {
        dserror("terminal resistance must be greater or equal to zero\n");
      }

      // calculate the 
      Rf = (R - dens*co/Ao)/ (R + dens*co/Ao);
      // calculate the backward charachteristic speed
      const double Wo = 4.0*sqrt(beta/(2.0*dens*sqrt(Ao)));
      Wb = -Rf*(Wf - Wo) - Wo;
    }// if(wk_type == "R")
    else if (wk_type == "RC") // an RC circuit with a peripheral Pressure (Pout)
    {
      dserror("So far, only the 3 element windkessel model is implimented\n");

      // ---------------------------------------------------------------
      // Read in the wind kessel parameters
      // ---------------------------------------------------------------
      // define the 2 element windkessel parameters
      double Pout, R, C;

      // Read in the periferal pressure of the wind kessel model
      if ((*curve)[0]>=0)
      {
        curvefac = DRT::Problem::Instance()->Curve((*curve)[0]).f(time);
        Pout = (*vals)[0]*curvefac;
      }
      else
      {
        Pout = (*vals)[2];
      }
      // read in the resistance value
      if ((*curve)[1]>=0)
      {
        curvefac = DRT::Problem::Instance()->Curve((*curve)[1]).f(time);
        R = (*vals)[1]*curvefac;
      }
      else
      {
        R = (*vals)[1];
      }
      // Read in the capacitance value
      if ((*curve)[2]>=0)
      {
        curvefac = DRT::Problem::Instance()->Curve((*curve)[2]).f(time);
        C = (*vals)[2]*curvefac;
      }
      else
      {
        C = (*vals)[2];
      }

      if (R <= 0.0 || C <= 0.0)
      {
        dserror("terminal resistance and capacitance must be always greater than zero\n");
      }

      // Calculate W2
      

    }//if (wk_type == "RC")
    else if (wk_type == "RCR") // The famous 3 element wind kessel model
    {
      // ---------------------------------------------------------------
      // Read in the wind kessel parameters
      // ---------------------------------------------------------------

      // define the 3 element windkessel parameters
      double Pout, R1, C, R2, Poutnm;

      // Read in the periferal pressure of the wind kessel model
      if ((*curve)[0]>=0)
      {
        curvefac = DRT::Problem::Instance()->Curve((*curve)[0]).f(time);
        Pout = (*vals)[0]*curvefac;
      }
      else
      {
        Pout = (*vals)[2];
      }
      // Read in Pout at time step n-1
      if ((*curve)[0]>=0)
      {
        double t;
        if (time <= dt)
          t = time;
        else 
          t = time -dt;
        curvefac = DRT::Problem::Instance()->Curve((*curve)[0]).f(t);
        Pout = (*vals)[0]*curvefac;
      }
      else
      {
        Pout = (*vals)[2];
      }
      // read in the source resistance value
      if ((*curve)[1]>=0)
      {
        curvefac = DRT::Problem::Instance()->Curve((*curve)[1]).f(time);
        R1 = (*vals)[1]*curvefac;
      }
      else
      {
        R1 = (*vals)[1];
      }
      // Read in the capacitance value
      if ((*curve)[2]>=0)
      {
        curvefac = DRT::Problem::Instance()->Curve((*curve)[2]).f(time);
        C = (*vals)[2]*curvefac;
      }
      else
      {
        C = (*vals)[2];
      }
      // read in the periferal resistance value
      if ((*curve)[3]>=0)
      {
        curvefac = DRT::Problem::Instance()->Curve((*curve)[3]).f(time);
        R2 = (*vals)[3]*curvefac;
      }
      else
      {
        R2 = (*vals)[3];
      }

      if (R1 < 0.0 || C <= 0.0 || R2 <= 0.0)
      {
        dserror("terminal resistances and capacitance must always be greater than zero\n");
      }

      // ---------------------------------------------------------------
      // Calculate Wb (backward characteristic wave speed)
      // ---------------------------------------------------------------

      // define inportant valriable
      double Pc, Qout, F, F_A, dFdA;

      // find Pc at time step n
      Pc   = P - Q*R1;

      // find Qout at time step n
      Qout = (Pc - Poutnm)/R2;

      // find Pc at n+1
      Pc   = Pc + dt/C*(Q - Qout);

      //solving the nonlinear equation F(A) = 0 using Newton-Raphson scheme

      F    = R1*Wf*A  - 4.0*R1*sqrt(beta/(2.0*dens*Ao)*sqrt(pow(A,5)))
                      - Pext - beta/Ao*(sqrt(A)-sqrt(Ao)) + Pc;
      int i = 0;
      F_A = pow(F*sqrt(Ao)/beta+1.0,2);
      while(fabs(F_A)>0.0000001)
      {
        dFdA = R1*Wf - 5.0*R1*sqrt(beta/(2.0*dens*Ao)*sqrt(A)) - 0.5*beta/(Ao*sqrt(A));
        A= A - F/dFdA;
        F = R1*Wf*A  - 4.0*R1*sqrt(beta/(2.0*dens*Ao)*sqrt(pow(A,5))) - Pext - beta/Ao*(sqrt(A)-sqrt(Ao)) + Pc;
        i++;
        F_A = pow(F*sqrt(Ao)/beta+1.0,2)-1.0;
        if(i>40)
          dserror("3 element windkessel Newton Raphson is not converging\n");
      }

      // finally find evaluate Wb
      Wb = Wf - 8.0*sqrt(beta/(2.0*dens*Ao)*sqrt(A));

    }// if (wk_type == "RCR")
    else if (wk_type == "RCRL") // four element windkessel model
    {

      dserror("So far, only the 3 element windkessel model is implimented\n");
      // ---------------------------------------------------------------
      // Read in the wind kessel parameters
      // ---------------------------------------------------------------

      // define the 4 element windkessel parameters
      double Pout, R1, C, R2, L;

      // Read in the periferal pressure of the wind kessel model
      if ((*curve)[0]>=0)
      {
        curvefac = DRT::Problem::Instance()->Curve((*curve)[0]).f(time);
        Pout = (*vals)[0]*curvefac;
      }
      else
      {
        Pout = (*vals)[2];
      }
      // read in the source resistance value
      if ((*curve)[1]>=0)
      {
        curvefac = DRT::Problem::Instance()->Curve((*curve)[1]).f(time);
        R1 = (*vals)[1]*curvefac;
      }
      else
      {
        R1 = (*vals)[1];
      }
      // Read in the capacitance value
      if ((*curve)[2]>=0)
      {
        curvefac = DRT::Problem::Instance()->Curve((*curve)[2]).f(time);
        C = (*vals)[2]*curvefac;
      }
      else
      {
        C = (*vals)[2];
      }
      // read in the periferal resistance value
      if ((*curve)[3]>=0)
      {
        curvefac = DRT::Problem::Instance()->Curve((*curve)[3]).f(time);
        R2 = (*vals)[3]*curvefac;
      }
      else
      {
        R2 = (*vals)[3];
      }
      // read in the inductance value
      if ((*curve)[4]>=0)
      {
        curvefac = DRT::Problem::Instance()->Curve((*curve)[4]).f(time);
        L = (*vals)[4]*curvefac;
      }
      else
      {
        L = (*vals)[4];
      }

      if (R1 <= 0.0 || C <= 0.0 || R2 <=0.0 || L <=0.0)
      {
        dserror("terminal resistance and capacitance must be always greater than zero\n");
      }
      // Calculate W2

    }//if (wk_type == "RCRL")
    else if (wk_type == "none")
    {
      dserror("So far, only the 3 element windkessel model is implimented\n");
    }
    else
    {
      dserror("\"%s\" is not supported type of windkessel model\n",wk_type.c_str());
    }

    // -----------------------------------------------------------------
    // return the calculated backward characteristic wave speed
    // -----------------------------------------------------------------
    params.set<double>("backward characteristic wave speed",Wb);

  }
  else
    dserror("so far windkessel BC supports only ExplicitWindkessel");

} //void ART::UTILS::SolveExplWindkesselBC


#endif // CCADISCRET
