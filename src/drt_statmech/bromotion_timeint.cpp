 /*!-------------------------------------------------------------------
\file BroMotion_timeInt.cpp

\brief This class provides a time integration for problems 
       regarding Brownian motion including mass terms
       and taking FSI computation results into account if present


<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*--------------------------------------------------------------------*/


#ifdef CCADISCRET

#include "bromotion_timeint.H"


/*-------------------------------------------------------------------*
 |  ctor (public)                                          umay 09/08|
 *-------------------------------------------------------------------*/
BroMotion_TimeInt::BroMotion_TimeInt(
    ParameterList&              params,
    DRT::Discretization&        dis,
    LINALG::Solver&             solver,
    IO::DiscretizationWriter&   output):
    StruGenAlpha(params,dis,solver,output)
{
  vector<DRT::Condition*> bromotioncond;
  discret_.GetCondition("BrownianMotion",bromotioncond);
  if (bromotioncond.size())
    bromotion_man_=rcp(new BroMotion_Manager(Discretization(),discret_));
  else
    dserror("apply Brownian motion conditions");
  
  dserror("not yet implemented");
      
}

    
   
/*-------------------------------------------------------------------*
 |  SetDefaults (public)                                   umay 09/08|
 *-------------------------------------------------------------------*/
void BroMotion_TimeInt::SetDefaults(ParameterList& params)
{
  params.set<string>("DYNAMICTYP"             ,"BrownianMotion");
  params.set<bool>  ("print to screen"        ,true);
  params.set<bool>  ("print to err"           ,false);
  params.set<FILE*> ("err file"               ,NULL);
  params.set<bool>  ("damping"                ,false);
  params.set<double>("damping factor K"       ,0.00001);
  params.set<double>("damping factor M"       ,0.00001);
  params.set<double>("beta"                   ,0.292);
#ifdef STRUGENALPHA_BE
  params.set<double>("delta"                  ,0.292);
#endif
  params.set<double>("gamma"                  ,0.581);
  params.set<double>("alpha m"                ,0.378);
  params.set<double>("alpha f"                ,0.459);
  params.set<double>("total time"             ,0.0);
  params.set<double>("delta time"             ,0.01);
  params.set<double>("max time"               ,0.02);
  params.set<int>   ("step"                   ,0);
  params.set<int>   ("nstep"                  ,5);
  params.set<int>   ("max iterations"         ,10);
  params.set<int>   ("num iterations"         ,-1);
  params.set<double>("tolerance displacements",1.0e-07);
  params.set<bool>  ("io structural disp"     ,false);
  params.set<int>   ("io disp every nstep"    ,10);
  params.set<string>("io structural stress"   ,"none");
  params.set<int>   ("io disp every nstep"    ,10);
  params.set<bool>  ("io structural strain"   ,false);
  params.set<int>   ("restart"                ,0);
  params.set<int>   ("write restart every"    ,0);
  params.set<bool>  ("contact"                ,false);
  // takes values "constant" consistent"
  params.set<string>("predictor"              ,"constant");
  // takes values "full newton" , "modified newton" , "nonlinear cg", "ptc"
  params.set<string>("equilibrium iteration"  ,"full newton");

  params.set<bool>  ("ADAPTCONV",false);
  params.set<double>("ADAPTCONV_BETTER",0.1);

  return;
}   



/*-------------------------------------------------------------------*
 |  integrate (public)                                     umay 09/08|
 *-------------------------------------------------------------------*/
void BroMotion_TimeInt::Integrate()
{
  ParameterList p;
  
  if (bromotion_man_!=null)
  {
    p.set("bromo_man", bromotion_man_);
    //bromotion_man_->EvaluatePotential(p,dism_,fint_,stiff_);
  }
}

#endif  // #ifdef CCADISCRET



