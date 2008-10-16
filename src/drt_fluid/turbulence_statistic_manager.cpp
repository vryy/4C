/*!

Manage the computation of averages for several 
canonical flows like channel flow, flow around a square 
cylinder, flow in a lid driven cavity etc.

The manager is intended to remove as much of the averaging
overhead as possible from the time integration method.

Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235

*/
#ifdef CCADISCRET

#include "turbulence_statistic_manager.H"

namespace FLD
{

  /*----------------------------------------------------------------------

    Standard Constructor for Genalpha time integration (public)

  ----------------------------------------------------------------------*/
  TurbulenceStatisticManager::TurbulenceStatisticManager(FluidGenAlphaIntegration& fluid)
    :
    dt_         (fluid.dt_       ),
    alphaM_     (fluid.alphaM_   ),
    alphaF_     (fluid.alphaF_   ),
    gamma_      (fluid.gamma_    ),
    discret_    (fluid.discret_  ),
    params_     (fluid.params_   ),
    alefluid_   (fluid.alefluid_ ),
    myaccnp_    (fluid.accnp_    ),
    myaccn_     (fluid.accn_     ),
    myaccam_    (fluid.accam_    ),
    myvelnp_    (fluid.velnp_    ),
    myveln_     (fluid.veln_     ),
    myvelaf_    (fluid.velaf_    ),
    myvedeaf_   (fluid.vedeaf_   ),
    myvedenp_   (fluid.vedenp_   ),
    myveden_    (fluid.veden_    ),
    mydispnp_   (fluid.dispnp_   ),
    mydispn_    (fluid.dispn_    ),
    mygridveln_ (fluid.gridveln_ ),
    mygridvelaf_(fluid.gridvelaf_),
    myforce_    (fluid.force_    )
  {
    // the flow parameter will control for which geometry the 
    // sampling is done
    if(fluid.special_flow_=="channel_flow_of_height_2")
    {
      flow_=channel_flow_of_height_2;

      // allocate one instance of the averaging procedure for 
      // the flow under consideration
      statistics_channel_=rcp(new TurbulenceStatistics(discret_,params_));

      // prepare time averaging for subscales and residual 
      {
        // activate the computation of subgrid dissipation, 
        // residuals etc
        subgrid_dissipation_=true;

        // get ordered layers of elements to average
        RefCountPtr<vector<double> > planecoords;
      
        if (planecoords == null)
        {
          int size=(statistics_channel_->ReturnNodePlaneCoords()).size();
        
          planecoords    = rcp( new vector<double>(size));
          (*planecoords) = statistics_channel_->ReturnNodePlaneCoords();
        }

        //--------------------------------------------------
        // local_incrres(_sq)        (in plane) averaged values of resM (^2)
        // local_incrsacc(_sq)       (in plane) averaged values of sacc (^2)
        // local_incrsvelaf(_sq)     (in plane) averaged values of svelaf (^2)
        // local_incrresC(_sq)       (in plane) averaged values of resC (^2)
        // local_incrspressacc(_sq)  (in plane) averaged values of spressacc (^2)
        // local_incrspressnp(_sq)   (in plane) averaged values of spressnp (^2)
        //--------------------------------------------------
        RefCountPtr<vector<double> > local_incrres          = rcp(new vector<double> (3*(planecoords->size()-1),0.0));
        RefCountPtr<vector<double> > local_incrres_sq       = rcp(new vector<double> (3*(planecoords->size()-1),0.0));
        RefCountPtr<vector<double> > local_incrsacc         = rcp(new vector<double> (3*(planecoords->size()-1),0.0));
        RefCountPtr<vector<double> > local_incrsacc_sq      = rcp(new vector<double> (3*(planecoords->size()-1),0.0));
        RefCountPtr<vector<double> > local_incrsvelaf       = rcp(new vector<double> (3*(planecoords->size()-1),0.0));
        RefCountPtr<vector<double> > local_incrsvelaf_sq    = rcp(new vector<double> (3*(planecoords->size()-1),0.0));
        RefCountPtr<vector<double> > local_incrresC         = rcp(new vector<double> ((planecoords->size()-1)  ,0.0));
        RefCountPtr<vector<double> > local_incrresC_sq      = rcp(new vector<double> ((planecoords->size()-1)  ,0.0));
        RefCountPtr<vector<double> > local_incrspressacc    = rcp(new vector<double> ((planecoords->size()-1)  ,0.0));
        RefCountPtr<vector<double> > local_incrspressacc_sq = rcp(new vector<double> ((planecoords->size()-1)  ,0.0));
        RefCountPtr<vector<double> > local_incrspressnp     = rcp(new vector<double> ((planecoords->size()-1)  ,0.0));
        RefCountPtr<vector<double> > local_incrspressnp_sq  = rcp(new vector<double> ((planecoords->size()-1)  ,0.0));

        // pass pointers to local sum vectors to the element
        eleparams_.set<RefCountPtr<vector<double> > >("planecoords_"    ,planecoords           );
        eleparams_.set<RefCountPtr<vector<double> > >("incrres"         ,local_incrres         );
        eleparams_.set<RefCountPtr<vector<double> > >("incrres_sq"      ,local_incrres_sq      );
        eleparams_.set<RefCountPtr<vector<double> > >("incrsacc"        ,local_incrsacc        );
        eleparams_.set<RefCountPtr<vector<double> > >("incrsacc_sq"     ,local_incrsacc_sq     );
        eleparams_.set<RefCountPtr<vector<double> > >("incrsvelaf"      ,local_incrsvelaf      );
        eleparams_.set<RefCountPtr<vector<double> > >("incrsvelaf_sq"   ,local_incrsvelaf_sq   );
        eleparams_.set<RefCountPtr<vector<double> > >("incrresC"        ,local_incrresC        );
        eleparams_.set<RefCountPtr<vector<double> > >("incrresC_sq"     ,local_incrresC_sq     );
        eleparams_.set<RefCountPtr<vector<double> > >("incrspressacc"   ,local_incrspressacc   );
        eleparams_.set<RefCountPtr<vector<double> > >("incrspressacc_sq",local_incrspressacc_sq);
        eleparams_.set<RefCountPtr<vector<double> > >("incrspressnp"    ,local_incrspressnp    );
        eleparams_.set<RefCountPtr<vector<double> > >("incrspressnp_sq" ,local_incrspressnp_sq );

        // action for elements
        eleparams_.set("action","time average for subscales and residual");

        // other parameters that might be needed by the elements
        {
          ParameterList& timelist = eleparams_.sublist("time integration parameters");
        
          timelist.set("alpha_M",alphaM_);
          timelist.set("alpha_F",alphaF_);
          timelist.set("gamma"  ,gamma_ );
          timelist.set("dt"     ,dt_    );
        }
      
        // parameters for stabilisation
        {
          eleparams_.sublist("STABILIZATION")    = params_.sublist("STABILIZATION");
        }
      }
    }
    else if(fluid.special_flow_=="lid_driven_cavity")
    {
      flow_=lid_driven_cavity;

      // allocate one instance of the averaging procedure for 
      // the flow under consideration
      statistics_ldc_    =rcp(new TurbulenceStatisticsLdc(discret_,params_));
    }
    else if(fluid.special_flow_=="square_cylinder")
    {
      flow_=square_cylinder;

      // allocate one instance of the averaging procedure for 
      // the flow under consideration
      statistics_sqc_    =rcp(new TurbulenceStatisticsSqc(discret_,params_));
    }
    else
    {
      flow_=no_special_flow;
    }

    // do the time integration independent setup
    Setup();

    return;
    
  }

  /*----------------------------------------------------------------------

    Standard Constructor for One-Step-Theta time integration (public)

  ----------------------------------------------------------------------*/
  TurbulenceStatisticManager::TurbulenceStatisticManager(FluidImplicitTimeInt& fluid)
    :
    dt_         (fluid.dta_         ),
    alphaM_     (0.0                ),
    alphaF_     (0.0                ),
    gamma_      (0.0                ),
    discret_    (fluid.discret_     ),
    params_     (fluid.params_      ),
    alefluid_   (fluid.alefluid_    ),
    myaccnp_    (null               ),
    myaccn_     (fluid.accn_        ),
    myaccam_    (null               ),
    myvelnp_    (fluid.velnp_       ),
    myveln_     (fluid.veln_        ),
    myvelaf_    (null               ),
    myvedeaf_   (null               ),
    myvedenp_   (fluid.vedenp_      ),
    myveden_    (fluid.veden_       ),
    mydispnp_   (fluid.dispnp_      ),
    mydispn_    (fluid.dispn_       ),
    mygridveln_ (fluid.gridv_       ),
    mygridvelaf_(null               ),
    myforce_    (fluid.trueresidual_)
  {

    // no subgrid dissipation computation is available for the 
    // one-steo-theta implementation
    subgrid_dissipation_=false;

    // the flow parameter will control for which geometry the 
    // sampling is done
    if(fluid.special_flow_=="channel_flow_of_height_2")
    {
      flow_=channel_flow_of_height_2;

      // allocate one instance of the averaging procedure for 
      // the flow under consideration
      statistics_channel_=rcp(new TurbulenceStatistics(discret_,params_));
    }
    else if(fluid.special_flow_=="lid_driven_cavity")
    {
      flow_=lid_driven_cavity;

      // allocate one instance of the averaging procedure for 
      // the flow under consideration
      statistics_ldc_    =rcp(new TurbulenceStatisticsLdc(discret_,params_));
    }
    else if(fluid.special_flow_=="square_cylinder")
    {
      flow_=square_cylinder;

      // allocate one instance of the averaging procedure for 
      // the flow under consideration
      statistics_sqc_    =rcp(new TurbulenceStatisticsSqc(discret_,params_));
    }
    else
    {
      flow_=no_special_flow;
    }

    // do the time integration independent setup
    Setup();

    return;
    
  }
  
  /*
    Destructor
  */
  TurbulenceStatisticManager::~TurbulenceStatisticManager()
  {

    return;
  }

  /*----------------------------------------------------------------------
  
    Time integration independent setup called by Constructor (private)

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::Setup()
  {
    
    ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));

    // check if we want to compute averages of Smagorinsky
    // constants, effective viscosities etc
    smagorinsky_=false;
    if (modelparams->get<string>("TURBULENCE_APPROACH","DNS_OR_RESVMM_LES")
        ==
        "CLASSICAL_LES")
    {
      if(modelparams->get<string>("PHYSICAL_MODEL","no_model")
         ==
         "Dynamic_Smagorinsky"
         ||
         params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
         ==
         "Smagorinsky_with_van_Driest_damping"
        )
      {
        smagorinsky_=true;    
      }
    }
    
    // parameters for sampling/dumping period
    if (flow_ != no_special_flow)
    {
      samstart_  = modelparams->get<int>("SAMPLING_START",1         );
      samstop_   = modelparams->get<int>("SAMPLING_STOP", 1000000000);
      dumperiod_ = modelparams->get<int>("DUMPING_PERIOD",1         );
    }

    if(discret_->Comm().MyPID()==0)
    {
      string hom_plane
        =
        modelparams->get<string>("CHANNEL_HOMPLANE","not specified");

      if (flow_ == channel_flow_of_height_2)
      {
        cout << "Additional output          : " ;
        cout << "Turbulence statistics are evaluated ";
        cout << "for a turbulent channel flow.\n";
        cout << "                             " ;
        cout << "The solution is averaged over the homogeneous ";
        cout << hom_plane;
        cout << " plane and over time.\n";
        cout << "\n";
        cout << "                             " ;
        cout << "Sampling period: steps " << samstart_ << " to ";
        cout << modelparams->get<int>("SAMPLING_STOP",1000000000) << ".\n";
        
        int dumperiod = modelparams->get<int>("DUMPING_PERIOD",1);


        if(dumperiod == 0)
        {
          cout << "                             " ;
          cout << "Using standalone records (i.e. start from 0 for a new record)\n";
        }
        else
        {
          cout << "                             " ;
          cout << "Volker-style incremental dumping is used (";
          cout << dumperiod << ")" << endl;
        }

        cout << endl;
        cout << endl;
      }
    }

    // allocate one instance of the averaging procedure for 
    // each special flow under consideration
    switch(flow_)
    {
    case channel_flow_of_height_2:
    {
      if(smagorinsky_)
      {
        // get ordered layers of elements to average
        RefCountPtr<vector<double> > planecoords;

        if (planecoords == null)
        {
          int size=(statistics_channel_->ReturnNodePlaneCoords()).size();

          planecoords    = rcp( new vector<double>(size));
          (*planecoords) = statistics_channel_->ReturnNodePlaneCoords();
        }
        modelparams->set<RefCountPtr<vector<double> > >("planecoords_"         ,planecoords          );

        // extended statistics (plane average of Cs, (Cs_delta)^2, visceff)
        // for dynamic Smagorinsky model
        RefCountPtr<vector<double> > local_Cs_sum          =  rcp(new vector<double> (planecoords->size()-1,0.0));
        RefCountPtr<vector<double> > local_Cs_delta_sq_sum =  rcp(new vector<double> (planecoords->size()-1,0.0));
        RefCountPtr<vector<double> > local_visceff_sum     =  rcp(new vector<double> (planecoords->size()-1,0.0));
          
        modelparams->set<RefCountPtr<vector<double> > >("local_Cs_sum"         ,local_Cs_sum         );
        modelparams->set<RefCountPtr<vector<double> > >("local_Cs_delta_sq_sum",local_Cs_delta_sq_sum);
        modelparams->set<RefCountPtr<vector<double> > >("local_visceff_sum"    ,local_visceff_sum    );
      }
      break;
    }
    default:
    {
      break;
    }
    }
    
    return;
  }


  /*----------------------------------------------------------------------

    Store values computed during the element call

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::StoreElementValues(int step)
  {

    // sampling takes place only in the sampling period
    if(step>=samstart_ && step<=samstop_ && flow_ != no_special_flow)
    {
      switch(flow_)
      {
      case channel_flow_of_height_2:
      {
 
        // add computed dynamic Smagorinsky quantities 
        // (effective viscosity etc. used during the computation)

        if(smagorinsky_)
        {
          ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));

          // get ordered layers of elements to average
          RefCountPtr<vector<double> > planecoords;

          if (planecoords == null)
          {
            int size=(statistics_channel_->ReturnNodePlaneCoords()).size();

            planecoords    = rcp( new vector<double>(size));
            (*planecoords) = statistics_channel_->ReturnNodePlaneCoords();
          }

          RefCountPtr<vector<double> > global_incr_Cs_sum;
          RefCountPtr<vector<double> > local_Cs_sum;
          global_incr_Cs_sum          = rcp(new vector<double> (planecoords->size()-1,0.0));
          local_Cs_sum                = modelparams->get<RefCountPtr<vector<double> > >("local_Cs_sum"         );

        
          RefCountPtr<vector<double> > global_incr_Cs_delta_sq_sum;
          RefCountPtr<vector<double> > local_Cs_delta_sq_sum;
          global_incr_Cs_delta_sq_sum = rcp(new vector<double> (planecoords->size()-1,0.0));
          local_Cs_delta_sq_sum       = modelparams->get<RefCountPtr<vector<double> > >("local_Cs_delta_sq_sum");
        
        
          RefCountPtr<vector<double> > global_incr_visceff_sum;
          RefCountPtr<vector<double> > local_visceff_sum;
          global_incr_visceff_sum     = rcp(new vector<double> (planecoords->size()-1,0.0));
          local_visceff_sum           = modelparams->get<RefCountPtr<vector<double> > >("local_visceff_sum"    );

          // now add all the stuff from the different processors
          discret_->Comm().SumAll(&((*local_Cs_sum               )[0]),
                                  &((*global_incr_Cs_sum         )[0]),
                                  local_Cs_sum->size());
          discret_->Comm().SumAll(&((*local_Cs_delta_sq_sum      )[0]),
                                  &((*global_incr_Cs_delta_sq_sum)[0]),
                                  local_Cs_delta_sq_sum->size());
          discret_->Comm().SumAll(&((*local_visceff_sum          )[0]),
                                  &((*global_incr_visceff_sum    )[0]),
                                  local_visceff_sum->size());
        
          statistics_channel_->ReplaceCsIncrement(global_incr_Cs_sum,
                                                  global_incr_Cs_delta_sq_sum,
                                                  global_incr_visceff_sum);
        
          // reinitialise to zero for next element call
          local_Cs_sum          =  rcp(new vector<double> (planecoords->size()-1,0.0));
          local_Cs_delta_sq_sum =  rcp(new vector<double> (planecoords->size()-1,0.0));
          local_visceff_sum     =  rcp(new vector<double> (planecoords->size()-1,0.0));
          
          modelparams->set<RefCountPtr<vector<double> > >("local_Cs_sum"         ,local_Cs_sum         );
          modelparams->set<RefCountPtr<vector<double> > >("local_Cs_delta_sq_sum",local_Cs_delta_sq_sum);
          modelparams->set<RefCountPtr<vector<double> > >("local_visceff_sum"    ,local_visceff_sum    );
      }
    
      break;
    }
    default:
    {
      // there are no values to be stored in these cases
      break;
    }
    }
  }

    return;
  }

  /*----------------------------------------------------------------------

    Include current quantities in the time averaging procedure

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::DoTimeSample(int step, double time)
  {

    // sampling takes place only in the sampling period
    if(step>=samstart_ && step<=samstop_ && flow_ != no_special_flow)
    {

      //--------------------------------------------------
      // calculate means, fluctuations etc of velocity, 
      // pressure, boundary forces etc.
      switch(flow_)
      {
      case channel_flow_of_height_2:
      {
        if(statistics_channel_==null)
        {
          dserror("need statistics_channel_ to do a time sample for a turbulent channel flow");
        }
        statistics_channel_->DoTimeSample(myvelnp_,*myforce_);
        break;
      }
      case lid_driven_cavity:
      {
        if(statistics_sqc_==null)
        {
          dserror("need statistics_ldc_ to do a time sample for a cavity flow");
        }
        statistics_ldc_->DoTimeSample(myvelnp_);

        break;
      }
      case square_cylinder:
      {
        if(statistics_sqc_==null)
        {
          dserror("need statistics_sqc_ to do a time sample for a flow around a square cylinder");
        }
        statistics_sqc_->DoTimeSample(myvelnp_);
 
        // computation of Lift&Drag statistics

        {
          RCP<map<int,vector<double> > > liftdragvals;

          FLD::UTILS::LiftDrag(*discret_,*myforce_,params_,liftdragvals);

          if((*liftdragvals).size()!=1)
          {
            dserror("expecting only one liftdrag label for the sampling of a flow around a square cylinder");
          }
          map<int,vector<double> >::iterator theonlyldval = (*liftdragvals).begin();

          statistics_sqc_->DoLiftDragTimeSample(((*theonlyldval).second)[0],
                                                ((*theonlyldval).second)[1]);
        }
        
        break;
      }
      default:
      {
        break;
      }
      }

      //--------------------------------------------------
      // do averaging of residuals, dissipation rates etc 
      // (all gausspoint-quantities)
      if(subgrid_dissipation_)
      {
        switch(flow_)
        {
        case channel_flow_of_height_2:
        {
          if(statistics_channel_==null)
          {
            dserror("need statistics_channel_ to do a time sample for a turbulent channel flow");
          }

          // set vector values needed by elements
          discret_->ClearState();
          discret_->SetState("u and p (n+1      ,trial)",myvelnp_);
          discret_->SetState("u and p (n+alpha_F,trial)",myvelaf_);
          discret_->SetState("acc     (n+alpha_M,trial)",myaccam_);

          if (alefluid_)
          {
            discret_->SetState("dispnp"    , mydispnp_   );
            discret_->SetState("gridvelaf" , mygridvelaf_);
          }

          {
            ParameterList& timelist = eleparams_.sublist("time integration parameters");
            timelist.set("time"   ,time);
          }

          // call loop over elements to compute means
          {
            discret_->Evaluate(eleparams_,null,null,null,null,null);
            discret_->ClearState();
          }

          RefCountPtr<vector<double> > local_incrres         =eleparams_.get<RefCountPtr<vector<double> > >("incrres"         );
          RefCountPtr<vector<double> > local_incrres_sq      =eleparams_.get<RefCountPtr<vector<double> > >("incrres_sq"      );
          RefCountPtr<vector<double> > local_incrsacc        =eleparams_.get<RefCountPtr<vector<double> > >("incrsacc"        );
          RefCountPtr<vector<double> > local_incrsacc_sq     =eleparams_.get<RefCountPtr<vector<double> > >("incrsacc_sq"     );
          RefCountPtr<vector<double> > local_incrsvelaf      =eleparams_.get<RefCountPtr<vector<double> > >("incrsvelaf"      );
          RefCountPtr<vector<double> > local_incrsvelaf_sq   =eleparams_.get<RefCountPtr<vector<double> > >("incrsvelaf_sq"   );
          RefCountPtr<vector<double> > local_incrresC        =eleparams_.get<RefCountPtr<vector<double> > >("incrresC"        );
          RefCountPtr<vector<double> > local_incrresC_sq     =eleparams_.get<RefCountPtr<vector<double> > >("incrresC_sq"     );
          RefCountPtr<vector<double> > local_incrspressacc   =eleparams_.get<RefCountPtr<vector<double> > >("incrspressacc"   );
          RefCountPtr<vector<double> > local_incrspressacc_sq=eleparams_.get<RefCountPtr<vector<double> > >("incrspressacc_sq");
          RefCountPtr<vector<double> > local_incrspressnp    =eleparams_.get<RefCountPtr<vector<double> > >("incrspressnp"    );
          RefCountPtr<vector<double> > local_incrspressnp_sq =eleparams_.get<RefCountPtr<vector<double> > >("incrspressnp_sq" );

          int presize = local_incrresC->size();
          int velsize = local_incrres ->size();


          //--------------------------------------------------
          // (in plane) averaged values of resM (^2)
          //--------------------------------------------------
          RefCountPtr<vector<double> > global_incrres;
          global_incrres=  rcp(new vector<double> (velsize,0.0));
      
          RefCountPtr<vector<double> > global_incrres_sq;
          global_incrres_sq=  rcp(new vector<double> (velsize,0.0));
      
          //--------------------------------------------------
          // (in plane) averaged values of sacc (^2)
          //--------------------------------------------------
          RefCountPtr<vector<double> > global_incrsacc;
          global_incrsacc=  rcp(new vector<double> (velsize,0.0));
      
          RefCountPtr<vector<double> > global_incrsacc_sq;
          global_incrsacc_sq=  rcp(new vector<double> (velsize,0.0));
      
          //--------------------------------------------------
          // (in plane) averaged values of svelaf (^2)
          //--------------------------------------------------
          RefCountPtr<vector<double> > global_incrsvelaf;
          global_incrsvelaf=  rcp(new vector<double> (velsize,0.0));

          RefCountPtr<vector<double> > global_incrsvelaf_sq;
          global_incrsvelaf_sq=  rcp(new vector<double> (velsize,0.0));

          //--------------------------------------------------
          // (in plane) averaged values of resC (^2)
          //--------------------------------------------------
          RefCountPtr<vector<double> > global_incrresC;
          global_incrresC=  rcp(new vector<double> (presize,0.0));

          RefCountPtr<vector<double> > global_incrresC_sq;
          global_incrresC_sq=  rcp(new vector<double> (presize,0.0));
      
          //--------------------------------------------------
          // (in plane) averaged values of spressacc (^2)
          //--------------------------------------------------
          RefCountPtr<vector<double> > global_incrspressacc;
          global_incrspressacc=  rcp(new vector<double> (presize,0.0));

          RefCountPtr<vector<double> > global_incrspressacc_sq;
          global_incrspressacc_sq=  rcp(new vector<double> (presize,0.0));

          //--------------------------------------------------
          // (in plane) averaged values of spressnp (^2)
          //--------------------------------------------------
          RefCountPtr<vector<double> > global_incrspressnp;
          global_incrspressnp=  rcp(new vector<double> (presize,0.0));

          RefCountPtr<vector<double> > global_incrspressnp_sq;
          global_incrspressnp_sq=  rcp(new vector<double> (presize,0.0));


          // compute global sums, momentum equation residuals
          discret_->Comm().SumAll(&((*local_incrres       )[0]),
                                  &((*global_incrres      )[0]),
                                  velsize);
          discret_->Comm().SumAll(&((*local_incrres_sq    )[0]),
                                  &((*global_incrres_sq   )[0]),
                                  velsize);
      
          discret_->Comm().SumAll(&((*local_incrsacc      )[0]),
                                  &((*global_incrsacc     )[0]),
                                  velsize);
          discret_->Comm().SumAll(&((*local_incrsacc_sq   )[0]),
                                  &((*global_incrsacc_sq  )[0]),
                                  velsize);
      
          discret_->Comm().SumAll(&((*local_incrsvelaf    )[0]),
                                  &((*global_incrsvelaf   )[0]),
                                  velsize);
          discret_->Comm().SumAll(&((*local_incrsvelaf_sq )[0]),
                                  &((*global_incrsvelaf_sq)[0]),
                                  velsize);
      
          // compute global sums, incompressibility residuals
          discret_->Comm().SumAll(&((*local_incrresC      )[0]),
                                  &((*global_incrresC     )[0]),
                                  presize);
          discret_->Comm().SumAll(&((*local_incrresC_sq   )[0]),
                                  &((*global_incrresC_sq  )[0]),
                                  presize);
      
          discret_->Comm().SumAll(&((*local_incrspressacc    )[0]),
                                  &((*global_incrspressacc   )[0]),
                                  presize);
          discret_->Comm().SumAll(&((*local_incrspressacc_sq )[0]),
                                  &((*global_incrspressacc_sq)[0]),
                                  presize);
      
          discret_->Comm().SumAll(&((*local_incrspressnp     )[0]),
                                  &((*global_incrspressnp    )[0]),
                                  presize);
          discret_->Comm().SumAll(&((*local_incrspressnp_sq  )[0]),
                                  &((*global_incrspressnp_sq )[0]),
                                  presize);

          statistics_channel_->AddToResAverage(global_incrres         ,
                                               global_incrres_sq      ,
                                               global_incrsacc        ,
                                               global_incrsacc_sq     ,
                                               global_incrsvelaf      ,
                                               global_incrsvelaf_sq   ,
                                               global_incrresC        ,
                                               global_incrresC_sq     ,
                                               global_incrspressacc   ,
                                               global_incrspressacc_sq,
                                               global_incrspressnp    ,
                                               global_incrspressnp_sq);
      
          // reset working arrays
          local_incrres         =  rcp(new vector<double> (velsize,0.0));
          local_incrres_sq      =  rcp(new vector<double> (velsize,0.0));
          local_incrsacc        =  rcp(new vector<double> (velsize,0.0));
          local_incrsacc_sq     =  rcp(new vector<double> (velsize,0.0));
          local_incrsvelaf      =  rcp(new vector<double> (velsize,0.0));
          local_incrsvelaf_sq   =  rcp(new vector<double> (velsize,0.0));
          local_incrresC        =  rcp(new vector<double> (presize,0.0));
          local_incrresC_sq     =  rcp(new vector<double> (presize,0.0));
          local_incrspressacc   =  rcp(new vector<double> (presize,0.0));
          local_incrspressacc_sq=  rcp(new vector<double> (presize,0.0));
          local_incrspressnp    =  rcp(new vector<double> (presize,0.0));
          local_incrspressnp_sq =  rcp(new vector<double> (presize,0.0));

          eleparams_.set<RefCountPtr<vector<double> > >("incrres"         ,local_incrres         );
          eleparams_.set<RefCountPtr<vector<double> > >("incrres_sq"      ,local_incrres_sq      );
          eleparams_.set<RefCountPtr<vector<double> > >("incrsacc"        ,local_incrsacc        );
          eleparams_.set<RefCountPtr<vector<double> > >("incrsacc_sq"     ,local_incrsacc_sq     );
          eleparams_.set<RefCountPtr<vector<double> > >("incrsvelaf"      ,local_incrsvelaf      );
          eleparams_.set<RefCountPtr<vector<double> > >("incrsvelaf_sq"   ,local_incrsvelaf_sq   );
          eleparams_.set<RefCountPtr<vector<double> > >("incrresC"        ,local_incrresC        );
          eleparams_.set<RefCountPtr<vector<double> > >("incrresC_sq"     ,local_incrresC_sq     );
          eleparams_.set<RefCountPtr<vector<double> > >("incrspressacc"   ,local_incrspressacc   );
          eleparams_.set<RefCountPtr<vector<double> > >("incrspressacc_sq",local_incrspressacc_sq);
          eleparams_.set<RefCountPtr<vector<double> > >("incrspressnp"    ,local_incrspressnp    );
          eleparams_.set<RefCountPtr<vector<double> > >("incrspressnp_sq" ,local_incrspressnp_sq );
      
          break;
        }
        default:
        {
          break;
        }
        }
      }
    } // end step in sampling period

    return;
  }
  
  /*----------------------------------------------------------------------

    Write (dump) the statistics to a file

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::DoOutput(int step)
  {

    // sampling takes place only in the sampling period
    if(step>=samstart_ && step<=samstop_ && flow_ != no_special_flow)
    {
      enum format {write_single_record   ,
                   write_multiple_records,
                   do_not_write          } outputformat=do_not_write;
    

      // sampling a la Volker --- single record is 
      // constantly updated
      if(dumperiod_!=0)
      {
        
        int samstep = step-samstart_+1;
        double dsamstep=samstep;
        double ddumperiod=dumperiod_;

        // dump every dumperiod steps
        if (fmod(dsamstep,ddumperiod)==0)
        {
          outputformat=write_single_record;
        }
      }

      // sampling a la Peter --- for each sampling period a
      // new record is written; they can be combined by a 
      // postprocessing script to a single long term sample
      // (allows restarts during sampling)
      if(dumperiod_==0)
      {

        int upres    =params_.get("write solution every", -1);
        int uprestart=params_.get("write restart every" , -1);
        
        // dump in combination with a restart/output
        if(step%upres == 0 || step%uprestart == 0)
        {
          outputformat=write_multiple_records;
        }
      }
      
      // do actual output (an time averaging)
      switch(flow_)
      {
      case channel_flow_of_height_2:
      {
        if(statistics_channel_==null)
        {
          dserror("need statistics_channel_ to do a time sample for a turbulent channel flow");
        }

        if(outputformat == write_multiple_records)
        {
          statistics_channel_->TimeAverageMeansAndOutputOfStatistics(step);
          statistics_channel_->ClearStatistics();
        }
        if(outputformat == write_single_record)
        {
          statistics_channel_->DumpStatistics(step);
        }
        break;
      }
      case lid_driven_cavity:
      {
        if(statistics_ldc_==null)
        {
          dserror("need statistics_ldc_ to do a time sample for a lid driven cavity");
        }

        if(outputformat == write_single_record)
        {
          statistics_ldc_->DumpStatistics(step);
        }        
        break;
      }
      case square_cylinder:
      {
        if(statistics_sqc_==null)
        {
          dserror("need statistics_ldc_ to do a time sample for a square cylinder flow");
        }

        if(outputformat == write_single_record)
        {
          statistics_sqc_->DumpStatistics(step);
        }        
        
        break;
      }
      default:
      {
        break;
      }
      }
    } // end step is in sampling period

    return;
  }
  
  /*----------------------------------------------------------------------
  
  Clear all statistics collected up to now

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::Reset()
  {
    switch(flow_)
    {
    case channel_flow_of_height_2:
    {
      if(statistics_channel_==null)
      {
        dserror("need statistics_channel_ to do a time sample for a turbulent channel flow");
      }
      statistics_channel_->ClearStatistics();

      break;
    }
    case lid_driven_cavity:
    {
      // nothing to do here
      break;
    }
    case square_cylinder:
    {
      // nothing to do here
      break;
    }
    default:
    {
      // nothing to do here
      break;
    }
    }

    return;
  }

} // end namespace FLD

#endif  // #ifdef CCADISCRET
