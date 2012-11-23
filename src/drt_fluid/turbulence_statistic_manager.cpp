/*!

Manage the computation of averages for several
canonical flows like channel flow, flow around a square
cylinder, flow in a lid driven cavity, flow over a backward-facing step etc.

The manager is intended to remove as much of the averaging
overhead as possible from the time integration method.

Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236

*/

#include "turbulence_statistic_manager.H"
#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_combust/combust_fluidimplicitintegration.H"
#include "../drt_fluid/fluid_utils.H" // for LiftDrag
#include "../drt_lib/drt_dofset_independent_pbc.H"
#include "../drt_io/io_pstream.H"
#include "../drt_fluid/turbulence_statistics_mean_general.H"
#include "../drt_fluid/turbulence_statistics_ccy.H"
#include "../drt_fluid/turbulence_statistics_cha.H"
#include "../drt_fluid/turbulence_statistics_bcf.H"
#include "../drt_fluid/turbulence_statistics_ldc.H"
#include "../drt_fluid/turbulence_statistics_bfs.H"
#include "../drt_fluid/turbulence_statistics_oracles.H"
#include "../drt_fluid/turbulence_statistics_sqc.H"

namespace FLD
{

  /*----------------------------------------------------------------------

    Standard Constructor for standard time integration (public)

  ----------------------------------------------------------------------*/
  TurbulenceStatisticManager::TurbulenceStatisticManager(FluidImplicitTimeInt& fluid)
    :
    dt_              (fluid.dta_           ),
    discret_         (fluid.discret_       ),
    params_          (fluid.params_        ),
    alefluid_        (fluid.alefluid_      ),
    myaccnp_         (fluid.accnp_         ),
    myaccn_          (fluid.accn_          ),
    myaccam_         (fluid.accam_         ),
    myvelnp_         (fluid.velnp_         ),
    myveln_          (fluid.veln_          ),
    myvelaf_         (fluid.velaf_         ),
    myhist_          (fluid.hist_          ),
    myscaaf_         (fluid.scaaf_         ),
    myscaam_         (fluid.scaam_         ),
    mydispnp_        (fluid.dispnp_        ),
    mydispn_         (fluid.dispn_         ),
    mygridvelaf_     (fluid.gridv_         ),
    myforce_         (fluid.trueresidual_  ),
    myfilteredvel_   (fluid.filteredvel_   ),
    myfilteredreystr_(fluid.filteredreystr_),
    myfsvelaf_       (fluid.fsvelaf_       ),
    myfsscaaf_       (fluid.fsscaaf_       ),
    flow_            (no_special_flow      ),
    withscatra_      (false                ),
    turbmodel_      (INPAR::FLUID::no_model),
    subgrid_dissipation_(false             ),
    inflow_(false                          ),
    statistics_general_mean_(Teuchos::null ),
    statistics_channel_(Teuchos::null      ),
    statistics_channel_multiphase_(Teuchos::null),
    statistics_ccy_(Teuchos::null          ),
    statistics_ldc_(Teuchos::null          ),
    statistics_bfs_(Teuchos::null          ),
    statistics_oracles_(Teuchos::null      ),
    statistics_sqc_(Teuchos::null          )
  {
    subgrid_dissipation_ = DRT::INPUT::IntegralValue<int>(params_->sublist("TURBULENCE MODEL"),"SUBGRID_DISSIPATION");
    // initialize
    withscatra_ = false;

    // toogle statistics output for turbulent inflow
    inflow_ = DRT::INPUT::IntegralValue<int>(params_->sublist("TURBULENT INFLOW"),"TURBULENTINFLOW")==true;
    
    // toogle output of mean velocity for paraview
    out_mean_ = DRT::INPUT::IntegralValue<int>(params_->sublist("TURBULENCE MODEL"),"OUTMEAN")==true;

    // the flow parameter will control for which geometry the
    // sampling is done
    if(fluid.special_flow_=="channel_flow_of_height_2")
    {
      flow_=channel_flow_of_height_2;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_channel_=Teuchos::rcp(new TurbulenceStatisticsCha(discret_            ,
                                                          alefluid_           ,
                                                          mydispnp_           ,
                                                          *params_             ,
                                                          subgrid_dissipation_));
    }
    else if(fluid.special_flow_=="loma_channel_flow_of_height_2")
    {
      flow_=loma_channel_flow_of_height_2;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_channel_=Teuchos::rcp(new TurbulenceStatisticsCha(discret_            ,
                                                          alefluid_           ,
                                                          mydispnp_           ,
                                                          *params_             ,
                                                          subgrid_dissipation_));
    }
    else if(fluid.special_flow_=="scatra_channel_flow_of_height_2")
    {
      flow_=scatra_channel_flow_of_height_2;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_channel_=Teuchos::rcp(new TurbulenceStatisticsCha(discret_            ,
                                                          alefluid_           ,
                                                          mydispnp_           ,
                                                          *params_             ,
                                                          subgrid_dissipation_));
    }
    else if(fluid.special_flow_=="lid_driven_cavity")
    {
      flow_=lid_driven_cavity;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_ldc_    =Teuchos::rcp(new TurbulenceStatisticsLdc(discret_,*params_));
    }
    else if(fluid.special_flow_=="loma_lid_driven_cavity")
    {
      flow_=loma_lid_driven_cavity;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_ldc_    =Teuchos::rcp(new TurbulenceStatisticsLdc(discret_,*params_));
    }
    else if(fluid.special_flow_=="backward_facing_step")
    {
      flow_=backward_facing_step;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_bfs_ = Teuchos::rcp(new TurbulenceStatisticsBfs(discret_,*params_,"geometry_DNS_incomp_flow"));

      // statistics manager for turbulent boundary layer not available
      if (inflow_)
       dserror("The backward-facing step based on the geometry the DNS requires a turbulent boundary layer inflow profile which is not supported, yet!");
    }
    else if(fluid.special_flow_=="loma_backward_facing_step")
    {
      flow_=loma_backward_facing_step;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_bfs_ = Teuchos::rcp(new TurbulenceStatisticsBfs(discret_,*params_,"geometry_LES_flow_with_heating"));

      // build statistics manager for inflow channel flow
      if (inflow_)
      {
        if(params_->sublist("TURBULENT INFLOW").get<string>("CANONICAL_INFLOW")=="channel_flow_of_height_2"
         or params_->sublist("TURBULENT INFLOW").get<string>("CANONICAL_INFLOW")=="loma_channel_flow_of_height_2"
         or params_->sublist("TURBULENT INFLOW").get<string>("CANONICAL_INFLOW")=="scatra_channel_flow_of_height_2")
        {
          // do not write any dissipation rates for inflow channels
          subgrid_dissipation_ = false;
          // allocate one instance of the averaging procedure for the flow under consideration
          statistics_channel_=Teuchos::rcp(new TurbulenceStatisticsCha(discret_,
                                                              alefluid_,
                                                              mydispnp_,
                                                              *params_,
                                                              subgrid_dissipation_));
        }
      }
    }
    else if(fluid.special_flow_=="combust_oracles")
    {
      flow_=combust_oracles;

      if(discret_->Comm().MyPID()==0)
        std::cout << "---  setting up turbulence statistics manager for ORACLES ..." << std::flush;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_oracles_ = Teuchos::rcp(new COMBUST::TurbulenceStatisticsORACLES(discret_,*params_,"geometry_ORACLES",false));

      // build statistics manager for inflow channel flow
      if (inflow_)
      {
        if(params_->sublist("TURBULENT INFLOW").get<string>("CANONICAL_INFLOW")=="channel_flow_of_height_2")
        {
          // do not write any dissipation rates for inflow channels
          subgrid_dissipation_ = false;
          // allocate one instance of the averaging procedure for the flow under consideration
          statistics_channel_=Teuchos::rcp(new TurbulenceStatisticsCha(discret_,
                                                              alefluid_,
                                                              mydispnp_,
                                                              *params_,
                                                              subgrid_dissipation_));
        }
      }

      if(discret_->Comm().MyPID()==0)
        std::cout << " done" << std::endl;
    }
    else if(fluid.special_flow_=="square_cylinder")
    {
      flow_=square_cylinder;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_sqc_    =Teuchos::rcp(new TurbulenceStatisticsSqc(discret_,*params_));
    }
    else if(fluid.special_flow_=="square_cylinder_nurbs")
    {
      flow_=square_cylinder_nurbs;

      // do the time integration independent setup
      Setup();
    }
    else if(fluid.special_flow_=="rotating_circular_cylinder_nurbs")
    {
      flow_=rotating_circular_cylinder_nurbs;
      const bool withscatra = false;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_ccy_=Teuchos::rcp(new TurbulenceStatisticsCcy(discret_            ,
                                                      alefluid_           ,
                                                      mydispnp_           ,
                                                      *params_             ,
                                                      withscatra));
    }
    else if(fluid.special_flow_=="rotating_circular_cylinder_nurbs_scatra")
    {
      flow_=rotating_circular_cylinder_nurbs_scatra;
      const bool withscatra = true;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_ccy_=Teuchos::rcp(new TurbulenceStatisticsCcy(discret_            ,
                                                      alefluid_           ,
                                                      mydispnp_           ,
                                                      *params_             ,
                                                      withscatra));
    }
    else if(fluid.special_flow_=="time_averaging")
    {
      flow_=time_averaging;

      // do the time integration independent setup
      Setup();
    }
    else
    {
      flow_=no_special_flow;

      // do the time integration independent setup
      Setup();
    }

    // allocate one instance of the flow independent averaging procedure
    // providing colorful output for paraview
    if (out_mean_)
    {
      Teuchos::ParameterList *  modelparams =&(params_->sublist("TURBULENCE MODEL"));

      string homdir = modelparams->get<string>("HOMDIR","not_specified");

      if(flow_==rotating_circular_cylinder_nurbs_scatra)
      {
        // additional averaging of scalar field
        statistics_general_mean_
        =Teuchos::rcp(new TurbulenceStatisticsGeneralMean(
            discret_,
            homdir,
            fluid.VelPresSplitter(),true));
      }
      else
      {
        statistics_general_mean_
        =Teuchos::rcp(new TurbulenceStatisticsGeneralMean(
            discret_,
            homdir,
            fluid.VelPresSplitter(),false));
      }    
    }
    else statistics_general_mean_=Teuchos::null;

    return;

  }


  /*----------------------------------------------------------------------

    Standard Constructor for combustion One-Step-Theta time integration (public)

  ----------------------------------------------------------------------*/
  TurbulenceStatisticManager::TurbulenceStatisticManager(CombustFluidImplicitTimeInt& timeint)
    :
    dt_              (timeint.dta_           ),
    discret_         (timeint.discret_     ),
    params_          (timeint.params_      ),
    alefluid_        (false                ),
    myaccnp_         (Teuchos::null        ), // size is not fixed as we deal with xfem problems
    myaccn_          (Teuchos::null        ), // size is not fixed
    myaccam_         (Teuchos::null        ), // size is not fixed
    myveln_          (Teuchos::null        ), // size is not fixed
    myvelaf_         (Teuchos::null        ), // size is not fixed
    myhist_          (Teuchos::null        ), // size is not fixed
    myscaaf_         (Teuchos::null        ),
    myscaam_         (Teuchos::null        ),
    mydispnp_        (Teuchos::null        ),
    mydispn_         (Teuchos::null        ),
    mygridvelaf_     (Teuchos::null        ),
    myfilteredvel_   (Teuchos::null        ),
    myfilteredreystr_(Teuchos::null        ),
    myfsvelaf_       (Teuchos::null        ),
    myfsscaaf_       (Teuchos::null        ),
    flow_            (no_special_flow      ),
    withscatra_      (false                ),
    turbmodel_      (INPAR::FLUID::no_model),
    subgrid_dissipation_(false             ),
    inflow_(false                          ),
    statistics_general_mean_(Teuchos::null ),
    statistics_channel_(Teuchos::null      ),
    statistics_channel_multiphase_(Teuchos::null),
    statistics_ccy_(Teuchos::null          ),
    statistics_ldc_(Teuchos::null          ),
    statistics_bfs_(Teuchos::null          ),
    statistics_oracles_(Teuchos::null      ),
    statistics_sqc_(Teuchos::null          )
  {

    // subgrid dissipation
    subgrid_dissipation_ = false;
    // boolean for statistics of transported scalar
    withscatra_ = false;
    // toogle statistics output for turbulent inflow
    inflow_ = DRT::INPUT::IntegralValue<int>(params_->sublist("TURBULENT INFLOW"),"TURBULENTINFLOW")==true;
    // toogle output of mean velocity for paraview
    out_mean_ = DRT::INPUT::IntegralValue<int>(params_->sublist("TURBULENCE MODEL"),"OUTMEAN")==true;

    // the flow parameter will control for which geometry the
    // sampling is done
    if(timeint.special_flow_=="bubbly_channel_flow")
    {
      flow_=bubbly_channel_flow;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_channel_multiphase_ = Teuchos::rcp(new COMBUST::TurbulenceStatisticsBcf(discret_, *params_ ));
    }
    else if(timeint.special_flow_=="combust_oracles")
    {
      flow_=combust_oracles;
      // statistics for transported scalar (G-function)
      withscatra_ = true;

      if(discret_->Comm().MyPID()==0)
        std::cout << "---  setting up turbulence statistics manager for ORACLES ..." << std::flush;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_oracles_ = Teuchos::rcp(new COMBUST::TurbulenceStatisticsORACLES(discret_,*params_,"geometry_ORACLES",withscatra_));

      // build statistics manager for inflow channel flow
      if (inflow_)
      {
        if(params_->sublist("TURBULENT INFLOW").get<string>("CANONICAL_INFLOW")=="channel_flow_of_height_2")
        {
          // do not write any dissipation rates for inflow channels
          subgrid_dissipation_ = false;
          // allocate one instance of the averaging procedure for the flow under consideration
          statistics_channel_=Teuchos::rcp(new TurbulenceStatisticsCha(discret_,
                                                              alefluid_,
                                                              mydispnp_,
                                                              *params_,
                                                              subgrid_dissipation_));
        }
      }

      if(discret_->Comm().MyPID()==0)
        std::cout << " done" << std::endl;
    }
    else
    {
      flow_=no_special_flow;

      // do the time integration independent setup
      Setup();
    }

    statistics_general_mean_ = Teuchos::null;
    // allocate one instance of the flow independent averaging procedure
    // providing colorful output for paraview
    if (out_mean_)
    {
      Teuchos::ParameterList *  modelparams =&(params_->sublist("TURBULENCE MODEL"));

      string homdir = modelparams->get<string>("HOMDIR","not_specified");

      statistics_general_mean_ = Teuchos::rcp(new TurbulenceStatisticsGeneralMean(
          discret_,
          timeint.standarddofset_,
          homdir,
          *timeint.velpressplitterForOutput_,
          withscatra_ // statistics for transported scalar
      ));
    }

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

    Teuchos::ParameterList *  modelparams =&(params_->sublist("TURBULENCE MODEL"));

    if (modelparams->get<string>("TURBULENCE_APPROACH","DNS_OR_RESVMM_LES")
        ==
        "CLASSICAL_LES")
    {
      // check if we want to compute averages of Smagorinsky
      // constants, effective viscosities etc
      if(modelparams->get<string>("PHYSICAL_MODEL","no_model")
         ==
         "Dynamic_Smagorinsky"
         ||
         modelparams->get<string>("PHYSICAL_MODEL","no_model")
         ==
         "Smagorinsky_with_van_Driest_damping"
         ||
         modelparams->get<string>("PHYSICAL_MODEL","no_model")
         ==
         "Smagorinsky"
        )
      {
        turbmodel_ = INPAR::FLUID::dynamic_smagorinsky;
      }
      // check if we want to compute averages of scale similarity
      // quantities (tau_SFS)
      else if(modelparams->get<string>("PHYSICAL_MODEL","no_model")
              ==
              "Scale_Similarity")
      {
        turbmodel_ = INPAR::FLUID::scale_similarity_basic;
      }
      // check if we want to compute averages of multifractal
      // quantities (N, B)
      else if(modelparams->get<string>("PHYSICAL_MODEL","no_model")
           ==
           "Multifractal_Subgrid_Scales")
      {
        turbmodel_ = INPAR::FLUID::multifractal_subgrid_scales;
      }
    }
    else
      turbmodel_ = INPAR::FLUID::no_model;

    // parameters for sampling/dumping period
    if (flow_ != no_special_flow)
    {
      samstart_  = modelparams->get<int>("SAMPLING_START",1         );
      samstop_   = modelparams->get<int>("SAMPLING_STOP", 1000000000);
      dumperiod_ = modelparams->get<int>("DUMPING_PERIOD",1         );
    }
    else
    {
      samstart_  =0;
      samstop_   =0;
      dumperiod_ =0;
    }


    if(discret_->Comm().MyPID()==0)
    {

      if (flow_ == channel_flow_of_height_2 or
          flow_ == loma_channel_flow_of_height_2 or
          flow_ == scatra_channel_flow_of_height_2 or
          flow_ == bubbly_channel_flow)
      {
        string homdir
          =
          modelparams->get<string>("HOMDIR","not_specified");

        if(homdir!="xy" && homdir!="xz" && homdir!="yz")
        {
          dserror("need two homogeneous directions to do averaging in plane channel flows\n");
        }

        cout << "Additional output          : " ;
        cout << "Turbulence statistics are evaluated ";
        cout << "for a turbulent channel flow.\n";
        cout << "                             " ;
        cout << "The solution is averaged over the homogeneous ";
        cout << homdir;
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

    if(discret_->Comm().MyPID()==0)
    {
      if (flow_ == combust_oracles)
      {
        samstart_  = modelparams->get<int>("SAMPLING_START",1);
        samstop_   = modelparams->get<int>("SAMPLING_STOP", 1000000000);
        dumperiod_ = 0; // used as switch for the multi-record statistic output

        string homdir = modelparams->get<string>("HOMDIR","not_specified");
        if(homdir!="not_specified")
          dserror("there is no homogeneous direction for the ORACLES problem\n");

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
      case loma_channel_flow_of_height_2:
      case scatra_channel_flow_of_height_2:
      {
        // add computed dynamic Smagorinsky quantities
        // (effective viscosity etc. used during the computation)
        if(turbmodel_ == INPAR::FLUID::dynamic_smagorinsky)
        {
          if(discret_->Comm().MyPID()==0)
          {
            cout << "\nSmagorinsky constant, effective viscosity, ... etc, ";
            cout << "all element-quantities \n";
            }
          statistics_channel_->AddDynamicSmagorinskyQuantities();
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

    Store values computed during the element call

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::StoreNodalValues(
       int                        step,
       const RCP<Epetra_Vector>   stress12)
  {
    // sampling takes place only in the sampling period
    if(step>=samstart_ && step<=samstop_ && flow_ != no_special_flow)
    {
      switch(flow_)
      {
      case channel_flow_of_height_2:
      case loma_channel_flow_of_height_2:
      case scatra_channel_flow_of_height_2:
      {
        // add computed subfilter stress
        if(turbmodel_ == INPAR::FLUID::scale_similarity_basic) statistics_channel_->AddSubfilterStresses(stress12);
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
  void TurbulenceStatisticManager::DoTimeSample(int          step,
                                                const double eosfac,
                                                const double thermpressaf,
                                                const double thermpressam,
                                                const double thermpressdtaf,
                                                const double thermpressdtam)
  {

    // store Smagorinsky statistics if used
    StoreElementValues(step);

    // sampling takes place only in the sampling period
    if(step>=samstart_ && step<=samstop_ && flow_ != no_special_flow)
    {
      double tcpu=Teuchos::Time::wallTime();

      //--------------------------------------------------
      // calculate means, fluctuations etc of velocity,
      // pressure, boundary forces etc.
      switch(flow_)
      {
      case channel_flow_of_height_2:
      {
        if(statistics_channel_==Teuchos::null)
          dserror("need statistics_channel_ to do a time sample for a turbulent channel flow");

        statistics_channel_->DoTimeSample(myvelnp_,myforce_);
        break;
      }
      case loma_channel_flow_of_height_2:
      {
        if(statistics_channel_==Teuchos::null)
          dserror("need statistics_channel_ to do a time sample for a turbulent channel flow at low Mach number");

        statistics_channel_->DoLomaTimeSample(myvelnp_,myscaaf_,myforce_,eosfac);
        break;
      }
      case scatra_channel_flow_of_height_2:
      {
        if(statistics_channel_==Teuchos::null)
          dserror("need statistics_channel_ to do a time sample for a turbulent passive scalar transport in channel");

        statistics_channel_->DoScatraTimeSample(myvelnp_,myscaaf_,myforce_);
        break;
      }
      case lid_driven_cavity:
      {
        if(statistics_ldc_==Teuchos::null)
          dserror("need statistics_ldc_ to do a time sample for a cavity flow");

        statistics_ldc_->DoTimeSample(myvelnp_);
        break;
      }
      case loma_lid_driven_cavity:
      {
        if(statistics_ldc_==Teuchos::null)
          dserror("need statistics_ldc_ to do a time sample for a cavity flow at low Mach number");

        statistics_ldc_->DoLomaTimeSample(myvelnp_,myscaaf_,*myforce_,eosfac);
        break;
      }
      case backward_facing_step:
      {
        if(statistics_bfs_==Teuchos::null)
          dserror("need statistics_bfs_ to do a time sample for a flow over a backward-facing step");

        statistics_bfs_->DoTimeSample(myvelnp_);
        break;
      }
      case loma_backward_facing_step:
      {
        if(statistics_bfs_==Teuchos::null)
          dserror("need statistics_bfs_ to do a time sample for a flow over a backward-facing step at low Mach number");

        if (DRT::INPUT::get<INPAR::FLUID::PhysicalType>(*params_, "Physical Type") == INPAR::FLUID::incompressible)
        {

          if (not withscatra_)
          {
            statistics_bfs_->DoTimeSample(myvelnp_);

            // do time sample for inflow channel flow
            if (inflow_)
            {
              if(params_->sublist("TURBULENT INFLOW").get<string>("CANONICAL_INFLOW")=="channel_flow_of_height_2")
                statistics_channel_->DoTimeSample(myvelnp_,myforce_);
              else
               dserror("channel_flow_of_height_2 expected!");
            }
          }
          else
          {
            statistics_bfs_->DoScatraTimeSample(myvelnp_,myscaaf_);

            // do time sample for inflow channel flow
            if (inflow_)
            {
              if(params_->sublist("TURBULENT INFLOW").get<string>("CANONICAL_INFLOW")=="scatra_channel_flow_of_height_2")
                statistics_channel_->DoScatraTimeSample(myvelnp_,myscaaf_,myforce_);
              else
                dserror("scatra_channel_flow_of_height_2 expected!");
            }
          }
        }
        else
        {
          statistics_bfs_->DoLomaTimeSample(myvelnp_,myscaaf_,eosfac);

          // do time sample for inflow channel flow
          if (inflow_)
          {
            if(params_->sublist("TURBULENT INFLOW").get<string>("CANONICAL_INFLOW")=="loma_channel_flow_of_height_2")
            {
              statistics_channel_->DoLomaTimeSample(myvelnp_,myscaaf_,myforce_,eosfac);
            }
          }
        }
        break;
      }
      case combust_oracles:
      {
        subgrid_dissipation_ = false;

        if(statistics_oracles_==Teuchos::null)
          dserror("need statistics_oracles_ to do a time sample for an ORACLES flow step");

        statistics_oracles_->DoTimeSample(myvelnp_,myforce_,Teuchos::null,Teuchos::null);

        // build statistics manager for inflow channel flow
        if (inflow_)
        {
          if(params_->sublist("TURBULENT INFLOW").get<string>("CANONICAL_INFLOW")=="channel_flow_of_height_2")
            statistics_channel_->DoTimeSample(myvelnp_,myforce_);
          else
            dserror("loma_channel_flow_of_height_2 expected!");
        }
        break;
      }
      case square_cylinder:
      {
        if(statistics_sqc_==Teuchos::null)
          dserror("need statistics_sqc_ to do a time sample for a flow around a square cylinder");

        statistics_sqc_->DoTimeSample(myvelnp_);

        // computation of Lift&Drag statistics
        {
          RCP<map<int,std::vector<double> > > liftdragvals;

          FLD::UTILS::LiftDrag(*discret_,*myforce_,*params_,liftdragvals);

          if((*liftdragvals).size()!=1)
          {
            dserror("expecting only one liftdrag label for the sampling of a flow around a square cylinder");
          }
          std::map<int,std::vector<double> >::iterator theonlyldval = (*liftdragvals).begin();

          statistics_sqc_->DoLiftDragTimeSample(((*theonlyldval).second)[0],
                                                ((*theonlyldval).second)[1]);
        }
        break;
      }
      case rotating_circular_cylinder_nurbs:
      {

        if(statistics_ccy_==Teuchos::null)
          dserror("need statistics_ccy_ to do a time sample for a flow in a rotating circular cylinder");

        statistics_ccy_->DoTimeSample(myvelnp_,Teuchos::null,Teuchos::null);
        break;
      }
      case rotating_circular_cylinder_nurbs_scatra:
      {

        if(statistics_ccy_==Teuchos::null)
          dserror("need statistics_ccy_ to do a time sample for a flow in a rotating circular cylinder");

        statistics_ccy_->DoTimeSample(myvelnp_,myscaaf_,myfullphinp_);
        break;
      }
      default:
      {
        break;
      }
      }

      if(discret_->Comm().MyPID()==0)
      {
        cout << "Computed statistics: mean values, fluctuations, boundary forces etc.             (";
        printf("%10.4E",Teuchos::Time::wallTime()-tcpu);
        cout << ")";
      }

      //--------------------------------------------------
      // do averaging of residuals, dissipation rates etc
      // (all gausspoint-quantities)
      if(subgrid_dissipation_)
      {
        tcpu=Teuchos::Time::wallTime();

        switch(flow_)
        {
        case channel_flow_of_height_2:
        case loma_channel_flow_of_height_2:
        case scatra_channel_flow_of_height_2:
        {
          if(statistics_channel_==Teuchos::null)
          {
            dserror("need statistics_channel_ to do a time sample for a turbulent channel flow");
          }
          //dserror("no subgrid dissipation");

          // set vector values needed by elements
          map<string,RCP<Epetra_Vector> > statevecs;
          map<string,RCP<Epetra_MultiVector> > statetenss;

          statevecs.insert(std::pair<string,RCP<Epetra_Vector> >("hist",myhist_));
          statevecs.insert(std::pair<string,RCP<Epetra_Vector> >("accam",myaccam_));
          statevecs.insert(std::pair<string,RCP<Epetra_Vector> >("scaaf",myscaaf_));
          statevecs.insert(std::pair<string,RCP<Epetra_Vector> >("scaam",myscaam_));

          if (alefluid_)
          {
            statevecs.insert(std::pair<string,RCP<Epetra_Vector> >("dispnp", mydispnp_   ));
            statevecs.insert(std::pair<string,RCP<Epetra_Vector> >("gridv" , mygridvelaf_));
          }

          if (DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo") == INPAR::FLUID::timeint_afgenalpha
            or DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo") == INPAR::FLUID::timeint_npgenalpha)
          {
            statevecs.insert(std::pair<string,RCP<Epetra_Vector> >("velaf",myvelaf_));
            if (DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo") == INPAR::FLUID::timeint_npgenalpha)
              statevecs.insert(std::pair<string,RCP<Epetra_Vector> >("velnp",myvelnp_));
          }
          else
            statevecs.insert(std::pair<string,RCP<Epetra_Vector> >("velaf",myvelnp_));

          if (params_->sublist("TURBULENCE MODEL").get<string>("FSSUGRVISC")!= "No"
              or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
          {
            statevecs.insert(std::pair<string,RCP<Epetra_Vector> >("fsvelaf",myfsvelaf_));
            if (myfsvelaf_==Teuchos::null)
              dserror ("Have not got fsvel!");
            statevecs.insert(std::pair<string,RCP<Epetra_Vector> >("fsscaaf",myfsscaaf_));
            if (myfsscaaf_==Teuchos::null)
                          dserror ("Have not got fssca!");
          }
          if (turbmodel_ == INPAR::FLUID::scale_similarity_basic)
          {
            statetenss.insert(std::pair<string,RCP<Epetra_MultiVector> >("filtered vel",myfilteredvel_));
            statetenss.insert(std::pair<string,RCP<Epetra_MultiVector> >("filtered reystr",myfilteredreystr_));
          }

          statistics_channel_->EvaluateResiduals(statevecs,statetenss,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam);

          break;
        }
        default:
        {
          break;
        }
        }

        if(discret_->Comm().MyPID()==0)
        {
          cout << "\nresiduals, dissipation rates etc, ";
          cout << "all gausspoint-quantities (";
          printf("%10.4E",Teuchos::Time::wallTime()-tcpu);
          cout << ")";
        }
      }
      if(discret_->Comm().MyPID()==0)
      {
        cout << "\n";
      }

      if(turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
      {
        switch(flow_)
        {
          case channel_flow_of_height_2:
          {
            // add parameters of multifractal subgrid-scales model
            statistics_channel_->AddModelParamsMultifractal(myvelaf_,myfsvelaf_,false);
            break;
          }
          case scatra_channel_flow_of_height_2:
          {
            // add parameters of multifractal subgrid-scales model
            statistics_channel_->AddModelParamsMultifractal(myvelaf_,myfsvelaf_,true);
            break;
          }
          default:
          {
            break;
          }
        }
      }

      // add vector(s) to general mean value computation
      // scatra vectors may be Teuchos::null
      if (statistics_general_mean_!=Teuchos::null)
        statistics_general_mean_->AddToCurrentTimeAverage(dt_,myvelnp_,myscaaf_,myfullphinp_);

    } // end step in sampling period

    return;
  }


  /*----------------------------------------------------------------------

    Include current quantities in the time averaging procedure

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::DoTimeSample(int                             step,
                                                Teuchos::RCP<Epetra_Vector>     velnp,
                                                Teuchos::RCP<Epetra_Vector>     force,
                                                Teuchos::RCP<Epetra_Vector>     phi,
                                                Teuchos::RCP<const DRT::DofSet> stddofset,
                                                Teuchos::RCP<const Epetra_Vector> discretmatchingvelnp /*= Teuchos::null */ // needed for 'bubbly_channel_flow'
)
  {
    // sampling takes place only in the sampling period
    if(step>=samstart_ && step<=samstop_ && flow_ != no_special_flow)
    {
      double tcpu=Teuchos::Time::wallTime();

      //--------------------------------------------------
      // calculate means, fluctuations etc of velocity,
      // pressure, boundary forces etc.
      switch(flow_)
      {
      case bubbly_channel_flow:
      {
        if(statistics_channel_multiphase_ == Teuchos::null)
          dserror("need statistics_channel_multiphase_ to do a time sample for a turbulent channel flow");

        if (velnp == Teuchos::null        or force == Teuchos::null
            or stddofset == Teuchos::null or discretmatchingvelnp == Teuchos::null
            or phi == Teuchos::null)
            dserror("The multi phase channel statistics need a current velnp, force, stddofset, discretmatchingvelnp, phinp.");

        statistics_channel_multiphase_->DoTimeSample(velnp, force, stddofset, discretmatchingvelnp, phi);
        break;
      }
      case combust_oracles:
      {
        subgrid_dissipation_ = false;

        if(statistics_oracles_==Teuchos::null)
          dserror("need statistics_oracles_ to do a time sample for an ORACLES flow step");

        statistics_oracles_->DoTimeSample(velnp,force,phi,stddofset);

        // build statistics manager for inflow channel flow
        if (inflow_)
        {
          if(params_->sublist("TURBULENT INFLOW").get<string>("CANONICAL_INFLOW")=="channel_flow_of_height_2")
          {
            statistics_channel_->DoTimeSample(velnp,force);
          }
        }
        break;
      }
      default:
      {
        dserror("called wrong DoTimeSample() for this kind of special flow");
        break;
      }
      }

      // add vector(s) to general mean value computation
      // scatra vectors may be Teuchos::null
      if (statistics_general_mean_!=Teuchos::null)
        statistics_general_mean_->AddToCurrentTimeAverage(dt_,velnp,myscaaf_,myfullphinp_);

      if(discret_->Comm().MyPID()==0)
      {
        cout << "                      taking time sample (";
        printf("%10.4E",Teuchos::Time::wallTime()-tcpu);
        cout << ")\n";
      }

    } // end step in sampling period

    return;
  }

  /*----------------------------------------------------------------------

    Write (dump) the statistics to a file

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::DoOutput(IO::DiscretizationWriter& output,
                                            int                       step,
                                            const double              inflow)
  {
    // sampling takes place only in the sampling period
    if(step>=samstart_ && step<=samstop_ && flow_ != no_special_flow)
    {
      enum format {write_single_record   ,
                   write_multiple_records,
                   do_not_write          } outputformat=do_not_write;
      bool output_inflow = false;

      // sampling a la Volker --- single record is constantly updated
      if(dumperiod_!=0)
      {
        int samstep = step-samstart_+1;

        // dump every dumperiod steps
        if (samstep%dumperiod_==0) outputformat=write_single_record;
      }

      // sampling a la Peter --- for each sampling period a
      // new record is written; they can be combined by a
      // postprocessing script to a single long term sample
      // (allows restarts during sampling)
      if(dumperiod_==0)
      {
        int upres    =params_->get<int>("write solution every");
        int uprestart=params_->get<int>("write restart every" );

        // dump in combination with a restart/output
        if((step%upres == 0 || ( uprestart > 0 && step%uprestart == 0) ) && step>samstart_)
          outputformat=write_multiple_records;
      }
      if(inflow_)
      {
        int upres    =params_->get<int>("write solution every");
        int uprestart=params_->get<int>("write restart every" );

        // dump in combination with a restart/output
        if((step%upres == 0 || ( uprestart > 0 && step%uprestart == 0) ) && step>samstart_)
          output_inflow=true;
      }

      if (discret_->Comm().MyPID()==0 && outputformat != do_not_write )
        std::cout << "---  statistics record: \n" << std::flush;

      // do actual output (time averaging)
      switch(flow_)
      {
      case channel_flow_of_height_2:
      {
        if(statistics_channel_==Teuchos::null)
          dserror("need statistics_channel_ to do a time sample for a turbulent channel flow");


        if(outputformat == write_multiple_records)
        {
          statistics_channel_->TimeAverageMeansAndOutputOfStatistics(step);
          statistics_channel_->ClearStatistics();
        }

        if(outputformat == write_single_record)
          statistics_channel_->DumpStatistics(step);
        break;
      }
      case loma_channel_flow_of_height_2:
      {
        if(statistics_channel_==Teuchos::null)
          dserror("need statistics_channel_ to do a time sample for a turbulent channel flow at low Mach number");

        if(outputformat == write_single_record)
          statistics_channel_->DumpLomaStatistics(step);
        break;
      }
      case scatra_channel_flow_of_height_2:
      {
        if(statistics_channel_==Teuchos::null)
          dserror("need statistics_channel_ to do a time sample for a turbulent channel flow at low Mach number");

        if(outputformat == write_single_record)
          statistics_channel_->DumpScatraStatistics(step);
        break;
      }
      case bubbly_channel_flow:
      {
        if(statistics_channel_multiphase_==Teuchos::null)
          dserror("need statistics_channel_multiphase_ to do a time sample for a turbulent channel flow");

        if(outputformat == write_multiple_records)
        {
          statistics_channel_multiphase_->TimeAverageMeansAndOutputOfStatistics(step);
          statistics_channel_multiphase_->ClearStatistics();
        }

        if(outputformat == write_single_record)
          statistics_channel_multiphase_->DumpStatistics(step);
        break;
      }
      case lid_driven_cavity:
      {
        if(statistics_ldc_==Teuchos::null)
          dserror("need statistics_ldc_ to do a time sample for a lid driven cavity");

        if(outputformat == write_single_record)
          statistics_ldc_->DumpStatistics(step);
        break;
      }
      case loma_lid_driven_cavity:
      {
        if(statistics_ldc_==Teuchos::null)
          dserror("need statistics_ldc_ to do a time sample for a lid driven cavity at low Mach number");

        if(outputformat == write_single_record)
          statistics_ldc_->DumpLomaStatistics(step);
        break;
      }
      case backward_facing_step:
      {
        if(statistics_bfs_==Teuchos::null)
          dserror("need statistics_bfs_ to do a time sample for a flow over a backward-facing step");

        if(outputformat == write_single_record)
          statistics_bfs_->DumpStatistics(step);

        //write statistics of inflow channel flow
        if (inflow_)
        {
          if(params_->sublist("TURBULENT INFLOW").get<string>("CANONICAL_INFLOW")=="channel_flow_of_height_2")
          {
            if(output_inflow)
            {
              statistics_channel_->TimeAverageMeansAndOutputOfStatistics(step);
              statistics_channel_->ClearStatistics();
            }
          }
        }
        break;
      }
      case loma_backward_facing_step:
      {
        if(statistics_bfs_==Teuchos::null)
          dserror("need statistics_bfs_ to do a time sample for a flow over a backward-facing step at low Mach number");

        if (DRT::INPUT::get<INPAR::FLUID::PhysicalType>(*params_, "Physical Type") == INPAR::FLUID::incompressible)
        {
          if (not withscatra_)
          {
            if(outputformat == write_single_record)
              statistics_bfs_->DumpStatistics(step);

            // write statistics of inflow channel flow
            if (inflow_)
            {
              if(params_->sublist("TURBULENT INFLOW").get<string>("CANONICAL_INFLOW")=="channel_flow_of_height_2")
              {
                if(output_inflow)
                {
                  statistics_channel_->TimeAverageMeansAndOutputOfStatistics(step);
                  statistics_channel_->ClearStatistics();
                }
              }
            }
          }
          else
          {
            if(outputformat == write_single_record)
              statistics_bfs_->DumpScatraStatistics(step);

            // write statistics of inflow channel flow
            if (inflow_)
            {
              if(params_->sublist("TURBULENT INFLOW").get<string>("CANONICAL_INFLOW")=="scatra_channel_flow_of_height_2")
               {
                if(outputformat == write_single_record)
                  statistics_channel_->DumpScatraStatistics(step);
               }
            }
         }
        }
        else // loma
        {
          if(outputformat == write_single_record)
            statistics_bfs_->DumpLomaStatistics(step);

          // write statistics of inflow channel flow
          if (inflow_)
          {
            if(params_->sublist("TURBULENT INFLOW").get<string>("CANONICAL_INFLOW")=="loma_channel_flow_of_height_2")
            {
              if(outputformat == write_single_record)
                statistics_channel_->DumpLomaStatistics(step);
            }
          }
        }
        break;
      }
      case combust_oracles:
      {
        if(statistics_oracles_==Teuchos::null)
          dserror("need statistics_oracles_ to do a time sample for an ORACLES flow step");
        if(outputformat == write_multiple_records)
        {
          statistics_oracles_->TimeAverageStatistics();
          statistics_oracles_->OutputStatistics(step);
          statistics_oracles_->ClearStatistics();
        }

        // build statistics manager for inflow channel flow
        if (inflow_)
        {
          if(params_->sublist("TURBULENT INFLOW").get<string>("CANONICAL_INFLOW")=="channel_flow_of_height_2")
          {
            if(outputformat == write_multiple_records)
            {
              statistics_channel_->TimeAverageMeansAndOutputOfStatistics(step);
              statistics_channel_->ClearStatistics();
            }
          }
        }
        break;
      }
      case square_cylinder:
      {
        if(statistics_sqc_==Teuchos::null)
          dserror("need statistics_sqc_ to do a time sample for a square cylinder flow");

        if (outputformat == write_single_record)
          statistics_sqc_->DumpStatistics(step);
        break;
      }
      case rotating_circular_cylinder_nurbs:
      case rotating_circular_cylinder_nurbs_scatra:
      {
        if(statistics_ccy_==Teuchos::null)
          dserror("need statistics_ccy_ to do a time sample for a flow in a rotating circular cylinder");

        statistics_ccy_->TimeAverageMeansAndOutputOfStatistics(step);
        statistics_ccy_->ClearStatistics();
        break;
      }
      default:
      {
        break;
      }
      }

//      if (discret_->Comm().MyPID()==0 && outputformat != do_not_write )
//        std::cout << "done" << std::endl;

      if(discret_->Comm().MyPID()==0 && outputformat != do_not_write)
      {
        cout << "XXXXXXXXXXXXXXXXXXXXX              ";
        cout << "wrote statistics record            ";
        cout << "XXXXXXXXXXXXXXXXXXXXX";
        cout << "\n\n";
      }


      // dump general mean value output in combination with a restart/output
      // don't write output if turbulent inflow or twophaseflow is computed
      if (!inflow and flow_ != bubbly_channel_flow)
      {
        int upres    =params_->get<int>("write solution every");
        int uprestart=params_->get<int>("write restart every" );

        if((step%upres == 0 || (uprestart > 0 && step%uprestart == 0)) &&  (statistics_general_mean_!=Teuchos::null))
          statistics_general_mean_->WriteOldAverageVec(output);
      }
    } // end step is in sampling period

    return;
  } // DoOutput


  /*----------------------------------------------------------------------

  Add results from scalar transport fields to statistics

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::AddScaTraResults(
      RCP<DRT::Discretization> scatradis,
      RCP<Epetra_Vector> phinp
  )
  {
    if(discret_->Comm().MyPID()==0)
    {
      IO::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
               << "TurbulenceStatisticManager: added access to ScaTra results\n"
               << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << IO::endl;
    }

    // store the relevant pointers to provide access
    scatradis_   = scatradis;
    myfullphinp_ = phinp;

    if (statistics_general_mean_!=Teuchos::null)
      statistics_general_mean_->AddScaTraResults(scatradis, phinp);

    if (statistics_ccy_!=Teuchos::null)
      statistics_ccy_->AddScaTraResults(scatradis, phinp);

    if(flow_==scatra_channel_flow_of_height_2)
    {
      // store scatra discretization (multifractal subgrid scales only)
      if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
        statistics_channel_->StoreScatraDiscret(scatradis_);
    }

    withscatra_ = true;

    return;
  }


  /*----------------------------------------------------------------------

  Write (dump) the scatra-specific mean fields to the result file

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::DoOutputForScaTra(
      IO::DiscretizationWriter& output,
      int                       step)
  {
    // sampling takes place only in the sampling period
    if(step>=samstart_ && step<=samstop_ && flow_ != no_special_flow)
    {
      // statistics for scatra fields was already written during DoOutput()
      // Thus, we have to care for the mean field only:

      // dump general mean value output for scatra results
      // in combination with a restart/output
      int upres    =params_->get("write solution every", -1);
      int uprestart=params_->get("write restart every" , -1);

      if((step%upres == 0 || step%uprestart == 0) && (statistics_general_mean_!=Teuchos::null))
        statistics_general_mean_->DoOutputForScaTra(output,step);
    }
    return;
  }


  /*----------------------------------------------------------------------

  Restart statistics collection

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::Restart(
    IO::DiscretizationReader& reader,
    int                       step
    )
  {

    if(statistics_general_mean_!=Teuchos::null)
    {
      if(samstart_<step && step<=samstop_)
      {
        if(discret_->Comm().MyPID()==0)
        {
          cout << "XXXXXXXXXXXXXXXXXXXXX              ";
          cout << "Read general mean values           ";
          cout << "XXXXXXXXXXXXXXXXXXXXX";
          cout << "\n\n";
        }

        statistics_general_mean_->ReadOldStatistics(reader);
      }
    }

    return;
  } // Restart


  /*----------------------------------------------------------------------

  Restart for scatra mean fields (statistics was restarted via Restart() )

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::RestartScaTra(
    IO::DiscretizationReader& scatrareader,
    int                       step
  )
  {
    // we have only to read in the mean field.
    // The rest of the restart was already done during the Restart() call
    if(statistics_general_mean_!=Teuchos::null)
    {
      if(samstart_<step && step<=samstop_)
      {
        if(discret_->Comm().MyPID()==0)
        {
          cout << "XXXXXXXXXXXXXXXXXXXXX        ";
          cout << "Read general mean values for ScaTra      ";
          cout << "XXXXXXXXXXXXXXXXXXXXX";
          cout << "\n\n";
        }

        statistics_general_mean_->ReadOldStatisticsScaTra(scatrareader);
      }
    }

    return;
  } // RestartScaTra


} // end namespace FLD

