/*!
\file turbulence_statistic_manager.cpp

\brief Manage the computation of averages for several
canonical flows like channel flow, flow around a square
cylinder, flow in a lid driven cavity, flow over a backward-facing step etc.

The manager is intended to remove as much of the averaging
overhead as possible from the time integration method.

\level 3

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*/

#include "turbulence_statistic_manager.H"
#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluid_timint_hdg.H"
#include "../drt_fluid/fluid_xwall.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_fluid/fluid_utils.H"  // for LiftDrag
#include "../drt_lib/drt_dofset_independent_pbc.H"
#include "../drt_io/io_pstream.H"
#include "../drt_fluid_turbulence/turbulence_statistics_mean_general.H"
#include "../drt_fluid_turbulence/turbulence_statistics_ccy.H"
#include "../drt_fluid_turbulence/turbulence_statistics_cha.H"
#include "../drt_fluid_turbulence/turbulence_statistics_ldc.H"
#include "../drt_fluid_turbulence/turbulence_statistics_bfs.H"
#include "../drt_fluid_turbulence/turbulence_statistics_bfda.H"
#include "../drt_fluid_turbulence/turbulence_statistics_sqc.H"
#include "../drt_fluid_turbulence/turbulence_statistics_hit.H"
#include "../drt_fluid_turbulence/turbulence_statistics_tgv.H"
#include "../drt_fluid_turbulence/turbulence_statistics_ph.H"
#include "../drt_lib/drt_globalproblem.H"

namespace FLD
{
  /*----------------------------------------------------------------------

    Standard Constructor for standard time integration (public)

  ----------------------------------------------------------------------*/
  TurbulenceStatisticManager::TurbulenceStatisticManager(FluidImplicitTimeInt& fluid)
      : dt_(fluid.dta_),
        discret_(fluid.discret_),
        params_(fluid.params_),
        alefluid_(fluid.alefluid_),
        myaccnp_(fluid.accnp_),
        myaccn_(fluid.accn_),
        myaccam_(fluid.accam_),
        myvelnp_(fluid.velnp_),
        myveln_(fluid.veln_),
        myvelaf_(fluid.velaf_),
        myhist_(fluid.hist_),
        myscaaf_(fluid.scaaf_),
        myscaam_(fluid.scaam_),
        mydispnp_(fluid.dispnp_),
        mydispn_(fluid.dispn_),
        mygridvelaf_(fluid.gridv_),
        myforce_(fluid.trueresidual_),
        myfsvelaf_(fluid.fsvelaf_),
        myfsscaaf_(fluid.fsscaaf_),
        myxwall_(fluid.xwall_),
        mystressmanager_(fluid.stressmanager_),
        myscatrandsvel_(-1),
        flow_(no_special_flow),
        withscatra_(false),
        turbmodel_(INPAR::FLUID::no_model),
        subgrid_dissipation_(false),
        inflow_(false),
        statistics_outfilename_(fluid.statistics_outfilename_),
        statistics_general_mean_(Teuchos::null),
        statistics_channel_(Teuchos::null),
        statistics_ccy_(Teuchos::null),
        statistics_ldc_(Teuchos::null),
        statistics_bfs_(Teuchos::null),
        statistics_ph_(Teuchos::null),
        statistics_bfda_(Teuchos::null),
        statistics_sqc_(Teuchos::null),
        statistics_hit_(Teuchos::null),
        statistics_tgv_(Teuchos::null)
  {
    subgrid_dissipation_ =
        DRT::INPUT::IntegralValue<int>(params_->sublist("TURBULENCE MODEL"), "SUBGRID_DISSIPATION");
    // initialize
    withscatra_ = false;

    // toogle statistics output for turbulent inflow
    inflow_ = DRT::INPUT::IntegralValue<int>(
                  params_->sublist("TURBULENT INFLOW"), "TURBULENTINFLOW") == true;

    // toogle output of mean velocity for paraview
    out_mean_ =
        DRT::INPUT::IntegralValue<int>(params_->sublist("TURBULENCE MODEL"), "OUTMEAN") == true;

    // the flow parameter will control for which geometry the
    // sampling is done
    if (fluid.special_flow_ == "channel_flow_of_height_2")
    {
      flow_ = channel_flow_of_height_2;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_channel_ = Teuchos::rcp(new TurbulenceStatisticsCha(discret_, alefluid_, mydispnp_,
          *params_, statistics_outfilename_, subgrid_dissipation_, myxwall_));
    }
    else if (fluid.special_flow_ == "loma_channel_flow_of_height_2")
    {
      flow_ = loma_channel_flow_of_height_2;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_channel_ = Teuchos::rcp(new TurbulenceStatisticsCha(discret_, alefluid_, mydispnp_,
          *params_, statistics_outfilename_, subgrid_dissipation_, Teuchos::null));
    }
    else if (fluid.special_flow_ == "scatra_channel_flow_of_height_2")
    {
      flow_ = scatra_channel_flow_of_height_2;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_channel_ = Teuchos::rcp(new TurbulenceStatisticsCha(discret_, alefluid_, mydispnp_,
          *params_, statistics_outfilename_, subgrid_dissipation_, Teuchos::null));
    }
    else if (fluid.special_flow_ == "decaying_homogeneous_isotropic_turbulence" or
             fluid.special_flow_ == "forced_homogeneous_isotropic_turbulence" or
             fluid.special_flow_ == "scatra_forced_homogeneous_isotropic_turbulence")
    {
      if (fluid.special_flow_ == "decaying_homogeneous_isotropic_turbulence")
        flow_ = decaying_homogeneous_isotropic_turbulence;
      else if (fluid.special_flow_ == "forced_homogeneous_isotropic_turbulence")
        flow_ = forced_homogeneous_isotropic_turbulence;
      else
        flow_ = scatra_forced_homogeneous_isotropic_turbulence;

      // do the time integration independent setup
      Setup();
      if (DRT::Problem::Instance()->SpatialApproximation() == "HDG")
      {
        TimIntHDG* hdgfluid = dynamic_cast<TimIntHDG*>(&fluid);
        if (hdgfluid == NULL) dserror("this should be a hdg time integer");

        // we want to use the interior velocity here
        myvelnp_ = hdgfluid->ReturnIntVelnp();

        // allocate one instance of the averaging procedure for
        // the flow under consideration
        if (flow_ == forced_homogeneous_isotropic_turbulence or
            flow_ == scatra_forced_homogeneous_isotropic_turbulence)
          statistics_hit_ = Teuchos::rcp(
              new TurbulenceStatisticsHitHDG(discret_, *params_, statistics_outfilename_, true));
        else
        {
          statistics_hit_ = Teuchos::null;
          dserror("decaying hit currently not implemented for HDG");
        }
      }
      else
      {
        // allocate one instance of the averaging procedure for
        // the flow under consideration
        if (flow_ == forced_homogeneous_isotropic_turbulence or
            flow_ == scatra_forced_homogeneous_isotropic_turbulence)
          statistics_hit_ = Teuchos::rcp(
              new TurbulenceStatisticsHit(discret_, *params_, statistics_outfilename_, true));
        else
          statistics_hit_ = Teuchos::rcp(
              new TurbulenceStatisticsHit(discret_, *params_, statistics_outfilename_, false));
      }
    }
    else if (fluid.special_flow_ == "taylor_green_vortex")
    {
      flow_ = taylor_green_vortex;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_tgv_ =
          Teuchos::rcp(new TurbulenceStatisticsTgv(discret_, *params_, statistics_outfilename_));
    }
    else if (fluid.special_flow_ == "lid_driven_cavity")
    {
      flow_ = lid_driven_cavity;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_ldc_ =
          Teuchos::rcp(new TurbulenceStatisticsLdc(discret_, *params_, statistics_outfilename_));
    }
    else if (fluid.special_flow_ == "loma_lid_driven_cavity")
    {
      flow_ = loma_lid_driven_cavity;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_ldc_ =
          Teuchos::rcp(new TurbulenceStatisticsLdc(discret_, *params_, statistics_outfilename_));
    }
    else if (fluid.special_flow_ == "backward_facing_step")
    {
      flow_ = backward_facing_step;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_bfs_ = Teuchos::rcp(new TurbulenceStatisticsBfs(
          discret_, *params_, statistics_outfilename_, "geometry_DNS_incomp_flow"));

      // statistics manager for turbulent boundary layer not available
      if (inflow_)
        dserror(
            "The backward-facing step based on the geometry the DNS requires a turbulent boundary "
            "layer inflow profile which is not supported, yet!");
    }
    else if (fluid.special_flow_ == "periodic_hill")
    {
      flow_ = periodic_hill;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_ph_ =
          Teuchos::rcp(new TurbulenceStatisticsPh(discret_, *params_, statistics_outfilename_));
    }
    else if (fluid.special_flow_ == "loma_backward_facing_step")
    {
      flow_ = loma_backward_facing_step;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_bfs_ = Teuchos::rcp(new TurbulenceStatisticsBfs(
          discret_, *params_, statistics_outfilename_, "geometry_LES_flow_with_heating"));

      // build statistics manager for inflow channel flow
      if (inflow_)
      {
        if (params_->sublist("TURBULENT INFLOW").get<std::string>("CANONICAL_INFLOW") ==
                "channel_flow_of_height_2" or
            params_->sublist("TURBULENT INFLOW").get<std::string>("CANONICAL_INFLOW") ==
                "loma_channel_flow_of_height_2" or
            params_->sublist("TURBULENT INFLOW").get<std::string>("CANONICAL_INFLOW") ==
                "scatra_channel_flow_of_height_2")
        {
          // do not write any dissipation rates for inflow channels
          subgrid_dissipation_ = false;
          // allocate one instance of the averaging procedure for the flow under consideration
          statistics_channel_ = Teuchos::rcp(new TurbulenceStatisticsCha(discret_, alefluid_,
              mydispnp_, *params_, statistics_outfilename_, subgrid_dissipation_, myxwall_));
        }
      }
    }
    else if (fluid.special_flow_ == "backward_facing_step2")
    {
      flow_ = backward_facing_step2;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_bfs_ = Teuchos::rcp(new TurbulenceStatisticsBfs(
          discret_, *params_, statistics_outfilename_, "geometry_EXP_vogel_eaton"));

      // build statistics manager for inflow channel flow
      if (inflow_)
      {
        if (params_->sublist("TURBULENT INFLOW").get<std::string>("CANONICAL_INFLOW") ==
                "channel_flow_of_height_2" or
            params_->sublist("TURBULENT INFLOW").get<std::string>("CANONICAL_INFLOW") ==
                "loma_channel_flow_of_height_2" or
            params_->sublist("TURBULENT INFLOW").get<std::string>("CANONICAL_INFLOW") ==
                "scatra_channel_flow_of_height_2")
        {
          // do not write any dissipation rates for inflow channels
          subgrid_dissipation_ = false;
          // allocate one instance of the averaging procedure for the flow under consideration
          statistics_channel_ = Teuchos::rcp(new TurbulenceStatisticsCha(discret_, alefluid_,
              mydispnp_, *params_, statistics_outfilename_, subgrid_dissipation_, myxwall_));
        }
      }
    }
    else if (fluid.special_flow_ == "square_cylinder")
    {
      flow_ = square_cylinder;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_sqc_ =
          Teuchos::rcp(new TurbulenceStatisticsSqc(discret_, *params_, statistics_outfilename_));
    }
    else if (fluid.special_flow_ == "square_cylinder_nurbs")
    {
      flow_ = square_cylinder_nurbs;

      // do the time integration independent setup
      Setup();
    }
    else if (fluid.special_flow_ == "rotating_circular_cylinder_nurbs")
    {
      flow_ = rotating_circular_cylinder_nurbs;
      const bool withscatra = false;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_ccy_ = Teuchos::rcp(new TurbulenceStatisticsCcy(
          discret_, alefluid_, mydispnp_, *params_, statistics_outfilename_, withscatra));
    }
    else if (fluid.special_flow_ == "rotating_circular_cylinder_nurbs_scatra")
    {
      flow_ = rotating_circular_cylinder_nurbs_scatra;
      const bool withscatra = true;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_ccy_ = Teuchos::rcp(new TurbulenceStatisticsCcy(
          discret_, alefluid_, mydispnp_, *params_, statistics_outfilename_, withscatra));
    }
    else if (fluid.special_flow_ == "time_averaging")
    {
      flow_ = time_averaging;

      // do the time integration independent setup
      Setup();
    }
    else if (fluid.special_flow_ == "blood_fda_flow")
    {
      flow_ = blood_fda_flow;

      // do the time integration independent setup
      Setup();

      // allocate one instance of the averaging procedure for
      // the flow under consideration
      statistics_bfda_ =
          Teuchos::rcp(new TurbulenceStatisticsBfda(discret_, *params_, statistics_outfilename_));
    }
    else
    {
      flow_ = no_special_flow;

      // do the time integration independent setup
      Setup();
    }

    // allocate one instance of the flow independent averaging procedure
    // providing colorful output for paraview
    if (out_mean_)
    {
      Teuchos::ParameterList* modelparams = &(params_->sublist("TURBULENCE MODEL"));

      std::string homdir = modelparams->get<std::string>("HOMDIR", "not_specified");

      if (flow_ == rotating_circular_cylinder_nurbs_scatra)
      {
        // additional averaging of scalar field
        statistics_general_mean_ = Teuchos::rcp(
            new TurbulenceStatisticsGeneralMean(discret_, homdir, *fluid.VelPresSplitter(), true));
      }
      else
      {
        statistics_general_mean_ = Teuchos::rcp(
            new TurbulenceStatisticsGeneralMean(discret_, homdir, *fluid.VelPresSplitter(), false));
      }
    }
    else
      statistics_general_mean_ = Teuchos::null;

    return;
  }


  /*
    Destructor
  */
  TurbulenceStatisticManager::~TurbulenceStatisticManager() { return; }

  /*----------------------------------------------------------------------

    Time integration independent setup called by Constructor (private)

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::Setup()
  {
    Teuchos::ParameterList* modelparams = &(params_->sublist("TURBULENCE MODEL"));

    if (modelparams->get<std::string>("TURBULENCE_APPROACH", "DNS_OR_RESVMM_LES") ==
        "CLASSICAL_LES")
    {
      // check if we want to compute averages of Smagorinsky
      // constants, effective viscosities etc
      if (modelparams->get<std::string>("PHYSICAL_MODEL", "no_model") == "Dynamic_Smagorinsky" ||
          modelparams->get<std::string>("PHYSICAL_MODEL", "no_model") ==
              "Smagorinsky_with_van_Driest_damping" ||
          modelparams->get<std::string>("PHYSICAL_MODEL", "no_model") == "Smagorinsky")
      {
        turbmodel_ = INPAR::FLUID::dynamic_smagorinsky;
      }
      // check if we want to compute averages of multifractal
      // quantities (N, B)
      else if (modelparams->get<std::string>("PHYSICAL_MODEL", "no_model") ==
               "Multifractal_Subgrid_Scales")
      {
        turbmodel_ = INPAR::FLUID::multifractal_subgrid_scales;
      }
      else if (modelparams->get<std::string>("PHYSICAL_MODEL", "no_model") == "Dynamic_Vreman")
      {
        turbmodel_ = INPAR::FLUID::dynamic_vreman;
        // some dummy values into the parameter list
        params_->set<double>("C_vreman", 0.0);
        params_->set<double>("C_vreman_theoretical", 0.0);
        params_->set<double>("Dt_vreman", 0.0);
      }
    }
    else
      turbmodel_ = INPAR::FLUID::no_model;

    // parameters for sampling/dumping period
    if (flow_ != no_special_flow)
    {
      samstart_ = modelparams->get<int>("SAMPLING_START", 1);
      samstop_ = modelparams->get<int>("SAMPLING_STOP", 1000000000);
      dumperiod_ = modelparams->get<int>("DUMPING_PERIOD", 1);
    }
    else
    {
      samstart_ = 0;
      samstop_ = 0;
      dumperiod_ = 0;
    }


    if (discret_->Comm().MyPID() == 0)
    {
      if (flow_ == channel_flow_of_height_2 or flow_ == loma_channel_flow_of_height_2 or
          flow_ == scatra_channel_flow_of_height_2 or flow_ == bubbly_channel_flow)
      {
        std::string homdir = modelparams->get<std::string>("HOMDIR", "not_specified");

        if (homdir != "xy" && homdir != "xz" && homdir != "yz")
        {
          dserror("need two homogeneous directions to do averaging in plane channel flows\n");
        }

        std::cout << "Additional output          : ";
        std::cout << "Turbulence statistics are evaluated ";
        std::cout << "for a turbulent channel flow.\n";
        std::cout << "                             ";
        std::cout << "The solution is averaged over the homogeneous ";
        std::cout << homdir;
        std::cout << " plane and over time.\n";
        std::cout << "\n";
        std::cout << "                             ";
        std::cout << "Sampling period: steps " << samstart_ << " to ";
        std::cout << modelparams->get<int>("SAMPLING_STOP", 1000000000) << ".\n";

        int dumperiod = modelparams->get<int>("DUMPING_PERIOD", 1);


        if (dumperiod == 0)
        {
          std::cout << "                             ";
          std::cout << "Using standalone records (i.e. start from 0 for a new record)\n";
        }
        else
        {
          std::cout << "                             ";
          std::cout << "Volker-style incremental dumping is used (";
          std::cout << dumperiod << ")" << std::endl;
        }

        std::cout << std::endl;
        std::cout << std::endl;
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
    if (step >= samstart_ && step <= samstop_ && flow_ != no_special_flow)
    {
      switch (flow_)
      {
        case channel_flow_of_height_2:
        case loma_channel_flow_of_height_2:
        case scatra_channel_flow_of_height_2:
        {
          // add computed dynamic Smagorinsky quantities
          // (effective viscosity etc. used during the computation)
          if (turbmodel_ == INPAR::FLUID::dynamic_smagorinsky)
          {
            if (discret_->Comm().MyPID() == 0)
            {
              std::cout << "\nSmagorinsky constant, effective viscosity, ... etc, ";
              std::cout << "all element-quantities \n";
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

    Include current quantities in the time averaging procedure

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::DoTimeSample(int step, const double eosfac,
      const double thermpressaf, const double thermpressam, const double thermpressdtaf,
      const double thermpressdtam)
  {
    // store Smagorinsky statistics if used
    StoreElementValues(step);

    // sampling takes place only in the sampling period
    if ((step >= samstart_ && step <= samstop_ &&
            flow_ != no_special_flow)  // usual case with statistical-stationary state
        or (step != 0 && flow_ == decaying_homogeneous_isotropic_turbulence))  // time-dependent!
    {
      double tcpu = Teuchos::Time::wallTime();

      //--------------------------------------------------
      // calculate means, fluctuations etc of velocity,
      // pressure, boundary forces etc.
      switch (flow_)
      {
        case channel_flow_of_height_2:
        {
          if (statistics_channel_ == Teuchos::null)
            dserror("need statistics_channel_ to do a time sample for a turbulent channel flow");

          statistics_channel_->DoTimeSample(myvelnp_, myforce_);
          break;
        }
        case loma_channel_flow_of_height_2:
        {
          if (statistics_channel_ == Teuchos::null)
            dserror(
                "need statistics_channel_ to do a time sample for a turbulent channel flow at low "
                "Mach number");

          statistics_channel_->DoLomaTimeSample(myvelnp_, myscaaf_, myforce_, eosfac);
          break;
        }
        case scatra_channel_flow_of_height_2:
        {
          if (statistics_channel_ == Teuchos::null)
            dserror(
                "need statistics_channel_ to do a time sample for a turbulent passive scalar "
                "transport in channel");

          statistics_channel_->DoScatraTimeSample(myvelnp_, myscaaf_, myforce_);
          break;
        }
        case decaying_homogeneous_isotropic_turbulence:
        case forced_homogeneous_isotropic_turbulence:
        {
          if (statistics_hit_ == Teuchos::null)
            dserror("need statistics_hit_ to do sampling for homogeneous isotropic turbulence");

          statistics_hit_->DoTimeSample(myvelnp_);
          break;
        }
        case scatra_forced_homogeneous_isotropic_turbulence:
        {
          if (statistics_hit_ == Teuchos::null)
            dserror("need statistics_hit_ to do sampling for homogeneous isotropic turbulence");

          statistics_hit_->DoScatraTimeSample(myvelnp_, myphinp_);
          break;
        }
        case lid_driven_cavity:
        {
          if (statistics_ldc_ == Teuchos::null)
            dserror("need statistics_ldc_ to do a time sample for a cavity flow");

          statistics_ldc_->DoTimeSample(myvelnp_);
          break;
        }
        case loma_lid_driven_cavity:
        {
          if (statistics_ldc_ == Teuchos::null)
            dserror(
                "need statistics_ldc_ to do a time sample for a cavity flow at low Mach number");

          statistics_ldc_->DoLomaTimeSample(myvelnp_, myscaaf_, *myforce_, eosfac);
          break;
        }
        case backward_facing_step:
        case backward_facing_step2:
        {
          if (statistics_bfs_ == Teuchos::null)
            dserror(
                "need statistics_bfs_ to do a time sample for a flow over a backward-facing step");

          statistics_bfs_->DoTimeSample(myvelnp_, mystressmanager_->GetStressesWOAgg(myforce_));

          // do time sample for inflow channel flow
          if (inflow_)
          {
            if (params_->sublist("TURBULENT INFLOW").get<std::string>("CANONICAL_INFLOW") ==
                "channel_flow_of_height_2")
              statistics_channel_->DoTimeSample(myvelnp_, myforce_);
            else
              dserror("channel_flow_of_height_2 expected!");
          }
          break;
        }
        case periodic_hill:
        {
          if (statistics_ph_ == Teuchos::null)
            dserror(
                "need statistics_ph_ to do a time sample for a flow over a backward-facing step");

          statistics_ph_->DoTimeSample(
              myvelnp_, mystressmanager_->GetWallShearStressesWOAgg(myforce_));
          break;
        }
        case loma_backward_facing_step:
        {
          if (statistics_bfs_ == Teuchos::null)
            dserror(
                "need statistics_bfs_ to do a time sample for a flow over a backward-facing step "
                "at low Mach number");

          if (DRT::INPUT::get<INPAR::FLUID::PhysicalType>(*params_, "Physical Type") ==
              INPAR::FLUID::incompressible)
          {
            if (not withscatra_)
            {
              statistics_bfs_->DoTimeSample(myvelnp_, mystressmanager_->GetStressesWOAgg(myforce_));

              // do time sample for inflow channel flow
              if (inflow_)
              {
                if (params_->sublist("TURBULENT INFLOW").get<std::string>("CANONICAL_INFLOW") ==
                    "channel_flow_of_height_2")
                  statistics_channel_->DoTimeSample(myvelnp_, myforce_);
                else
                  dserror("channel_flow_of_height_2 expected!");
              }
            }
            else
            {
              statistics_bfs_->DoScatraTimeSample(myvelnp_, myscaaf_);

              // do time sample for inflow channel flow
              if (inflow_)
              {
                if (params_->sublist("TURBULENT INFLOW").get<std::string>("CANONICAL_INFLOW") ==
                    "scatra_channel_flow_of_height_2")
                  dserror(
                      "Use channel_flow_of_height_2 instead of scatra_channel_flow_of_height_2!");
                // to get DoScatraTimeSample() running also for inflow problems pointswise
                // evaluation as available for channel flow is required
                // statistics_channel_->DoScatraTimeSample(myvelnp_,myscaaf_,myforce_);
                // if the inflow channel is not subject to scalar transport using the functions for
                // usual channel flow is ok
                else if (params_->sublist("TURBULENT INFLOW")
                             .get<std::string>("CANONICAL_INFLOW") == "channel_flow_of_height_2")
                  statistics_channel_->DoTimeSample(myvelnp_, myforce_);
                else
                  dserror("channel_flow_of_height_2 expected!");
              }
            }
          }
          else
          {
            statistics_bfs_->DoLomaTimeSample(myvelnp_, myscaaf_, eosfac);

            // do time sample for inflow channel flow
            if (inflow_)
            {
              if (params_->sublist("TURBULENT INFLOW").get<std::string>("CANONICAL_INFLOW") ==
                  "loma_channel_flow_of_height_2")
              {
                std::cout << "Warning: no statistics for loma inflow channel!" << std::endl;
                // to get DoLomaTimeSample() running also for inflow problems pointswise evaluation
                // as available for channel flow is required
                // statistics_channel_->DoLomaTimeSample(myvelnp_,myscaaf_,myforce_,eosfac);
                // statistics_channel_->DoTimeSample(myvelnp_,myforce_) is not an option,
                // since it requires visc and dens
                // you may use the statistics of the inlet section instead
              }
            }
          }
          break;
        }
        case blood_fda_flow:
        {
          if (statistics_bfda_ == Teuchos::null)
            dserror("need statistics_bfda_ to do a time sample for a blood fda flow");

          statistics_bfda_->DoTimeSample(myvelnp_);
          break;
        }
        case square_cylinder:
        {
          if (statistics_sqc_ == Teuchos::null)
            dserror("need statistics_sqc_ to do a time sample for a flow around a square cylinder");

          statistics_sqc_->DoTimeSample(myvelnp_);

          // computation of Lift&Drag statistics, if required
          if (params_->get<bool>("LIFTDRAG"))
          {
            Teuchos::RCP<std::map<int, std::vector<double>>> liftdragvals;

            // spatial dimension of problem
            const int ndim = params_->get<int>("number of velocity degrees of freedom");

            FLD::UTILS::LiftDrag(discret_, myforce_, mydispnp_, ndim, liftdragvals, alefluid_);

            if ((*liftdragvals).size() != 1)
            {
              dserror(
                  "expecting only one liftdrag label for the sampling of a flow around a square "
                  "cylinder");
            }
            std::map<int, std::vector<double>>::iterator theonlyldval = (*liftdragvals).begin();

            statistics_sqc_->DoLiftDragTimeSample(
                ((*theonlyldval).second)[0], ((*theonlyldval).second)[1]);
          }
          break;
        }
        case rotating_circular_cylinder_nurbs:
        {
          if (statistics_ccy_ == Teuchos::null)
            dserror(
                "need statistics_ccy_ to do a time sample for a flow in a rotating circular "
                "cylinder");

          statistics_ccy_->DoTimeSample(myvelnp_, Teuchos::null, Teuchos::null);
          break;
        }
        case rotating_circular_cylinder_nurbs_scatra:
        {
          if (statistics_ccy_ == Teuchos::null)
            dserror(
                "need statistics_ccy_ to do a time sample for a flow in a rotating circular "
                "cylinder");

          statistics_ccy_->DoTimeSample(myvelnp_, myscaaf_, myphinp_);
          break;
        }
        default:
        {
          break;
        }
      }

      if (discret_->Comm().MyPID() == 0)
      {
        std::cout
            << "Computed statistics: mean values, fluctuations, boundary forces etc.             (";
        printf("%10.4E", Teuchos::Time::wallTime() - tcpu);
        std::cout << ")";
      }

      //--------------------------------------------------
      // do averaging of residuals, dissipation rates etc
      // (all gausspoint-quantities)
      if (subgrid_dissipation_)
      {
        tcpu = Teuchos::Time::wallTime();

        switch (flow_)
        {
          case channel_flow_of_height_2:
          case loma_channel_flow_of_height_2:
          case scatra_channel_flow_of_height_2:
          case taylor_green_vortex:
          {
            if (statistics_channel_ == Teuchos::null and statistics_tgv_ == Teuchos::null)
              dserror("No dissipation rates for this flow type!");

            // set vector values needed by elements
            std::map<std::string, Teuchos::RCP<Epetra_Vector>> statevecs;
            std::map<std::string, Teuchos::RCP<Epetra_MultiVector>> statetenss;
            std::map<std::string, Teuchos::RCP<Epetra_Vector>> scatrastatevecs;
            std::map<std::string, Teuchos::RCP<Epetra_MultiVector>> scatrafieldvecs;

            statevecs.insert(std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("hist", myhist_));
            statevecs.insert(
                std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("accam", myaccam_));
            statevecs.insert(
                std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("scaaf", myscaaf_));
            statevecs.insert(
                std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("scaam", myscaam_));

            if (alefluid_)
            {
              statevecs.insert(
                  std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("dispnp", mydispnp_));
              statevecs.insert(
                  std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("gridv", mygridvelaf_));
              if (scatradis_ != Teuchos::null) dserror("Not supported!");
            }

            if (DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo") ==
                    INPAR::FLUID::timeint_afgenalpha or
                DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo") ==
                    INPAR::FLUID::timeint_npgenalpha)
            {
              statevecs.insert(
                  std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("velaf", myvelaf_));
              if (DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo") ==
                  INPAR::FLUID::timeint_npgenalpha)
                statevecs.insert(
                    std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("velnp", myvelnp_));

              // additional scatra vectors
              scatrastatevecs.insert(
                  std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("phinp", myphiaf_));
              scatrastatevecs.insert(
                  std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("phiam", myphiam_));
              scatrastatevecs.insert(
                  std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("hist", myphidtam_));
            }
            else
            {
              statevecs.insert(
                  std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("velaf", myvelnp_));
              statevecs.insert(
                  std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("velaf", myvelnp_));

              // additional scatra vectors
              scatrastatevecs.insert(
                  std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("phinp", myphinp_));
              scatrastatevecs.insert(
                  std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("hist", myscatrahist_));
            }

            if (params_->sublist("TURBULENCE MODEL").get<std::string>("FSSUGRVISC") != "No" or
                turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
            {
              statevecs.insert(
                  std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("fsvelaf", myfsvelaf_));
              if (myfsvelaf_ == Teuchos::null) dserror("Have not got fsvel!");

              if (DRT::INPUT::get<INPAR::FLUID::PhysicalType>(*params_, "Physical Type") ==
                      INPAR::FLUID::loma and
                  turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
              {
                statevecs.insert(
                    std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("fsscaaf", myfsscaaf_));
                if (myfsscaaf_ == Teuchos::null) dserror("Have not got fssca!");
              }

              // additional scatra vectors
              if (withscatra_)
              {
                scatrastatevecs.insert(
                    std::pair<std::string, Teuchos::RCP<Epetra_Vector>>("fsphinp", myfsphi_));
                if (myfsphi_ == Teuchos::null) dserror("Have not got fsphi!");
              }
            }

            switch (flow_)
            {
              case channel_flow_of_height_2:
              case loma_channel_flow_of_height_2:
              case scatra_channel_flow_of_height_2:
              {
                statistics_channel_->EvaluateResiduals(statevecs, statetenss, thermpressaf,
                    thermpressam, thermpressdtaf, thermpressdtam, scatrastatevecs, scatrafieldvecs,
                    myscatrandsvel_);
                break;
              }
              case taylor_green_vortex:
              {
                statistics_tgv_->EvaluateResiduals(statevecs, statetenss, thermpressaf,
                    thermpressam, thermpressdtaf, thermpressdtam);
                if (step == 0)  // sorry, this function is not called for step=0
                {
                  statistics_tgv_->DumpStatistics(0);
                  statistics_tgv_->ClearStatistics();
                }
                break;
              }
              default:
                dserror("Dissipation not supported for this flow type!");
                break;
            }

            break;
          }
          default:
          {
            break;
          }
        }

        if (discret_->Comm().MyPID() == 0)
        {
          std::cout << "\nresiduals, dissipation rates etc, ";
          std::cout << "all gausspoint-quantities (";
          printf("%10.4E", Teuchos::Time::wallTime() - tcpu);
          std::cout << ")";
        }
      }
      if (discret_->Comm().MyPID() == 0)
      {
        std::cout << "\n";
      }

      if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales and inflow_ == false and
          myxwall_ == Teuchos::null)
      {
        switch (flow_)
        {
          case channel_flow_of_height_2:
          {
            // add parameters of multifractal subgrid-scales model
            if (DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo") ==
                    INPAR::FLUID::timeint_afgenalpha or
                DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo") ==
                    INPAR::FLUID::timeint_npgenalpha)
              statistics_channel_->AddModelParamsMultifractal(myvelaf_, myfsvelaf_, false);
            else
              statistics_channel_->AddModelParamsMultifractal(myvelnp_, myfsvelaf_, false);
            break;
          }
          case scatra_channel_flow_of_height_2:
          {
            // add parameters of multifractal subgrid-scales model
            if (DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo") ==
                    INPAR::FLUID::timeint_afgenalpha or
                DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo") ==
                    INPAR::FLUID::timeint_npgenalpha)
              statistics_channel_->AddModelParamsMultifractal(myvelaf_, myfsvelaf_, true);
            else
              statistics_channel_->AddModelParamsMultifractal(myvelnp_, myfsvelaf_, false);
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
      if (statistics_general_mean_ != Teuchos::null)
        statistics_general_mean_->AddToCurrentTimeAverage(dt_, myvelnp_, myscaaf_, myphinp_);

    }  // end step in sampling period

    // for homogeneous isotropic turbulence, the initial field is averaged
    // to get the amount of energy at the beginning which depends on the
    // resolution
    if (flow_ == decaying_homogeneous_isotropic_turbulence and step == 0)
    {
      if (discret_->Comm().MyPID() == 0)
      {
        std::cout << "XXXXXXXXXXXXXXXXXXXXX              ";
        std::cout << "calculate initial energy spectrum  ";
        std::cout << "XXXXXXXXXXXXXXXXXXXXX";
        std::cout << "\n\n";
      }

      statistics_hit_->DoTimeSample(myvelnp_);
      statistics_hit_->DumpStatistics(0);
      statistics_hit_->ClearStatistics();

      if (discret_->Comm().MyPID() == 0)
      {
        std::cout << "XXXXXXXXXXXXXXXXXXXXX              ";
        std::cout << "wrote statistics record            ";
        std::cout << "XXXXXXXXXXXXXXXXXXXXX";
        std::cout << "\n\n";
      }
    }

    return;
  }


  /*----------------------------------------------------------------------

    Include current quantities in the time averaging procedure

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::DoTimeSample(int step, Teuchos::RCP<Epetra_Vector> velnp,
      Teuchos::RCP<Epetra_Vector> force, Teuchos::RCP<Epetra_Vector> phi,
      Teuchos::RCP<const DRT::DofSet> stddofset)
  {
    // sampling takes place only in the sampling period
    if (step >= samstart_ && step <= samstop_ && flow_ != no_special_flow)
    {
      double tcpu = Teuchos::Time::wallTime();

      //--------------------------------------------------
      // calculate means, fluctuations etc of velocity,
      // pressure, boundary forces etc.
      switch (flow_)
      {
        default:
        {
          dserror("called wrong DoTimeSample() for this kind of special flow");
          break;
        }
      }

      // add vector(s) to general mean value computation
      // scatra vectors may be Teuchos::null
      if (statistics_general_mean_ != Teuchos::null)
        statistics_general_mean_->AddToCurrentTimeAverage(dt_, velnp, myscaaf_, myphinp_);

      if (discret_->Comm().MyPID() == 0)
      {
        std::cout << "                      taking time sample (";
        printf("%10.4E", Teuchos::Time::wallTime() - tcpu);
        std::cout << ")\n";
      }

    }  // end step in sampling period

    return;
  }


  /*----------------------------------------------------------------------

    get current velnp pointer from fluid
    necessary for meshtying                                    bk 02/14
  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::GetCurrentVelnp(Teuchos::RCP<Epetra_Vector> velnp)
  {
    myvelnp_ = velnp;
    return;
  }

  /*----------------------------------------------------------------------

    Write (dump) the statistics to a file

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::DoOutput(
      IO::DiscretizationWriter& output, int step, const bool inflow)
  {
    // sampling takes place only in the sampling period
    if (step >= samstart_ && step <= samstop_ && flow_ != no_special_flow)
    {
      enum format
      {
        write_single_record,
        write_multiple_records,
        do_not_write
      } outputformat = do_not_write;
      bool output_inflow = false;

      // sampling a la Volker --- single record is constantly updated
      if (dumperiod_ != 0)
      {
        int samstep = step - samstart_ + 1;

        // dump every dumperiod steps
        if (samstep % dumperiod_ == 0) outputformat = write_single_record;
      }

      // sampling a la Peter --- for each sampling period a
      // new record is written; they can be combined by a
      // postprocessing script to a single long term sample
      // (allows restarts during sampling)
      if (dumperiod_ == 0)
      {
        int upres = params_->get<int>("write solution every");
        int uprestart = params_->get<int>("write restart every");

        // dump in combination with a restart/output
        if ((step % upres == 0 || (uprestart > 0 && step % uprestart == 0)) && step > samstart_)
          outputformat = write_multiple_records;
      }
      if (inflow_)
      {
        int upres = params_->get<int>("write solution every");
        int uprestart = params_->get<int>("write restart every");

        // dump in combination with a restart/output
        if ((step % upres == 0 || (uprestart > 0 && step % uprestart == 0)) && step > samstart_)
          output_inflow = true;
      }

      if (discret_->Comm().MyPID() == 0 && outputformat != do_not_write)
        std::cout << "---  statistics record: \n" << std::flush;

      // do actual output (time averaging)
      switch (flow_)
      {
        case channel_flow_of_height_2:
        {
          if (statistics_channel_ == Teuchos::null)
            dserror("need statistics_channel_ to do a time sample for a turbulent channel flow");


          if (outputformat == write_multiple_records)
          {
            statistics_channel_->TimeAverageMeansAndOutputOfStatistics(step);
            statistics_channel_->ClearStatistics();
          }

          if (outputformat == write_single_record) statistics_channel_->DumpStatistics(step);
          break;
        }
        case loma_channel_flow_of_height_2:
        {
          if (statistics_channel_ == Teuchos::null)
            dserror(
                "need statistics_channel_ to do a time sample for a turbulent channel flow at low "
                "Mach number");

          if (outputformat == write_single_record) statistics_channel_->DumpLomaStatistics(step);
          break;
        }
        case scatra_channel_flow_of_height_2:
        {
          if (statistics_channel_ == Teuchos::null)
            dserror(
                "need statistics_channel_ to do a time sample for a turbulent channel flow at low "
                "Mach number");

          if (outputformat == write_single_record) statistics_channel_->DumpScatraStatistics(step);
          break;
        }
        case decaying_homogeneous_isotropic_turbulence:
        case forced_homogeneous_isotropic_turbulence:
        {
          if (statistics_hit_ == Teuchos::null)
            dserror("need statistics_hit_ to do sampling for homogeneous isotropic turbulence");

          if (flow_ == forced_homogeneous_isotropic_turbulence)
          {
            // write statistics only during sampling period
            if (outputformat == write_single_record)
              statistics_hit_->DumpStatistics(step);
            else if (outputformat == write_multiple_records)
            {
              statistics_hit_->DumpStatistics(step, true);
              statistics_hit_->ClearStatistics();
            }
          }
          else
          {
            // write statistics for every time step,
            // since there is not any statistical-stationary state
            statistics_hit_->DumpStatistics(step);
            statistics_hit_->ClearStatistics();
          }
          break;
        }
        case scatra_forced_homogeneous_isotropic_turbulence:
        {
          if (statistics_hit_ == Teuchos::null)
            dserror("need statistics_hit_ to do sampling for homogeneous isotropic turbulence");

          // write statistics only during sampling period
          if (outputformat == write_single_record)
            statistics_hit_->DumpScatraStatistics(step);
          else if (outputformat == write_multiple_records)
          {
            statistics_hit_->DumpScatraStatistics(step, true);
            statistics_hit_->ClearScatraStatistics();
          }
          break;
        }
        case taylor_green_vortex:
        {
          // write statistics for every time step,
          // since there is not any statistical-stationary state
          statistics_tgv_->DumpStatistics(step);
          statistics_tgv_->ClearStatistics();
          break;
        }
        case lid_driven_cavity:
        {
          if (statistics_ldc_ == Teuchos::null)
            dserror("need statistics_ldc_ to do a time sample for a lid driven cavity");

          if (outputformat == write_single_record) statistics_ldc_->DumpStatistics(step);

          if (outputformat != do_not_write) statistics_ldc_->WriteRestart(output);
          break;
        }
        case loma_lid_driven_cavity:
        {
          if (statistics_ldc_ == Teuchos::null)
            dserror(
                "need statistics_ldc_ to do a time sample for a lid driven cavity at low Mach "
                "number");

          if (outputformat == write_single_record) statistics_ldc_->DumpLomaStatistics(step);
          break;
        }
        case backward_facing_step:
        case backward_facing_step2:
        {
          if (statistics_bfs_ == Teuchos::null)
            dserror(
                "need statistics_bfs_ to do a time sample for a flow over a backward-facing step");

          if (outputformat == write_single_record) statistics_bfs_->DumpStatistics(step);

          // write statistics of inflow channel flow
          if (inflow_)
          {
            if (params_->sublist("TURBULENT INFLOW").get<std::string>("CANONICAL_INFLOW") ==
                "channel_flow_of_height_2")
            {
              if (output_inflow)
              {
                statistics_channel_->TimeAverageMeansAndOutputOfStatistics(step);
                statistics_channel_->ClearStatistics();
              }
            }
          }
          break;
        }
        case periodic_hill:
        {
          if (statistics_ph_ == Teuchos::null)
            dserror("need statistics_ph_ to do a time sample for a flow over a periodic hill");

          if (outputformat == write_single_record) statistics_ph_->DumpStatistics(step);

          break;
        }
        case loma_backward_facing_step:
        {
          if (statistics_bfs_ == Teuchos::null)
            dserror(
                "need statistics_bfs_ to do a time sample for a flow over a backward-facing step "
                "at low Mach number");

          if (DRT::INPUT::get<INPAR::FLUID::PhysicalType>(*params_, "Physical Type") ==
              INPAR::FLUID::incompressible)
          {
            if (not withscatra_)
            {
              if (outputformat == write_single_record) statistics_bfs_->DumpStatistics(step);

              // write statistics of inflow channel flow
              if (inflow_)
              {
                if (params_->sublist("TURBULENT INFLOW").get<std::string>("CANONICAL_INFLOW") ==
                    "channel_flow_of_height_2")
                {
                  if (output_inflow)
                  {
                    statistics_channel_->TimeAverageMeansAndOutputOfStatistics(step);
                    statistics_channel_->ClearStatistics();
                  }
                }
              }
            }
            else
            {
              if (outputformat == write_single_record) statistics_bfs_->DumpScatraStatistics(step);

              // write statistics of inflow channel flow
              if (inflow_)
              {
                if (params_->sublist("TURBULENT INFLOW").get<std::string>("CANONICAL_INFLOW") ==
                    "scatra_channel_flow_of_height_2")
                {
                  if (outputformat == write_single_record)
                    statistics_channel_->DumpScatraStatistics(step);
                }
                else if (params_->sublist("TURBULENT INFLOW")
                             .get<std::string>("CANONICAL_INFLOW") == "channel_flow_of_height_2")
                {
                  if (output_inflow)
                  {
                    statistics_channel_->TimeAverageMeansAndOutputOfStatistics(step);
                    statistics_channel_->ClearStatistics();
                  }
                }
              }
            }
          }
          else  // loma
          {
            if (outputformat == write_single_record) statistics_bfs_->DumpLomaStatistics(step);

            // write statistics of inflow channel flow
            if (inflow_)
            {
              if (params_->sublist("TURBULENT INFLOW").get<std::string>("CANONICAL_INFLOW") ==
                  "loma_channel_flow_of_height_2")
              {
                std::cout << "Warning: no statistics for loma inflow channel!" << std::endl;
                //              if(outputformat == write_single_record)
                //                statistics_channel_->DumpLomaStatistics(step);
              }
            }
          }
          break;
        }
        case blood_fda_flow:
        {
          if (statistics_bfda_ == Teuchos::null)
            dserror("need statistics_bfda_ to do a time sample for a blood fda flow");

          if (outputformat == write_single_record) statistics_bfda_->DumpStatistics(step);

          break;
        }
        case square_cylinder:
        {
          if (statistics_sqc_ == Teuchos::null)
            dserror("need statistics_sqc_ to do a time sample for a square cylinder flow");

          if (outputformat == write_single_record) statistics_sqc_->DumpStatistics(step);
          break;
        }
        case rotating_circular_cylinder_nurbs:
        case rotating_circular_cylinder_nurbs_scatra:
        {
          if (statistics_ccy_ == Teuchos::null)
            dserror(
                "need statistics_ccy_ to do a time sample for a flow in a rotating circular "
                "cylinder");

          statistics_ccy_->TimeAverageMeansAndOutputOfStatistics(step);
          statistics_ccy_->ClearStatistics();
          break;
        }
        default:
        {
          break;
        }
      }

      if (discret_->Comm().MyPID() == 0 && outputformat != do_not_write)
      {
        std::cout << "XXXXXXXXXXXXXXXXXXXXX              ";
        std::cout << "wrote statistics record            ";
        std::cout << "XXXXXXXXXXXXXXXXXXXXX";
        std::cout << "\n\n";
      }


      // dump general mean value output in combination with a restart/output
      // don't write output if turbulent inflow or twophaseflow is computed
      if (!inflow and flow_ != bubbly_channel_flow)
      {
        int upres = params_->get<int>("write solution every");
        int uprestart = params_->get<int>("write restart every");

        if ((step % upres == 0 || (uprestart > 0 && step % uprestart == 0)) &&
            (statistics_general_mean_ != Teuchos::null))
          statistics_general_mean_->WriteOldAverageVec(output);
      }
    }  // end step is in sampling period

    if (discret_->Comm().MyPID() == 0 and turbmodel_ == INPAR::FLUID::dynamic_vreman)
    {
      std::string fnamevreman(statistics_outfilename_);

      fnamevreman.append(".vremanconstant");
      double Cv;
      double Cv_theo;
      double Dt = 0.0;
      std::cout << __LINE__ << std::endl;
      Cv = params_->get<double>("C_vreman", 0.0);
      Cv_theo = params_->get<double>("C_vreman_theoretical");
      if (withscatra_) Dt = params_->get<double>("Dt_vreman");
      std::ofstream fvreman;
      if (step <= 1)
      {
        fvreman.open(fnamevreman.c_str(), std::ofstream::trunc);
        fvreman << "time step C_v         C_v clipped D_T (required for subgrid diffusivity)\n";
      }
      else
        fvreman.open(fnamevreman.c_str(), std::ofstream::app);
      fvreman.width(10);
      fvreman << step;
      fvreman.width(12);
      fvreman << Cv_theo;
      fvreman.width(12);
      fvreman << Cv;
      if (withscatra_)
      {
        fvreman.width(12);
        fvreman << Dt * Cv;  // to make it comparable to the the original paper
      }
      fvreman << "\n";
      fvreman.flush();
      fvreman.close();
    }

    return;
  }  // DoOutput


  /*----------------------------------------------------------------------

  Provide access to scalar transport field

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::AddScaTraField(
      Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_timeint)
  {
    if (discret_->Comm().MyPID() == 0)
    {
      IO::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
               << IO::endl;
      IO::cout << "TurbulenceStatisticManager: provided access to ScaTra time integration"
               << IO::endl;
      IO::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
               << IO::endl;
    }

    // store the relevant pointers to provide access
    scatradis_ = scatra_timeint->Discretization();
    scatraparams_ = scatra_timeint->ScatraParameterList();
    // and sublists from extraparams
    scatraextraparams_ = scatra_timeint->ScatraExtraParameterList();
    // remark: this is not a good idea, since the sublists are copied and modifications
    //         during simulation are never seen in the statistics manager; moreover
    //         this would add the following sublists also to the parameter list in the
    //         scatra time integration
    // scatraparams_->sublist("TURBULENCE MODEL") =
    // scatra_timeint->ScatraExtraParameterList()->sublist("TURBULENCE MODEL");
    // scatraparams_->sublist("SUBGRID VISCOSITY") =
    // scatra_timeint->ScatraExtraParameterList()->sublist("SUBGRID VISCOSITY");
    // scatraparams_->sublist("MULTIFRACTAL SUBGRID SCALES") =
    // scatra_timeint->ScatraExtraParameterList()->sublist("MULTIFRACTAL SUBGRID SCALES");
    // scatraparams_->sublist("LOMA") = scatra_timeint->ScatraExtraParameterList()->sublist("LOMA");
    scatratimeparams_ = scatra_timeint->ScatraTimeParameterList();
    // required vectors
    // remark: Although some of these field are already set for the fluid,
    //         we set set them here once more. They are required for integration
    //         on the element level of the scatra field and have to be set in the
    //         specific form using MultiVectors (cf. scatra time integration). If
    //         these vectors are not taken form the scatra time integration we
    //         would have to transfer them to the scatra dofs here!
    myphinp_ = scatra_timeint->Phinp();
    myphiaf_ = scatra_timeint->Phiaf();
    myphiam_ = scatra_timeint->Phiam();
    myscatrandsvel_ = scatra_timeint->NdsVel();
    myscatrahist_ = scatra_timeint->Hist();
    myphidtam_ = scatra_timeint->Phidtam();
    myfsphi_ = scatra_timeint->FsPhi();

    if (statistics_general_mean_ != Teuchos::null)
      statistics_general_mean_->AddScaTraResults(scatradis_, myphinp_);

    if (statistics_ccy_ != Teuchos::null) statistics_ccy_->AddScaTraResults(scatradis_, myphinp_);

    if (flow_ == scatra_channel_flow_of_height_2 or flow_ == loma_channel_flow_of_height_2)
    {
      statistics_channel_->StoreScatraDiscretAndParams(
          scatradis_, scatraparams_, scatraextraparams_, scatratimeparams_);
    }

    if (flow_ == scatra_forced_homogeneous_isotropic_turbulence)
      statistics_hit_->StoreScatraDiscret(scatradis_);

    withscatra_ = true;

    return;
  }


  /*----------------------------------------------------------------------

  Write (dump) the scatra-specific mean fields to the result file

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::DoOutputForScaTra(IO::DiscretizationWriter& output, int step)
  {
    // sampling takes place only in the sampling period
    if (step >= samstart_ && step <= samstop_ && flow_ != no_special_flow)
    {
      // statistics for scatra fields was already written during DoOutput()
      // Thus, we have to care for the mean field only:

      // dump general mean value output for scatra results
      // in combination with a restart/output
      int upres = params_->get("write solution every", -1);
      int uprestart = params_->get("write restart every", -1);

      if ((step % upres == 0 || step % uprestart == 0) &&
          (statistics_general_mean_ != Teuchos::null))
        statistics_general_mean_->DoOutputForScaTra(output, step);
    }
    return;
  }


  /*----------------------------------------------------------------------

  Restart statistics collection

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::ReadRestart(IO::DiscretizationReader& reader, int step)
  {
    if (samstart_ < step && step <= samstop_)
    {
      if (statistics_general_mean_ != Teuchos::null)
      {
        if (discret_->Comm().MyPID() == 0)
        {
          std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXX          ";
          std::cout << "Read general mean values           ";
          std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXX";
          std::cout << "\n\n";
        }

        statistics_general_mean_->ReadOldStatistics(reader);
      }

      if (statistics_ldc_ != Teuchos::null)
      {
        if (discret_->Comm().MyPID() == 0)
        {
          std::cout << "XXXXXXXXXXXXXXXXXXXXX              ";
          std::cout << "Read ldc statistics                ";
          std::cout << "XXXXXXXXXXXXXXXXXXXXX";
          std::cout << "\n\n";
        }

        statistics_ldc_->ReadRestart(reader);
      }
    }

    return;
  }  // Restart


  /*----------------------------------------------------------------------

  Restart for scatra mean fields (statistics was restarted via Restart() )

  ----------------------------------------------------------------------*/
  void TurbulenceStatisticManager::ReadRestartScaTra(
      IO::DiscretizationReader& scatrareader, int step)
  {
    // we have only to read in the mean field.
    // The rest of the restart was already done during the Restart() call
    if (statistics_general_mean_ != Teuchos::null)
    {
      if (samstart_ < step && step <= samstop_)
      {
        if (discret_->Comm().MyPID() == 0)
        {
          std::cout << "XXXXXXXXXXXXXXXXXXXXX        ";
          std::cout << "Read general mean values for ScaTra      ";
          std::cout << "XXXXXXXXXXXXXXXXXXXXX";
          std::cout << "\n\n";
        }

        statistics_general_mean_->ReadOldStatisticsScaTra(scatrareader);
      }
    }

    return;
  }  // RestartScaTra


}  // end namespace FLD
