/*!----------------------------------------------------------------------
\file fluid_projectionmethod.cpp
\brief Control routine for fluid solver with pressure correction,

<pre>
Created on: Jun 10, 2009
Maintainer: Tobias Wiesner
wiesner@lnm.mw.tum.de
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "fluid_projectionmethod.H"
#include "time_integration_scheme.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/linalg_ana.H"
#include "../drt_lib/drt_nodematchingoctree.H"
#include "../drt_lib/drt_function.H"
#include "fluid_utils.H"
#include <Epetra_Vector.h>

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                                |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::FluidProjectionMethod::FluidProjectionMethod(RefCountPtr<DRT::Discretization> actdis,
        LINALG::Solver& solver,
        LINALG::Solver& psolver,
        ParameterList& params,
        IO::DiscretizationWriter& output,
        bool alefluid)
:
    discret_(actdis),
    solver_(solver),
    psolver_(psolver),
    params_ (params),
    output_(output),
    alefluid_(alefluid),
    time_(0.0),
    step_(0),
    extrapolationpredictor_(params.get("do explicit predictor",true)),
    uprestart_(params.get("write restart every", -1)),
    upres_(params.get("write solution every", -1)),
    surfacesplitter_(NULL)
    {
    // get processor ID from the communicator
    myrank_ = discret_->Comm().MyPID();

    // get basic parameters first

    // type of time-integration
    timealgo_ = params_.get<FLUID_TIMEINTTYPE>("time int algo");
    // time-step size
    dtp_ = dta_ = params_.get<double>("time step size");
    // maximum number of timesteps
    stepmax_ = params_.get<int> ("max number timesteps");
    // maximum simulation time
    maxtime_ = params_.get<double>("total time");
    // parameter theta for time-integration schemes
    theta_ = params_.get<double>("theta");
    // af-generalized-alpha parameters: gamma_ = 0.5 + alphaM_ - alphaF_
    // (may be reset below when starting algorithm is used)
    alphaM_   = params_.get<double>("alpha_M");
    alphaF_   = params_.get<double>("alpha_F");
    gamma_    = params_.get<double>("gamma");

    // number of steps for starting algorithm
    numstasteps_ = params_.get<int> ("number of start steps");
    // starting algorithm only for af-generalized-alpha so far
    // -> check for time-integration scheme and reasonability of number of steps
    startalgo_ = false;
    if (numstasteps_ > 0)
    {
        if (timealgo_ != timeint_afgenalpha)
            dserror("no starting algorithm supported for schemes other than af-gen-alpha");
        else startalgo_= true;
        if (numstasteps_>stepmax_)
            dserror("more steps for starting algorithm than steps overall");
    }

    // parameter for linearization scheme (fixed-point-like or Newton)
    newton_ = params_.get<string>("Linearisation");

    // ensure that degrees of freedom in the discretization have been set
    if(!discret_->Filled() || !actdis->HaveDofs()) discret_->FillComplete();

    // get vector layout velocity dofs <-> pressure dofs
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    const int numdim = params_.get<int>("number of velocity degrees of freedom");
    FLD::UTILS::SetupFluidSplit(*discret_,numdim,velpressplitter_);

    // matrices for pressure correction method
    pressmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*PressureRowMap(),24,false,true));
    massmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*VelocityRowMap(),24,false,true));
    gradop_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap));
    gradopwithoutbc_= Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap));

    // Vectors for projection method
    veltilde_ = LINALG::CreateVector(*VelocityRowMap(),true);
    lmassinvvec_ = LINALG::CreateVector(*VelocityRowMap(),true);
    phi_ = LINALG::CreateVector(*PressureRowMap(),true);

    // Vectors passed to the element
    // -----------------------------
    // velocity/pressure at time n+1
    velnp_        = LINALG::CreateVector(*dofrowmap,true);

    // velocity/pressure at time n and n-1
    veln_         = LINALG::CreateVector(*dofrowmap,true);
    velnm_        = LINALG::CreateVector(*dofrowmap,true);

    // acceleration at time n+1, n and n-1
    accnp_        = LINALG::CreateVector(*dofrowmap,true);
    accn_         = LINALG::CreateVector(*dofrowmap,true);
    accnm_        = LINALG::CreateVector(*dofrowmap,true);

    // vectors only required for af-generalized-alpha scheme
    if (timealgo_==timeint_afgenalpha)
    {
        // velocity/pressure at time n+alpha_F
        velaf_        = LINALG::CreateVector(*dofrowmap,true);

        // acceleration/(density time derivative) at time n+alpha_M
        accam_        = LINALG::CreateVector(*dofrowmap,true);

        // velocity/density at time n+alpha_M and n+alpha_F
        //vedeam_       = LINALG::CreateVector(*dofrowmap,true);
        //vedeaf_       = LINALG::CreateVector(*dofrowmap,true);
    }

    // velocity/density at time n+1
    vedenp_       = LINALG::CreateVector(*dofrowmap,true);

    // history vector
    hist_           = LINALG::CreateVector(*dofrowmap,true);

    if (alefluid_)
    {
        dispnp_       = LINALG::CreateVector(*dofrowmap,true);
        dispn_        = LINALG::CreateVector(*dofrowmap,true);
        dispnm_       = LINALG::CreateVector(*dofrowmap,true);
        gridv_        = LINALG::CreateVector(*dofrowmap,true);
    }

    // Vectors associated to boundary conditions
    // -----------------------------------------

    // a vector of zeros to be used to enforce zero dirichlet boundary conditions
    zeros_   = LINALG::CreateVector(*dofrowmap,true);

    // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
    dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
    {
        ParameterList eleparams;
        // other parameters needed by the elements
        eleparams.set("total time",time_);
        discret_->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null,
                Teuchos::null, dbcmaps_);
        zeros_->PutScalar(0.0); // just in case of change
    }

    // the vector containing body and surface forces
    neumann_loads_= LINALG::CreateVector(*dofrowmap,true);

    // residual vector (only needed for partitioned FSI)
    trueresidual_ = LINALG::CreateVector(*dofrowmap,true);

    // get constant density variable for incompressible flow
    // set vedenp-vector values to 1.0 for incompressible flow
    // set density variable to 1.0 for low-Mach-number flow
    {
        ParameterList eleparams;
        eleparams.set("action","get_density");
        discret_->Evaluate(eleparams,null,null,null,null,null);
        density_ = eleparams.get("density", 1.0);
        if (density_ <= 0.0) dserror("received illegal density value");
        vedenp_->PutScalar(1.0);
    }
    }

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
| Start the time integration. Allows                                   |
|                                                                      |
|  o starting steps with different algorithms                          |
|  o the "standard" time integration                                   |
|                                                           gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidProjectionMethod::Integrate()
{
    // output of stabilization details
    if (myrank_==0)
    {
        cout << "Solving strategy         : Projection method (pressure correction)\n";
        cout << "\n";
    }

    // distinguish stationary and instationary case
    if (timealgo_==timeint_stationary) dserror("no stationary problems supported");
    else TimeLoop();

    // print the results of time measurements
    TimeMonitor::summarize();

    return;
} // FluidProjectionMethod::Integrate


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
| contains the time loop                                    gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidProjectionMethod::TimeLoop()
{
    // time measurement: time loop
    TEUCHOS_FUNC_TIME_MONITOR(" + time loop");

    // how do we want to solve or fluid equations?
    const int dyntype = params_.get<int>("type of nonlinear solve");


    while (step_<stepmax_ and time_<maxtime_)
    {
        PrepareTimeStep();
        // -------------------------------------------------------------------
        //                       output to screen
        // -------------------------------------------------------------------
        if (myrank_==0)
        {
            switch (timealgo_)
            {
            case timeint_one_step_theta:
                printf("TIME: %11.4E/%11.4E  DT = %11.4E   One-Step-Theta    STEP = %4d/%4d \n",
                        time_,maxtime_,dta_,step_,stepmax_);
                break;
            case timeint_afgenalpha:
                printf("TIME: %11.4E/%11.4E  DT = %11.4E  Generalized-Alpha  STEP = %4d/%4d \n",
                        time_,maxtime_,dta_,step_,stepmax_);
                break;
            case timeint_bdf2:
                printf("TIME: %11.4E/%11.4E  DT = %11.4E       BDF2          STEP = %4d/%4d \n",
                        time_,maxtime_,dta_,step_,stepmax_);
                break;
            default:
                dserror("parameter out of range: IOP\n");
            } /* end of switch(timealgo) */
        }

        switch (dyntype)
        {
        case 0:
            // -----------------------------------------------------------------
            //                     solve nonlinear equation
            // -----------------------------------------------------------------
            ProjectionSolve();
            break;
        default:
            dserror("Type of dynamics unknown!!");
        }

        // -------------------------------------------------------------------
        //                         update solution
        //        current solution becomes old solution of next timestep
        // -------------------------------------------------------------------
        TimeUpdate();

        // -------------------------------------------------------------------
        //                         output of solution
        // -------------------------------------------------------------------
        Output();

        // -------------------------------------------------------------------
        // evaluate error for test flows with analytical solutions
        // -------------------------------------------------------------------
        //EvaluateErrorComparedToAnalyticalSol(); TODO

        // -------------------------------------------------------------------
        //                       update time step sizes
        // -------------------------------------------------------------------
        dtp_ = dta_;

        // -------------------------------------------------------------------
        //                    stop criterium for timeloop
        // -------------------------------------------------------------------
    }
} // FluidProjectionMethod::TimeLoop


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | setup the variables to do a new time step                 u.kue 06/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidProjectionMethod::PrepareTimeStep()
{
    // -------------------------------------------------------------------
    //              set time dependent parameters
    // -------------------------------------------------------------------
    step_ += 1;
    time_ += dta_;

    // for BDF2, theta is set by the time-step sizes, 2/3 for const. dt
    if (timealgo_==timeint_bdf2) theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);

    // NOTE: we don't want any predictor or RHS for pressure part
    // we bypass it by storing the pressure part in a temporary variable
    // and write back later
    Teuchos::RCP<Epetra_Vector> pretemp = LINALG::CreateVector(*PressureRowMap(),true);
    velpressplitter_.ExtractCondVector(velnp_,pretemp);


    // -------------------------------------------------------------------
    // set part(s) of the rhs vector(s) belonging to the old timestep
    // for low-Mach-number flow: distinguish momentum and continuity part
    // (continuity part only meaningful for low-Mach-number flow)
    // For low-Mach-number flow, all velocity values are multiplied by the
    // respective density values.
    //
    // Stationary/af-generalized-alpha:
    //
    //               mom: hist_ = 0.0
    //              (con: hist_ = 0.0)
    //
    // One-step-Theta:
    //
    //               mom: hist_ = veln_  + dt*(1-Theta)*accn_
    //              (con: hist_ = densn_ + dt*(1-Theta)*densdtn_)
    //
    //
    // BDF2: for constant time step:
    //
    //               mom: hist_ = 4/3 veln_  - 1/3 velnm_
    //              (con: hist_ = 4/3 densn_ - 1/3 densnm_)
    //
    // -------------------------------------------------------------------
    TIMEINT_THETA_BDF2::SetOldPartOfRighthandside(
            veln_, velnm_, accn_,
            timealgo_, dta_, theta_,
            hist_);

    // -------------------------------------------------------------------
    //                     do explicit predictor step
    //
    //                      +-                                      -+
    //                      | /     dta \          dta  veln_-velnm_ |
    // velnp_ = veln_ + dta | | 1 + --- | accn_ - ----- ------------ |
    //                      | \     dtp /          dtp     dtp       |
    //                      +-                                      -+
    //
    // -------------------------------------------------------------------
    //
    // We cannot have a predictor in case of monolithic FSI here. There needs to
    // be a way to turn this off.
    if (extrapolationpredictor_)
    {
        if (step_>1)
        {
            TIMEINT_THETA_BDF2::ExplicitPredictor(
                    veln_, velnm_, accn_,
                    timealgo_, dta_, dtp_,
                    velnp_);
        }
    }

    // restore pressure. write back values from temporary vector to velnp_
    velpressplitter_.InsertCondVector(pretemp,velnp_);

    // -------------------------------------------------------------------
    //  evaluate Dirichlet and Neumann boundary conditions
    //  (the latter here only if no Neumann inflow)
    // -------------------------------------------------------------------
    {
        ParameterList eleparams;

        // total time required for Dirichlet conditions
        eleparams.set("total time",time_);

        // set vector values needed by elements
        discret_->ClearState();
        discret_->SetState("velnp",velnp_);

        // predicted dirichlet values
        // velnp then also holds prescribed new dirichlet values
        discret_->EvaluateDirichlet(eleparams,velnp_,null,null,null);
        discret_->ClearState();

        // set all parameters and states required for Neumann conditions
        if (timealgo_==timeint_afgenalpha)
        {
            eleparams.set("total time",time_-(1-alphaF_)*dta_);
            eleparams.set("thsl",1.0);
        }
        else
        {
            eleparams.set("total time",time_);
            eleparams.set("thsl",theta_*dta_);
        }
        // For generalized-alpha, the following is an approximation,
        // since actually vedeaf would be required, which is not yet
        // available, though.
        discret_->SetState("vedenp",vedenp_);
        neumann_loads_->PutScalar(0.0);

        // evaluate Neumann conditions
        discret_->EvaluateNeumann(eleparams,*neumann_loads_);
        discret_->ClearState();
    }

    // -------------------------------------------------------------------
    //  For af-generalized-alpha time-integration scheme:
    //  set "pseudo-theta", calculate initial accelerations according to
    //  prescribed Dirichlet values for generalized-alpha time
    //  integration.
    // -------------------------------------------------------------------
    if (timealgo_==timeint_afgenalpha)
    {
        // starting algorithm
        if (startalgo_)
        {
            // use backward-Euler-type parameter combination
            if (step_<=numstasteps_)
            {
                alphaM_ = 1.0;
                alphaF_ = 1.0;
                gamma_  = 1.0;
            }
            else
            {
                alphaM_ = params_.get<double>("alpha_M");
                alphaF_ = params_.get<double>("alpha_F");
                gamma_  = params_.get<double>("gamma");
                startalgo_ = false;
            }
        }

        // set "pseudo-theta" for af-generalized-alpha scheme
        theta_ = alphaF_*gamma_/alphaM_;

        // --------------------------------------------------
        // adjust accnp according to Dirichlet values of velnp
        //
        //                                  n+1     n
        //                               vel   - vel
        //       n+1      n  gamma-1.0      (0)
        //    acc    = acc * --------- + ------------
        //       (0)           gamma      gamma * dt
        //
        // in case of conservative form: velocity*density
        accnp_->Update(1.0,*velnp_,-1.0,*veln_,0.0);
        accnp_->Update((gamma_-1.0)/gamma_,*accn_,1.0/(gamma_*dta_));

        // -------------------------------------------------------------------
        // compute values at intermediate time steps for generalized-alpha
        // -------------------------------------------------------------------
        GenAlphaIntermediateValues();
    }
}

void FLD::FluidProjectionMethod::GenAlphaIntermediateValues()
{
    //       n+alphaM                n+1                      n
    //    acc         = alpha_M * acc     + (1-alpha_M) *  acc
    //       (i)                     (i)
    accam_->Update((alphaM_),*accnp_,(1.0-alphaM_),*accn_,0.0);

    //       n+alphaF              n+1                   n
    //      u         = alpha_F * u     + (1-alpha_F) * u
    //       (i)                   (i)
    velaf_->Update((alphaF_),*velnp_,(1.0-alphaF_),*veln_,0.0);
} // FluidProjectionMethod::GenAlphaIntermediateValues

/*-------------------------------------------------------
 * ------------------------------------------------------*/
void FLD::FluidProjectionMethod::TimeUpdate()
{
    // compute accelerations
    TIMEINT_THETA_BDF2::CalculateAcceleration(
            velnp_, veln_, velnm_, accn_,
            timealgo_, step_, theta_, dta_, dtp_,
            accnp_);

    // update old acceleration
    accn_->Update(1.0,*accnp_,0.0);

    // velocities/pressures of this step become most recent
    // velocities/pressures of the last step
    velnm_->Update(1.0,*veln_ ,0.0);
    veln_ ->Update(1.0,*velnp_,0.0);

    if (alefluid_)
    {
        dispnm_->Update(1.0,*dispn_,0.0);
        dispn_ ->Update(1.0,*dispnp_,0.0);
    }

    return;
}// FluidProjectionMethod::TimeUpdate

/*-------------------------------------------------------
 * ------------------------------------------------------*/
void FLD::FluidProjectionMethod::Output()
{
    // output of solution
    if (step_%upres_ == 0)
    {
        // step number and time
        output_.NewStep(step_,time_);

        // velocity/pressure vector
        output_.WriteVector("velnp",velnp_);

        // (hydrodynamic) pressure
        Teuchos::RCP<Epetra_Vector> pressure = velpressplitter_.ExtractCondVector(velnp_);
        pressure->Scale(density_);
        output_.WriteVector("pressure", pressure);

        //output_.WriteVector("residual", trueresidual_);
        if (alefluid_) output_.WriteVector("dispnp", dispnp_);

        // write domain decomposition for visualization (only once!)
        if (step_==upres_) output_.WriteElementData();

        if (uprestart_ != 0 && step_%uprestart_ == 0) //add restart data
        {
            // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
            output_.WriteVector("accnp",accnp_);
            output_.WriteVector("accn", accn_);
            output_.WriteVector("veln", veln_);
            output_.WriteVector("velnm",velnm_);

            if (alefluid_)
            {
                output_.WriteVector("dispn", dispn_);
                output_.WriteVector("dispnm",dispnm_);
            }
        }
    }
    // write restart also when uprestart_ is not a integer multiple of upres_
    else if (uprestart_ != 0 && step_%uprestart_ == 0)
    {
        // step number and time
        output_.NewStep(step_,time_);

        // velocity/pressure vector
        output_.WriteVector("velnp",velnp_);

        //output_.WriteVector("residual", trueresidual_);
        if (alefluid_)
        {
            output_.WriteVector("dispnp", dispnp_);
            output_.WriteVector("dispn", dispn_);
            output_.WriteVector("dispnm",dispnm_);
        }

        // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
        output_.WriteVector("accnp",accnp_);
        output_.WriteVector("accn", accn_);
        output_.WriteVector("veln", veln_);
        output_.WriteVector("velnm",velnm_);
    }

    return;
} // FluidProjectionMethod::Output

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the pressure correction method
 | ready for partitioned FSI
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidProjectionMethod::ProjectionSolve()
{
#ifdef DEBUG
    // grobe Zeitmessung
    const double tcpu=ds_cputime();
#endif

    // Weiche ob semiimplicit oder implicit solve
    const bool implicitSolve  = params_.get<bool>("PROJ_IMPLICIT",true);
    if(implicitSolve)
        ProjectionSolveImplicit();
    else
        ProjectionSolveSemiImplicit();

#ifdef DEBUG
    const double tcpu2=ds_cputime();
    //(*log_) << tcpu2-tcpu << endl;
    cout << "comput. time projection solve: " << tcpu2-tcpu << " s" << endl;
    //log_->flush();
#endif
}

/*-------------------------------------------------------
 * ------------------------------------------------------*/
void FLD::FluidProjectionMethod::ProjectionSolveImplicit()
{
    // wipe out global increment variables (perhaps these variables would better be local)
    veltilde_->PutScalar(0.0);
    phi_->PutScalar(0.0);

    // get mass matrix, inverse lumped mass matrix and gradient operator
    // for non-ALE problem only once per time-step
    if(step_ == 1 || alefluid_==true) getMatrices();

    // step1: solve simplified impulse equation (without div-free restriction)
    // input: velnp (with pressure solution for time n), ALE specific info
    // returns veltilde (and velnp with veltilde as velocity part and "old" pressure part (time n as before))
    SolveImpulseEqn();

    // step2: solve PPE for projection
    // input: veltilde
    // returns: phi (pressure increment)
    ProjectionStep();

    // step3: update
    // update velocity (velocity solution in velnp -> velocity output)
    // update pressure (pressure solution in pressure part of velnp is used for SolveImpulseEqnEx in next time step!)
    VelUpdate();

    // we suppose ALE case = FSI case
    if(alefluid_==true) CalcResidual();
}

/*-------------------------------------------------------
 * ------------------------------------------------------*/
void FLD::FluidProjectionMethod::ProjectionSolveSemiImplicit()
{
    // wipe out global increment variables (perhaps these variables would better be local)
    veltilde_->PutScalar(0.0);
    phi_->PutScalar(0.0);

    // get mass matrix, inverse lumped mass matrix and gradient operator
    // for non-ALE problem only once per time-step
    if(step_ == 1 || alefluid_==true) getMatrices();

    // step1: solve simplified impulse equation (without div-free restriction)
    // input: velnp (with pressure solution for time n), ALE specific info
    // returns veltilde (and velnp with veltilde as velocity part and "old" pressure part (time n as before))
    if(step_ == 1)	SolveImpulseEqn();
    else			SolveImpulseEqnSemi();

    // step2: solve PPE for projection
    // input: veltilde
    // returns: phi (pressure increment)
    ProjectionStep();

    // step3: update
    // update velocity (velocity solution in velnp -> velocity output)
    // update pressure (pressure solution in pressure part of velnp is used for SolveImpulseEqnEx in next time step!)
    VelUpdate();

    // we suppose ALE case = FSI case
    if(alefluid_==true) CalcResidual();
}

/*!
* calculates and stores massmat, gradop and lmassinv
* G^T * (M_L^{-1}*G)
*/
void FLD::FluidProjectionMethod::getMatrices()
{
    // time measurement
    TEUCHOS_FUNC_TIME_MONITOR(" ~ getMatrices");

    // element call for massmat and gradop

    // define local variables (filled by element call)
    RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > blockmat = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(velpressplitter_,velpressplitter_,72,false,true));
    RCP<Epetra_Vector> lmassmat = LINALG::CreateVector(*(discret_->DofRowMap()),true);  // vector for lumped mass matrix (global)

    // zero out
    blockmat->Zero();

    // set action for element call
    params_.set("action","calc_gradop_and_massmatrix");

    // element call
    discret_->ClearState();
    if(alefluid_) discret_->SetState("dispnp", dispnp_);
    discret_->Evaluate(params_,blockmat,lmassmat);
    discret_->ClearState();

    // set dirichlet BC entries to one
    Teuchos::RCP<Epetra_Vector> dirichones = LINALG::CreateVector(*(dbcmaps_->CondMap()),false);
    dirichones->PutScalar(1.0);
    dbcmaps_->InsertCondVector(dirichones, lmassmat);

    // build inverse lumped mass "matrix"
    lmassinvvec_ = velpressplitter_.ExtractOtherVector(lmassmat);
    lmassinvvec_->Reciprocal(*lmassinvvec_);

    // complete matrices
    blockmat->Complete();

    gradopwithoutbc_ = rcp(new LINALG::SparseMatrix(blockmat->Matrix(0,1)));
    blockmat->ApplyDirichlet(*(dbcmaps_->CondMap()));	// dirichlet bc

    // store matrices in global variables
    massmat_ = rcp(new LINALG::SparseMatrix(blockmat->Matrix(0,0),View));
    gradop_ = rcp(new LINALG::SparseMatrix(blockmat->Matrix(0,1),View));

    // calculate pressure matrix

    // multiply M_L^{-1} * G
    RCP<LINALG::SparseMatrix> er1 = Teuchos::rcp(new LINALG::SparseMatrix(*VelocityRowMap(),24,false,true));
    er1->Add(*gradop_,false,1.0,0.0);
    er1->Complete(gradop_->DomainMap(),gradop_->RangeMap());
    er1->LeftScale(*lmassinvvec_);

    // Multiply G^T * (M_L^{-1}*G)
    pressmat_ = LINALG::Multiply(*gradop_,true,*er1,false,true);
}	// end FLD::FluidProjectionMethod::getMatrices()

/*!*****************************************************
* \brief
* Solve simplified impulse equation (impulse equation without algebraic div-free restriction)
* use iterative fixed point solver with explicit treatment of pressure term M
* nonlinear solver
* dot{u} + Ku + N(u)u = f - M*M_L^{-1}*G p(n)
* Input: velnp (full vector but only velocity part for timestep n+1 from predictor or elsewhere)
* ALE-specific terms
* Output: solution veltilde (intermediate velocity, stored in velnp)
*/
void FLD::FluidProjectionMethod::SolveImpulseEqn()
{
    // time measurement: SolveImpulseEqn
    TEUCHOS_FUNC_TIME_MONITOR(" ~ SolveImpulseEqn");

    // define local variables
    Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > blockmat = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(velpressplitter_,velpressplitter_,72,false,true));
    Teuchos::RCP<Epetra_Vector> imprhs = LINALG::CreateVector(*(discret_->DofRowMap()),true);

    // loacl increment variable for velocity
    Teuchos::RCP<Epetra_Vector> incvel = LINALG::CreateVector(*VelocityRowMap(),true);

    // vectors for calculation of RHS
    Teuchos::RCP<Epetra_Vector> MMlGpn = LINALG::CreateVector(*VelocityRowMap(),true);
    Teuchos::RCP<Epetra_Vector> temp = LINALG::CreateVector(*VelocityRowMap(),true);

    // ---------------------------------------------- nonlinear iteration
    // ------------------------------- stop nonlinear iteration when both
    //                                 increment-norms are below this bound
    const double  ittol     =params_.get<double>("tolerance for nonlin iter");

    int               itnum = 0;
    int               itemax = 0;
    bool              stopnonliniter = false;

    itemax  = params_.get<int>   ("max nonlin iter steps");
    //------------------------------ turn adaptive solver tolerance on/off
    const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
    const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);

    // iteration solver loop for nonlinear system
    while (stopnonliniter==false)
    {
        itnum++;

        //------------------------------------------------------------
        // prepare element call
        //------------------------------------------------------------
        // zero out iterative variables
        incvel->Scale(0.0);
        MMlGpn->Scale(0.0);
        temp->Scale(0.0);
        imprhs->Scale(0.0);
        blockmat->Zero();

        // add Neumann loads
        imprhs->Update(1.0,*neumann_loads_,0.0);

        // create the parameters for the discretization
        ParameterList eleparams;
        eleparams.set("action","calc_impulseeqn_implicit");
        eleparams.set("total time",time_);
        eleparams.set("thsl",theta_*dta_);
        eleparams.set("dt",dta_);
        eleparams.set("fs subgrid viscosity","no");
        eleparams.set("Linearisation",newton_);
        eleparams.set("low-Mach-number solver","no");
        eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");
        eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

        // set vector values needed by elements
        discret_->ClearState();
        discret_->SetState("velnp",velnp_);
        discret_->SetState("vedenp",vedenp_);
        discret_->SetState("hist"  ,hist_ );
        if (alefluid_)
        {
            discret_->SetState("dispnp", dispnp_);
            discret_->SetState("gridv", gridv_);
        }

        //------------------------------------------------------------
        // element call
        //------------------------------------------------------------
        // convergence check at itemax is skipped for speedup if
        // CONVCHECK is set to L_2_norm_without_residual_at_itemax
        if ((itnum != itemax)
                ||
                (params_.get<string>("CONVCHECK","L_2_norm")
                        !=
                                "L_2_norm_without_residual_at_itemax"))
        {
            // call standard loop over elements
            discret_->Evaluate(eleparams,blockmat,imprhs);

            discret_->ClearState();

            // finalize the complete matrix
            blockmat->Complete();
        }


        //------------------------------------------------------------
        // set residual DOFs for dirichlet BCs
        //------------------------------------------------------------

        // blank residual DOFs which are on Dirichlet BC
        // We can do this because the values at the dirichlet positions
        // are not used anyway.
        // We could avoid this though, if velrowmap_ and prerowmap_ would
        // not include the dirichlet values as well. But it is expensive
        // to avoid that.
        {
            Teuchos::RCP<Epetra_Vector> dirichzeros = LINALG::CreateVector(*(dbcmaps_->CondMap()),true);
            dbcmaps_->InsertCondVector(dirichzeros, imprhs);
        }

        //-------------------------------------------------
        // adapt RHS
        //-------------------------------------------------

        // compute M * M_L^{-1} * [G*p(n)] -> temp

        // use gradop from getMatrices

        // calculate (G*p(n)) -> MMlGpn
        RCP<Epetra_Vector> velnp2 = velpressplitter_.ExtractCondVector(velnp_);
        gradop_->Apply(*velnp2,*MMlGpn);

        // Multiply lmassinv M_L11{-1}*(G*p(n)) -> temp
        temp->Multiply(1.0,*lmassinvvec_,*MMlGpn,0.0);  // 09.01 VZ

        // get M * [M_L^{-1}*[G*p(n)]] -> MMlGpn
        massmat_->Apply(*temp,*MMlGpn);

        // velocity-part of imprhs for RHS
        RCP<Epetra_Vector> rhs = velpressplitter_.ExtractOtherVector(imprhs);

        // update RHS: RHS - theta*dt * [M * [M_L^{-1}*[G*p(n)]]]
        rhs->Update(-theta_*dta_,*MMlGpn,1.0);

        //------------------------------------------------------------
        // prepare exit conditions
        //------------------------------------------------------------
        double incvelnorm_L2;
        double velnorm_L2;
        double vresnorm;

        // Normen bestimmen
        rhs->Norm2(&vresnorm);  // residual norm RHS (velocity-part)
        incvel->Norm2(&incvelnorm_L2);  // norm of velocity increment
        Teuchos::RCP<Epetra_Vector> onlyvel = velpressplitter_.ExtractOtherVector(velnp_);  // voller Vektor
        velpressplitter_.ExtractOtherVector(velnp_,onlyvel);
        onlyvel->Norm2(&velnorm_L2);  // Norm velnp_ (velocity-part)

        // care for the case that nothing really happens in the velocity field
        if (velnorm_L2 < 1e-5)
        {
            velnorm_L2 = 1.0;
        }

        //------------------------------------------------------------
        // convergence check
        //------------------------------------------------------------
        // this is the convergence check
        // We always require at least one solve. Otherwise the
        // perturbation at the FSI interface might get by unnoticed.
        if (vresnorm <= ittol and
                incvelnorm_L2/velnorm_L2 <= ittol)
        {
            if(itnum==1)
            {
                if(myrank_==0)
                {
                    printf("+------------+-------------------+--------------+\n");
                    printf("| step/max   | tol       [norm]  |   vel res    |\n");
                    printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   |",
                            itnum,itemax,ittol,vresnorm);
                    printf("\n");
                }
            }
            else
            {
                stopnonliniter=true;
                // Screenoutput
                if (myrank_ == 0)
                {
                    printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   |",
                            itnum,itemax,ittol,vresnorm);
                    printf(")\n");
                    printf("+------------+-------------------+--------------+\n");

                    FILE* errfile = params_.get<FILE*>("err file",NULL);
                    if (errfile!=NULL)
                    {
                        fprintf(errfile,"fluid solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E\n",
                                itnum,itemax,ittol,vresnorm);
                    }
                }
                break;
            }
        }
        else // if not yet converged
            if (myrank_ == 0)
            {
                if(itnum==1)
                {
                    printf("+------------+-------------------+--------------+\n");
                    printf("| step/max   | tol       [norm]  |   vel res    |\n");
                }
                printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   |",
                        itnum,itemax,ittol,vresnorm);

                printf("\n");
            }

        // warn if itemax is reached without convergence, but proceed to
        // next timestep...
        if ((itnum == itemax) and (vresnorm > ittol or
                incvelnorm_L2/velnorm_L2 > ittol))
        {
            stopnonliniter=true;
            if (myrank_ == 0)
            {
                printf("+---------------------------------------------------------------+\n");
                printf("|            >>>>>> not converged in itemax steps!              |\n");
                printf("+---------------------------------------------------------------+\n");

                FILE* errfile = params_.get<FILE*>("err file",NULL);
                if (errfile!=NULL)
                {
                    fprintf(errfile,"fluid unconverged solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
                            itnum,itemax,ittol,vresnorm,0.0,
                            incvelnorm_L2/velnorm_L2,0.0);
                }
            }
            break;
        }


        //------------------------------------------------------------
        // apply dirichlet BCs in GLS
        //------------------------------------------------------------
        //--------- Apply dirichlet boundary conditions to system of equations
        //          residual discplacements are supposed to be zero at
        //          boundary conditions
        incvel->PutScalar(0.0);

        //Teuchos::RCP<LINALG::SparseMatrix> bop = rcp(new LINALG::SparseMatrix(blockmat->Matrix(0,0))); // INEFFIZIENT
        Teuchos::RCP<LINALG::SparseMatrix> bop = rcp(&(blockmat->Matrix(0,0)),false);
        {
            // time measurement: application of dbc
            TEUCHOS_FUNC_TIME_MONITOR("      + apply DBC");

            Teuchos::RCP<Epetra_Vector> onlyvelzeros = velpressplitter_.ExtractOtherVector(zeros_);
            Teuchos::RCP<Epetra_Vector> dirichones = LINALG::CreateVector(*(dbcmaps_->CondMap()),false);
            dirichones->PutScalar(1.0);
            Teuchos::RCP<Epetra_Vector> dirichtoggle = LINALG::CreateVector(*(discret_->DofRowMap()),true);
            dbcmaps_->InsertCondVector(dirichones, dirichtoggle);
            Teuchos::RCP<Epetra_Vector> onlyveldirichtoggle = velpressplitter_.ExtractOtherVector(dirichtoggle);
            LINALG::ApplyDirichlettoSystem(bop,incvel,rhs,onlyvelzeros,onlyveldirichtoggle);
        }

        //------------------------------------------------------------
        // solve GLS
        //------------------------------------------------------------
        // do adaptive linear solver tolerance (not in first solve)
        if (isadapttol && itnum>1)
        {
            double currresidual = vresnorm;
            currresidual = max(currresidual,incvelnorm_L2/velnorm_L2);
            solver_.AdaptTolerance(ittol,currresidual,adaptolbetter);
        }
        solver_.Solve(bop->EpetraOperator(),incvel,rhs,true,itnum==1);
        solver_.ResetTolerance();

        //------------------------------------------------------------
        // incemental Update
        //------------------------------------------------------------
        velpressplitter_.AddOtherVector(incvel,velnp_);
    } // <- end iterative solver

    veltilde_ = velpressplitter_.ExtractOtherVector(velnp_);
}

/*!*****************************************************
* \brief
* Solve simplified impulse equation (impulse equation without algebraic div-free restriction)
* use iterative fixed point solver with explicit treatment of pressure term M
* dot{u} + Ku + N(u)u = f - M*M_L^{-1}*G p(n)
* semi-implicit treatment of convection!
* Input: velnp (full vector but only velocity part for timestep n+1 from predictor or elsewhere)
* ALE-specific terms
* Output: solution veltilde (intermediate velocity, stored in velnp)
*/
void FLD::FluidProjectionMethod::SolveImpulseEqnSemi()
{
    // semi-implicit treatment of convection N(u_n)\tilde{u}
    // -> solve linear equation

    // time measurement: solveImpulseEqnSemi
    TEUCHOS_FUNC_TIME_MONITOR(" ~ SolveImpulseEqnSemi");

    // local variables
    Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > blockmat = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(velpressplitter_,velpressplitter_,72,false,true));
    Teuchos::RCP<Epetra_Vector> imprhs = LINALG::CreateVector(*(Discretization()->DofRowMap()),true);

    // local increment variable for velocity
    Teuchos::RCP<Epetra_Vector> incvel = LINALG::CreateVector(*VelocityRowMap(),true);

    // vector variables for calculation of RHS
    Teuchos::RCP<Epetra_Vector> MMlGpn = LINALG::CreateVector(*VelocityRowMap(),true);
    Teuchos::RCP<Epetra_Vector> temp = LINALG::CreateVector(*VelocityRowMap(),true);

    //------------------------------------------------------------
    // prepare element call
    //------------------------------------------------------------
    incvel->Scale(0.0);
    MMlGpn->Scale(0.0);
    temp->Scale(0.0);
    imprhs->Scale(0.0);
    blockmat->Zero();

    // add Neumann loads
    imprhs->Update(1.0,*neumann_loads_,0.0);

    // create the parameters for the discretization
    ParameterList eleparams;
    eleparams.set("action","calc_impulseeqn_semiimplicit");
    eleparams.set("total time",time_);
    eleparams.set("thsl",theta_*dta_);
    eleparams.set("dt",dta_);
    eleparams.set("fs subgrid viscosity","no");
    eleparams.set("Linearisation",newton_);
    eleparams.set("low-Mach-number solver","no");
    eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");
    eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("velnp",velnp_);
    discret_->SetState("veln",veln_);   // alte Daten
    discret_->SetState("vedenp",vedenp_);
    discret_->SetState("hist"  ,hist_ );
    if (alefluid_)
    {
        discret_->SetState("dispnp", dispnp_);
        discret_->SetState("gridv", gridv_);  // brauch ich gridv wirklich?
    }

    //------------------------------------------------------------
    // element call
    //------------------------------------------------------------
    {
        // call standard loop over elements
        discret_->Evaluate(eleparams,blockmat,imprhs);

        discret_->ClearState();

        // finalize the complete matrix
        blockmat->Complete();
    }

    //------------------------------------------------------------
    // set residual DOFs for dirichlet BCs
    //------------------------------------------------------------

    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the dirichlet positions
    // are not used anyway.
    // We could avoid this though, if velrowmap_ and prerowmap_ would
    // not include the dirichlet values as well. But it is expensive
    // to avoid that.
    {
        Teuchos::RCP<Epetra_Vector> dirichzeros = LINALG::CreateVector(*(dbcmaps_->CondMap()),true);
        dbcmaps_->InsertCondVector(dirichzeros, imprhs);
    }

    //-------------------------------------------------
    // adapt RHS
    //-------------------------------------------------

    // compute M * M_L^{-1} * [G*p(n)] -> temp

    // use gradop from getMatrices

    // calculate (G*p(n)) -> MMlGpn
    RCP<Epetra_Vector> velnp2 = velpressplitter_.ExtractCondVector(velnp_);
    gradop_->Apply(*velnp2,*MMlGpn);

    // Multiply lmassinv M_L11{-1}*(G*p(n)) -> temp
    temp->Multiply(1.0,*lmassinvvec_,*MMlGpn,0.0);    // 09.01 VZ

    // get M * [M_L^{-1}*[G*p(n)]] -> MMlGpn
    massmat_->Apply(*temp,*MMlGpn);

    // velocity-part of imprhs for RHS
    RCP<Epetra_Vector> rhs = velpressplitter_.ExtractOtherVector(imprhs);

    // update RHS: RHS - theta*dt * [M * [M_L^{-1}*[G*p(n)]]]
    rhs->Update(-theta_*dta_,*MMlGpn,1.0);

    //------------------------------------------------------------
    // apply dirichlet BCs in GLS
    //------------------------------------------------------------
    //--------- Apply dirichlet boundary conditions to system of equations
    //          residual discplacements are supposed to be zero at
    //          boundary conditions
    incvel->PutScalar(0.0);

    Teuchos::RCP<LINALG::SparseMatrix> bop = rcp(&(blockmat->Matrix(0,0)),false);
    {
        // time measurement: application of dbc
        TEUCHOS_FUNC_TIME_MONITOR("      + apply DBC");

        Teuchos::RCP<Epetra_Vector> onlyvelzeros = velpressplitter_.ExtractOtherVector(zeros_);
        Teuchos::RCP<Epetra_Vector> dirichones = LINALG::CreateVector(*(dbcmaps_->CondMap()),false);
        dirichones->PutScalar(1.0);
        Teuchos::RCP<Epetra_Vector> dirichtoggle = LINALG::CreateVector(*(discret_->DofRowMap()),true);
        dbcmaps_->InsertCondVector(dirichones, dirichtoggle);
        Teuchos::RCP<Epetra_Vector> onlyveldirichtoggle = velpressplitter_.ExtractOtherVector(dirichtoggle);
        LINALG::ApplyDirichlettoSystem(bop,incvel,rhs,onlyvelzeros,onlyveldirichtoggle);
    }

    //------------------------------------------------------------
    // solve GLS
    //------------------------------------------------------------
    solver_.Solve(bop->EpetraOperator(),incvel,rhs,true,true);

    //------------------------------------------------------------
    // incemental Update
    //------------------------------------------------------------
    velpressplitter_.AddOtherVector(incvel,velnp_);

    //------------------------------------------------------------
    // calculate resnorm for printout
    //------------------------------------------------------------
    double vresnorm;

    // Normen bestimmen
    rhs->Norm2(&vresnorm);  // residual norm RHS (velocity-part)

    if(myrank_==0)
    {
        printf("+-------------------------------------------+\n| vresnorm  | %10.3E   |\n+-------------------------------------------+",vresnorm);
        printf("\n");
    }

    veltilde_ = velpressplitter_.ExtractOtherVector(velnp_);
}

/*!
* function for projection step solve A phi = G^T veltilde with
* A = pressmat from getMatrices G = gradop
* veltilde: intermediate velocity (solution of simplified impulse equation)
* input: velnp, ALE information
* output: phi (pressure increment)
*/
void FLD::FluidProjectionMethod::ProjectionStep()
{
    // time measurement: ProjectionStep
    TEUCHOS_FUNC_TIME_MONITOR(" ~ ProjectionStep");

    // calculate G^T veltilde
    RCP<Epetra_Vector> veltemp = velpressplitter_.ExtractOtherVector(velnp_);
    RCP<Epetra_Vector> prhs = LINALG::CreateVector(*PressureRowMap(),true);
    gradopwithoutbc_->Multiply(true,*veltemp,*prhs);
    // solve pressure poisson equation
    psolver_.Solve(pressmat_->EpetraMatrix(),phi_,prhs,true,true);
}

/*!
* Update step for projection method update velocity part:
* velnp = veltilde - M_L^{-1} G
* phi update
* pressure part: velnp = velnp + phi/(theta*dta)
* input: phi, ALE information
* output:
* velnp (velocity part of velocity at n+1)
* velnp (pressure part of pressure at n+1)
*/
void FLD::FluidProjectionMethod::VelUpdate()
{
    // time measurement: VelUpdate
    TEUCHOS_FUNC_TIME_MONITOR(" ~ VelUpdate");

    // calculate G*phi with gradop

    // calculate (G*phi) -> tempphi
    Teuchos::RCP<Epetra_Vector> tempphi= LINALG::CreateVector(*VelocityRowMap(),true);
    gradop_->Apply(*phi_,*tempphi);

    // Multiply lmassinv M_L11{-1}*(G*p(n)) -> tempphi2
    Teuchos::RCP<Epetra_Vector> tempphi2= LINALG::CreateVector(*VelocityRowMap(),true);
    tempphi2->Multiply(1.0,*lmassinvvec_,*tempphi,0.0); // 09.01 VZ

    // velocity update, provide end solution for time step n+1 in velnp
    veltilde_->Update(-1.0,*tempphi2,1.0);
    velpressplitter_.InsertOtherVector(veltilde_,velnp_);

    // pressure update
    Teuchos::RCP<Epetra_Vector> temp3= LINALG::CreateVector(*PressureRowMap(),true);
    temp3->Update(1/(dta_*theta_),*phi_,0.0);

    // update pressure part:
    // velnp_ = velnp_ + phi/(theta dt)
    // because velnp_ (pressure part) hasn't changed until now it is pressure part of velnp for time n and we can just add the new values
    velpressplitter_.AddCondVector(temp3,velnp_);
}

/*!
* calculates residual only (needed for partitioned FSI)
*/
void FLD::FluidProjectionMethod::CalcResidual()
{
    // element call for calculation of full residual
    // needed for interface-forces (partitioned methods for FSI use trueresidual_ for extracting interface-forces)

    Teuchos::RCP<Epetra_Vector> residual = LINALG::CreateVector(*(discret_->DofRowMap()),true);
    residual->PutScalar(0.0);

    // add Neumann loads
    residual->Update(1.0,*neumann_loads_,0.0);  // brauch ich das?

    // create the parameters for the discretization
    ParameterList eleparams;
    eleparams.set("action","calc_fluid_residual");
    eleparams.set("total time",time_);
    eleparams.set("thsl",theta_*dta_);
    eleparams.set("dt",dta_);
    eleparams.set("fs subgrid viscosity","no");
    eleparams.set("Linearisation",newton_);
    eleparams.set("low-Mach-number solver","no");
    eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");
    eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("velnp",velnp_);
    discret_->SetState("vedenp",vedenp_);
    discret_->SetState("hist"  ,hist_ );
    if (alefluid_)
    {
        discret_->SetState("dispnp", dispnp_);
        discret_->SetState("gridv", gridv_);
    }

    //------------------------------------------------------------
    // element call
    //------------------------------------------------------------

    // call standard loop over elements
    discret_->Evaluate(eleparams,Teuchos::null,residual);

    discret_->ClearState();

    // true residual: full vector
    trueresidual_->Update(ResidualScaling(),*residual,0.0);
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
|                                                           chfoe 01/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidProjectionMethod::UpdateGridv()
{
    // get order of accuracy of grid velocity determination
    // from input file data
    const int order  = params_.get<int>("order gridvel");

    switch (order)
    {
    case 1:
        /* get gridvelocity from BE time discretisation of mesh motion:
        -> cheap
        -> easy
        -> limits FSI algorithm to first order accuracy in time

                 x^n+1 - x^n
             uG = -----------
                    Delta t                        */
        gridv_->Update(1/dta_, *dispnp_, -1/dta_, *dispn_, 0.0);
        break;
    case 2:
        /* get gridvelocity from BDF2 time discretisation of mesh motion:
          -> requires one more previous mesh position or displacemnt
          -> somewhat more complicated
          -> allows second order accuracy for the overall flow solution  */
        gridv_->Update(1.5/dta_, *dispnp_, -2.0, *dispn_, 0.0);
        gridv_->Update(0.5, *dispnm_, 1.0);
        break;
    }
}

void FLD::FluidProjectionMethod::ReadRestart(int step)
{
    IO::DiscretizationReader reader(discret_,step);
    time_ = reader.ReadDouble("time");
    step_ = reader.ReadInt("step");

    reader.ReadVector(velnp_,"velnp");
    reader.ReadVector(veln_, "veln");
    reader.ReadVector(velnm_,"velnm");
    reader.ReadVector(accn_ ,"accn");

    if (alefluid_)
    {
        reader.ReadVector(dispnp_,"dispnp");
        reader.ReadVector(dispn_ , "dispn");
        reader.ReadVector(dispnm_,"dispnm");
    }
    // also read impedance bc information if required
    // Note: this method acts only if there is an impedance BC
    //impedancebc_->ReadRestart(reader);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::FluidProjectionMethod::AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd)
{
    std::vector<Teuchos::RCP<const Epetra_Map> > condmaps;
    condmaps.push_back(maptoadd);
    condmaps.push_back(dbcmaps_->CondMap());
    Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);
    *dbcmaps_ = LINALG::MapExtractor(*(discret_->DofRowMap()), condmerged);
    return;
}

#endif /* CCADISCRET */

