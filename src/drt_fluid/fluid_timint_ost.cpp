/*----------------------------------------------------------------------*/
/*!
\file fluid_timint_ost.cpp
\brief One-Step-Theta time-integration scheme

<pre>
Maintainers: Ursula Rasthofer & Martin Kronbichler
             {rasthofer,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_timint_ost.H"
#include "../drt_fluid_turbulence/scale_sep_gmo.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_turbulence/dyn_smag.H"
#include "../drt_fluid_turbulence/dyn_vreman.H"
#include "../drt_fluid_turbulence/boxfilter.H"

#include "../linalg/linalg_solver.H"



/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntOneStepTheta::TimIntOneStepTheta(
    const Teuchos::RCP<DRT::Discretization>&      actdis,
    const Teuchos::RCP<LINALG::Solver>&           solver,
    const Teuchos::RCP<Teuchos::ParameterList>&   params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output,
    bool                                          alefluid /*= false*/)
: FluidImplicitTimeInt(actdis,solver,params,output,alefluid)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  FLD::FluidImplicitTimeInt::Init();

  //check, if starting algorithm is desired
  if (numstasteps_ > 0)
    dserror("no starting algorithm supported for schemes other than af-gen-alpha");

  SetElementTimeParameter();

  CompleteGeneralInit();

  return;
}

/*----------------------------------------------------------------------*
| Destructor dtor (public)                                     bk 11/13 |
*-----------------------------------------------------------------------*/
FLD::TimIntOneStepTheta::~TimIntOneStepTheta()
{
  return;
}

/*----------------------------------------------------------------------*
| Print information about current time step to screen          bk 11/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::PrintTimeStepInfo()
{
  if (myrank_==0)
  {
    printf("TIME: %11.4E/%11.4E  DT = %11.4E   One-Step-Theta (%0.2f)   STEP = %4d/%4d \n",
          time_,maxtime_,dta_,theta_,step_,stepmax_);
  }
  return;
}

/*----------------------------------------------------------------------*
| set old part of right hand side                              bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::SetOldPartOfRighthandside()
{
  /*
     for low-Mach-number flow: distinguish momentum and continuity part
     (continuity part only meaningful for low-Mach-number flow)

     One-step-Theta:

                   mom: hist_ = veln_  + dt*(1-Theta)*accn_
                  (con: hist_ = densn_ + dt*(1-Theta)*densdtn_)
  */

      hist_->Update(1.0, *veln_, dta_*(1.0-theta_), *accn_, 0.0);


  return;
}

/*----------------------------------------------------------------------*
| set integration-scheme-specific state                        bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::SetStateTimInt()
{

  discret_->SetState("velaf",velnp_);

  return;
}

/*----------------------------------------------------------------------*
| calculate acceleration                                       bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::CalculateAcceleration(
    const Teuchos::RCP<const Epetra_Vector>    velnp,
    const Teuchos::RCP<const Epetra_Vector>    veln,
    const Teuchos::RCP<const Epetra_Vector>    velnm,
    const Teuchos::RCP<const Epetra_Vector>    accn,
    const Teuchos::RCP<Epetra_Vector>          accnp
)
{
  /*

  Following formulations are for n+1; acceleration values, however, are
  directly stored in vectors at time n (velocity has not yet been updated).

  One-step-Theta:

   acc(n+1) = (vel(n+1)-vel(n)) / (Theta * dt(n)) - (1/Theta -1) * acc(n)

  */

  const double fact1 = 1.0/(theta_*dta_);
  const double fact2 =-1.0/theta_ +1.0;   /* = -1/Theta + 1 */

  accnp->Update( fact1,*velnp,0.0);
  accnp->Update(-fact1,*veln ,1.0);
  accnp->Update( fact2,*accn,1.0);

  return;
}

/*----------------------------------------------------------------------*
| set gamma                                                    bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::SetGamma(Teuchos::ParameterList& eleparams)
{

  eleparams.set("gamma"  ,theta_);
  return;
}

/*----------------------------------------------------------------------*
| scale separation                                             bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::Sep_Multiply()
{

  Sep_->Multiply(false,*velnp_,*fsvelaf_);
  return;
}

/*----------------------------------------------------------------------*
 | paraview output of filtered velocity                  rasthofer 02/11|
 *----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::OutputofFilteredVel(
     Teuchos::RCP<Epetra_Vector> outvec,
     Teuchos::RCP<Epetra_Vector> fsoutvec)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  Teuchos::RCP<Epetra_Vector> row_finescaleveltmp;
  row_finescaleveltmp = Teuchos::rcp(new Epetra_Vector(*dofrowmap,true));

  // get fine scale velocity
  if (scale_sep_ == INPAR::FLUID::algebraic_multigrid_operator)
    Sep_->Multiply(false,*velnp_,*row_finescaleveltmp);
  else
    ScaleSepGMO_->ApplyScaleSeparation(velnp_,row_finescaleveltmp);
  // get filtered or coarse scale velocity
  outvec->Update(1.0,*velnp_,-1.0,*row_finescaleveltmp,0.0);

  fsoutvec->Update(1.0,*row_finescaleveltmp,0.0);

  return;
}

// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::TimIntOneStepTheta::SetElementTimeParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action",FLD::set_time_parameter);

  //set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);

  // set general element parameters
  eleparams.set("dt",dta_);
  eleparams.set("theta",theta_);
  eleparams.set("omtheta",1.0-theta_);

  // set scheme-specific element parameters and vector values
  eleparams.set("total time",time_);


  // call standard loop over elements
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  return;
}

/*----------------------------------------------------------------------*
| return time integration factor                               bk 12/13 |
*-----------------------------------------------------------------------*/
const double FLD::TimIntOneStepTheta::TimIntParam() const
{
  double retval = 0.0;
    // this is the interpolation weight for quantities from last time step
    retval = 0.0;
  return retval;
}


