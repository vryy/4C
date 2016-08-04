/*!----------------------------------------------------------------------
\file scatra_timint_ost_endoexocytosis.cpp

\brief One-Step-Theta time-integration scheme for endo-/exocytosis model

\level 2

\maintainer Andreas Rauch

*----------------------------------------------------------------------*/

#include "scatra_timint_ost_endoexocytosis.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 |  Constructor (public)                                    rauch 08/16 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntOneStepThetaEndoExocytosis::TimIntOneStepThetaEndoExocytosis(
    Teuchos::RCP<DRT::Discretization>        actdis,        //!< discretization
    Teuchos::RCP<LINALG::Solver>             solver,        //!< linear solver
    Teuchos::RCP<Teuchos::ParameterList>     params,        //!< parameter list
    Teuchos::RCP<Teuchos::ParameterList>     extraparams,   //!< supplementary parameter list
    Teuchos::RCP<IO::DiscretizationWriter>   output,        //!< output writer
    const int                                probnum        //!< global problem number
    )
:  ScaTraTimIntImpl(actdis,solver,params,extraparams,output),
   TimIntOneStepTheta(actdis,solver,params,extraparams,output),
   endoexo_delay_(-1),
   internalization_steps_(0),
   exo_surface_area_(0.0),
   delta_phi_(0.0),
   Delta_phi_(0.0),
   mean_scalars_(0,0.0),
   internalization_vec_(0),
   source_(Teuchos::null)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting
  // this has to be done before all state vectors are initialized

  return;
}


/*------------------------------------------------------------------------*
 |  initialize time integration                               rauch 08/16 |
 *------------------------------------------------------------------------*/
void SCATRA::TimIntOneStepThetaEndoExocytosis::Init()
{
  // initialize base class
  TimIntOneStepTheta::Init();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices: local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // create exocytotic source vector
  source_ = LINALG::CreateVector(*dofrowmap,true);

  // get delay between internalization and begin of exocytotic transport at exo. surface
  endoexo_delay_ =
      DRT::Problem::Instance()->CellMigrationParams().sublist("ENDOEXOCYTOSIS MODULE")
                                                     .get<double>("ENDOEXO_DELAY");

  int check = (int)(endoexo_delay_/dta_);
  if(double(check)!=(endoexo_delay_/dta_))
    dserror("ENDOEXO_DELAY needs to be a multiple of TIMESTEP !");

  // initialize internalization vector
  internalization_vec_.resize((int)(endoexo_delay_/dta_),0.0);
  internalization_steps_ = (int)(internalization_vec_.size());

  // initialize vector in which the deltas (n -> n+1) are stored
  // entry 0 contains the difference at t-T
  Delta_phi_.resize(internalization_steps_,0.0);
}


/*------------------------------------------------------------------------*
 |  manipulations of scatra algorithm prior to Solve()        rauch 08/16 |
 *------------------------------------------------------------------------*/
void SCATRA::TimIntOneStepThetaEndoExocytosis::PreSolve()
{
  // temp variables
  const int dof_id_internalized_scalar = discret_->GetCondition("ScaTraCellExt")->GetInt("ScalarID");

  // calculate the new difference in molecule numbers from n -> n+1
  if(step_>internalization_steps_)
  {
    // Member step_ is incremented before Solve. Thus, it starts with 1 for the first time step.
    int entry =((step_-1)-internalization_steps_)%(internalization_steps_);
    delta_phi_ = Delta_phi_[entry];
  }
  else
    delta_phi_ = 0.0;

  /////////////////////////////////////////////////////////////////////////////////////
  // Calc area of exocytotic surface.
  // The surface needs to have the condition 'ScaTraCellExt'.
  // At this surface the biomolecules internalized at t-endoexo_delay_
  // are segregated.

  // declare parameter list
  Teuchos::ParameterList eleparams;

  // reinitialize surface area
  exo_surface_area_=0.0;

  // set action for surface area calculation
  eleparams.set<int>("action",SCATRA::bd_calc_boundary_integral);

  // create result vector
  Teuchos::RCP<Epetra_SerialDenseVector> area = Teuchos::rcp(new Epetra_SerialDenseVector(1));

  // evaluate over surface "ScaTraCellExt"
  discret_->EvaluateScalars(eleparams,area,"ScaTraCellExt");

  // extract surface are from result vector
  exo_surface_area_ = (*area)[0];

  // safety check
  if (std::abs(exo_surface_area_) < 1E-15)
    dserror("Exocytosis surface area is close to zero!\n"
            "Something went wrong.\n"
            "Check definition of 'ScaTraCellExt' conditioned surface.");

  /////////////////////////////////////////////////////////////////////////////////////
  // Calc exocytotic source term at externalization surface.
  // The surface needs to have the condition 'ScaTraCellExt'.
  // At this surface the biomolecules internalized at t-endoexo_delay_
  // are segregated.

  // reinitialize source vector
  source_->Scale(0.0);

  // set action for source evaluation
  eleparams.set<int>("action",SCATRA::bd_integrate_weighted_scalar);

  // provide other data to evaluate routine
  eleparams.set<double>("scalar",(delta_phi_/dta_));
  eleparams.set<double>("user defined prefac",-(1.0/exo_surface_area_));
  eleparams.set<int>("ScalarID",0);

  // evaluate the source
  discret_->EvaluateCondition(eleparams,source_,"ScaTraCellExt");

  double sourcenorm=-1234.0;
  source_->Norm2(&sourcenorm);


  if(discret_->Comm().MyPID()==0)
  {
    std::cout<<"########################################################"<<std::endl;
    std::cout<<"# Externalized "<<std::setprecision(7)<<delta_phi_<<" molecules of species "<< dof_id_internalized_scalar<<std::endl;
    std::cout<<"# L2-Norm of source vector = "<<std::setprecision(7)<<sourcenorm<<std::endl;
    std::cout<<"########################################################"<<std::endl;
  }

  // write source_ to rhs
  AddContributionToRHS(source_);

  return;
}


/*------------------------------------------------------------------------*
 |  manipulations of scatra algorithm after call to Solve()   rauch 08/16 |
 *------------------------------------------------------------------------*/
void SCATRA::TimIntOneStepThetaEndoExocytosis::PostSolve()
{
  // temp variables
  int dof_id_internalized_scalar = discret_->GetCondition("ScaTraCellExt")->GetInt("ScalarID");
  int numscal = NumScalInCondition(*(discret_->GetCondition("ScaTraCellInt")));

  /////////////////////////////////////////////////////////////////////////////////////
  // clear state and set new state phinp
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);

  // set action for elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",SCATRA::calc_total_and_mean_scalars);
  eleparams.set<bool>("inverting",false);

  // provide displacement field in case of ALE
  if (isale_)
    eleparams.set<int>("ndsdisp",nds_disp_);

  // evaluate integrals of concentrations and domain
  Teuchos::RCP<Epetra_SerialDenseVector> scalars =
      Teuchos::rcp(new Epetra_SerialDenseVector(numscal+1));
  discret_->EvaluateScalars(eleparams, scalars, "ScaTraCellInt");

  // clear all states
  discret_->ClearState();

  /////////////////////////////////////////////////////////////////////////////////////
  // calculate mean concentrations
  mean_scalars_.resize(numscal);
  if (std::abs((*scalars)[numscal]) < 1E-15)
    dserror("Domain has zero volume!\n"
            "Something went wrong.\n"
            "Check ---DESIGN TOTAL AND MEAN SCALAR VOL CONDITIONS in your .dat file!");

  // loop over all scalars and store mean scalars
  for(int k=0; k<numscal; ++k)
  {
    mean_scalars_[k] = (*scalars)[k];
  }

  // entry associated with current time step 'step_'
  int entry = -1234;
  if(step_>internalization_steps_)
  {
    // Member step_ is incremented before Solve. Thus, it starts with 1 for the first time step.
    entry =((step_-1)-internalization_steps_)%(internalization_steps_);
  }
  else
  {
    // Member step_ is incremented before Solve. Thus, it starts with 1 for the first time step.
    entry =((step_-1))%(internalization_steps_);
  }

  // in case we write to the first vector entry, the previous value is the last vector entry
  double previousvectorentry=-1234.0;
  if (entry==0)
    previousvectorentry=internalization_vec_[internalization_steps_-1];
  else
    previousvectorentry=internalization_vec_[entry-1];

  // save current amount of internalized scalar at correct position in vector belonging to time-delay
  internalization_vec_[entry] = mean_scalars_[dof_id_internalized_scalar];

  // Do in first time step
  if(step_ == 1)
  {
    // The first delta equals the amount of internalized biomolecules, since
    // initially we assume zero internalized biomolcules.
    Delta_phi_[0]=internalization_vec_[0];
  }
  else
  {
    // calc the difference between the internalized biomolecules in the current and the previous step
    Delta_phi_[entry] = ( internalization_vec_[entry] - previousvectorentry );
  }

  if(discret_->Comm().MyPID()==0)
  {
    std::cout<<"########################################################"<<std::endl;
    std::cout<<"# Internalized "<<std::setprecision(7)<<Delta_phi_[entry]<<" molecules of species "<< dof_id_internalized_scalar<<std::endl;
    std::cout<<"########################################################"<<std::endl;
  }

  return;

}


/*----------------------------------------------------------------------*
 | Destructor dtor                                 (public) rauch 08/16 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntOneStepThetaEndoExocytosis::~TimIntOneStepThetaEndoExocytosis()
{
  return;
}
