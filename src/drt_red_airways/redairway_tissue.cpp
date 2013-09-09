/*!----------------------------------------------------------------------
\file redairway_tissue.cpp
\brief Control routine for coupled reduced airways and continuum tissue models

<pre>
Maintainer: Lena Yoshihara
            yoshihara@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*----------------------------------------------------------------------*/


#include <stdio.h>
//#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "redairway_tissue.H"
#include "airwayimplicitintegration.H"
#include "../drt_adapter/ad_str_redairway.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"

/*----------------------------------------------------------------------*
 |  Constructor (public)                                 yoshihara 09/12|
 *----------------------------------------------------------------------*/
AIRWAY::RedAirwayTissue::RedAirwayTissue(const Epetra_Comm& comm,
                                         const Teuchos::ParameterList& timeparams)
  : ADAPTER::AlgorithmBase(comm,timeparams)
{
  // before setting up the structure time integrator, manipulate coupling conditions -> turn them
  // into neumann orthopressure conditions

  std::vector<DRT::Condition*> surfneumcond;
  std::vector<int> tmp;
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
  if (structdis == Teuchos::null)
    dserror("no structure discretization available");

  // first get all Neumann conditions on structure
  structdis->GetCondition("SurfaceNeumann",surfneumcond);
  unsigned int numneumcond = surfneumcond.size();
  if (numneumcond == 0) dserror("no Neumann conditions on structure");

  // now filter those Neumann conditions that are due to the coupling
  std::vector<DRT::Condition*> coupcond;
  for (unsigned int i = 0; i < numneumcond; ++i)
  {
    DRT::Condition* actcond = surfneumcond[i];
    if (actcond->Type() == DRT::Condition::RedAirwayTissue)
      coupcond.push_back(actcond);
  }
  unsigned int numcond = coupcond.size();
  if (numcond == 0) dserror("no coupling conditions found");

  for (unsigned int i = 0; i < numcond; ++i)
  {
    DRT::Condition* cond = coupcond[i];
    cond->Add("type","neum_orthopressure");
    std::vector<int> onoff(6,0);
    onoff[0] = 1;
    cond->Add("onoff",onoff);
    std::vector<double> val(6,0.0);
    cond->Add("val",val);

    int condID = (coupcond[i])->GetInt("coupling id");
    tmp.push_back(condID);
  }

  std::vector<DRT::Condition*> nodecond;
  Teuchos::RCP<DRT::Discretization> redairwaydis = DRT::Problem::Instance()->GetDis("red_airway");
  if (redairwaydis == Teuchos::null)
    dserror("no redairway discretization available");

  // first get all redairway prescribed conditions on structure
  redairwaydis->GetCondition("RedAirwayPrescribedCond",nodecond);
  unsigned int numnodecond = nodecond.size();
  if (numnodecond == 0) dserror("no redairway prescribed conditions on redairway discretization");

  // now filter those node conditions that are due to the coupling
  std::vector<DRT::Condition*> nodecoupcond;
  for (unsigned int i = 0; i < numnodecond; ++i)
  {
    DRT::Condition* actcond = nodecond[i];
    if (actcond->Type() == DRT::Condition::RedAirwayNodeTissue)
      nodecoupcond.push_back(actcond);
  }
  unsigned int numnodecoupcond = nodecoupcond.size();
  if (numnodecoupcond == 0) dserror("no coupling conditions found");

  for (unsigned int i = 0; i < numnodecoupcond; ++i)
  {
    DRT::Condition* cond = nodecoupcond[i];
    cond->Add("boundarycond","flow");
    std::vector<double> val(1,0.0);
    cond->Add("val",val);
  }


  Epetra_Map redundantmap(tmp.size(),tmp.size(),&tmp[0],0,comm);
  couppres_ip_ = Teuchos::rcp(new Epetra_Vector(redundantmap, true));
  couppres_im_ = Teuchos::rcp(new Epetra_Vector(redundantmap, true));
  coupflux_ip_ = Teuchos::rcp(new Epetra_Vector(redundantmap, true));
  coupflux_im_ = Teuchos::rcp(new Epetra_Vector(redundantmap, true));
  coupvol_ip_ = Teuchos::rcp(new Epetra_Vector(redundantmap, true));
  coupvol_im_ = Teuchos::rcp(new Epetra_Vector(redundantmap, true));

  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure = Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(sdyn, const_cast<Teuchos::ParameterList&>(sdyn), structdis));
  structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::StructureRedAirway>(structure->StructureFieldrcp());

  SetupRedAirways();
  const Teuchos::ParameterList& rawdyn   = DRT::Problem::Instance()->ReducedDAirwayDynamicParams();

  // check certain parameters
  if (sdyn.get<double>("TIMESTEP") != timeparams.get<double>("TIMESTEP") or
      sdyn.get<int>("NUMSTEP") != timeparams.get<int>("NUMSTEP") or
      sdyn.get<double>("MAXTIME") != timeparams.get<double>("MAXTIME") or
      rawdyn.get<double>("TIMESTEP") != timeparams.get<double>("TIMESTEP") or
      rawdyn.get<int>("NUMSTEP") != timeparams.get<int>("NUMSTEP"))
    dserror("Parameter(s) for time integrators inconsistent");


  // get coupling parameters
  const Teuchos::ParameterList& rawtisdyn   = DRT::Problem::Instance()->RedAirwayTissueDynamicParams();
  // get max iterations
  itermax_ = rawtisdyn.get<int>("MAXITER");

  // get tolerance for pressure
  tolp_ = rawtisdyn.get<double>("CONVTOL_P");

  // get tolerance for flux
  tolq_ = rawtisdyn.get<double>("CONVTOL_Q");;
  
  // get normal direction 
  // -> if normal == 1.0 : the pressure will be implimented from inside the element to the outside
  // -> if normal ==-1.0 : the pressure will be implimented from outside the element to the inside
  normal_ = rawtisdyn.get<double>("NORMAL");
  // determine initial volume
  structure_->InitVol();
}


/*----------------------------------------------------------------------*
 |  Read restart                                         yoshihara 09/12|
 *----------------------------------------------------------------------*/
void AIRWAY::RedAirwayTissue::ReadRestart(int step)
{
  dserror("not implemented yet");
}


/*----------------------------------------------------------------------*
 |  Integrate                                            yoshihara 09/12|
 *----------------------------------------------------------------------*/
void AIRWAY::RedAirwayTissue::Integrate()
{
  while (NotFinished())
  {
    int iter=0;
    IncrementTimeAndStep();
    redairways_->PrepareTimeStep();

    do
    {
      DoRedAirwayStep();
      DoStructureStep();
      iter++;
    }while (NotConverged(iter)&&iter<itermax_);

    if (iter>= itermax_)
    {
      dserror("FIELD ITERATION NOT CONVERGED IN %d STEPS AT TIME T=%f",itermax_,Time() );
    }

    UpdateAndOutput();
  }
}


/*----------------------------------------------------------------------*
 |  Output of one iteration between fields               yoshihara 09/12|
 *----------------------------------------------------------------------*/
void AIRWAY::RedAirwayTissue::OutputIteration(Teuchos::RCP<Epetra_Vector> pres_inc, Teuchos::RCP<Epetra_Vector> scaled_pres_inc,
    Teuchos::RCP<Epetra_Vector> flux_inc, Teuchos::RCP<Epetra_Vector> scaled_flux_inc, int iter)
{
  if (couppres_ip_->Comm().MyPID() == 0)
  {
    printf("\nFIELD ITERATION: %i / %i\n", iter,itermax_);
    printf(" Tolerances:                                                                        %4.2e      %4.2e\n", tolp_, tolq_);
    printf(" Volume ID      Vol            P             Q            dP            dQ           dP_scal       dQ_scal\n");
    for (int i=0; i<couppres_ip_->Map().NumMyElements(); ++i)
    {
      printf("     %d       %4.3e     %4.3e     %4.3e     %4.3e     %4.3e     %4.3e     %4.3e\n",
          couppres_ip_->Map().GID(i), (*coupvol_ip_)[i], (*couppres_ip_)[i],  (*coupflux_ip_)[i], (*pres_inc)[i],
          (*flux_inc)[i], (*scaled_pres_inc)[i], (*scaled_flux_inc)[i]);
    }
    printf("----------------------------------------------------------------------------------------------------------\n\n");
  }
}


/*----------------------------------------------------------------------*
 |  Update and output                                    yoshihara 09/12|
 *----------------------------------------------------------------------*/
void AIRWAY::RedAirwayTissue::UpdateAndOutput()
{
  structure_->PrepareOutput();
  structure_->Update();
  structure_->Output();

  redairways_->TimeUpdate();
  redairways_->Output();
}


/*----------------------------------------------------------------------*
 |  Integrate airways time step and calculate pressures  yoshihara 09/12|
 *----------------------------------------------------------------------*/
void AIRWAY::RedAirwayTissue::DoRedAirwayStep()
{
  // scale with -1 (redairway convention: outflow is negative)
  coupflux_ip_->Scale(-1.0);

  redairways_->SetAirwayFluxFromTissue(coupflux_ip_);
  redairways_->IntegrateStep();
  redairways_->ExtractPressure(couppres_ip_);
  couppres_ip_->Update(0.0,*couppres_ip_,normal_);
}

/*----------------------------------------------------------------------*
 |  Integrate structure time step and calculate fluxes   yoshihara 09/12|
 *----------------------------------------------------------------------*/
void AIRWAY::RedAirwayTissue::DoStructureStep()
{
  structure_->SetPressure(couppres_ip_);
  structure_->IntegrateStep();
  structure_->CalcFlux(coupflux_ip_,coupvol_ip_,Dt());
 
}


/*----------------------------------------------------------------------*
 |  Check for convergence between fields                 yoshihara 09/12|
 *----------------------------------------------------------------------*/
/* Note: This has to be done for each field on its own, not as a Norm2
 * over all fields.
 */
bool AIRWAY::RedAirwayTissue::NotConverged(int iter)
{
  Teuchos::RCP<Epetra_Vector> pres_inc = Teuchos::rcp(new Epetra_Vector(*couppres_ip_));
  Teuchos::RCP<Epetra_Vector> scaled_pres_inc = Teuchos::rcp(new Epetra_Vector(*couppres_ip_));
  Teuchos::RCP<Epetra_Vector> flux_inc = Teuchos::rcp(new Epetra_Vector(*coupflux_ip_));
  Teuchos::RCP<Epetra_Vector> scaled_flux_inc = Teuchos::rcp(new Epetra_Vector(*coupflux_ip_));

  /*elegant solution not fully finished...
  pres_inc->Update(1.0,*couppres_ip_, -1.0, *couppres_im_, 0.0);
  scaled_pres_inc->ReciprocalMultiply(1.0,*couppres_ip_,*pres_inc, 0.0);

  flux_inc->Update(1.0,*coupflux_ip_, -1.0, *coupflux_im_, 0.0);
  scaled_flux_inc->ReciprocalMultiply(1.0,*coupflux_ip_,*flux_inc, 0.0);
  */

  //Calculate Pressure Norm
  for (int i=0; i<couppres_ip_->Map().NumMyElements(); ++i)
  {
        //Calculate pressure increment
        (*pres_inc)[i] = abs((*couppres_ip_)[i] - (*couppres_im_)[i]);

        //Calculate scaled pressure increment
        if(abs((*couppres_ip_)[i]) > 1e-05)
          (*scaled_pres_inc)[i] = (*pres_inc)[i] / abs((*couppres_ip_)[i]);
        else
          (*scaled_pres_inc)[i] = (*pres_inc)[i];
   }

  //Calculate Flux Norm
  for (int i=0; i<coupflux_ip_->Map().NumMyElements(); ++i)
  {
          //Calculate flux increment
          (*flux_inc)[i] = abs((*coupflux_ip_)[i] - (*coupflux_im_)[i]);

          //Calculate scaled flux increment
          if(abs((*coupflux_ip_)[i]) > 1e-05)
            (*scaled_flux_inc)[i] = (*flux_inc)[i] / abs((*coupflux_ip_)[i]);
          else
            (*scaled_flux_inc)[i] = (*flux_inc)[i];
  }

  //Output
  OutputIteration(pres_inc, scaled_pres_inc, flux_inc, scaled_flux_inc, iter);

  //Update values
  couppres_im_->Update(1.0,*couppres_ip_,0.0);
  coupflux_im_->Update(1.0,*coupflux_ip_,0.0);
  coupvol_im_->Update(1.0,*coupvol_ip_,0.0);

  double pres_max, flux_max;
  scaled_pres_inc->NormInf(&pres_max);
  scaled_flux_inc->NormInf(&flux_max);

  if (pres_max < tolp_ and flux_max < tolq_ and iter > 1)
    return false;

  return true;
}


void AIRWAY::RedAirwayTissue::SetupRedAirways()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis("red_airway");

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled())
  {
    actdis->FillComplete();
  }

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RCP<IO::DiscretizationWriter>  output = Teuchos::rcp( new IO::DiscretizationWriter(actdis),false);
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  // const Teuchos::ParameterList& probsize = DRT::Problem::Instance()->ProblemSizeParams();
  //  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& rawdyn   = DRT::Problem::Instance()->ReducedDAirwayDynamicParams();

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // get the solver number
  const int linsolvernumber = rawdyn.get<int>("LINEAR_SOLVER");
  // check if the present solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined. Please set LINEAR_SOLVER in REDUCED DIMENSIONAL AIRWAYS DYNAMIC to a valid number!");
  RCP<LINALG::Solver> solver = Teuchos::rcp( new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                                                   actdis->Comm(),
                                                   DRT::Problem::Instance()->ErrorFile()->Handle()),
                                    false);
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  ParameterList airwaystimeparams;

  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  const int ndim = DRT::Problem::Instance()->NDim();
  airwaystimeparams.set<int>              ("number of degrees of freedom" ,1*ndim);

  // -------------------------------------------------- time integration
  // the default time step size
  airwaystimeparams.set<double>           ("time step size"           ,rawdyn.get<double>("TIMESTEP"));
  // maximum number of timesteps
  airwaystimeparams.set<int>              ("max number timesteps"     ,rawdyn.get<int>("NUMSTEP"));

  // ----------------------------------------------- restart and output
  // restart
  airwaystimeparams.set                  ("write restart every"       ,rawdyn.get<int>("RESTARTEVRY"));
  // solution output
  airwaystimeparams.set                  ("write solution every"      ,rawdyn.get<int>("UPRES"));

  // ----------------------------------------------- solver parameters
  // solver type
  airwaystimeparams.set                  ("solver type"             ,rawdyn.get<std::string>("SOLVERTYPE"));
  // tolerance
  airwaystimeparams.set                  ("tolerance"               ,rawdyn.get<double>("TOLERANCE"));
  // maximum number of iterations
  airwaystimeparams.set                  ("maximum iteration steps" ,rawdyn.get<int>("MAXITERATIONS"));

  airwaystimeparams.set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

  //------------------------------------------------------------------
  // create all vectors and variables associated with the time
  // integration (call the constructor);
  // the only parameter from the list required here is the number of
  // velocity degrees of freedom
  //------------------------------------------------------------------
  redairways_ = Teuchos::rcp(new AIRWAY::RedAirwayImplicitTimeInt(actdis,*solver,airwaystimeparams,*output));

  redairways_->SetupForCoupling();
}


