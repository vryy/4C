/*!----------------------------------------------------------------------
\file redairway_tissue.cpp
\brief Control routine for coupled reduced airways and continuum tissue models
\level 2
<pre>
\maintainer Lena Yoshihara
            yoshihara@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*----------------------------------------------------------------------*/


#include <stdio.h>

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
  //Before setting up the structure time integrator, manipulate coupling conditions -> turn them
  //into neumann orthopressure conditions
  std::vector<DRT::Condition*> surfneumcond;
  std::vector<int> tmp;
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
  if (structdis == Teuchos::null)
    dserror("no structure discretization available");

  //First get all Neumann conditions on structure
  structdis->GetCondition("SurfaceNeumann",surfneumcond);
  unsigned int numneumcond = surfneumcond.size();
  if (numneumcond == 0) dserror("no Neumann conditions on structure");

  //Now filter those Neumann conditions that are due to the coupling
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

  //First get all redairway prescribed conditions on structure
  redairwaydis->GetCondition("RedAirwayPrescribedCond",nodecond);
  unsigned int numnodecond = nodecond.size();
  if (numnodecond == 0) dserror("no redairway prescribed conditions on redairway discretization");

  //Now filter those node conditions that are due to the coupling
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
  couppres_ip_tilde_ = Teuchos::rcp(new Epetra_Vector(redundantmap, true));
  couppres_im_ = Teuchos::rcp(new Epetra_Vector(redundantmap, true));
  couppres_im_tilde_ = Teuchos::rcp(new Epetra_Vector(redundantmap, true));
  couppres_il_ = Teuchos::rcp(new Epetra_Vector(redundantmap, true));
  omega_np_ = Teuchos::rcp(new Epetra_Vector(redundantmap, true));
  coupflux_ip_ = Teuchos::rcp(new Epetra_Vector(redundantmap, true));
  coupflux_im_ = Teuchos::rcp(new Epetra_Vector(redundantmap, true));
  coupvol_ip_ = Teuchos::rcp(new Epetra_Vector(redundantmap, true));
  coupvol_im_ = Teuchos::rcp(new Epetra_Vector(redundantmap, true));

  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure = Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(sdyn, const_cast<Teuchos::ParameterList&>(sdyn), structdis));
  structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::StructureRedAirway>(structure->StructureField());
  structure_->Setup();

  SetupRedAirways();
  const Teuchos::ParameterList& rawdyn   = DRT::Problem::Instance()->ReducedDAirwayDynamicParams();

  //Check Time integration parameters
  if (sdyn.get<double>("TIMESTEP") != timeparams.get<double>("TIMESTEP") or
      sdyn.get<int>("NUMSTEP") != timeparams.get<int>("NUMSTEP") or
      sdyn.get<double>("MAXTIME") != timeparams.get<double>("MAXTIME") or
      rawdyn.get<double>("TIMESTEP") != timeparams.get<double>("TIMESTEP") or
      rawdyn.get<int>("NUMSTEP") != timeparams.get<int>("NUMSTEP"))
    dserror("Parameter(s) for time integrators inconsistent");

  //Check Time integration parameters
  if (sdyn.get<int>("RESTARTEVRY") != rawdyn.get<int>("RESTARTEVRY"))
    dserror("Parameters for restart inconsistent");

  //Get coupling parameters
  const Teuchos::ParameterList& rawtisdyn   = DRT::Problem::Instance()->RedAirwayTissueDynamicParams();
  //Get max iterations
  itermax_ = rawtisdyn.get<int>("MAXITER");

  //Get tolerance for pressure
  tolp_ = rawtisdyn.get<double>("CONVTOL_P");

  //Get tolerance for flux
  tolq_ = rawtisdyn.get<double>("CONVTOL_Q");

  //Dynamic relaxation type
  relaxtype_ = DRT::INPUT::IntegralValue<INPAR::ARTNET::RELAXTYPE_3D_0D>(rawtisdyn,"RELAXTYPE");

  //Get normal direction
  // -> if normal == 1.0 : the pressure will be implemented from inside the element to the outside
  // -> if normal ==-1.0 : the pressure will be implemented from outside the element to the inside
  normal_ = rawtisdyn.get<double>("NORMAL");

  //Determine initial volume
  structure_->InitVol();

} //end of constructor


/*----------------------------------------------------------------------*
 |  Read restart                                              roth 10/13|
 *----------------------------------------------------------------------*/
void AIRWAY::RedAirwayTissue::ReadRestart(const int step)
{
  structure_->ReadRestart(step);
  structure_->InitVol();
  redairways_->ReadRestart(step);

  //Read the coupling variables at restart
  redairways_->GetRestartReader(step)->ReadVector(couppres_im_,"couppres_im");
  redairways_->GetRestartReader(step)->ReadVector(coupflux_im_,"coupflux_im");
  redairways_->GetRestartReader(step)->ReadVector(coupflux_ip_,"coupflux_im");
  redairways_->GetRestartReader(step)->ReadVector(coupvol_im_,"coupvol_im");

  //Set timestep and time after restart
  double restartTime = redairways_->Time();
  SetTimeStep(restartTime, step);

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
      RelaxPressure(iter);
      DoStructureStep();
      iter++;
    }while (NotConverged(iter)&&iter<itermax_);

    if ((iter>= itermax_) && (couppres_ip_->Comm().MyPID() == 0))
    {
      dserror("FIELD ITERATION NOT CONVERGED IN %d STEPS AT TIME T=%f",itermax_,Time() );
    }

    UpdateAndOutput();
  }

}


/*----------------------------------------------------------------------*
 |  Integrate airways time step and calculate pressures  yoshihara 09/12|
 *----------------------------------------------------------------------*/
void AIRWAY::RedAirwayTissue::DoRedAirwayStep()
{
  //Scale with -1 (redairway convention: outflow is negative)
  coupflux_ip_->Scale(-1.0);

  redairways_->SetAirwayFluxFromTissue(coupflux_ip_);
  redairways_->IntegrateStep();
  redairways_->ExtractPressure(couppres_ip_tilde_);
  couppres_ip_tilde_->Update(0.0,*couppres_ip_tilde_,normal_);

}

/*----------------------------------------------------------------------*
 |  Relax the coupling pressure as described in Kuettler et al.:        |
 |  "Fixed-point-fluid-structure interaction solvers with dynamic       |
 |   relaxation" (2008)                                       roth 10/13|
 *----------------------------------------------------------------------*/
void AIRWAY::RedAirwayTissue::RelaxPressure(int iter)
{

  switch(relaxtype_)
  {
    case INPAR::ARTNET::norelaxation:
    {
      printf("No dynamic relaxation \n");
      couppres_ip_->Update(1.0,*couppres_ip_tilde_,0.0);
    }
    break;

    case INPAR::ARTNET::fixedrelaxation:
    {
      // Fixed \omega_i
      omega_ = 0.85;
      printf("Fixed Relaxation Parameter: %f \n",omega_);

      //Dynamic Relaxation formula (35)
      // p^{n+1}_{i+1} = \omega_i \tilde{p}^{n+1}_{i+1} + (1-\omega_i) p^{n+1}_{i}
      // where \tilde{p}^{n+1}_{i+1} = couppres_ip_
      //                 p^{n+1}_{i} = couppres_im_
      couppres_ip_->Update(omega_, *couppres_ip_tilde_, (1.0-omega_), *couppres_im_, 0.0);
    }
    break;

    case INPAR::ARTNET::Aitken:
    {
      // \omega_{i+1} = omega_np_
      // Kuettler (2008), Remark 10: Two previous steps are required to calculate the \omega,
      // thus the relaxation factor is fixed for the first two steps to \omega = 1.0
      if (iter < 2)
      {
        omega_np_->PutScalar(1.0);
      }
      else if (iter >= 2)
      {
        //Relaxation factor formula (41)
        // \omega_{i+1} = (p_{i} - p_{i+1}) / (p_{i} - \tilde{p}_{i+1} - p_{i+1} + \tilde{p}_{i+2})
        //      where p_{i} = couppres_il_
        //          p_{i+1} = couppres_im_
        //  \tilde{p}_{i+1} = couppres_im_tilde_
        //  \tilde{p}_{i+2} = couppres_ip_tilde_
        omega_np_->Update(1.0,*couppres_il_,-1.0, *couppres_im_,0.0);

        Teuchos::RCP<Epetra_Vector> denominator = Teuchos::rcp(new Epetra_Vector(*omega_np_));
        denominator->Update(-1.0,*couppres_im_tilde_,+1.0,*couppres_ip_tilde_,1.0);

        omega_np_->ReciprocalMultiply(1.0,*denominator,*omega_np_,0.0);

        //Safety check for \omega_i+1
        for (int i=0; i<couppres_ip_->Map().NumMyElements(); ++i)
        {
          if ((*omega_np_)[i] < 0.0)
            (*omega_np_)[i] = 0.0;
          else if ((*omega_np_)[i] > 1.0)
            (*omega_np_)[i] = 1.0;
        }
      }

      //Aitken Relaxation formula (35)
      // p^{n+1}_{i+1} = \omega_i \tilde{p}^{n+1}_{i+1} + (1-\omega_i) p^{n+1}_{i}
      for (int i=0; i<couppres_ip_->Map().NumMyElements(); ++i)
      {
        (*couppres_ip_)[i] = (*omega_np_)[i]*(*couppres_ip_tilde_)[i] + (1-(*omega_np_)[i])*(*couppres_im_)[i];
      }

      //Print relaxation factor \omega_np_
      if(couppres_ip_->Comm().MyPID() == 0)
      {
        printf("Aitken Relaxation: \n");
        for (int i=0; i<couppres_ip_->Map().NumMyElements(); ++i)
          std::cout << "omega_np_["<< i << "]: " << (*omega_np_)[i] << std::endl;
        std::cout <<  std::endl;
      }
    }
    break;

    // not implemented yet...
    case INPAR::ARTNET::SD:
    {
      dserror("Currently only two types of RELAXTYPE possible: norelaxation, fixedrelaxation, Aitken");
    }
    break;

    default:
    {
      dserror("Currently only two types of RELAXTYPE possible: norelaxation, fixedrelaxation, Aitken");
    }
    break;
  }

  /*/Print coupling pressure for validation of dynamic relaxation
  std::cout << "\nSingle couppres_ip_: " << std::endl;
  for (int i=0; i<couppres_ip_->Map().NumMyElements(); ++i)
  {
    printf(" Time: %f couppres_ip_: %6.3e \n",Time(), (*couppres_ip_)[i]);
  }*/

}


/*----------------------------------------------------------------------*
 |  Integrate structure time step and calculate fluxes   yoshihara 09/12|
 *----------------------------------------------------------------------*/
void AIRWAY::RedAirwayTissue::DoStructureStep()
{
  structure_->SetPressure(couppres_ip_);
  structure_->PrepareTimeStep();
  structure_->Solve();
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
  couppres_il_->Update(1.0,*couppres_im_,0.0);
  couppres_im_->Update(1.0,*couppres_ip_,0.0);
  couppres_im_tilde_->Update(1.0,*couppres_ip_tilde_,0.0);
  coupflux_im_->Update(1.0,*coupflux_ip_,0.0);
  coupvol_im_->Update(1.0,*coupvol_ip_,0.0);

  double pres_max, flux_max;
  scaled_pres_inc->NormInf(&pres_max);
  scaled_flux_inc->NormInf(&flux_max);

  if (pres_max < tolp_ and flux_max < tolq_ and iter > 1)
    return false;

  return true;
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

  //In case of restart write all coupling variables to restart file
  if (redairways_->Step()%uprestart_ == 0)
  {
    redairways_->GetOutputWriter().WriteVector("couppres_im",couppres_im_);
    redairways_->GetOutputWriter().WriteVector("coupflux_im",coupflux_im_);
    redairways_->GetOutputWriter().WriteVector("coupvol_im",coupvol_im_);
  }
}


void AIRWAY::RedAirwayTissue::SetupRedAirways()
{
  //Access the discretization
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis("red_airway");

  //Set degrees of freedom in the discretization
  if (!actdis->Filled())
  {
    actdis->FillComplete();
  }

  //Context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter>  output = actdis->Writer();
  output->WriteMesh(0,0.0);

  //Set some pointers and variables
  const Teuchos::ParameterList& rawdyn   = DRT::Problem::Instance()->ReducedDAirwayDynamicParams();

  //Create a solver
  // Get the solver number
  const int linsolvernumber = rawdyn.get<int>("LINEAR_SOLVER");
  // Check if the present solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined. Please set LINEAR_SOLVER in REDUCED DIMENSIONAL AIRWAYS DYNAMIC to a valid number!");
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp( new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                                                   actdis->Comm(),
                                                   DRT::Problem::Instance()->ErrorFile()->Handle()),
                                    false);
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  //Set parameters in list required for all schemes
  Teuchos::ParameterList airwaystimeparams;

  //Number of degrees of freedom
  const int ndim = DRT::Problem::Instance()->NDim();
  airwaystimeparams.set<int>              ("number of degrees of freedom" ,1*ndim);

  //Time integration
  // Default time step size
  airwaystimeparams.set<double>           ("time step size"           ,rawdyn.get<double>("TIMESTEP"));
  // Maximum number of timesteps
  airwaystimeparams.set<int>              ("max number timesteps"     ,rawdyn.get<int>("NUMSTEP"));

  // Restart and output
  // Restart
  airwaystimeparams.set                  ("write restart every"       ,rawdyn.get<int>("RESTARTEVRY"));
  // Solution output
  airwaystimeparams.set                  ("write solution every"      ,rawdyn.get<int>("RESULTSEVRY"));

  // Solver parameters
  // Solver type
  airwaystimeparams.set                  ("solver type"             ,rawdyn.get<std::string>("SOLVERTYPE"));
  // Tolerance
  airwaystimeparams.set                  ("tolerance"               ,rawdyn.get<double>("TOLERANCE"));
  // Maximum number of iterations
  airwaystimeparams.set                  ("maximum iteration steps" ,rawdyn.get<int>("MAXITERATIONS"));
  // solve scatra flag
  if (rawdyn.get<std::string>("SOLVESCATRA")=="yes")
    airwaystimeparams.set                  ("SolveScatra" ,true);
  else
    airwaystimeparams.set                  ("SolveScatra" ,false);
  //Adjust acini volume with pre-stress condition
  if (rawdyn.get<std::string>("CALCV0PRESTRESS")=="yes")
  {
    airwaystimeparams.set                  ("CalcV0PreStress" ,true);
    airwaystimeparams.set                  ("transpulmpress"          ,rawdyn.get<double>("TRANSPULMPRESS"));
  }
  else
    airwaystimeparams.set                  ("CalcV0PreStress" ,false);


  airwaystimeparams.set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

  //Get restart timestep
  uprestart_ = rawdyn.get<int>("RESTARTEVRY");

  //------------------------------------------------------------------
  // create all vectors and variables associated with the time
  // integration (call the constructor);
  // the only parameter from the list required here is the number of
  // velocity degrees of freedom
  //------------------------------------------------------------------
  redairways_ = Teuchos::rcp(new AIRWAY::RedAirwayImplicitTimeInt(actdis,*solver,airwaystimeparams,*output));

  redairways_->SetupForCoupling();
}
