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

  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure = Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(sdyn));
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
    IncrementTimeAndStep();
    redairways_->PrepareTimeStep();

    do
    {
      DoRedAirwayStep();
      DoStructureStep();
    }while (NotConverged());

    UpdateAndOutput();
  }
}


/*----------------------------------------------------------------------*
 |  Output of one iteration between fields               yoshihara 09/12|
 *----------------------------------------------------------------------*/
void AIRWAY::RedAirwayTissue::OutputIteration(double pres_inc_norm, double flux_inc_norm)
{
  if (couppres_ip_->Comm().MyPID() == 0)
  {
    cout << "-------------------------  FIELD ITERATION ---------------------------" << endl;
    for (int i=0; i<couppres_ip_->Map().NumMyElements(); ++i)
    {
      cout << "\t time:\t" << Time() << "\t ID:\t" << couppres_ip_->Map().GID(i) << "\t P:\t" <<  (*couppres_ip_)[i]
           << "\t Q:\t" <<  (*coupflux_ip_)[i] << "\t DP2:\t" << pres_inc_norm << "\t DQ2:\t" << flux_inc_norm << endl;
    }
    cout << "---------------------------------------------------------------------" << endl;
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
}

/*----------------------------------------------------------------------*
 |  Integrate structure time step and calculate fluxes   yoshihara 09/12|
 *----------------------------------------------------------------------*/
void AIRWAY::RedAirwayTissue::DoStructureStep()
{
  structure_->SetPressure(couppres_ip_);
  structure_->IntegrateStep();
  structure_->CalcFlux(coupflux_ip_,Dt());
}


/*----------------------------------------------------------------------*
 |  Check for convergence between fields                 yoshihara 09/12|
 *----------------------------------------------------------------------*/
bool AIRWAY::RedAirwayTissue::NotConverged()
{
  Teuchos::RCP<Epetra_Vector> pres_inc = Teuchos::rcp(new Epetra_Vector(*couppres_ip_));
  pres_inc->Update(-1.0,*couppres_im_,1.0);
  double pres_inc_norm;
  pres_inc->Norm2(&pres_inc_norm);

  Teuchos::RCP<Epetra_Vector> flux_inc = Teuchos::rcp(new Epetra_Vector(*coupflux_ip_));
  flux_inc->Update(-1.0,*coupflux_im_,1.0);
  double flux_inc_norm;
  flux_inc->Norm2(&flux_inc_norm);

  double tol = 1e-6;

  OutputIteration(pres_inc_norm, flux_inc_norm);

  couppres_im_->Update(1.0,*couppres_ip_,0.0);
  coupflux_im_->Update(1.0,*coupflux_ip_,0.0);
//  couppres_ip_->PutScalar(0.0);
//  coupflux_ip_->PutScalar(0.0);

  if (pres_inc_norm < tol and flux_inc_norm < tol)
    return false;

  return true;
}


void AIRWAY::RedAirwayTissue::SetupRedAirways()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
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
  RCP<IO::DiscretizationWriter>  output = rcp( new IO::DiscretizationWriter(actdis),false);
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
  RCP<LINALG::Solver> solver = rcp( new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
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
  airwaystimeparams.set                  ("solver type"             ,rawdyn.get<string>("SOLVERTYPE"));
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


