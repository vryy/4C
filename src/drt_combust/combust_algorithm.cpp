/*!----------------------------------------------------------------------*
\file combust_algorithm.cpp

\brief base combustion algorithm

	detailed description in header file combust_algorithm.H
	
<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "combust_algorithm.H"
#include "combust_utils.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 * constructor                                              henke 06/08 * 
 *----------------------------------------------------------------------*/
COMBUST::Algorithm::Algorithm(Epetra_Comm& comm)
:  FluidBaseAlgorithm(DRT::Problem::Instance()->CombustionDynamicParams(),false),
   ScaTraBaseAlgorithm(DRT::Problem::Instance()->CombustionDynamicParams()),
   comm_(comm),
   step_(0),
   time_(0.0)
{
  /* Der constructor sollte den gesamten Algorithmus initialisiern.
   * Das heisst:
   * o G-Funktionsvektor (t und ig+1) auf initial value setzen
   * o Geschwindigkeitsvektor (t und iu+1) auf initial value setzen
   * o alle Zähler auf 0 setzen (step_(0), f_giter_(0), g_iter_(0), f_iter_(0))
   * o alle Normen und Grenzwerte auf 0 setzen
  */ 

  // taking time loop control parameters out of fluid dynamics section
  const Teuchos::ParameterList& combustdyn = DRT::Problem::Instance()->CombustionDynamicParams();
  // maximum simulation time
  timemax_=combustdyn.get<double>("MAXTIME");
  // maximum number of timesteps
  stepmax_ = combustdyn.get<int>("NUMSTEP");
  // time step size
  dt_ = combustdyn.get<double>("TIMESTEP");
  
  // counter FGI iterations
  fgiter_ = 0;
  // maximum number of Fluid - G-function iterations
  fgitermax_ = combustdyn.get<int>("ITEMAX");
 
  // set options for reinitialization; they should be read fro the parameter list
  signeddistfunc = 1;
  reinitialize_action = 1; // reinitialize_action = signeddistfunc;

  // velocitynp_ ist nicht initialisiert!
}

/*----------------------------------------------------------------------*
 * destructor                                               henke 06/08 * 
 *----------------------------------------------------------------------*/
COMBUST::Algorithm::~Algorithm()
{
}

/*----------------------------------------------------------------------*
 * public: time loop of combustion algorithm                henke 06/08 *
 *----------------------------------------------------------------------*/
void COMBUST::Algorithm::TimeLoop()
{
    // solve con-dif equation without coupling to Navier-Stokes
    // ConDifField().SetVelocityField(1,1);
    // ConDifField().Integrate();

  // time loop
  while (NotFinished())
  {
    // prepare next time step
    PrepareTimeStep();
    
    //create integration cells
    CreateIntegrationCells();

    // reinitialize G-function 
    ReinitializeGfunc();
    
    // Fluid-G-function-Interaction loop
    while (NotConvergedFGI())
    {
    	// prepare Fluid-G-function iteration
    	PrepareFGIteration();
    	
    	// solve nonlinear Navier-Stokes system
    	DoFluidField();

    	// solve linear G-function equation
    	DoGfuncField();
    	
    	// update field vectors
    	UpdateFGIteration();
    	
    } // Fluid-G-function-Interaction loop
    
    // update all field solvers
    UpdateTimeStep();

    // write output to screen and files
    Output();

  } // time loop

return;
} // TimeLoop()

/*---------------------------------------------------------------------*
 *  protected: reinitialize G-function                     henke 06/08 *
 *---------------------------------------------------------------------*/
void COMBUST::Algorithm::ReinitializeGfunc()
{
	/* Here, the G-function is reinitialized, because we suspect that the
	 * signed distance property has been lost due to inaccuracies. There
	 * are various options to reinitialize the G-function.
	 * For the time being we use an exact, but expensive, procedure to ensure
	 * this property (SignedDistFunc()).                    henke 06/08 */
	switch (reinitialize_action)
	{
		case 1:
			SignedDistFunc();
			break;
		default:
			dserror ("Unknown option to reinitialize the G-function");
	}
	return;
}

/*---------------------------------------------------------------------*
 *  protected: build signed distance function              henke 06/08 *
 *---------------------------------------------------------------------*/
void COMBUST::Algorithm::SignedDistFunc()
{
	/* This member function constructs a G-function field that meets the 
	 * signed distance property. The algorithm assigns the value of the 
	 * distance between each node and the surface defined by G=0 as a scalar 
	 * value to every node in the G-function discretization.*/
	return;
}

/*----------------------------------------------------------------------*
 *  protected: FGI iteration converged?                     henke 06/08 *
 *----------------------------------------------------------------------*/
bool COMBUST::Algorithm::NotConvergedFGI()
{
	bool notconverged = false;
	if (fgiter_ < fgitermax_ and true) // (fgiter <= fgitermax and ComputeGfuncNorm() < maxepsg and ComputeFluidNorm() < maxepsf)
	{
		notconverged = true;
	}
	return notconverged; // return notconverged;
}

/*----------------------------------------------------------------------*
 * protected: prepare time step for combustion algorithm    henke 06/08 *
 *----------------------------------------------------------------------*/
void COMBUST::Algorithm::PrepareTimeStep()
{
  step_ += 1;
  time_ += dt_;
  fgiter_ = 0;
  // fgnormgfunc = large value, determined in Input file
  fgvelnormL2_ = 1.0;
  // fgnormfluid = large value
  fggfuncnormL2_ = 1.0;
    
  if (Comm().MyPID()==0)
  {
	//cout<<"---------------------------------------  time step  ------------------------------------------\n";
	printf("---------------------------------------  time step %2d ----------------------------------------\n",step_);
    printf("TIME: %11.4E/%11.4E  DT = %11.4E STEP = %4d/%4d \n",time_,timemax_,dt_,step_,stepmax_);
  }
  
  FluidField().PrepareTimeStep();
  ScaTraField().PrepareTimeStep();
  
  /* Ich wuerde hier gerne vergleichen, ob die Zeitschritte des Fluids und des ConDif
   * identisch sind. Es existiert eine Variable time_ in jedem feld */
  return;
}

/*----------------------------------------------------------------------*
 * protected: prepare time step for combustion algorithm    henke 06/08 *
 *----------------------------------------------------------------------*/
void COMBUST::Algorithm::PrepareFGIteration()
{
	fgiter_ += 1;
	if (Comm().MyPID()==0)
	{
	  //cout<<"\n---------------------------------------  FGI loop  -------------------------------------------\n";
	  printf("\n---------------------------------------  FGI loop: iteration number: %2d ----------------------\n",fgiter_);
	}	
}

/*----------------------------------------------------------------------*
 * protected: perform a fluid time integration step         henke 06/08 *
 *----------------------------------------------------------------------*/
void COMBUST::Algorithm::DoFluidField()
{
  // determine location of flame interface
  // FlameInterface(Gfunction);	
  
  if (Comm().MyPID()==0)
  {
    cout<<"\n---------------------------------------  FLUID SOLVER  ---------------------------------------\n";
  }
  
  // solve nonlinear Navier-Stokes system
  FluidField().NonlinearSolve();
  return;
}

/*----------------------------------------------------------------------*
 * protected: perform a G-function time integration step    henke 06/08 *
 *----------------------------------------------------------------------*/
void COMBUST::Algorithm::DoGfuncField()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n---------------------------------------  G-FUNCTION SOLVER  ----------------------------------\n";
  }
  
  /* Hier muss eigentlich nur die Geschwindigkeit bei G=0 (Flammeninterface gelesen 
     werden). Tatsächlich kann vorerst aber das gesamte Navier-Stokes Fluid Feld als Input
     genommen werden. Die G-function sucht sich dann an jedem Knoten die Geschwindigkeit
     als Konvektions-Diffusions Geschwindigkeit.*/

  // exract new velocity from fluid field as provided by the Navier-Stokes solver
  velocitynp_=FluidField().ExtractVelocityPart(FluidField().Velnp());

  // assign the fluid velocity to the G-function field as convective velocity
  ScaTraField().SetVelocityField(2,velocitynp_);

  // solve convection-diffusion equation
  ScaTraField().Solve();
  /* Hier muss später eine nichtlineare G-Funktion gelöst werden. Diese function existiert aber noch 
   * nicht im condifimplicitintegration.cpp */
  //ConDifField().NonlinearSolve();
  return;
}

/*----------------------------------------------------------------------*
 * protected: update                                        henke 06/08 *
 *----------------------------------------------------------------------*/
void COMBUST::Algorithm::UpdateFGIteration()
{
	// update the Fluid and the FGI vector at the end of the FGI loop
	return;
}

/*----------------------------------------------------------------------*
 * protected: update                                        henke 06/08 *
 *----------------------------------------------------------------------*/
void COMBUST::Algorithm::UpdateTimeStep()
{
  FluidField().Update();
  ScaTraField().Update();
  return;
}

/*----------------------------------------------------------------------*
 * protected: output                                        henke 06/08 *
 *----------------------------------------------------------------------*/
void COMBUST::Algorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  FluidField().Output();
  FluidField().LiftDrag();
  ScaTraField().Output();

  // debug IO
#if 0
  // print out convective velocity vector
  cout<<*velocitynp_<<endl;
#endif

  return;
}

#endif // #ifdef CCADISCRET
