/*!-----------------------------------------------------------------------------------------------*
\file combust_algorithm.cpp

\brief base combustion algorithm

    detailed description in header file combust_algorithm.H

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "combust_algorithm.H"
#include "combust_utils.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 * constructor                                              henke 06/08 * 
 *----------------------------------------------------------------------*/
COMBUST::Algorithm::Algorithm(Epetra_Comm& comm, Teuchos::ParameterList& combustdyn)
: ScaTraFluidCouplingAlgorithm(comm, combustdyn),
// initialize member variables
  fgiter_(0),
  fgitermax_(combustdyn.get<int>("ITEMAX")),
  reinitializationaction_(compute_signeddistancefunction) // combustdyn.get<ReinitializeAction>("REINITIALIZEGFUNC")
{
  /* Der constructor sollte den gesamten Algorithmus initialisieren.
   * Hier müssen alle Variablen, die den Einzelfeldern übergeordnet sind, initialisiert werden.
   * 
   * das heisst:
   * o G-Funktionsvektor (t und ig+1) auf initial value setzen
   * o Geschwindigkeitsvektor (t und iu+1) auf initial value setzen
   * o alle Zähler auf 0 setzen (step_(0), f_giter_(0), g_iter_(0), f_iter_(0))
   * o alle Normen und Grenzwerte auf 0 setzen
   * 
   * Zusammenfassend kann an sagen, dass alles was in der Verbrennungsrechnung vor der Zeitschleife 
   * passieren soll, hier passieren muss, weil die combust dyn gleich die Zeitschleife ruft.
  */ 

  /*----------------------------------------------------------------------------------------------*
   * initialize all data structures needed for the combustion algorithm
   * 
   * - initialize G-function by scalar function field
   * - capture the flame front and create interface geometry (triangulation)
   * - determine initial enrichment (DofManager wird bereits mit dem Element d.h. Diskretisierung angelegt)
   * - ...
   *----------------------------------------------------------------------------------------------*/

  // brauch ich das hier? Das ist ein relikt von früher    henke 10/08
  ReinitializeGfunc(initialize);
  // CaptureFlameFront();
  // DetermineEnrichment();
  
}

/*------------------------------------------------------------------------------------------------*
 | destructor                                                                         henke 06/08 | 
 *------------------------------------------------------------------------------------------------*/
COMBUST::Algorithm::~Algorithm()
{
}

/*------------------------------------------------------------------------------------------------*
 | public: time loop of algorithm for dynamic combustion problem                      henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::TimeLoop()
{

/*  
    if (Comm().MyPID()==0)
	{
		cout<<"\n Combustion Algorithmus Timeloop \n";
		cout<<"\n time step in timeloop: "<< Step() << &endl;
	}
*/

  // time loop
  while (NotFinished())
  {
    // prepare next time step
    PrepareTimeStep();

    //create integration cells
    CreateIntegrationCells();

    // reinitialize G-function 
    ReinitializeGfunc(reinitializationaction_);
    
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

/*------------------------------------------------------------------------------------------------*
 | public: algorithm for static combustion problem                                    henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::SolveStationaryProblem()
{
  dserror("Der Algorithmus für statische Verbrennungproblem kann noch nix!");
  return;
}

/*------------------------------------------------------------------------------------------------*
 |  protected: reinitialize G-function                                                henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::ReinitializeGfunc(ReinitializationAction action)
{
  /* Here, the G-function is reinitialized, because we suspect that the signed distance property has
   * been lost due to inaccuracies. There are various options to reinitialize the G-function.
   * For the time being we use an exact, but expensive, procedure to ensure this property which 
   * computes the distance for every point. (SignedDistFunc)                          henke 06/08 */
  switch (action)
  {
    case initialize:
      InitializeGFunc();
      break;
    case compute_signeddistancefunction:
      SignedDistFunc();
      break;
    default:
      dserror ("Unknown option to reinitialize the G-function");
  }
  return;
}

/*------------------------------------------------------------------------------------------------*
 |  protected: initialize G-function by scalar function field                         henke 09/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::InitializeGFunc()
{
  // This function reads a CURVE from the input file and initializes the G-function with its values
  return;
}

/*------------------------------------------------------------------------------------------------*
 |  protected: build signed distance function                                         henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::SignedDistFunc()
{
  /* This member function constructs a G-function field that meets the 
   * signed distance property. The algorithm assigns the value of the 
   * distance between each node and the surface defined by G=0 as a scalar 
   * value to every node in the G-function discretization.*/
  return;
}

/*------------------------------------------------------------------------------------------------*
 |  protected: FGI iteration converged?                                               henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
bool COMBUST::Algorithm::NotConvergedFGI()
{
  // if (fgiter <= fgitermax and ComputeGfuncNorm() < maxepsg and ComputeFluidNorm() < maxepsf)
  return (fgiter_ < fgitermax_ and true);
}

/*------------------------------------------------------------------------------------------------*
 | protected: prepare time step for combustion algorithm                              henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::PrepareTimeStep()
{
  IncrementTimeAndStep();
  fgiter_ = 0;
  // fgnormgfunc = large value, determined in Input file
  fgvelnormL2_ = 1.0;
  // fgnormfluid = large value
  fggfuncnormL2_ = 1.0;
    
  if (Comm().MyPID()==0)
  {
    //cout<<"---------------------------------------  time step  ------------------------------------------\n";
    printf("----------------------Combustion-------  time step %2d ----------------------------------------\n",Step());
    printf("TIME: %11.4E/%11.4E  DT = %11.4E STEP = %4d/%4d \n",Time(),MaxTime(),Dt(),Step(),NStep());
  }
  
  FluidField().PrepareTimeStep();

  // transfer the initial(!!) convective velocity
  //(fluid initial field was set inside the constructor of fluid base class)
  if (Step()==1) ScaTraField().SetVelocityField(2,ConvectiveVelocity());

  // prepare time step (+ initialize one-step-theta scheme correctly with 
  // velocity given above)
  ScaTraField().PrepareTimeStep();
  
  // synchronicity check between combust algorithm and base algorithms
  if (FluidField().Time() != Time())
    dserror("Time in Fluid differs from time in combustion algorithm");
  if (ScaTraField().Time() != Time())
    dserror("Time in ScaTra differs from time in combustion algorithm");

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: prepare time step for combustion algorithm                              henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::PrepareFGIteration()
{
  fgiter_ += 1;
  if (Comm().MyPID()==0)
  {
    //cout<<"\n---------------------------------------  FGI loop  -------------------------------------------\n";
    printf("\n---------------------------------------  FGI loop: iteration number: %2d ----------------------\n",fgiter_);
  }
}

/*------------------------------------------------------------------------------------------------*
 * protected: perform a fluid time integration step                                   henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::DoFluidField()
{
  // determine location of flame interface
  // FlameInterface(Gfunction);	
  
  if (Comm().MyPID()==0)
  {
    cout<<"\n---------------------------------------  FLUID SOLVER  ---------------------------------------\n";
  }
  
  // CaptureFlameFront();
  // DetermineEnrichment();

  // solve nonlinear Navier-Stokes system
  FluidField().NonlinearSolve();
  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: perform a G-function time integration step                              henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
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
  GetCurrentFluidVelocity();

  // assign the fluid velocity to the G-function field as convective velocity
  ScaTraField().SetVelocityField(2,ConvectiveVelocity());

  // solve convection-diffusion equation
  ScaTraField().Solve();
  /* Hier muss später eine nichtlineare G-Funktion gelöst werden. Diese function existiert aber noch 
   * nicht in einem ScaTra Zeitintegrationsverfahren*/
  //ScaTraField().NonlinearSolve();
  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: update                                                                  henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::UpdateFGIteration()
{
	// update the Fluid and the FGI vector at the end of the FGI loop
	return;
}

/*------------------------------------------------------------------------------------------------*
 * protected: update                                                                  henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::UpdateTimeStep()
{
  FluidField().Update();
  ScaTraField().Update();
  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: output                                                                  henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
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
