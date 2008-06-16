/*!----------------------------------------------------------------------*
\file combust_algorithm.H

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
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 * constructor                                              henke 06/08 * 
 *----------------------------------------------------------------------*/
COMBUST::Algorithm::Algorithm(Epetra_Comm& comm)
:  FluidBaseAlgorithm(DRT::Problem::Instance()->FluidDynamicParams(),false),
   ConDifBaseAlgorithm(DRT::Problem::Instance()->FluidDynamicParams()),
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
  const Teuchos::ParameterList& fluiddyn = DRT::Problem::Instance()->FluidDynamicParams();
  // maximum simulation time
  maxtime_=fluiddyn.get<double>("MAXTIME");
  // maximum number of timesteps
  numstep_ = fluiddyn.get<int>("NUMSTEP");
  // time step size
  dt_ = fluiddyn.get<double>("TIMESTEP");
  
  // counter FGI iterations
  fgiter = 0;
  maxfgiter = 1; // Dieser Parameter muss aus der Parameterliste ausgelesen werden
  // counter Fluid iterations
  fiter = 0;
  // counter G-function iterations
  giter = 0;

  // get RCP to actual velocity field (time n+1)
  velocitynp_=FluidField().ExtractVelocityPart(FluidField().Velnp());
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

    // reinitialize G-function 
    ReinitializeGfunc();
    
    // Fluid-G-function-Interaction loop
    while (NotConverged())
    {
    	// prepare Fluid-G-function iteration
    	PrepareFGIteration();
    	
    	// solve nonlinear Navier-Stokes system
    	DoFluidStep();

    	// solve G-function equation
    	DoGfuncStep();
    	
    	// update field vectors
    	UpdateFGI();
    	
    } // Fluid-G-function-Interaction loop
    
    // update all field solvers
    UpdateTime();

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
	/* Hier soll die G-Funktion reinitialisiert werden.
	   Dafür gibt es mehrere Möglichkeiten. Ein switch über alle Optionen 
	   ist hier nötig.
	   Zunächst wird einfach die signed-distance function gerufen.
	*/
	return;
}

/*----------------------------------------------------------------------*
 *  protected: FGI iteration converged?                     henke 06/08 *
 *----------------------------------------------------------------------*/
bool COMBUST::Algorithm::NotConverged()
{
	bool notconverged = false;
	if (fgiter <= maxfgiter and true) // (fgiter <= maxfgiter and ComputeGfuncNorm() < maxepsg and ComputeFluidNorm() < maxepsf)
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
  fgiter = 0;
  // fgnormgfunc = large value, determined in Input file
  // fgnormfluid = large value
  if (Comm().MyPID()==0)
  {
    cout<<"\n******************\n   FLUID SOLVER  \n******************\n";
  }

  FluidField().PrepareTimeStep();
  //ConDifField().PrepareTimeStep();
  return;
}

/*----------------------------------------------------------------------*
 * protected: prepare time step for combustion algorithm    henke 06/08 *
 *----------------------------------------------------------------------*/
void COMBUST::Algorithm::PrepareFGIteration()
{
	fgiter += 1;
	if (Comm().MyPID()==0)
	{
	  cout<<"\n******************\n   FGI   \n******************\n";
	}	
}

/*----------------------------------------------------------------------*
 * protected: perform a fluid time integration step         henke 06/08 *
 *----------------------------------------------------------------------*/
void COMBUST::Algorithm::DoFluidStep()
{
  // determine location of flame interface
  // FlameInterface(Gfunction);	
	
  // solve nonlinear Navier-Stokes system
  FluidField().NonlinearSolve();
  return;
}

/*----------------------------------------------------------------------*
 * protected: perform a G-function time integration step    henke 06/08 *
 *----------------------------------------------------------------------*/
void COMBUST::Algorithm::DoGfuncStep()
{
  
  if (Comm().MyPID()==0)
  {
    cout<<"\n****************\n G-FUNCTION SOLVER \n****************\n";
  }
  
  /* Hier muss eigentlich nur die Geschwindigkeit bei G=0 (Flammeninterface gelesen 
     werden). Tatsächlich kann vorerst aber das gesamte Fluid Feld als Input
     genommen werden. Die G-function sucht sich dann an jedem Knoten die Geschwindigkeit
     als Konvektions-Diffusions Geschwindigkeit.*/

  // get new velocity from Navier-Stokes solver
  velocitynp_=FluidField().ExtractVelocityPart(FluidField().Velnp());

  // prepare time step
  ConDifField().PrepareTimeStep();
  // transfer convective velocity
  ConDifField().SetVelocityField(2,velocitynp_);
  //ConDifField().SetVelocityField(1,1);

  // solve convection-diffusion equation
  ConDifField().Solve();
  return;
}

/*----------------------------------------------------------------------*
 * protected: update                                        henke 06/08 *
 *----------------------------------------------------------------------*/
void COMBUST::Algorithm::UpdateFGI()
{
	// update the Fluid and the FGI vector at the end of the FGI loop
	return;
}

/*----------------------------------------------------------------------*
 * protected: update                                        henke 06/08 *
 *----------------------------------------------------------------------*/
void COMBUST::Algorithm::UpdateTime()
{
  FluidField().Update();
  ConDifField().Update();
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
  ConDifField().Output();

  // debug IO
#if 0
  // print out convective velocity vector
  cout<<*velocitynp_<<endl;
#endif

  return;
}

#endif // #ifdef CCADISCRET
