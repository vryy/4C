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
#include "combust_flamefront.H"
#include "combust_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "combust_fluidimplicitintegration.H"

/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
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
   * 
   * scalar transport velocity field has been initialized in ScaTraFluidCouplingAlgorithm()
  */
COMBUST::Algorithm::Algorithm(Epetra_Comm& comm, const Teuchos::ParameterList& combustdyn)
: ScaTraFluidCouplingAlgorithm(comm, combustdyn),
// initialize member variables
  fgiter_(0),
  fgitermax_(combustdyn.get<int>("ITEMAX")),
//  fgvelnormL2_(?),
//  fggfuncnormL2_(?),
/* hier müssen noch mehr Parameter eingelesen werden, nämlich genau die, die für die FGI-Schleife
   nötig sind. Die für die Schleifen über die Einzelfelder existieren in den Zeitintegrationsschemata. */
//  reinitializationaction_(combustdyn.sublist("COMBUSTION GFUNCTION").get<INPAR::COMBUST::ReInitialActionGfunc>("REINITIALIZATION")),
  reinitializationaction_(Teuchos::getIntegralValue<INPAR::COMBUST::ReInitialActionGfunc>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REINITIALIZATION")),
  combustdyn_(combustdyn),
  interfacehandle_(Teuchos::null),
  flamefront_(Teuchos::null)
  {

  // get pointers to the discretizations from the time integration scheme of each field
  /* remark: fluiddis cannot be of type "const Teuchos::RCP<const DRT::Dis...>", because parent
   * class. InterfaceHandle only accepts "const Teuchos::RCP<DRT::Dis...>"            henke 01/09 */
  const Teuchos::RCP<DRT::Discretization> fluiddis = FluidField().Discretization();
  const Teuchos::RCP<DRT::Discretization> gfuncdis = ScaTraField().Discretization();

  /*----------------------------------------------------------------------------------------------*
   * initialize all data structures needed for the combustion algorithm
   *
   * - capture the flame front and create interface geometry (triangulation)
   * - determine initial enrichment (DofManager wird bereits mit dem Element d.h. Diskretisierung angelegt)
   * - ...
   *----------------------------------------------------------------------------------------------*/
  // construct initial flame front
  flamefront_ = rcp(new COMBUST::FlameFront(fluiddis,gfuncdis));
  flamefront_->ProcessFlameFront(combustdyn_,ScaTraField().Phinp());

  // construct interfacehandle using initial flame front
  interfacehandle_ = rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis,flamefront_));

  std::cout << "Combustion Algorithm constructor done \n" << endl;
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

  if (Comm().MyPID()==0)
  {
    std::cout << "Combustion Algorithm Timeloop starting \n" << endl;
  }

  // time loop
  while (NotFinished())
  {
    // prepare next time step
    PrepareTimeStep();

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
  dserror("Der Algorithmus für statische Verbrennungprobleme kann noch nix!");
  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: reinitialize G-function                                                 henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::ReinitializeGfunc(INPAR::COMBUST::ReInitialActionGfunc action)
{
  /* Here, the G-function is reinitialized, because we suspect that the signed distance property has
   * been lost due to numerical inaccuracies. There are various options to reinitialize the
   * G-function. For the time being we use an exact, but expensive, procedure to ensure this
   * property which computes the distance for every point. (SignedDistFunc)           henke 06/08 */
  std::cout << "This is the ReinitializeGfunc() function \n" << endl;
  switch (action)
  {
    case INPAR::COMBUST::reinitialize_by_function:
      // read a FUNCTION from the input file and reinitialize the G-function with its values
      ScaTraField().SetInitialField(1,combustdyn_.sublist("COMBUSTION GFUNCTION").get<int>("REINITFUNCNO"));
      break;
    case INPAR::COMBUST::compute_signeddistancefunction:
      SignedDistFunc();
      break;
    default:
      dserror ("Unknown option to reinitialize the G-function");
  }
  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: build signed distance function                                          henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::SignedDistFunc()
{
  /* This member function constructs a G-function field that meets the signed distance property. The
   * algorithm assigns the value of the distance between each node and the surface defined by G=0 as
   * a scalar value to every node in the G-function discretization.*/

  // get a pointer to the G-function discretization
  Teuchos::RCP<DRT::Discretization> gfuncdis = ScaTraField().Discretization();

  const Epetra_Map* dofrowmap = gfuncdis->DofRowMap();

  // loop all nodes on the processor
  for(int lnodeid=0; lnodeid < gfuncdis->NumMyRowNodes(); lnodeid++)
  {
    // get the processor local node
    DRT::Node*  lnode      = gfuncdis->lRowNode(lnodeid);
    // the set of degrees of freedom associated with the node
    vector<int> nodedofset = gfuncdis->Dof(lnode); // this should not be a vector, but a scalar

    int numdofs = nodedofset.size(); // this should be 1 (scalar field!)
    if (numdofs != 1)
      dserror("There are more than 1 dof for a node in a scalar field!");

    for (int k=0;k< numdofs;++k) // Ich lasse die Schleife erstmal drin -> mehr als eine level set Function möglich
    {
      const int dofgid = nodedofset[k];
      int doflid = dofrowmap->LID(dofgid);
      // compute value of signed distance function for this doflid
      // evaluate component k of spatial function
      // remark: if there is only 1 level set function, there should be only 1 component (scalar field!)
//      double distance = GEO::computeDistance(lnode->X(),flamefront_,k);
      double distance = 2.3; // hack to compile the code: Complete nonsense!
/*
 *  klären ob es eine solche Funktion schon gibt und ob sie im GEO namespace sein soll
 *
 *  Der Input für "computeDistance" sind:
 *  - Koordinaten für einen Knoten
 *  - Beschreibung aller Flächen, aus denen sich das Interface zusammensetzt
 *  - "k" wird wahrscheinlich sehr lange nicht gebraucht werden.
 */
        // assign new values of signed distance function to the G-function field
        ScaTraField().Phin()->ReplaceMyValues(1,&distance,&doflid);
//        ScaTraField().Phin()->ReplaceMyValues(1,&distance,&doflid);
    }
  }
  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: FGI iteration converged?                                                henke 06/08 |
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

  // prepare time step
  /* remark: initial velocity field has been transfered to scalar transport field in constructor of 
   * ScaTraFluidCouplingAlgorithm (initialvelset_ == true). Time integration schemes, such as 
   * the one-step-theta scheme, are thus initialized correctly.
   */
  ScaTraField().PrepareTimeStep();

  // synchronicity check between combust algorithm and base algorithms
  if (FluidField().Time() != Time())
    dserror("Time in Fluid time integration differs from time in combustion algorithm");
  if (ScaTraField().Time() != Time())
    dserror("Time in ScaTra time integration  differs from time in combustion algorithm");

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
 | protected: perform a fluid time integration step                                   henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::DoFluidField()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n---------------------------------------  FLUID SOLVER  ---------------------------------------\n";
  }

  // export interface information to the fluid time integration
  FluidField().ImportInterface(interfacehandle_);

  // solve nonlinear Navier-Stokes equations
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

  /* Hier muss eigentlich nur die Geschwindigkeit bei G=0 (Flammeninterface gelesen werden).
   * Tatsächlich kann vorerst aber das gesamte Navier-Stokes Fluid Feld als Input genommen werden.
   * Die G-function sucht sich dann an jedem Knoten die Geschwindigkeit als Konvektions-Diffusions
   * Geschwindigkeit.
   */

  // assign the fluid velocity to the G-function field as convective velocity
  ScaTraField().SetVelocityField(
      FluidField().Velnp(),
      FluidField().SgVelVisc(),
      FluidField().Discretization()
  );

  // solve nonlinear convection-diffusion equation
  ScaTraField().NonlinearSolve();
  // solve linear convection-diffusion equation
//  ScaTraField().Solve();

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: update                                                                  henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::UpdateFGIteration()
{
  // update flame front according to evolved G-function field
  flamefront_->ProcessFlameFront(combustdyn_,ScaTraField().Phinp());

  // update interfacehandle according to updated flame front
  interfacehandle_->UpdateInterfaceHandle();

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
