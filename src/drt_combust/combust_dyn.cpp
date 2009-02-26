/*!---------------------------------------------------------------------*
\file combust_dyn.cpp

\brief control routine for dynamic combustion analysis

	The approach to simulation of premixed combustion followed by 
	the Emmy-Noether-Group at LNM consists of solving a coupled 
	two-field problem. One field is the fluid and the second field 
	is the level-set function indicating the location of the flame 
	interface. The level-set function is a scalar transport equation 
	of convection (and diffusion)-type. The scalar is called "G" and
	therefore the level-set function is often called the "G-function".
	
	Implementing the interaction between both fields with the FEM,
	results in solving both (non-linear) systems of equations separately 
	and coupling them via exchanging information at the flame interface. 
	Two discretizations are created in BACI, one of type "fluid" and 
	another one of type "condif".
	
	"combust_dyn" is the upper most control routine for the dynamic 
	combustion analysis in the DRT::COMBUST module.
	It controls the following tasks:
	
	o create the G-function discretization from the fluid discretization
	o call a combustion algorithm to solve the multi-field problem
	o check the results

	remark: in this function, the notation is changed from 	"condif"-
	style to "gfunc"-style . That is, even though we are still dealing 
	with a convection-diffusion field using condif elements and so on,
	for convenience the "condif" discretization is from now on named 
	after the G-funciton: 
	
	disnumcdf -> disnumgff
	condifdis -> gfuncdis
	...
	

\author henke
\date 06/08

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifdef PARALLEL
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "../drt_lib/drt_globalproblem.H"
#include "combust_dyn.H"
#include "combust_utils.H"
#include "combust_algorithm.H"
#include "../drt_scatra/scatra_utils.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 | main  control routine for dynamic combustion analysis    henke 06/08 |
 *----------------------------------------------------------------------*/
void combust_dyn()
{
  // create a communicator
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  /*----------------------------------------------------------------------------------------------*
   * print COMBUST-Logo on screen
   *----------------------------------------------------------------------------------------------*/
  if (comm.MyPID()==0) COMBUST::printlogo();

  /*----------------------------------------------------------------------------------------------*
   * create G-function discretization and set degrees of freedom in fluid and gfunc discretization
   *----------------------------------------------------------------------------------------------*/
  // get discretization ids
  int disnumff = genprob.numff; // discretization number fluid; typically 0
  int disnumgff = genprob.numscatra; // discretization number G-function; typically 1
  /* remark: precisely here the ScaTra discretization (genprob.numscatra) is named according to its 
   * physical meaning in the combustion context, (namely the G-function (disnumgff)).
   * -> Here ScaTra becomes Combustion! */
  
  // access fluid discretization
  RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(disnumff,0);
  if (!fluiddis->Filled()) fluiddis->FillComplete();
  if (fluiddis->NumGlobalNodes()==0)
    dserror("No fluid discretization found!");

  // access G-function discretization (it should be empty)
  RCP<DRT::Discretization> gfuncdis = DRT::Problem::Instance()->Dis(disnumgff,0);
  if (!gfuncdis->Filled()) gfuncdis->FillComplete();

  // create G-function discretization (fill with condif elements)
  if (gfuncdis->NumGlobalNodes()==0)
  {
    Epetra_Time time(comm);
    std::map<string,string> conditions_to_copy;
    conditions_to_copy.insert(pair<string,string>("TransportDirichlet","Dirichlet"));
    //conditions_to_copy.insert("FluidStressCalc","FluxCalculation"); // a hack

    // access the scalar transport parameter list
    const Teuchos::ParameterList& scatracontrol = DRT::Problem::Instance()->ScalarTransportDynamicParams();
    const int matid = scatracontrol.get<int>("MATID");

    SCATRA::CreateScaTraDiscretization(fluiddis,gfuncdis,conditions_to_copy,matid,false);

    if (comm.MyPID()==0)
      cout<<"Created G-function discretization from fluid discretization in...."
          <<time.ElapsedTime() << " secs\n\n";
  }
  else
    dserror("G-function discretization is not empty. Fluid and G-function already present!");

  /*----------------------------------------------------------------------------------------------*
   * create a combustion algorithm
   *----------------------------------------------------------------------------------------------*/
  // get the combustion parameter list
  Teuchos::ParameterList combustdyn = DRT::Problem::Instance()->CombustionDynamicParams();
  // create a COMBUST::Algorithm instance
  Teuchos::RCP<COMBUST::Algorithm> combust_ = Teuchos::rcp(new COMBUST::Algorithm(comm, combustdyn));

  /*----------------------------------------------------------------------------------------------* 
   * restart stuff; what exactly does this do?; comment has to be completed
   *----------------------------------------------------------------------------------------------*/ 
  if (genprob.restart)
  {
    // read the restart information, set vectors and variables
    combust_->ReadRestart(genprob.restart);
  }
  /*----------------------------------------------------------------------------------------------*
   * call one of the available time integration schemes
   *----------------------------------------------------------------------------------------------*/
  FLUID_TIMEINTTYPE timeintscheme = Teuchos::getIntegralValue<FLUID_TIMEINTTYPE>(combustdyn,"TIMEINTEGR");

  if (timeintscheme == timeint_one_step_theta)
  {
    // solve a dynamic combustion problem
    combust_->TimeLoop();
    /* Hier kann auch z.B. Genalpha mit combust->TimeLoop() gerufen werden, weil der combustion 
     * Algorithmus ja schon weiss welche Zeitintegration er hat. Es muss dann eine Klasse
     * "GenalphaTimeInt" existieren, die eine Funktion TimeLoop() hat. Dann muss allerdings auch 
     * das ADAPTER::FluidCombust ein entsprechendes Object GenalphaTimeInt haben. Momentan hat ein
     * FluidCombust immer ein CombustFluidImplicitTimeInt! */
  }
  else if (timeintscheme == timeint_stationary)
  {
    // solve a static combustion problem
    combust_->SolveStationaryProblem();
  }
  else
  {
    // error: impossible, every time integration scheme is static or dynamic
    dserror("the combustion module can not handle this type of time integration scheme");
  }

  /*----------------------------------------------------------------------------------------------*
   * validate the results
   *----------------------------------------------------------------------------------------------*/
  // summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform the result test
  DRT::ResultTestManager testmanager(comm);
  testmanager.AddFieldTest(combust_->FluidField().CreateFieldTest());
  testmanager.AddFieldTest(combust_->CreateScaTraFieldTest());
  testmanager.TestAll();

  return;

} // combust_dyn()

#endif  // #ifdef CCADISCRET
