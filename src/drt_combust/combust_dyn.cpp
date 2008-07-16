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

// übernommen von Axels xfluid_dyn_nln_drt.cpp
// die grünen Ergänzungen sind bei Georg zusätzlich eingebunden
// Dinge die Georg hat, aber Axel nicht, sind mit "? nötig" markiert

#include <ctime>   // ? nötig
#include <cstdlib> // ? nötig
#include <iostream>
#include <string>  // ? nötig

#include <Teuchos_Time.hpp>  // ? nötig
#include <Teuchos_TimeMonitor.hpp>
//#include <Teuchos_StandardParameterEntryValidators.hpp> // von Axel!

#include <Epetra_Time.h>  // ? nötig

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
#include "combust_create_gfunction.H"
//#include "xfluidimplicitintegration.H"
//#include "../drt_lib/drt_resulttest.H"
//#include "xfluidresulttest.H"
//
//#include "../drt_lib/drt_validparameters.H" // ? nötig
//#include "../drt_lib/drt_condition_utils.H" // ? nötig


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

  // print COMBUST-Logo on screen
  if (comm.MyPID()==0) COMBUST::printlogo();

  // -------------------------------------------------------------------
  // prepare the discretizations
  // -------------------------------------------------------------------

  // get discretization ids
  int disnumff = genprob.numff; // discretization number fluid; typically 0
  int disnumgff = genprob.numscatra; // discretization number G-function; typically 1
  // remark: precisely here the convection-diffusion discretization 
  // (genprob.numscatra) is called and named according to its meaning 
  // in the combustion context, namely the G-function (disnumgff)
  
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
    COMBUST::CreateGfuncDiscretization(disnumff,disnumgff);
    if (comm.MyPID()==0)
    cout<<"Created G-function discretization from fluid discretization in...."
    <<time.ElapsedTime() << " secs\n\n";
  }
  else
    dserror("G-function discretization is not empty. Fluid and G-function already present. This is not supported.");

  
/*  const int fmyrank = fluiddis->Comm().MyPID();
 *  std::cout << "FluidProc: " << fmyrank << endl;
 * flush(cout);
 * const int smyrank = soliddis->Comm().MyPID();
 * std::cout << "SolidProc: " << smyrank << endl;
 * flush(cout);
 * //cout << *soliddis;
 */
  
 // -------------------------------------------------------------------
 // set degrees of freedom in the discretization
 // -------------------------------------------------------------------
//  if (!fluiddis->Filled()) fluiddis->FillComplete();
//  if (!gfuncdis->Filled()) gfuncdis->FillComplete();
  
  // create an COMBUST::Algorithm instance
  /* Muss das nicht eigentlich combust_ heissen? Wegen Konvention   henke 07/08*/
    Teuchos::RCP<COMBUST::Algorithm> combust = Teuchos::rcp(new COMBUST::Algorithm(comm));

    if (genprob.restart)
    {
      // read the restart information, set vectors and variables
      //combust->ReadRestart(genprob.restart);
      dserror("restart not yet available");
      exit(1);
    }
    
    // solve the whole combustion problem
    combust->TimeLoop();

    // summarize the performance measurements
    Teuchos::TimeMonitor::summarize();

    // perform the result test
    DRT::ResultTestManager testmanager(comm);
    testmanager.AddFieldTest(combust->FluidField().CreateFieldTest());
    testmanager.AddFieldTest(combust->CreateScaTraFieldTest());
    testmanager.TestAll();

    return;

} // combust_dyn()  

#endif  // #ifdef CCADISCRET
