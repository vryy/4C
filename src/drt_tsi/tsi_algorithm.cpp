/*----------------------------------------------------------------------*/
/*!
\file tsi_algorithm.cpp

\brief  Basis of all TSI algorithms that perform a coupling between the linear momentum equation
        and the heat conduction equation

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 12/09 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*
 | headers                                                   dano 12/09 |
 *----------------------------------------------------------------------*/
#include "tsi_algorithm.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.
/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 12/09 |
 *----------------------------------------------------------------------*/
TSI::Algorithm::Algorithm(
  Epetra_Comm& comm
  )
  : AlgorithmBase(comm,DRT::Problem::Instance()->TSIDynamicParams()),
    StructureBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams()),
    ThermoBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams()),
    tempincnp_(rcp(new Epetra_Vector(*(ThermoField().Tempnp()))))
{
  // build the proxy dofsets
  // build a proxy of the temperature discretization for the structure field
  Teuchos::RCP<DRT::DofSet> thermodofset =
    ThermoField().Discretization()->GetDofSetProxy();
  if (StructureField().Discretization()->AddDofSet(thermodofset)!=1)
    dserror("unexpected dof sets in structure field");
  // build a proxy of the structure discretization for the temperature field
  Teuchos::RCP<DRT::DofSet> structdofset =
    StructureField().Discretization()->GetDofSetProxy();
  if (ThermoField().Discretization()->AddDofSet(structdofset)!=1)
    dserror("unexpected dof sets in thermo field");

//  // now check if the two dofmaps are available and then bye bye
//  cout << "structure dofmap" << endl;
//  cout << *StructureField().DofRowMap(0) << endl;
//  cout << "thermo dofmap" << endl;
//  cout << *StructureField().DofRowMap(1) << endl;
//  cout << "thermo dofmap" << endl;
//  cout << *ThermoField().DofRowMap(0) << endl;
//  cout << "structure dofmap" << endl;
//  cout << *ThermoField().DofRowMap(1) << endl;
//  exit(0);
}

/*----------------------------------------------------------------------*
 | destructor (public)                                      dano 12/09 |
 *----------------------------------------------------------------------*/
TSI::Algorithm::~Algorithm()
{
}

/*----------------------------------------------------------------------*
 | (public)                                                  dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::ReadRestart(int step)
{
//  //fsi-like part 10.12.09
//  StructureField().ReadRestart(step);
//  double time = ThermoField().ReadRestart(step);
//  SetTimeStep(time,step);

  StructureField().ReadRestart(step);
  ThermoField().ReadRestart(step);
//  //20.01.2010 passt das so? entspricht das derselben Zeit??
//  double time = ThermoField().GetTime();
//  SetTimeStep(time,step);

//  // Georg´s version 10.12.09
//  ThermoField().ReadRestart(step);
//  StructureField().ReadRestart(step);
//  SetTimeStep(StructureField().Time(),step);
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TSI::Algorithm::TimeLoop()
{
  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn
    = DRT::Problem::Instance()->TSIDynamicParams();
  // Get the parameters for the ConvergenceCheck
  itmax_ = tsidyn.get<int>("ITEMAX"); // default: = 1
  ittol_ = tsidyn.get<double>("CONVTOL"); // default: =1e-6

  // time loop
  while (NotFinished())
  {
    // iterate between the two fields
    int  itnum = 0;  // itnum in Scatra
    bool stopnonliniter = false;

    // Start OUTER ITERATION
    while (stopnonliniter == false)
    {
      itnum ++;

      // store temperature from first solution for convergence check
      tempincnp_->Update(1.0,*ThermoField().Tempnp(),0.0);

      // prepare next time step
      PrepareTimeStep();

      // solve structure system
      DoStructureStep();

      // solve thermo system
      DoThermoStep();

      // update all single field solvers
      Update();

      // write output to screen and files
      Output();

      // check convergence of temperature field for "partitioned scheme" 10.02
      stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);
    } // end OUTER ITERATION

  } // time loop
}

/*----------------------------------------------------------------------*
 | (protected)                                               dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::PrepareTimeStep()
{
  IncrementTimeAndStep();

  PrintHeader();

  StructureField().PrepareTimeStep();
  ThermoField().PrepareTimeStep();
}

/*----------------------------------------------------------------------*
 | (protected)                                               dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::DoStructureStep()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"************************\n";
    cout<<"    STRUCTURE SOLVER    \n";
    cout<<"************************\n";
  }

  // solve structure system
  /// Do the nonlinear solve for the time step. All boundary conditions have
  /// been set.
  StructureField().Solve();
  return;
}


/*----------------------------------------------------------------------*
 | (protected)                                               dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::DoThermoStep()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"*********************\n";
    cout<<"    THERMO SOLVER    \n";
    cout<<"*********************\n";
  }

  /// solve temperature system
  /// Do the nonlinear solve for the time step. All boundary conditions have
  /// been set.
  ThermoField().Solve();
  return;
}


/*----------------------------------------------------------------------*
 | (protected)                                               dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::Update()
{
  StructureField().Update();
  ThermoField().Update();
  return;
}

/*----------------------------------------------------------------------*
 | (protected)                                               dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  StructureField().Output();
  ThermoField().Output();
}

/*----------------------------------------------------------------------*
 | convergence check for the temperature field               dano 02/10 |
 | originally by vg 01/09                                               |
 *----------------------------------------------------------------------*/
bool TSI::Algorithm::ConvergenceCheck(
  int itnum,
  const int itmax,
  const double ittol
  )
{
  // convergence check based on the temperature increment
  bool stopnonliniter = false;

  //    | temperature increment |_2
  //  -------------------------------- < Tolerance
  //     | temperature_n+1 |_2

  // Variables to save different L2 - Norms
  // define L2-norm of incremental temperature and temperature
  // here: only the temperature field is checked for convergence!!!
  double tempincnorm_L2(0.0); // pot
  double tempnorm_L2(0.0);

  // build the current temperature increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  tempincnp_->Update(1.0,*(ThermoField().Tempnp()),-1.0);

  // build the L2-norm of the temperature increment and the temperature
  tempincnp_->Norm2(&tempincnorm_L2);
  ThermoField().Tempnp()->Norm2(&tempnorm_L2);

  // care for the case that there is (almost) zero temperature
  // (usually not required for temperature)
  // if (tempnorm_L2 < 1e-6) tempnorm_L2 = 1.0;

  // Print the incremental based convergence check to the screen
  if (Comm().MyPID() == 0)
  {
    cout<<"\n";
    cout<<"****************************\n";
    cout<<"    OUTER ITERATION STEP    \n";
    cout<<"****************************\n";
    printf("+------------+-------------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|- temp-inc   -|\n");
    printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   |",
         itnum,itmax,ittol,tempincnorm_L2/tempnorm_L2);
    printf("\n");
    printf("+------------+-------------------+--------------+\n");
  }

  // Converged
  if (tempincnorm_L2/tempnorm_L2 <= ittol)
  {
    stopnonliniter = true;
    if (Comm().MyPID() == 0)
    {
      printf("\n");
      printf("|  >>>>>> converged after iteration %3d/%3d ! |\n", itnum,itmax);
      printf("+-----------------------------------------------+\n");
    }
  }

  // warn if itemax is reached without convergence, but proceed to next
  // timestep
  if ((itnum == itmax) and (tempincnorm_L2/tempnorm_L2 > ittol))
  {
    stopnonliniter = true;
    if ((Comm().MyPID() == 0))
    {
      printf("|     >>>>>> not converged in itemax steps!     |\n");
      printf("+-----------------------------------------------+\n");
      printf("\n");
      printf("\n");
    }
  }

  return stopnonliniter;
}  // TSI::Algorithm::ConvergenceCheck

/*----------------------------------------------------------------------*/
#endif  // CCADISCRET
