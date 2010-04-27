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
  Teuchos::RCP<DRT::DofSet> thermodofset
    = ThermoField().Discretization()->GetDofSetProxy();
  if (StructureField().Discretization()->AddDofSet(thermodofset)!=1)
    dserror("unexpected dof sets in structure field");
  // build a proxy of the structure discretization for the temperature field
  Teuchos::RCP<DRT::DofSet> structdofset
    = StructureField().Discretization()->GetDofSetProxy();
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

  // now build the matrices again and consider the temperature dependence
  StructureField().TSIMatrix();

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
  // 14.04.10
  ThermoField().ReadRestart(step);
  StructureField().ReadRestart(step);
  SetTimeStep(ThermoField().GetTime(),step);

  //  //fsi-like part 10.12.09
  //  StructureField().ReadRestart(step);
  //  double time = ThermoField().ReadRestart(step);
  //  SetTimeStep(time,step);

  //  StructureField().ReadRestart(step);
  //  ThermoField().ReadRestart(step);
  //  //20.01.2010 passt das so? entspricht das derselben Zeit??
  //  double time = ThermoField().GetTime();
  //  SetTimeStep(time,step);

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
  itmax_ = tsidyn.get<int>("ITEMAX"); // default: =1
  ittol_ = tsidyn.get<double>("CONVTOL"); // default: =1e-6

  // get an idea of the temperatures (like in partitioned FSI)
  itempn_ = ThermoField().ExtractTemperatures();

  // ==================================================================

  // time loop
  while (NotFinished())
  {
    // prepare next time step ==> muss ausserhalb der Newton-Schleife sein!! 21.04.10 ULI!!!
    PrepareTimeStep();

    // get active nodes from structural contact simulation
    RCP<MORTAR::ManagerBase> cmtman = StructureField().ContactManager();

    // tsi with or without contact
    if (cmtman == Teuchos::null)
    {

      // Begin Nonlinear Solver / Outer Iteration *******************************

      // iterate between the two fields
      int  itnum = 0;  // itnum in Scatra
      bool stopnonliniter = false;

      // Outer Iteration loop starts
      if (Comm().MyPID()==0)
      {
        cout<<"\n";
        cout<<"**************************************************************\n";
        cout<<"      OUTER ITERATION LOOP \n";
        printf("      Time Step %3d/%3d \n",ThermoField().GetTimeStep(), ThermoField().GetTimeNumStep());
        cout<<"**************************************************************\n";
      }

      // Start OUTER ITERATION
      while (stopnonliniter == false)
      {
        itnum ++;

        // store temperature from first solution for convergence check
        tempincnp_->Update(1.0,*ThermoField().Tempnp(),0.0);

        const Teuchos::RCP<Epetra_Vector> itemp = DoThermoStep();

        // solve structure system
        DoStructureStep(itemp);

        // check convergence of temperature field for "partitioned scheme"
        stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);

      } // end OUTER ITERATION

      // End Nonlinear Solver **************************************

      // ==================================================================

      // update all single field solvers
      Update();

      // extract final temperatures,
      // since we did update, this is very easy to extract
      itempn_ = ThermoField().ExtractTemperatures();

      // write output to screen and files
      Output();
    }

    else
    {
      // not yet implemented
      dserror("TSI with CONTACT not yet implemented");

      // prepare next time step
      PrepareTimeStep();

      // solve thermo system
      const Teuchos::RCP<Epetra_Vector> itemp = DoThermoStep();

      // solve structure system
      DoStructureStep(itemp);

      // update all single field solvers
      Update();

      // extract final temperatures,
      // since we did update, this is very easy to extract
      itempn_ = ThermoField().ExtractTemperatures();

      // write output to screen and files
      Output();
    }
  } // time loop
}


/*----------------------------------------------------------------------*
 | prepare time step (protected)                             dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::PrepareTimeStep()
{
  IncrementTimeAndStep();

  PrintHeader();

  // predict
  // 14.04.10 check again if PrepareTimeStep() is needed here!
  StructureField().PrepareTimeStep();
  ThermoField().PrepareTimeStep();
}


/*----------------------------------------------------------------------*
 | Solve the structure system (protected)                    dano 03/10 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::DoStructureStep(Teuchos::RCP<Epetra_Vector> itemp)
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"************************\n";
    cout<<"    STRUCTURE SOLVER    \n";
    cout<<"************************\n";
  }
  // call the current temperatures
  StructureField().ApplyTemperatures(itemp);

  // call the predictor here, because the temperature is considered
  StructureField().PrepareTimeStep();

  // solve structure system
  /// Do the nonlinear solve for the time step. All boundary conditions have
  /// been set.
  StructureField().Solve();
  return;
}


/*----------------------------------------------------------------------*
 | Solve the thermo system (protected)                       dano 03/10 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> TSI::Algorithm::DoThermoStep(
  //Teuchos::RCP<Epetra_Vector> idisp // 25.3.10 necessary if coupling displacemnents to temperature
  )
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
  // now extract the current temperatures and pass it to the structure
  return ThermoField().ExtractTemperatures();
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

//  // 14.04.10 NAN
//  cout << "\ntempincnorm_L2\n" << tempincnorm_L2 << endl;
//  cout << "\ntempnorm_L2\n" << tempnorm_L2 << endl;

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
      printf("|  Outer Iteration loop converged after iteration %3d/%3d ! |\n", itnum,itmax);
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
