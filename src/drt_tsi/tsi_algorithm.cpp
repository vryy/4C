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
    // TSICoupling does not exist at the moment 10.12.09
 //   ThermoBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams(),"TSICoupling")
    ThermoBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams())
{
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
  //fsi-like part 10.12.09
  StructureField().ReadRestart(step);
//  double time = ThermoField().ReadRestart(step);
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
    // time loop
    while (NotFinished())
    {
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

// 10.12.09 does not work at the moment until there will be multiple dofsets in BACI
///*----------------------------------------------------------------------*
// | (protected)                                               dano 12/09 |
// *----------------------------------------------------------------------*/
//Teuchos::RCP<Epetra_Vector> TSI::Algorithm::StructToThermo(Teuchos::RCP<Epetra_Vector> iv) const
//{
//  return coupst_.MasterToSlave(iv);
//}
//
///*----------------------------------------------------------------------*
// | (protected)                                               dano 12/09 |
// *----------------------------------------------------------------------*/
//Teuchos::RCP<Epetra_Vector> TSI::Algorithm::ThermoToStruct(Teuchos::RCP<Epetra_Vector> iv) const
//{
//  return coupst_.SlaveToMaster(iv);
//}
//
///*----------------------------------------------------------------------*
// | (protected)                                               dano 12/09 |
// *----------------------------------------------------------------------*/
//Teuchos::RCP<Epetra_Vector> TSI::Algorithm::StructToThermo(
//  Teuchos::RCP<const Epetra_Vector> iv
//  ) const
//{
//  return coupst_.MasterToSlave(iv);
//}
//
///*----------------------------------------------------------------------*
// | (protected)                                               dano 12/09 |
// *----------------------------------------------------------------------*/
//Teuchos::RCP<Epetra_Vector> TSI::Algorithm::ThermoToStruct(
//  Teuchos::RCP<const Epetra_Vector> iv
//  ) const
//{
//  return coupst_.SlaveToMaster(iv);
//}

/*----------------------------------------------------------------------*/
#endif  // CCADISCRET
