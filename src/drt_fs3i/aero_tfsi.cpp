/*----------------------------------------------------------------------*/
/*!
\file aero_tfsi.cpp

\brief The aero code INCA is coupled with the monolithic TSI system
       based on communication with MPI

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-152537
</pre>
*/


/*----------------------------------------------------------------------*
 | headers                                                  ghamm 12/11 |
 *----------------------------------------------------------------------*/
#include "aero_tfsi.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_adapter/adapter_thermo.H"
#include "../drt_tsi/tsi_monolithic.H"
#include "../drt_tsi/tsi_partitioned.H"
#include "../drt_tsi/tsi_utils.H"
#include "aero_tfsi_serv.H"



/*----------------------------------------------------------------------*
 | AeroTFSI                                                 ghamm 12/11 |
 *----------------------------------------------------------------------*/
FS3I::AeroTFSI::AeroTFSI(
  const Epetra_Comm& lcomm
  ) :
  FS3I_Base(),
  lcomm_(lcomm),
  tsi_(Teuchos::null)
{
  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn = DRT::Problem::Instance()->TSIDynamicParams();
  // decide if monolithic or partitioned coupling
  coupling_ = DRT::INPUT::IntegralValue<INPAR::TSI::SolutionSchemeOverFields>(tsidyn,"COUPALGO");

  // check if INCA is called first when starting the coupled simulation
  int worldrank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldrank);

  if(worldrank == lcomm.MyPID())
    dserror("ERROR: INCA must!! be always started first in a coupled simulation");

  // first: typecast the Epetra_Comm to Epetra_MpiComm
  const Epetra_MpiComm& epetrampicomm = dynamic_cast<const Epetra_MpiComm&>(lcomm);
  if (!(&epetrampicomm))
    dserror("ERROR: casting Epetra_Comm -> Epetra_MpiComm failed");
  const MPI_Comm localcomm = epetrampicomm.GetMpiComm();

  // leaders for each code are hard coded as follows
  localBACIleader_ = 0;
  INCAleader_ = 0;

  // intercommunicator is created; tag is important here because the same tag is used in INCA
  int tag = 0;
  int remoteleader = INCAleader_;
  MPI_Intercomm_create(localcomm, 0, MPI_COMM_WORLD, remoteleader, tag, &intercomm_);

  MPI_Barrier(intercomm_);

  // setup of the discretizations, including clone strategy
  TSI::UTILS::SetupTSI(lcomm);

  switch(coupling_)
  {
  case INPAR::TSI::Monolithic:
  {
    // create an TSI::Monolithic instance
    const Teuchos::ParameterList& sdynparams = DRT::Problem::Instance()->StructuralDynamicParams();
    tsi_ = Teuchos::rcp(new TSI::Monolithic(lcomm,sdynparams));
    break;
  }
  case INPAR::TSI::IterStagg:
  {
    // create an TSI::Partitioned instance
    tsi_ = Teuchos::rcp(new TSI::Partitioned(lcomm));
    break;
  }
  default:
  {
    dserror("For TFSI coupling only Monolithic or IterStagg is available.");
    break;
  }
  }

  // setup of the helper class
  aerocoupling_ = Teuchos::rcp(new FS3I::UTILS::AeroCouplingUtils(tsi_->StructureField()->Discretization(),
                            tsi_->ThermoField()->Discretization()));

}


/*----------------------------------------------------------------------*
 | time loop of the coupled tfsi system                     ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::Timeloop()
{
  // in case of restart, initial interface data must be send at this position
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    tsi_->ReadRestart(restart);

    // interface data has to be send to INCA in case of restart at the very beginning
    // of the time loop
    Teuchos::RCP<Epetra_Vector> idispnRestart = aerocoupling_->StrExtractInterfaceVal(tsi_->StructureField()->ExtractDispn());
    Teuchos::RCP<Epetra_Vector> ivelnRestart = aerocoupling_->StrExtractInterfaceVal(tsi_->StructureField()->ExtractVeln());
    Teuchos::RCP<Epetra_Vector> ithermoloadRestart = aerocoupling_->ThrExtractInterfaceVal(tsi_->ThermoField()->ExtractTempn());

    vector<double> aerosenddataRestart;
    aerocoupling_->PackData(idispnRestart, ithermoloadRestart, aerosenddataRestart);
    SendAeroData(aerosenddataRestart);
  }

  vector<double> timestep(1);
  // get the first time step from INCA; timen_ and dt has to be set correctly
  GetTimeStep(timestep);
  SetInitialTimeStepAndTime(timestep);

  // time loop
  while (tsi_->NotFinished())
  {
    // receive data from INCA and make it available on all BACI procs
    vector<double> aerodata;
    ReceiveAeroData(aerodata);

    //fill data from INCA into suitable variables
    map<int, map<int, LINALG::Matrix<3,1> > > aerocoords;
    map<int, map<int, LINALG::Matrix<4,1> > > aeroforces;
    SplitData(aerodata, aerocoords, aeroforces);

    //get new vectors for mapping the fluid interface data
    Teuchos::RCP<Epetra_Vector> iforce = LINALG::CreateVector(*(aerocoupling_->GetInterfaceStructDis()->DofRowMap()), true);
    Teuchos::RCP<Epetra_Vector> ithermoload = LINALG::CreateVector(*(aerocoupling_->GetInterfaceThermoDis()->DofRowMap()), true);

    // current displacement of the interface
    Teuchos::RCP<Epetra_Vector> idispn = aerocoupling_->StrExtractInterfaceVal(tsi_->StructureField()->ExtractDispn());

    aerocoupling_->ProjectForceOnStruct(idispn, aerocoords, aeroforces, iforce, ithermoload);

    ApplyInterfaceData(iforce, ithermoload);

    // counter and print header; predict solution of both fields
    tsi_->PrepareTimeStep();

    // TSI time step is solved
    SolveTSIstep();

    // calculate stresses, strains, energies
    tsi_->PrepareOutput();

    // communication with AERO-code to send interface position and temperature
    // skip last sending procedure in the simulation because INCA has already finished
    if(tsi_->NotFinished())
    {
      // gather data and send it to INCA
      idispn->PutScalar(0.0);
      ithermoload->PutScalar(0.0);
      // reuse of idispn; data is at n+1
      idispn = aerocoupling_->StrExtractInterfaceVal(tsi_->StructureField()->ExtractDispnp());
      // calculate velocity of the structure currently not needed
      // extract interface temperatures
      ithermoload = aerocoupling_->ThrExtractInterfaceVal(tsi_->ThermoField()->ExtractTempnp());

      vector<double> aerosenddata;
      aerocoupling_->PackData(idispn, ithermoload, aerosenddata);

      SendAeroData(aerosenddata);

      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      // INCA is working at the moment
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      // get the time step from INCA for the next time step in BACI
      GetTimeStep(timestep);

      // set the time step received from INCA
      SetTimeStep(timestep);
    }

    // update all single field solvers
    tsi_->Update();

    // write output to screen and files
    tsi_->Output();

  }  // NotFinished

}  // TimeLoop


/*----------------------------------------------------------------------*
 | setup of the TSI system                                  ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::SetupSystem()
{
  tsi_->SetupSystem();

  return;
}


/*----------------------------------------------------------------------*
 | Apply fluid interface loads to the TSI system           ghamm 12/11  |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::ApplyInterfaceData(
  Teuchos::RCP<Epetra_Vector> iforce,
  Teuchos::RCP<Epetra_Vector> ithermoload
  )
{
  // apply structural interface tractions to the structural field
  Teuchos::RCP<Epetra_Vector> strfifc = LINALG::CreateVector(*tsi_->StructureField()->Discretization()->DofRowMap(), true);
  aerocoupling_->StrApplyInterfaceVal(iforce, strfifc);
  tsi_->StructureField()->SetForceInterface(strfifc);
  tsi_->StructureField()->PreparePartitionStep();

  // apply thermal interface heat flux to the thermal field
  Teuchos::RCP<Epetra_Vector> thrfifc = LINALG::CreateVector(*tsi_->ThermoField()->Discretization()->DofRowMap(), true);
  aerocoupling_->ThrApplyInterfaceVal(ithermoload, thrfifc);
  tsi_->ThermoField()->SetForceInterface(thrfifc);

  return;
}


/*----------------------------------------------------------------------*
 | Solve the current TSI time step                         ghamm 02/12  |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::SolveTSIstep()
{
  switch(coupling_)
  {
  case INPAR::TSI::Monolithic:
  {
    Teuchos::rcp_dynamic_cast<TSI::Monolithic>(tsi_,true)->NewtonFull();
    break;
  }
  case INPAR::TSI::IterStagg:
  {
    Teuchos::rcp_dynamic_cast<TSI::Partitioned>(tsi_,true)->TimeLoopFull();
    break;
  }
  default:
  {
    dserror("For TFSI coupling only Monolithic or IterStagg is available.");
    break;
  }
  }

  return;
}


/*----------------------------------------------------------------------*
 | Receive time step from INCA via MPI                      ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::GetTimeStep(
  vector<double>& timestep
  )
{
  MPI_Barrier(intercomm_);

  // get time step from INCA
  if(lcomm_.MyPID() == 0)
  {
    MPI_Status status;
    int tag_timestep = 3000;
    MPI_Recv(&timestep[0], 1, MPI_DOUBLE, INCAleader_, tag_timestep, intercomm_, &status);
  }
  // broadcast time step to all processors in BACI
  lcomm_.Broadcast(&timestep[0], 1 , localBACIleader_);

  return;
}


/*----------------------------------------------------------------------*
 | apply time step from INCA to BACI                        ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::SetTimeStep(
  vector<double>& timestepsize
  )
{
  tsi_->StructureField()->SetTimeStepSize(timestepsize[0]);

  tsi_->ThermoField()->SetTimeStepSize(timestepsize[0]);

  tsi_->SetDt(timestepsize[0]);

  return;
}


/*----------------------------------------------------------------------*
 | apply first time step and adapt time                     ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::SetInitialTimeStepAndTime(
  vector<double>& timestepsize
  )
{
  tsi_->StructureField()->SetInitialTimeStepAndTime(timestepsize[0]);

  tsi_->ThermoField()->SetInitialTimeStepAndTime(timestepsize[0]);

  tsi_->SetDt(timestepsize[0]);

  return;
}


/*----------------------------------------------------------------------*
 | communication with AERO-code to receive the interface forces         |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::ReceiveAeroData(
  vector<double>& aerodata
  )
{
  // ==================================================================
  // intercommunicator is used to receive data on proc 0
  int lengthRecv = 0;

  if(lcomm_.MyPID() == 0)
  {
    MPI_Status status;
    //first: receive length of data from aero code
    int tag = 4000;
    MPI_Recv(&lengthRecv, 1, MPI_INT, INCAleader_, tag, intercomm_, &status);
    if(lengthRecv == 0)
      dserror("Length of data from INCA is not received correctly!!");

    //second: receive data
    tag = 5000;
    aerodata.resize(lengthRecv);
    MPI_Recv(&aerodata[0], lengthRecv, MPI_DOUBLE, INCAleader_, tag, intercomm_, &status);
  }

  MPI_Barrier(intercomm_);

  // ==================================================================
  // broadcast received data in BACI to all processors
  // again: length first, then data
  lcomm_.Broadcast(&lengthRecv, 1 , localBACIleader_);

  aerodata.resize(lengthRecv);

  //aero forces have to be redundant on all processors
  lcomm_.Broadcast(&aerodata[0], lengthRecv , localBACIleader_);

  lcomm_.Barrier();

  return;
}


/*----------------------------------------------------------------------*
 | cast data from INCA in a different format                ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::SplitData(
  vector<double> aerodata,
  map<int, map<int, LINALG::Matrix<3,1> > >& aerocoords,
  map<int, map<int, LINALG::Matrix<4,1> > >& aeroforces
  )
{
  int length = aerodata.size();
  for(int out=0; out<(length/4); out++)
  {
    LINALG::Matrix<3,1> tmp1;
    LINALG::Matrix<4,1> tmp2;
    for(int in=0; in<3; in++)
    {
      tmp1(in)=aerodata[0 + out*4 + in];
    }
    // note: currently heat fluxes are transferred but forces not yet --> zeros
    for(int in=0; in<3; in++)
    {
      tmp2(in)=0.0;
    }
    tmp2(3)=aerodata[3 + out*4];

    aerocoords[0][out] = tmp1;
    aeroforces[0][out] = tmp2;
  }

  return;
}


/*----------------------------------------------------------------------*
 | communication with AERO-code to send interface data      ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::SendAeroData(
  vector<double>& aerosenddata
  )
{
  // ==================================================================
  // intercommunicator is used to send data to INCA proc 0
  MPI_Barrier(intercomm_);

  if(lcomm_.MyPID() == 0)
  {
    int lengthSend = aerosenddata.size();
    //first: send length of data from aero code
    int tag = 6000;
    MPI_Send(&lengthSend, 1, MPI_INT, INCAleader_, tag, intercomm_);

    //second: send data
    tag = 7000;
    MPI_Send(&aerosenddata[0], lengthSend, MPI_DOUBLE, INCAleader_, tag, intercomm_);
  }

  MPI_Barrier(intercomm_);

  return;
}


/*----------------------------------------------------------------------*
| read restart information for given time step              ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::ReadRestart()
{
  // read the restart information, set vectors and variables ---
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    tsi_->ReadRestart(restart);
  }

  return;
}


/*----------------------------------------------------------------------*
| single fields are tested                                  ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem::Instance()->AddFieldTest(tsi_->StructureField()->CreateFieldTest());
  DRT::Problem::Instance()->AddFieldTest(tsi_->ThermoField()->CreateFieldTest());

  DRT::Problem::Instance()->TestAll(comm);

  return;
}
