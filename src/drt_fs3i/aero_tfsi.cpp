/*----------------------------------------------------------------------*/
/*!
\file aero_tfsi.cpp

\brief The aero code INCA is coupled with the monolithic TSI system
       based on communication with MPI

\level 3

\maintainer Georg Hammerl

*/


/*----------------------------------------------------------------------*
 | headers                                                  ghamm 12/11 |
 *----------------------------------------------------------------------*/
#include "aero_tfsi.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_adapter/adapter_thermo.H"
#include "../drt_tsi/tsi_monolithic.H"
#include "../drt_tsi/tsi_utils.H"
#include "aero_tfsi_serv.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include <Teuchos_Time.hpp>


/*----------------------------------------------------------------------*
 | AeroTFSI                                                 ghamm 12/11 |
 *----------------------------------------------------------------------*/
FS3I::AeroTFSI::AeroTFSI(
  const Epetra_Comm& lcomm
  ) :
  FS3I_Base(),
  lcomm_(lcomm),
  mpilcomm_(dynamic_cast<const Epetra_MpiComm&>(lcomm).GetMpiComm()),
  tsi_(Teuchos::null),
  structure_(Teuchos::null),
  thermo_(Teuchos::null),
  stopflag_(0),
  time_(0.0),
  diskoutput_(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->IOParams(),"OUTPUT_BIN"))
{
  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn = DRT::Problem::Instance()->TSIDynamicParams();
  // coupling strategy for INCA and BACI
  tfsi_coupling_ = DRT::INPUT::IntegralValue<INPAR::TSI::BaciIncaCoupling>(tsidyn,"TFSI_COUPALGO");

  // leaders for each code are hard coded as follows
  localBACIleader_ = 0;
  INCAleader_ = 0;

  if(tfsi_coupling_ != INPAR::TSI::NoIncaFSI)
  {
    // check if INCA is called first when starting the coupled simulation
    int worldrank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldrank);

    if(worldrank == lcomm.MyPID())
      dserror("ERROR: INCA must!! be always started first in a coupled simulation");

    // intercommunicator is created; tag is important here because the same tag is used in INCA
    int tag = 0;
    int remoteleader = INCAleader_;
    MPI_Intercomm_create(mpilcomm_, 0, MPI_COMM_WORLD, remoteleader, tag, &intercomm_);

    MPI_Barrier(intercomm_);
  }

  // setup of the integrator on Baci side and the coupling helper class
  switch(tfsi_coupling_)
  {
  case INPAR::TSI::TFSI:
  {
    // setup of the discretizations, including clone strategy
    TSI::UTILS::SetupTSI(lcomm);

    probdyn_ = Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->TSIDynamicParams()));
    // create a TSI::Monolithic instance
    const Teuchos::ParameterList& sdynparams = DRT::Problem::Instance()->StructuralDynamicParams();
    tsi_ = Teuchos::rcp(new TSI::Monolithic(lcomm,sdynparams));

    // setup of the helper class
    aerocoupling_ = Teuchos::rcp(new FS3I::UTILS::AeroCouplingUtils(tsi_->StructureField()->Discretization(),
                                                                    tsi_->ThermoField()->Discretization()));
    break;
  }
  case INPAR::TSI::FSI:
  case INPAR::TSI::NoIncaFSI:
  {
    // access the structural discretization
    Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");

    // access structural dynamic params list which will be possibly modified while creating the time integrator
    probdyn_ = Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->StructuralDynamicParams()));

    Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure =
        Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(*probdyn_, *probdyn_, structdis));
    structure_ = structure->StructureField();

    // setup of the helper class
    aerocoupling_ = Teuchos::rcp(new FS3I::UTILS::AeroCouplingUtils(structure_->Discretization(), true));
    break;
  }
  case INPAR::TSI::ConjHeatTransfer:
  {
    // access the thermal discretization
    Teuchos::RCP<DRT::Discretization> thermodis = DRT::Problem::Instance()->GetDis("thermo");

    // access thermal dynamic params list
    probdyn_ = Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->ThermalDynamicParams()));

    Teuchos::RCP<ADAPTER::ThermoBaseAlgorithm> thermo =
        Teuchos::rcp(new ADAPTER::ThermoBaseAlgorithm(*probdyn_, thermodis));
    thermo_ = thermo->ThermoFieldrcp();

    // setup of the helper class
    aerocoupling_ = Teuchos::rcp(new FS3I::UTILS::AeroCouplingUtils(thermo_->Discretization(), false));
    break;
  }
  default:
    dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
    break;
  }

  PrintCouplingStrategy();
}


/*----------------------------------------------------------------------*
 | time loop of the coupled system                          ghamm 06/14 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::Timeloop()
{
  switch(tfsi_coupling_)
  {
  case INPAR::TSI::TFSI:
  {
    TimeloopT<TSI::Monolithic>(tsi_);
    break;
  }
  case INPAR::TSI::FSI:
  case INPAR::TSI::NoIncaFSI:
  {
    TimeloopT<ADAPTER::Structure>(structure_);
    break;
  }
  case INPAR::TSI::ConjHeatTransfer:
  {
    TimeloopT<ADAPTER::Thermo>(thermo_);
    break;
  }
  default:
    dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
    break;
  }

  return;
}


/*----------------------------------------------------------------------*
 | time loop of the coupled tfsi system                     ghamm 12/11 |
 *----------------------------------------------------------------------*/
template<class A>
void FS3I::AeroTFSI::TimeloopT(Teuchos::RCP<A> solver)
{
  if(DRT::Problem::Instance()->Restart() == 0 && diskoutput_ == true)
  {
    FILE* outFile;
    outFile = fopen("interfaceDisp.txt", "w");
    fclose(outFile);
    outFile = fopen("interfaceFlux.txt", "w");
    fclose(outFile);
    outFile = fopen("coupling_status_baci.dat", "w");
    fprintf(outFile, "#   itstep          time\n");
    fclose(outFile);
  }

  double time = Teuchos::Time::wallTime();
  double CommTime = 0.0;
  double INCATime = 0.0;
  double BACITime = 0.0;
  double BACICouplingTime = 0.0;

  // it is expected to have identical temperatures in the interface at the beginning
  // --> structural data can always be sent (and applied) to INCA
  if(tfsi_coupling_ != INPAR::TSI::NoIncaFSI)
  {
    for(int interf=0; interf<aerocoupling_->NumInterfaces(); interf++)
    {
      // sending initial data to INCA
      Teuchos::RCP<const Epetra_Vector> idispn = Teuchos::null;
      Teuchos::RCP<const Epetra_Vector> iveln = Teuchos::null;
      Teuchos::RCP<const Epetra_Vector> itempn = Teuchos::null;
      GetInterfaceDatan(interf, idispn, iveln, itempn);

      std::vector<double> aerosenddata;
      aerocoupling_->TransferStructValuesToFluid(interf, idispn, iveln, itempn, aerosenddata, false);

      ComputeTiming(time, BACICouplingTime);

      SendAeroData(aerosenddata);

      ComputeTiming(time, CommTime);
    }
  }

  // time need to be set correctly in case of pure structural or thermal problems
  ResetInitialTimen();
  time_ = solver->Time();

  // time loop
  while (solver->NotFinished() and stopflag_ == 0)
  {
    if(lcomm_.MyPID() == 0 && diskoutput_ == true)
    {
      FILE* outFile;
      outFile = fopen("interfaceDisp.txt", "a");
      fprintf(outFile, "%.8e  ", (solver->Time()/aerocoupling_->TimeScaling()));
      fclose(outFile);
      outFile = fopen("interfaceFlux.txt", "a");
      fprintf(outFile, "%.12e  ", solver->Time()/aerocoupling_->TimeScaling());
      fclose(outFile);
    }
    double flux_in = 0.0;
    double flux_struct = 0.0;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // INCA is working at the moment
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // get the time step from INCA and set it in the TSI system
    GetAndSetDt(solver);

    ComputeTiming(time, INCATime);

    // stopflag_ is received from INCA which can have the following meanings
    // 0: continue
    // 1: shutdown without computing the time step, just write restart output
    // 2: shutdown after computing last time step and output
    INCAfinshed(solver);

    if(stopflag_ == 1)
    {
      // write output to screen and files and leave time loop
      Output(solver);
      break;
    }

    // full vectors of structural and thermal field to be filled and applied to the solid problem
    Teuchos::RCP<Epetra_Vector> strfifc = Teuchos::null;
    Teuchos::RCP<Epetra_Vector> thrfifc =Teuchos::null;
    GetFullVectors(strfifc, thrfifc);

    std::vector<LINALG::Matrix<3,1> > aerocoords;
    std::vector<LINALG::Matrix<4,1> > aeroforces;

    ComputeTiming(time, BACITime);

    for(int interf=0; interf<aerocoupling_->NumInterfaces(); interf++)
    {
      std::vector<double> aerodata;
      // receive data from INCA
      const int nodeoffset = ReceiveAeroData(aerodata);

      ComputeTiming(time, CommTime);

      // fill data from INCA into suitable variables
      double flux_serial = SplitData(aerodata, aerocoords, aeroforces);

      // get new vectors for mapping the fluid interface data
      Teuchos::RCP<Epetra_Vector> iforce =Teuchos::null;
      if(aerocoupling_->MechCoupling())
        iforce = LINALG::CreateVector(*(aerocoupling_->GetInterfaceStructDis(interf)->DofRowMap()), true);
      Teuchos::RCP<Epetra_Vector> ithermoload = Teuchos::null;
      if(aerocoupling_->ThermoCoupling())
        ithermoload = LINALG::CreateVector(*(aerocoupling_->GetInterfaceThermoDis(interf)->DofRowMap()), true);

      // current displacement of the interface
      Teuchos::RCP<Epetra_Vector> idispn = Teuchos::null;
      if(aerocoupling_->MechCoupling())
      {
        if(tfsi_coupling_ == INPAR::TSI::TFSI)
          idispn = aerocoupling_->StrExtractInterfaceVal(interf, tsi_->StructureField()->Dispn());
        else
          idispn = aerocoupling_->StrExtractInterfaceVal(interf, structure_->Dispn());
      }

      // setup mortar coupling and transfer loads
      aerocoupling_->BuildMortarCoupling(interf, idispn, aerocoords, nodeoffset);
      aerocoupling_->TransferFluidLoadsToStruct(interf, aeroforces, iforce, ithermoload);

      // apply structural interface tractions to the structural field
      aerocoupling_->StrInsertInterfaceVal(interf, iforce, strfifc);
      // apply thermal interface heat flux to the thermal field
      aerocoupling_->ThrInsertInterfaceVal(interf, ithermoload, thrfifc);

      // incoming flux is summed over all procs
      double flux_global = 0.0;
      lcomm_.SumAll(&flux_serial, &flux_global, 1);
      flux_in += flux_global;

      // transferred flux is summed over all procs
      flux_serial = 0.0;
      flux_global = 0.0;
      if(aerocoupling_->ThermoCoupling())
      {
        for(int k=0; k<ithermoload->MyLength(); k++)
          flux_serial += (*ithermoload)[k];
      }
      lcomm_.SumAll(&flux_serial, &flux_global, 1);
      flux_struct += flux_global;

      ComputeTiming(time, BACICouplingTime);
    }

    if(lcomm_.MyPID() == 0 && diskoutput_ == true)
    {
      FILE* outFile = fopen("interfaceFlux.txt", "a");
      fprintf(outFile, "%.12e  ", flux_in);
      fprintf(outFile, "%.12e  ", flux_struct);
      if(flux_struct != 0.0)
        fprintf(outFile, "%.12e\n", flux_in/flux_struct);
      else
        fprintf(outFile, "\n");
      fclose(outFile);
    }

    // apply filled vectors to the TSI problem
    ApplyInterfaceData(strfifc, thrfifc);

    // counter and print header; predict solution of both fields
    solver->PrepareTimeStep();

    // store time for further use
    time_ = solver->Time();

    // write file to disk which is used for coupled restarts
    WriteStatusFile(solver);

    // TSI time step is solved
    solver->Solve();

    // calculate stresses, strains, energies for the structure
    solver->PrepareOutput();

    ComputeTiming(time, BACITime);

    // communication with AERO-code to send interface position and temperature
    // skip last sending procedure in the simulation if INCA has already finished
    if(stopflag_ != 2 and solver->NotFinished())
    {
      for(int interf=0; interf<aerocoupling_->NumInterfaces(); interf++)
      {
        // gather data and send it to INCA
        Teuchos::RCP<const Epetra_Vector> idispnp = Teuchos::null;
        Teuchos::RCP<const Epetra_Vector> ivelnp = Teuchos::null;
        Teuchos::RCP<const Epetra_Vector> itempnp = Teuchos::null;
        GetInterfaceDatanp(interf, idispnp, ivelnp, itempnp);

        std::vector<double> aerosenddata;
        aerocoupling_->TransferStructValuesToFluid(interf, idispnp, ivelnp, itempnp, aerosenddata);

        ComputeTiming(time, BACICouplingTime);

        SendAeroData(aerosenddata);

        ComputeTiming(time, CommTime);
      }
    }

    // update all single field solvers
    solver->Update();

    // write output to screen and files
    Output(solver);

  }  // NotFinished

  ComputeTiming(time, BACITime);

  FinishTimeMonitoring(CommTime, INCATime, BACITime, BACICouplingTime);

}  // TimeLoop


/*----------------------------------------------------------------------*
 | setup of the TSI system                                  ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::SetupSystem()
{
  if(tfsi_coupling_ == INPAR::TSI::TFSI)
  {
    tsi_->SetupSystem();
  }

  return;
}


/*----------------------------------------------------------------------*
 | Apply fluid interface loads to the TSI system           ghamm 12/11  |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::ApplyInterfaceData(
  Teuchos::RCP<Epetra_Vector> strfifc,
  Teuchos::RCP<Epetra_Vector> thrfifc
  )
{
  switch(tfsi_coupling_)
  {
  case INPAR::TSI::TFSI:
  {
    // apply structural interface tractions to the structural field
    tsi_->StructureField()->SetForceInterface(strfifc);

    // apply thermal interface heat flux to the thermal field
    tsi_->ThermoField()->SetForceInterface(thrfifc);
    break;
  }
  case INPAR::TSI::FSI:
  case INPAR::TSI::NoIncaFSI:
  {
    // apply structural interface tractions to the structural field
    structure_->SetForceInterface(strfifc);
    break;
  }
  case INPAR::TSI::ConjHeatTransfer:
  {
    // apply thermal interface heat flux to the thermal field
    thermo_->SetForceInterface(thrfifc);
    break;
  }
  default:
    dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
    break;
  }

  return;
}


/*----------------------------------------------------------------------*
 | Receive time step from INCA via MPI and set it in Baci   ghamm 12/11 |
 *----------------------------------------------------------------------*/
template<class A>
void FS3I::AeroTFSI::GetAndSetDt(Teuchos::RCP<A> solver)
{
  // this is just a dummy time for debug reasons when INCA is not available
  double timen = (solver->Step()+1)*0.001;
  if(tfsi_coupling_ != INPAR::TSI::NoIncaFSI)
  {
    // get time step from INCA
    if(lcomm_.MyPID() == localBACIleader_)
    {
      MPI_Status status;
      int tag_timestep = 3000;
      MPI_Recv(&timen, 1, MPI_DOUBLE, INCAleader_, tag_timestep, intercomm_, &status);
      timen *= aerocoupling_->TimeScaling();
    }
  }
  // broadcast time step to all processors in BACI
  lcomm_.Broadcast(&timen, 1 , localBACIleader_);

  const double dt = timen - time_;

  // set time step in underlying field; timen_ is overwritten in PrepareTimeStep()
  if(tfsi_coupling_ == INPAR::TSI::TFSI)
  {
    tsi_->StructureField()->SetDt(dt);
    tsi_->ThermoField()->SetDt(dt);
  }
  solver->SetDt(dt);

  return;
}


/*----------------------------------------------------------------------*
 | Correct initial time in single fields                    ghamm 06/14 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::ResetInitialTimen()
{
  // time need to be set correctly in case of pure structural or thermal problems
  switch(tfsi_coupling_)
  {
  case INPAR::TSI::TFSI:
  {
    // do nothing
    break;
  }
  case INPAR::TSI::FSI:
  case INPAR::TSI::NoIncaFSI:
  {
    structure_->SetTimen( structure_->Time() - structure_->Dt() );
    break;
  }
  case INPAR::TSI::ConjHeatTransfer:
  {
    thermo_->SetTimen( thermo_->Time() - thermo_->Dt() );
    break;
  }
  default:
    dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
    break;
  }

  return;
}


/*----------------------------------------------------------------------*
 | Receive flag from INCA via MPI whether to stop           ghamm 10/13 |
 *----------------------------------------------------------------------*/
template<class A>
void FS3I::AeroTFSI::INCAfinshed(
  Teuchos::RCP<A> solver
  )
{
  if(tfsi_coupling_ != INPAR::TSI::NoIncaFSI)
  {
    // get stop tag
    if(lcomm_.MyPID() == localBACIleader_)
    {
      MPI_Status status;
      int tag_stopflag = 2000;
      // note here: on INCA side an MPI_LOGICAL is used
      MPI_Recv(&stopflag_, 1, MPI_INT, INCAleader_, tag_stopflag, intercomm_, &status);

      if(stopflag_ != 0)
      {
        FILE* outFile = fopen("result_status_baci.inp", "w");
        // assumption that restart is always written when simulation is stopped by INCA
        fprintf(outFile, "%s\t%d", DRT::Problem::Instance()->OutputControlFile()->FileName().c_str(), solver->Step());
        fclose(outFile);
      }
    }
  }
  // broadcast stop flag to all processors in BACI
  lcomm_.Broadcast(&stopflag_, 1 , localBACIleader_);

  return;
}


/*----------------------------------------------------------------------*
 | write status file                                        ghamm 12/11 |
 *----------------------------------------------------------------------*/
template<class A>
void FS3I::AeroTFSI::WriteStatusFile(Teuchos::RCP<A> solver)
{
  // write output to file which is used for a coupled restart
  if(lcomm_.MyPID() == 0 && diskoutput_ == true)
  {
    FILE* outFile = fopen("coupling_status_baci.dat", "a");
    fprintf(outFile, "%10d", solver->Step());
    fprintf(outFile, " %20.14f\n", (solver->Time()/aerocoupling_->TimeScaling()));
    fclose(outFile);
  }

  return;
}


/*----------------------------------------------------------------------*
 | communication with AERO-code to receive the interface forces         |
 *----------------------------------------------------------------------*/
int FS3I::AeroTFSI::ReceiveAeroData(
  std::vector<double>& aerodata
  )
{
  int lengthRecv = 0;
  std::vector<double> receivebuf;
  if(tfsi_coupling_ != INPAR::TSI::NoIncaFSI)
  {
    // ==================================================================
    // intercommunicator is used to receive data on proc 0
    if(lcomm_.MyPID() == localBACIleader_)
    {
      MPI_Status status;
      //first: receive length of data from aero code
      int tag = 4000;
      MPI_Recv(&lengthRecv, 1, MPI_INT, INCAleader_, tag, intercomm_, &status);
      if(lengthRecv == 0)
        dserror("Length of data from INCA is not received correctly!!");

      //second: receive data
      tag = 5000;
      receivebuf.resize(lengthRecv);
      MPI_Recv(&receivebuf[0], lengthRecv, MPI_DOUBLE, INCAleader_, tag, intercomm_, &status);
    }
  }
  else
  {
    // for testing without INCA --> patch test
    // loading of a flat plate coupled at x=[0.1, 0.3], y=0 and z=[0,0.05]
    // with a constant load in y-direction and represented with 8 triangles
    lengthRecv = 8*13; // 8 tris with 13 doubles each
    receivebuf.resize(lengthRecv);

    //tri 1
    receivebuf[0] = 0.1;  receivebuf[1] = 0.0;  receivebuf[2] = 0.0;
    receivebuf[3] = 0.15; receivebuf[4] = 0.0;  receivebuf[5] = 0.0;
    receivebuf[6] = 0.1;  receivebuf[7] = 0.0;  receivebuf[8] = 0.05;

    receivebuf[9] = 0.0;
    receivebuf[10] = 1.0e-3 * 1.2;
    receivebuf[11] = 0.0;

    receivebuf[12] = 0.0; //1.0e-3 * 0.01;

    //tri 2
    receivebuf[13+0] = 0.1;  receivebuf[13+1] = 0.0;  receivebuf[13+2] = 0.05;
    receivebuf[13+3] = 0.15; receivebuf[13+4] = 0.0;  receivebuf[13+5] = 0.0;
    receivebuf[13+6] = 0.15; receivebuf[13+7] = 0.0;  receivebuf[13+8] = 0.05;

    receivebuf[13+9] = 0.0;
    receivebuf[13+10] = 1.0e-3 * 1.2;
    receivebuf[13+11] = 0.0;

    receivebuf[13+12] = 0.0; //1.0e-3 * 0.01;

    //tri 3
    receivebuf[0+26] = 0.15; receivebuf[1+26] = 0.0;  receivebuf[2+26] = 0.0;
    receivebuf[3+26] = 0.2;  receivebuf[4+26] = 0.0;  receivebuf[5+26] = 0.0;
    receivebuf[6+26] = 0.15; receivebuf[7+26] = 0.0;  receivebuf[8+26] = 0.05;

    receivebuf[9+26] = 0.0;
    receivebuf[10+26] = 1.0e-3 * 1.2;
    receivebuf[11+26] = 0.0;

    receivebuf[12+26] = 0.0; //1.0e-3 * 0.01;

    //tri 4
    receivebuf[13+0+26] = 0.15; receivebuf[13+1+26] = 0.0;  receivebuf[13+2+26] = 0.05;
    receivebuf[13+3+26] = 0.2;  receivebuf[13+4+26] = 0.0;  receivebuf[13+5+26] = 0.0;
    receivebuf[13+6+26] = 0.2;  receivebuf[13+7+26] = 0.0;  receivebuf[13+8+26] = 0.05;

    receivebuf[13+9+26] = 0.0;
    receivebuf[13+10+26] = 1.0e-3 * 1.2;
    receivebuf[13+11+26] = 0.0;

    receivebuf[13+12+26] = 0.0; //1.0e-3 * 0.01;

    //tri 5
    receivebuf[0+52] = 0.2;  receivebuf[1+52] = 0.0;  receivebuf[2+52] = 0.0;
    receivebuf[3+52] = 0.25; receivebuf[4+52] = 0.0;  receivebuf[5+52] = 0.0;
    receivebuf[6+52] = 0.2;  receivebuf[7+52] = 0.0;  receivebuf[8+52] = 0.05;

    receivebuf[9+52] = 0.0;
    receivebuf[10+52] = 1.0e-3 * 1.2;
    receivebuf[11+52] = 0.0;

    receivebuf[12+52] = 0.0; //1.0e-3 * 0.01;

    //tri 6
    receivebuf[13+0+52] = 0.2;  receivebuf[13+1+52] = 0.0;  receivebuf[13+2+52] = 0.05;
    receivebuf[13+3+52] = 0.25; receivebuf[13+4+52] = 0.0;  receivebuf[13+5+52] = 0.0;
    receivebuf[13+6+52] = 0.25; receivebuf[13+7+52] = 0.0;  receivebuf[13+8+52] = 0.05;

    receivebuf[13+9+52] = 0.0;
    receivebuf[13+10+52] = 1.0e-3 * 1.2;
    receivebuf[13+11+52] = 0.0;

    receivebuf[13+12+52] = 0.0; //1.0e-3 * 0.01;

    //tri 7
    receivebuf[0+78] = 0.25; receivebuf[1+78] = 0.0;  receivebuf[2+78] = 0.0;
    receivebuf[3+78] = 0.3;  receivebuf[4+78] = 0.0;  receivebuf[5+78] = 0.0;
    receivebuf[6+78] = 0.25; receivebuf[7+78] = 0.0;  receivebuf[8+78] = 0.05;

    receivebuf[9+78] = 0.0;
    receivebuf[10+78] = 1.0e-3 * 1.2;
    receivebuf[11+78] = 0.0;

    receivebuf[12+78] = 0.0; //1.0e-3 * 0.01;

    //tri 8
    receivebuf[13+0+78] = 0.25; receivebuf[13+1+78] = 0.0;  receivebuf[13+2+78] = 0.05;
    receivebuf[13+3+78] = 0.3;  receivebuf[13+4+78] = 0.0;  receivebuf[13+5+78] = 0.0;
    receivebuf[13+6+78] = 0.3;  receivebuf[13+7+78] = 0.0;  receivebuf[13+8+78] = 0.05;

    receivebuf[13+9+78] = 0.0;
    receivebuf[13+10+78] = 1.0e-3 * 1.2;
    receivebuf[13+11+78] = 0.0;

    receivebuf[13+12+78] = 0.0; //1.0e-3 * 0.01;
  }

  // ==================================================================
  // scatter received data in BACI to all processors
  lcomm_.Broadcast(&lengthRecv, 1 , localBACIleader_);

  if(lengthRecv%13 != 0)
    dserror("received data from INCA does not match the expected layout");

  // one data set (=tri) consists of 13 doubles (9*position + 3*forces + 1*heatflux)
  int numtris = lengthRecv / 13;
  int numproc = lcomm_.NumProc();

  // equal distribution of the tris for the beginning
  int triperproc = int(floor(numtris/numproc));
  std::vector<int> distribution(numproc, triperproc);
  // remaining tris are added to the procs in the beginning
  int remainingtris = numtris - triperproc*numproc;
  if(remainingtris >= numproc)
    dserror("distribution failed");
  for(int k=0; k<remainingtris; ++k)
  {
    distribution[k]++;
  }

  // compute offset in send buffer for data scattering
  int sumdoubles = 0;
  int displs[numproc];
  for(int k=0; k<numproc; ++k)
  {
    displs[k] = sumdoubles;
    // scaling with 13 in order to obtain number of doubles instead of number of tris
    distribution[k] *= 13;
    sumdoubles += distribution[k];
  }

  int nummydoubles = distribution[lcomm_.MyPID()];
  aerodata.resize(nummydoubles);

  // aero data is scattered
  MPI_Scatterv(&receivebuf[0], &distribution[0], &displs[0], MPI_DOUBLE, &aerodata[0], nummydoubles, MPI_DOUBLE, localBACIleader_, mpilcomm_);

  // get nodal distribution over procs for offsets in mortar (3 nodes per tri and 13 doubles per tri)
  for(int k=0; k<numproc; ++k)
  {
    distribution[k] = displs[k] * 3 / 13;
  }

  lcomm_.Barrier();

  // returns node offset from proc to proc for mortar
  return distribution[lcomm_.MyPID()];
}


/*----------------------------------------------------------------------*
 | cast data from INCA in a different format                ghamm 12/11 |
 *----------------------------------------------------------------------*/
double FS3I::AeroTFSI::SplitData(
  std::vector<double>& aerodata,
  std::vector<LINALG::Matrix<3,1> >& aerocoords,
  std::vector<LINALG::Matrix<4,1> >& aeroloads
  )
{
  double flux_in = 0.0;

  const double lengthscaling = aerocoupling_->LengthScaling();
  const double timescaling = aerocoupling_->TimeScaling();
  const double invtimescaling = 1.0 / timescaling;
  // heat flux scaling
  const double heatscalingfac = lengthscaling*lengthscaling*invtimescaling*invtimescaling*invtimescaling;
  // force scaling
  const double forcescalingfac = lengthscaling*invtimescaling*invtimescaling;

  // incoming data from INCA (x1y1z1 x2y2z2 x3y3z3 fx fy fz hf)
  // due to the Finite Volume scheme on fluid side a constant load per tri is a valid assumption
  size_t numtris = aerodata.size()/13;
  aerocoords.resize(numtris*3);
  aeroloads.resize(numtris);
  const size_t dim = 3;

  LINALG::Matrix<3,1> tmp1;
  LINALG::Matrix<4,1> tmp2;
  int nodecounter = 0;
  // loop over all tris
  for(size_t tri=0; tri<numtris; ++tri)
  {
    size_t numnode = 3;
    // loop over nodes of this tri --> count inverse in order to flip normal (pointing outward of fluid)
    for(int node=numnode-1; node>=0; --node)
    {
      // geometry has to be scaled with lengthscaling
      for(size_t d=0; d<dim; ++d)
      {
        // pull out fluid position
        tmp1(d) = aerodata[tri*13 + node*dim + d] * lengthscaling;
      }
      aerocoords[nodecounter] = tmp1;

      nodecounter++;
    }

    // pull out fluid forces
    for(size_t d=0; d<dim; ++d)
    {
      tmp2(d) = aerodata[tri*13 + 9 + d] * forcescalingfac;
    }

    // pull out fluid heat flux
    tmp2(3) = aerodata[tri*13 + 12] * heatscalingfac;

    aeroloads[tri] = tmp2;
    flux_in += tmp2(3);
  }

  return flux_in;
}


/*----------------------------------------------------------------------*
 | communication with AERO-code to send interface data      ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::SendAeroData(
  std::vector<double>& aerosenddata
  )
{
  if(tfsi_coupling_ != INPAR::TSI::NoIncaFSI)
  {
    // ==================================================================
    // intercommunicator is used to send data to INCA proc 0
    if(lcomm_.MyPID() == localBACIleader_)
    {
      int lengthSend = aerosenddata.size();
      //first: send length of data from aero code
      int tag = 6000;
      MPI_Send(&lengthSend, 1, MPI_INT, INCAleader_, tag, intercomm_);

      //second: send data
      tag = 7000;
      MPI_Send(&aerosenddata[0], lengthSend, MPI_DOUBLE, INCAleader_, tag, intercomm_);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | write output of TSI system and additional restart info   ghamm 12/13 |
 *----------------------------------------------------------------------*/
template<class A>
void FS3I::AeroTFSI::Output(Teuchos::RCP<A> solver)
{
  // call the problem dyn parameter list
  int uprestart = probdyn_->get<int>("RESTARTEVRY");
  if(stopflag_ != 2)
  {
    // do either normal output or forced output when stopflag == 1
    solver->Output(stopflag_);
  }
  else
  {
    // write standard output first
    solver->Output(false);
    // add missing restart information if necessary
    solver->Output(true);
  }

  // write latest step to file for proper restarting
  if(((uprestart != 0) and (solver->Step()%uprestart == 0)) or stopflag_!=0)
  {
    if(lcomm_.MyPID() == 0 && diskoutput_ == true)
    {
      FILE *outFile;
      outFile = fopen("result_status_baci.inp", "w");
      // StepOld() is written to file and thus it is necessary here to decrement step
      fprintf(outFile, "%s\t%d", DRT::Problem::Instance()->OutputControlFile()->FileName().c_str(), (solver->Step()-1));
      fclose(outFile);
    }
  }

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
    switch(tfsi_coupling_)
    {
    case INPAR::TSI::TFSI:
    {
      tsi_->ReadRestart(restart);
      break;
    }
    case INPAR::TSI::FSI:
    case INPAR::TSI::NoIncaFSI:
    {
      structure_->ReadRestart(restart);
      break;
    }
    case INPAR::TSI::ConjHeatTransfer:
    {
      thermo_->ReadRestart(restart);
      break;
    }
    default:
      dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
      break;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
| single fields are tested                                  ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::TestResults(const Epetra_Comm& comm)
{
  switch(tfsi_coupling_)
  {
  case INPAR::TSI::TFSI:
  {
    DRT::Problem::Instance()->AddFieldTest(tsi_->StructureField()->CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(tsi_->ThermoField()->CreateFieldTest());
    break;
  }
  case INPAR::TSI::FSI:
  case INPAR::TSI::NoIncaFSI:
  {
    DRT::Problem::Instance()->AddFieldTest(structure_->CreateFieldTest());
    break;
  }
  case INPAR::TSI::ConjHeatTransfer:
  {
    DRT::Problem::Instance()->AddFieldTest(thermo_->CreateFieldTest());
    break;
  }
  default:
    dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
    break;
  }

  DRT::Problem::Instance()->TestAll(comm);

  return;
}


/*----------------------------------------------------------------------*
| get full vectors of structure                             ghamm 06/14 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::GetFullVectors(Teuchos::RCP<Epetra_Vector>& strfifc, Teuchos::RCP<Epetra_Vector>& thrfifc)
{
  switch(tfsi_coupling_)
  {
  case INPAR::TSI::TFSI:
  {
    strfifc = LINALG::CreateVector(*tsi_->StructureField()->Discretization()->DofRowMap(), true);
    thrfifc = LINALG::CreateVector(*tsi_->ThermoField()->Discretization()->DofRowMap(), true);
    break;
  }
  case INPAR::TSI::FSI:
  case INPAR::TSI::NoIncaFSI:
  {
    strfifc = LINALG::CreateVector(*structure_->Discretization()->DofRowMap(), true);
    break;
  }
  case INPAR::TSI::ConjHeatTransfer:
  {
    thrfifc = LINALG::CreateVector(*thermo_->Discretization()->DofRowMap(), true);
    break;
  }
  default:
    dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
    break;
  }

  return;
}


/*----------------------------------------------------------------------*
| get vectors with interface data of the solid at t_{n+1}   ghamm 06/14 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::GetInterfaceDatanp(
  const int interf,
  Teuchos::RCP<const Epetra_Vector>& idispnp,
  Teuchos::RCP<const Epetra_Vector>& ivelnp,
  Teuchos::RCP<const Epetra_Vector>& itempnp
  )
{
  switch(tfsi_coupling_)
  {
  case INPAR::TSI::TFSI:
  {
    idispnp = aerocoupling_->StrExtractInterfaceVal(interf, tsi_->StructureField()->Dispnp());
    ivelnp = aerocoupling_->StrExtractInterfaceVal(interf, tsi_->StructureField()->Velnp());
    itempnp = aerocoupling_->ThrExtractInterfaceVal(interf, tsi_->ThermoField()->Tempnp());
    break;
  }
  case INPAR::TSI::FSI:
  case INPAR::TSI::NoIncaFSI:
  {
    idispnp = aerocoupling_->StrExtractInterfaceVal(interf, structure_->Dispnp());
    ivelnp = aerocoupling_->StrExtractInterfaceVal(interf, structure_->Velnp());
    break;
  }
  case INPAR::TSI::ConjHeatTransfer:
  {
    itempnp = aerocoupling_->ThrExtractInterfaceVal(interf, thermo_->Tempnp());
    break;
  }
  default:
    dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
    break;
  }

  return;
}


/*----------------------------------------------------------------------*
| get vectors with interface data of the solid at t_{n}     ghamm 06/14 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::GetInterfaceDatan(
  const int interf,
  Teuchos::RCP<const Epetra_Vector>& idispn,
  Teuchos::RCP<const Epetra_Vector>& iveln,
  Teuchos::RCP<const Epetra_Vector>& itempn
  )
{
  switch(tfsi_coupling_)
  {
  case INPAR::TSI::TFSI:
  {
    idispn = aerocoupling_->StrExtractInterfaceVal(interf, tsi_->StructureField()->Dispn());
    iveln = aerocoupling_->StrExtractInterfaceVal(interf, tsi_->StructureField()->Veln());
    itempn = aerocoupling_->ThrExtractInterfaceVal(interf, tsi_->ThermoField()->Tempn());
    break;
  }
  case INPAR::TSI::FSI:
  case INPAR::TSI::NoIncaFSI:
  {
    idispn= aerocoupling_->StrExtractInterfaceVal(interf, structure_->Dispn());
    iveln = aerocoupling_->StrExtractInterfaceVal(interf, structure_->Veln());
    break;
  }
  case INPAR::TSI::ConjHeatTransfer:
  {
    itempn = aerocoupling_->ThrExtractInterfaceVal(interf, thermo_->Tempn());
    break;
  }
  default:
    dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
    break;
  }

  return;
}


/*----------------------------------------------------------------------*
| print chosen coupling strategy for TFSI to screen         ghamm 02/13 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::PrintCouplingStrategy()
{
  if(lcomm_.MyPID() == 0)
  {
    switch(tfsi_coupling_)
    {
    case INPAR::TSI::TFSI:
    {
      IO::cout << "\nThermo-Fluid-Structure-Interaction with mortar coupling using dual shape functions for Lagrange multiplier. \n" << IO::endl;
      break;
    }
    case INPAR::TSI::FSI:
    case INPAR::TSI::NoIncaFSI:
    {
      IO::cout << "\nFluid-Structure-Interaction with mortar coupling using dual shape functions for Lagrange multiplier. \n" << IO::endl;
      break;
    }
    case INPAR::TSI::ConjHeatTransfer:
    {
      IO::cout << "\nConjugate Heat Transfer problem with mortar coupling using dual shape functions for Lagrange multiplier. \n" << IO::endl;
      break;
    }
    default:
      dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
      break;
    }

    IO::cout << "\n using the following scaling parameters for length: " << aerocoupling_->LengthScaling()
        << " and time: " << aerocoupling_->TimeScaling() << "\n" << IO::endl;
  }

  return;
}

/*----------------------------------------------------------------------*
| compute wall clock timings                                ghamm 06/14 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::ComputeTiming(double& time, double& summation)
{
  // compute passed time
  double t_end = Teuchos::Time::wallTime()-time;
  summation += t_end;

  // reset start point
  time = Teuchos::Time::wallTime();

  return;
}


/*----------------------------------------------------------------------*
| summarize time monitoring and print it to screen          ghamm 06/14 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::FinishTimeMonitoring(
  double CommTime,
  double INCATime,
  double BACITime,
  double BACICouplingTime
  )
{
  // brief (and not totally correct) time monitoring at the end of the simulation
  double ParCommTime = 0.0;
  double ParINCATime = 0.0;
  double ParBACITime = 0.0;
  double ParBACICouplingTime = 0.0;
  lcomm_.MaxAll(&CommTime,&ParCommTime,1);
  lcomm_.MaxAll(&INCATime,&ParINCATime,1);
  lcomm_.MaxAll(&BACITime,&ParBACITime,1);
  lcomm_.MaxAll(&BACICouplingTime,&ParBACICouplingTime,1);
  if(lcomm_.MyPID()==0)
  {
    std::cout<<"Brief (and rough) time monitoring for a coupled simulation:" << std::endl;
    std::cout <<"Communication time: "<< ParCommTime << std::endl;
    std::cout <<"INCA solution time: "<< ParINCATime << std::endl;
    std::cout <<"BACI solution time: "<< ParBACITime << std::endl;
    std::cout <<"BACI coupling time: "<< ParBACICouplingTime << std::endl;
  }

  return;
}

