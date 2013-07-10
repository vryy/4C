/*----------------------------------------------------------------------*/
/*!
\file aero_tfsi.cpp

\brief The aero code INCA is coupled with the monolithic TSI system
       based on communication with MPI

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
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
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_io/io_pstream.H"
#include <Teuchos_Time.hpp>

// define flag for debug reasons
#define INCA_COUPLING

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
  tfsi_coupling_ = DRT::INPUT::IntegralValue<INPAR::TSI::BaciIncaCoupling>(tsidyn,"TFSI_COUPALGO");
  additional_boundary_layer_ = DRT::INPUT::IntegralValue<INPAR::TSI::BaciIncaCoupling>(tsidyn,"TFSI_MORTAR_ADDITIONAL_BOUNDLAYER");
  fluid_forces_transfer_ = DRT::INPUT::IntegralValue<INPAR::TSI::BaciIncaCoupling>(tsidyn,"TFSI_FLUID_FORCES_TRANSFER");
  PrintCouplingStrategy();

#ifdef INCA_COUPLING
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
#endif

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

  switch(tfsi_coupling_)
  {
  case INPAR::TSI::conforming :
  case INPAR::TSI::proj_RBFI :
  {
    if(additional_boundary_layer_ == true)
      dserror("There is no additional boundary layer available for non-mortar coupled TFSI problems");
    break;

  }
  default:
    break;
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
  FILE *outFile;
  outFile = fopen("interfaceTemp.txt", "w");
  fclose(outFile);
  outFile = fopen("interfaceFlux.txt", "w");
  fclose(outFile);

  double t_start = 0.0;
  double t_end = 0.0;
  double CommTime = 0.0;
  double INCATime = 0.0;
  double BACITime = 0.0;

  // it is expected to have identical temperatures in the interface at the beginning
  // --> structural data can always be sent (and applied) to INCA
#ifdef INCA_COUPLING
  for(int interf=0; interf<aerocoupling_->NumInterfaces(); interf++)
  {
    std::map<int, LINALG::Matrix<3,1> > aerocoords;
    std::map<int, LINALG::Matrix<4,1> > aeroforces;
    int lengthphysicaldomain = 0;
    if( tfsi_coupling_ == INPAR::TSI::mortar_mortar_dual or
        tfsi_coupling_ == INPAR::TSI::proj_mortar_dual)
    {
      t_start = Teuchos::Time::wallTime();

      std::vector<double> aerodata;
      // receive data from INCA for physical domain
      ReceiveAeroData(aerodata);

      //fill data from INCA into suitable variables
      SplitData(aerodata, aerocoords, aeroforces,0);

      lengthphysicaldomain = aerocoords.size();

      aerodata.clear();
      // receive data from INCA for boundary layer around physical domain
      ReceiveAeroData(aerodata);

      SplitData(aerodata, aerocoords, aeroforces,lengthphysicaldomain);

      t_end = Teuchos::Time::wallTime()-t_start;
      CommTime += t_end;

      t_start = Teuchos::Time::wallTime();

      // current displacement of the interface
      Teuchos::RCP<Epetra_Vector> idispn = aerocoupling_->StrExtractInterfaceVal(interf, tsi_->StructureField()->ExtractDispn());

      aerocoupling_->BuildMortarCoupling(interf, idispn, aerocoords);

      t_end = Teuchos::Time::wallTime()-t_start;
      BACITime += t_end;
    }

    t_start = Teuchos::Time::wallTime();

    // sending initial data to INCA
    Teuchos::RCP<Epetra_Vector> idispnStart = aerocoupling_->StrExtractInterfaceVal(interf, tsi_->StructureField()->ExtractDispn());
//    Teuchos::RCP<Epetra_Vector> ivelnRestart = aerocoupling_->StrExtractInterfaceVal(interf, tsi_->StructureField()->ExtractVeln());
    Teuchos::RCP<Epetra_Vector> ithermoloadStart = aerocoupling_->ThrExtractInterfaceVal(interf, tsi_->ThermoField()->ExtractTempn());

    std::vector<double> aerosenddata;
    switch(tfsi_coupling_)
    {
    case INPAR::TSI::conforming :
    {
      aerocoupling_->PackData(interf, idispnStart, ithermoloadStart, aerosenddata, false);
      break;
    }
    case INPAR::TSI::mortar_mortar_dual :
    case INPAR::TSI::proj_mortar_dual :
    {
      aerocoupling_->TransferStructValuesToFluidDual(interf, lengthphysicaldomain, aerocoords, idispnStart,  ithermoloadStart, aerosenddata);
      break;
    }
    case INPAR::TSI::proj_RBFI :
    {
      aerocoupling_->PackDataRBFI(interf, idispnStart, ithermoloadStart, aerosenddata, false);
      // additional send which will be removed 15.05.2013
      SendAeroData(aerosenddata);
      break;
    }
    default:
      dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
      break;
    }

    t_end = Teuchos::Time::wallTime()-t_start;
    BACITime += t_end;

    t_start = Teuchos::Time::wallTime();

    SendAeroData(aerosenddata);

    t_end = Teuchos::Time::wallTime()-t_start;
    CommTime += t_end;
  }
#endif

  std::vector<double> timestep(1);
  // get the first time step from INCA; timen_ and dt has to be set correctly
  GetTimeStep(timestep);
  SetInitialTimeStepAndTime(timestep);

  // time loop
  while (tsi_->NotFinished())
  {
    if(lcomm_.MyPID() == 0)
    {
      FILE *outFile;
      outFile = fopen("interfaceTemp.txt", "a");
      fprintf(outFile, "%.8e  ", tsi_->Time());
      fclose(outFile);
      outFile = fopen("interfaceFlux.txt", "a");
      fprintf(outFile, "%.12e  ", tsi_->Time());
      fclose(outFile);
    }
    double flux_in = 0.0;
    double flux_struct = 0.0;

    // full vectors of structural and thermal field to be filled and applied to the TSI problem
    Teuchos::RCP<Epetra_Vector> strfifc = LINALG::CreateVector(*tsi_->StructureField()->Discretization()->DofRowMap(), true);
    Teuchos::RCP<Epetra_Vector> thrfifc = LINALG::CreateVector(*tsi_->ThermoField()->Discretization()->DofRowMap(), true);
    std::map<int, LINALG::Matrix<3,1> > aerocoords;
    std::map<int, LINALG::Matrix<4,1> > aeroforces;
    int lengthphysicaldomain[aerocoupling_->NumInterfaces()];
    for(int interf=0; interf<aerocoupling_->NumInterfaces(); interf++)
    {
      t_start = Teuchos::Time::wallTime();

      std::vector<double> aerodata;
      // receive data from INCA
      ReceiveAeroData(aerodata);

      //fill data from INCA into suitable variables
      flux_in += SplitData(aerodata, aerocoords, aeroforces, 0);

      // receive data from INCA for boundary layer around physical domain
      // only geometry is important, forces are dropped
      if( tfsi_coupling_ == INPAR::TSI::mortar_mortar_dual or
          tfsi_coupling_ == INPAR::TSI::proj_mortar_dual)
      {
        lengthphysicaldomain[interf] = aerocoords.size();
        aerodata.clear();
        ReceiveAeroData(aerodata);
        SplitData(aerodata, aerocoords, aeroforces,lengthphysicaldomain[interf]);
      }
      t_end = Teuchos::Time::wallTime()-t_start;
      CommTime += t_end;

      t_start = Teuchos::Time::wallTime();

      //get new vectors for mapping the fluid interface data
      Teuchos::RCP<Epetra_Vector> iforce = LINALG::CreateVector(*(aerocoupling_->GetInterfaceStructDis(interf)->DofRowMap()), true);
      Teuchos::RCP<Epetra_Vector> ithermoload = LINALG::CreateVector(*(aerocoupling_->GetInterfaceThermoDis(interf)->DofRowMap()), true);

      // current displacement of the interface
      Teuchos::RCP<Epetra_Vector> idispn = aerocoupling_->StrExtractInterfaceVal(interf, tsi_->StructureField()->ExtractDispn());

      switch(tfsi_coupling_)
      {
      case INPAR::TSI::proj_mortar_dual :
      {
        aerocoupling_->BuildMortarCoupling(interf, idispn, aerocoords);
        // no break here!
      }
      case INPAR::TSI::conforming :
      case INPAR::TSI::proj_RBFI :
      {
        aerocoupling_->ProjectForceOnStruct(interf, idispn, aerocoords, aeroforces, iforce, ithermoload);
        break;
      }
      case INPAR::TSI::mortar_mortar_dual :
      {
        aerocoupling_->BuildMortarCoupling(interf, idispn, aerocoords);
        aerocoupling_->TransferFluidLoadsToStructDual(interf, aeroforces, iforce, ithermoload);
        break;
      }
      default:
        dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
        break;
      }

      // apply structural interface tractions to the structural field
      aerocoupling_->StrApplyInterfaceVal(interf, iforce, strfifc);
      // apply thermal interface heat flux to the thermal field
      aerocoupling_->ThrApplyInterfaceVal(interf, ithermoload, thrfifc);

      double flux_serial = 0.0;
      double flux_global = 0.0;
      for(int k=0; k<ithermoload->MyLength(); k++)
        flux_serial += (*ithermoload)[k];
      lcomm_.SumAll(&flux_serial, &flux_global, 1);
      flux_struct += flux_global;

      t_end = Teuchos::Time::wallTime()-t_start;
      BACITime += t_end;
    }

    if(lcomm_.MyPID() == 0)
    {
      outFile = fopen("interfaceFlux.txt", "a");
      fprintf(outFile, "%.12e  ", flux_in);
      fprintf(outFile, "%.12e  ", flux_struct);
      fprintf(outFile, "%.12e\n", flux_in/flux_struct);
      fclose(outFile);
    }

    t_start = Teuchos::Time::wallTime();

    // apply filled vectors to the TSI problem
    ApplyInterfaceData(strfifc, thrfifc);

    // counter and print header; predict solution of both fields
    tsi_->PrepareTimeStep();

    // TSI time step is solved
    SolveTSIstep();

    // calculate stresses, strains, energies
    tsi_->PrepareOutput();

    t_end = Teuchos::Time::wallTime()-t_start;
    BACITime += t_end;

    // communication with AERO-code to send interface position and temperature
    // skip last sending procedure in the simulation because INCA has already finished
    if(tsi_->NotFinished())
    {
      for(int interf=0; interf<aerocoupling_->NumInterfaces(); interf++)
      {
        t_start = Teuchos::Time::wallTime();

        // gather data and send it to INCA
        Teuchos::RCP<Epetra_Vector> idispnp = aerocoupling_->StrExtractInterfaceVal(interf,tsi_->StructureField()->ExtractDispnp());
        // calculate velocity of the structure currently not needed
        Teuchos::RCP<Epetra_Vector> itemp = aerocoupling_->ThrExtractInterfaceVal(interf, tsi_->ThermoField()->ExtractTempnp());

        std::vector<double> aerosenddata;
        switch(tfsi_coupling_)
        {
        case INPAR::TSI::conforming :
        {
          aerocoupling_->PackData(interf, idispnp, itemp, aerosenddata);
          break;
        }
        case INPAR::TSI::mortar_mortar_dual :
        case INPAR::TSI::proj_mortar_dual :
        {
          aerocoupling_->TransferStructValuesToFluidDual(interf, lengthphysicaldomain[interf], aerocoords, idispnp, itemp, aerosenddata);
          break;
        }
        case INPAR::TSI::proj_RBFI :
        {
          aerocoupling_->PackDataRBFI(interf, idispnp, itemp, aerosenddata, true);
          break;
        }
        default:
          dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
          break;
        }

        t_end = Teuchos::Time::wallTime()-t_start;
        BACITime += t_end;

        t_start = Teuchos::Time::wallTime();

        SendAeroData(aerosenddata);

        t_end = Teuchos::Time::wallTime()-t_start;
        CommTime += t_end;
      }

      t_start = Teuchos::Time::wallTime();


      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      // INCA is working at the moment
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      // get the time step from INCA for the next time step in BACI
      GetTimeStep(timestep);

      t_end = Teuchos::Time::wallTime()-t_start;
      INCATime += t_end;

      // set the time step received from INCA
      SetTimeStep(timestep);
    }

    t_start = Teuchos::Time::wallTime();

    // update all single field solvers
    tsi_->Update();

    // write output to screen and files
    tsi_->Output();

    t_end = Teuchos::Time::wallTime()-t_start;
    BACITime += t_end;

  }  // NotFinished


  // brief (and not totally correct) time monitoring at the end of the simulation
  double ParCommTime = 0.0;
  double ParINCATime = 0.0;
  double ParBACITime = 0.0;
  lcomm_.MaxAll(&CommTime,&ParCommTime,1);
  lcomm_.MaxAll(&INCATime,&ParINCATime,1);
  lcomm_.MaxAll(&BACITime,&ParBACITime,1);
  if(lcomm_.MyPID()==0)
  {
    std::cout<<"Brief (and rough) time monitoring for a coupled simulation:" << std::endl;
    std::cout <<"Communication time: "<< ParCommTime << std::endl;
    std::cout <<"INCA solution time: "<< ParINCATime << std::endl;
    std::cout <<"BACI solution time: "<< ParBACITime << std::endl;
  }

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
  Teuchos::RCP<Epetra_Vector> strfifc,
  Teuchos::RCP<Epetra_Vector> thrfifc
  )
{
  // apply structural interface tractions to the structural field
  tsi_->StructureField()->SetForceInterface(strfifc);
  tsi_->StructureField()->PreparePartitionStep();

  // apply thermal interface heat flux to the thermal field
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
  std::vector<double>& timestep
  )
{
#ifdef INCA_COUPLING
  // get time step from INCA
  if(lcomm_.MyPID() == 0)
  {
    MPI_Status status;
    int tag_timestep = 3000;
    MPI_Recv(&timestep[0], 1, MPI_DOUBLE, INCAleader_, tag_timestep, intercomm_, &status);
  }
  // broadcast time step to all processors in BACI
  lcomm_.Broadcast(&timestep[0], 1 , localBACIleader_);
#endif

  return;
}


/*----------------------------------------------------------------------*
 | apply time step from INCA to BACI                        ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::SetTimeStep(
  std::vector<double>& timestepsize
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
  std::vector<double>& timestepsize
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
  std::vector<double>& aerodata
  )
{

#ifdef INCA_COUPLING
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

    FILE *outFile;
    outFile = fopen("interfaceTemp.txt", "a");
    fprintf(outFile, "%.10e ", aerodata[699]);
    fclose(outFile);

  }

  // ==================================================================
  // broadcast received data in BACI to all processors
  // again: length first, then data
  lcomm_.Broadcast(&lengthRecv, 1 , localBACIleader_);

  aerodata.resize(lengthRecv);

  //aero forces have to be redundant on all processors
  lcomm_.Broadcast(&aerodata[0], lengthRecv , localBACIleader_);

  lcomm_.Barrier();
#endif

  return;
}


/*----------------------------------------------------------------------*
 | cast data from INCA in a different format                ghamm 12/11 |
 *----------------------------------------------------------------------*/
double FS3I::AeroTFSI::SplitData(
  std::vector<double>& aerodata,
  std::map<int, LINALG::Matrix<3,1> >& aerocoords,
  std::map<int, LINALG::Matrix<4,1> >& aeroforces,
  int startingvalue
  )
{
  // without additional boundary layer leave here
  if(startingvalue != 0 and additional_boundary_layer_ == false)
    return 0.0;

  double flux_in = 0.0;

  // incoming data from INCA (xxxxyyyyzzzzffff) splitted into xyzxyzxyzxyz & ffff
  // or in case of forces that are included:
  // incoming data from INCA (xxxx yyyy zzzz hhhh fxfxfxfx fyfyfyfy fzfzfzfz) splitted into xyzxyzxyzxyz & hfxfyfz.hfxfyfz.hfxfyfz
  size_t length = aerodata.size();
  size_t numfluidpoints;
  if(fluid_forces_transfer_)
    numfluidpoints = length/7;
  else
    numfluidpoints = length/4;

  for(size_t out=startingvalue; out<(numfluidpoints+startingvalue); out++)
  {
    LINALG::Matrix<3,1> tmp1;
    LINALG::Matrix<4,1> tmp2(true);

    for(size_t in=0; in<3; in++)
    {
      tmp1(in)=aerodata[in*numfluidpoints + out - startingvalue];
    }
    aerocoords[out] = tmp1;

    // only geometry is important in case of the additional boundary layer --> no fluxes/forces
    if(startingvalue == 0)
    {
      // pull out heat flux
      tmp2(0)=aerodata[3*numfluidpoints + out];

      // pull out forces (if available)
      if(fluid_forces_transfer_)
      {
        for(size_t in=1; in<4; in++)
        {
          tmp2(in)=aerodata[(in+3)*numfluidpoints + out];
        }
      }

      aeroforces[out] = tmp2;
      flux_in += tmp2(0);
    }
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
#ifdef INCA_COUPLING
  // ==================================================================
  // intercommunicator is used to send data to INCA proc 0
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

#endif

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


/*----------------------------------------------------------------------*
| print chosen coupling strategy for TFSI to screen         ghamm 02/13 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::PrintCouplingStrategy()
{
  switch(tfsi_coupling_)
  {
  case INPAR::TSI::conforming :
  {
    IO::cout << "\nTFSI coupling: conforming meshes. \n" << IO::endl;
    break;
  }
  case INPAR::TSI::mortar_mortar_dual :
  {
    IO::cout << "\nTFSI coupling: mortar coupling with dual shape functions for Lagrange multiplier. \n" << IO::endl;
    break;
  }
  case INPAR::TSI::proj_mortar_dual:
  {
    IO::cout << "\nTFSI coupling: conservative projection from fluid to structure and mortar coupling \n"
        "with dual shape functions for Lagrange multiplier for structure to fluid. \n" << IO::endl;
    break;
  }
  case INPAR::TSI::proj_RBFI:
  {
    IO::cout << "\nTFSI coupling: conservative projection from fluid to structure and radial basis function \n"
        "interpolation for structure to fluid. \n" << IO::endl;
    break;
  }
  default:
    dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
    break;
  }
}
