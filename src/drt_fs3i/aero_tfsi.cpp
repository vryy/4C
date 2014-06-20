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
#include "../drt_tsi/tsi_utils.H"
#include "aero_tfsi_serv.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_io/io_control.H"
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
  mpilcomm_(dynamic_cast<const Epetra_MpiComm&>(lcomm).GetMpiComm()),
  tsi_(Teuchos::null),
  stopflag_(0)
{
  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn = DRT::Problem::Instance()->TSIDynamicParams();
  // coupling strategy for INCA and BACI
  tfsi_coupling_ = DRT::INPUT::IntegralValue<INPAR::TSI::BaciIncaCoupling>(tsidyn,"TFSI_COUPALGO");

  // leaders for each code are hard coded as follows
  localBACIleader_ = 0;
  INCAleader_ = 0;

#ifdef INCA_COUPLING
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
#endif

  // setup of the discretizations, including clone strategy
  TSI::UTILS::SetupTSI(lcomm);

  // create a TSI::Monolithic instance
  const Teuchos::ParameterList& sdynparams = DRT::Problem::Instance()->StructuralDynamicParams();
  tsi_ = Teuchos::rcp(new TSI::Monolithic(lcomm,sdynparams));

  // setup of the helper class
  aerocoupling_ = Teuchos::rcp(new FS3I::UTILS::AeroCouplingUtils(tsi_->StructureField()->Discretization(),
                                                                  tsi_->ThermoField()->Discretization()));

  PrintCouplingStrategy();

}


/*----------------------------------------------------------------------*
 | time loop of the coupled tfsi system                     ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::Timeloop()
{
  if(DRT::Problem::Instance()->Restart() == 0)
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

  double t_start = 0.0;
  double t_end = 0.0;
  double CommTime = 0.0;
  double INCATime = 0.0;
  double BACITime = 0.0;
  double BACICouplingTime = 0.0;

  // it is expected to have identical temperatures in the interface at the beginning
  // --> structural data can always be sent (and applied) to INCA
#ifdef INCA_COUPLING
  for(int interf=0; interf<aerocoupling_->NumInterfaces(); interf++)
  {
    t_start = Teuchos::Time::wallTime();

    // sending initial data to INCA
    Teuchos::RCP<const Epetra_Vector> idispnStart = aerocoupling_->StrExtractInterfaceVal(interf, tsi_->StructureField()->Dispn());
    Teuchos::RCP<const Epetra_Vector> ivelnStart = aerocoupling_->StrExtractInterfaceVal(interf, tsi_->StructureField()->Veln());
    Teuchos::RCP<const Epetra_Vector> ithermoloadStart = aerocoupling_->ThrExtractInterfaceVal(interf, tsi_->ThermoField()->ExtractTempn());

    std::vector<double> aerosenddata;
    switch(tfsi_coupling_)
    {
    case INPAR::TSI::conforming :
    {
      aerocoupling_->TransferStructValuesToFluidConforming(interf, idispnStart, ithermoloadStart, aerosenddata, false);
      break;
    }
    case INPAR::TSI::mortar_mortar_dual :
    {
      aerocoupling_->TransferStructValuesToFluid(interf, idispnStart, ivelnStart, ithermoloadStart, aerosenddata, false);
      break;
    }
    default:
      dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
      break;
    }

    t_end = Teuchos::Time::wallTime()-t_start;
    BACICouplingTime += t_end;

    t_start = Teuchos::Time::wallTime();

    SendAeroData(aerosenddata);

    t_end = Teuchos::Time::wallTime()-t_start;
    CommTime += t_end;
  }
#endif
  t_start = Teuchos::Time::wallTime();

  // time loop
  while (tsi_->NotFinished() and stopflag_ == 0)
  {
    if(lcomm_.MyPID() == 0)
    {
      FILE* outFile;
      outFile = fopen("interfaceDisp.txt", "a");
      fprintf(outFile, "%.8e  ", (tsi_->Time()/aerocoupling_->TimeScaling()));
      fclose(outFile);
      outFile = fopen("interfaceFlux.txt", "a");
      fprintf(outFile, "%.12e  ", tsi_->Time()/aerocoupling_->TimeScaling());
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
    GetAndSetDt();

    t_end = Teuchos::Time::wallTime()-t_start;
    INCATime += t_end;
    t_start = Teuchos::Time::wallTime();

    // stopflag_ is received from INCA which can have the following meanings
    // 0: continue
    // 1: shutdown without computing the time step, just write restart output
    // 2: shutdown after computing last time step and output
    INCAfinshed();

    if(stopflag_ == 1)
    {
      // write output to screen and files and leave time loop
      Output();
      break;
    }

    // full vectors of structural and thermal field to be filled and applied to the TSI problem
    Teuchos::RCP<Epetra_Vector> strfifc = LINALG::CreateVector(*tsi_->StructureField()->Discretization()->DofRowMap(), true);
    Teuchos::RCP<Epetra_Vector> thrfifc = LINALG::CreateVector(*tsi_->ThermoField()->Discretization()->DofRowMap(), true);
    std::vector<LINALG::Matrix<3,1> > aerocoords;
    std::vector<LINALG::Matrix<4,1> > aeroforces;
    t_end = Teuchos::Time::wallTime()-t_start;
    BACITime += t_end;
    for(int interf=0; interf<aerocoupling_->NumInterfaces(); interf++)
    {
      t_start = Teuchos::Time::wallTime();

      std::vector<double> aerodata;
      // receive data from INCA
      int nodeoffset = ReceiveAeroData(aerodata);

      t_end = Teuchos::Time::wallTime()-t_start;
      CommTime += t_end;
      t_start = Teuchos::Time::wallTime();

      // fill data from INCA into suitable variables
      double flux_serial = SplitData(aerodata, aerocoords, aeroforces);

      // get new vectors for mapping the fluid interface data
      Teuchos::RCP<Epetra_Vector> iforce = LINALG::CreateVector(*(aerocoupling_->GetInterfaceStructDis(interf)->DofRowMap()), true);
      Teuchos::RCP<Epetra_Vector> ithermoload = LINALG::CreateVector(*(aerocoupling_->GetInterfaceThermoDis(interf)->DofRowMap()), true);

      // current displacement of the interface
      Teuchos::RCP<Epetra_Vector> idispn = aerocoupling_->StrExtractInterfaceVal(interf, tsi_->StructureField()->Dispn());

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

      // transfered flux is summed over all procs
      flux_serial = 0.0;
      flux_global = 0.0;
      for(int k=0; k<ithermoload->MyLength(); k++)
        flux_serial += (*ithermoload)[k];
      lcomm_.SumAll(&flux_serial, &flux_global, 1);
      flux_struct += flux_global;

      t_end = Teuchos::Time::wallTime()-t_start;
      BACICouplingTime += t_end;
    }

    t_start = Teuchos::Time::wallTime();

    if(lcomm_.MyPID() == 0)
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
    PrepareTimeStep();

    // TSI time step is solved
    tsi_->NewtonFull();

    // calculate stresses, strains, energies
    tsi_->PrepareOutput();

    t_end = Teuchos::Time::wallTime()-t_start;
    BACITime += t_end;

    // communication with AERO-code to send interface position and temperature
    // skip last sending procedure in the simulation if INCA has already finished
    if(stopflag_ != 2 and tsi_->NotFinished())
    {
      for(int interf=0; interf<aerocoupling_->NumInterfaces(); interf++)
      {
        t_start = Teuchos::Time::wallTime();

        // gather data and send it to INCA
        Teuchos::RCP<Epetra_Vector> idispnp = aerocoupling_->StrExtractInterfaceVal(interf,tsi_->StructureField()->Dispnp());
        Teuchos::RCP<Epetra_Vector> ivelnp = aerocoupling_->StrExtractInterfaceVal(interf,tsi_->StructureField()->Velnp());
        Teuchos::RCP<Epetra_Vector> itempnp = aerocoupling_->ThrExtractInterfaceVal(interf, tsi_->ThermoField()->ExtractTempnp());

        std::vector<double> aerosenddata;
        switch(tfsi_coupling_)
        {
        case INPAR::TSI::conforming :
        {
          aerocoupling_->TransferStructValuesToFluidConforming(interf, idispnp, itempnp, aerosenddata);
          break;
        }
        case INPAR::TSI::mortar_mortar_dual :
        {
          aerocoupling_->TransferStructValuesToFluid(interf, idispnp, ivelnp, itempnp, aerosenddata);
          break;
        }
        default:
          dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
          break;
        }

        t_end = Teuchos::Time::wallTime()-t_start;
        BACICouplingTime += t_end;
        t_start = Teuchos::Time::wallTime();

        SendAeroData(aerosenddata);

        t_end = Teuchos::Time::wallTime()-t_start;
        CommTime += t_end;
      }
    }

    t_start = Teuchos::Time::wallTime();

    // update all single field solvers
    tsi_->Update();

    // write output to screen and files
    Output();

  }  // NotFinished

  t_end = Teuchos::Time::wallTime()-t_start;
  BACITime += t_end;

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

  // apply thermal interface heat flux to the thermal field
  tsi_->ThermoField()->SetForceInterface(thrfifc);

  return;
}


/*----------------------------------------------------------------------*
 | Receive time step from INCA via MPI and set it in Baci   ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::GetAndSetDt()
{
  // this is just a dummy time for debug reasons when INCA is not available
  double timen = (tsi_->Step()+1)*0.001;
#ifdef INCA_COUPLING
  // get time step from INCA
  if(lcomm_.MyPID() == localBACIleader_)
  {
    MPI_Status status;
    int tag_timestep = 3000;
    MPI_Recv(&timen, 1, MPI_DOUBLE, INCAleader_, tag_timestep, intercomm_, &status);
    timen *= aerocoupling_->TimeScaling();
  }
#endif
  // broadcast time step to all processors in BACI
  lcomm_.Broadcast(&timen, 1 , localBACIleader_);

  const double dt = timen - tsi_->Time();

  // set time step in TSI; timen_ is overwritten in PrepareTimeStep()
  SetDt(dt);


  return;
}


/*----------------------------------------------------------------------*
 | Receive flag from INCA via MPI whether to stop           ghamm 10/13 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::INCAfinshed()
{
#ifdef INCA_COUPLING
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
      fprintf(outFile, "%s\t%d", DRT::Problem::Instance()->OutputControlFile()->FileName().c_str(), tsi_->Step());
      fclose(outFile);
    }
  }
#endif
  // broadcast stop flag to all processors in BACI
  lcomm_.Broadcast(&stopflag_, 1 , localBACIleader_);

  return;
}


/*----------------------------------------------------------------------*
 | apply time step from INCA to BACI                        ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::SetDt(
  const double dt
  )
{
  tsi_->StructureField()->SetDt(dt);

  tsi_->ThermoField()->SetDt(dt);

  tsi_->SetDt(dt);

  return;
}


/*----------------------------------------------------------------------*
 | prepare time step, increment counters                    ghamm 12/11 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::PrepareTimeStep()
{
  tsi_->PrepareTimeStep();

  if(lcomm_.MyPID() == 0)
  {
    FILE* outFile = fopen("coupling_status_baci.dat", "a");
    fprintf(outFile, "%10d", tsi_->Step());
    fprintf(outFile, " %20.14f\n", (tsi_->Time()/aerocoupling_->TimeScaling()));
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
#ifdef INCA_COUPLING
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
#else
  // for testing with a flat plate
  lengthRecv = 26;
  receivebuf.resize(lengthRecv);

  //tri 1
  receivebuf[0] = 0.1;
  receivebuf[1] = 0.0;
  receivebuf[2] = 0.0;
  receivebuf[3] = 0.3;
  receivebuf[4] = 0.0;
  receivebuf[5] = 0.0;
  receivebuf[6] = 0.1;
  receivebuf[7] = 0.0;
  receivebuf[8] = 0.008;

  receivebuf[9] = 0.0;
  receivebuf[10] = 1.0e-3 * 1.2*3.0;
  receivebuf[11] = 0.0;

  receivebuf[12] = 1.0e-3 * 0.01*3.0;

  //tri 2
  receivebuf[13+0] = 0.1;
  receivebuf[13+1] = 0.0;
  receivebuf[13+2] = 0.008;
  receivebuf[13+3] = 0.3;
  receivebuf[13+4] = 0.0;
  receivebuf[13+5] = 0.0;
  receivebuf[13+6] = 0.3;
  receivebuf[13+7] = 0.0;
  receivebuf[13+8] = 0.008;

  receivebuf[13+9] = 0.0;
  receivebuf[13+10] = 1.0e-3 * 1.2*3.0;
  receivebuf[13+11] = 0.0;

  receivebuf[13+12] = 1.0e-3 * 0.01*3.0;
#endif

  // ==================================================================
  // scatter received data in BACI to all processors
  lcomm_.Broadcast(&lengthRecv, 1 , localBACIleader_);

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

  // returns node offset for mortar
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

  double lengthscaling = aerocoupling_->LengthScaling();
  double timescaling = aerocoupling_->TimeScaling();
  // heat flux scaling
  double heatscalingfac = lengthscaling*lengthscaling/timescaling/timescaling/timescaling;
  // force scaling
  double forcescalingfac = lengthscaling/timescaling/timescaling;

  // incoming data from INCA (x1y1z1 x2y2z2 x3y3z3 fx fy fz hf)
  // split into coordinates and loads per node
  // loads are equally divided onto the three nodes of a tri because constant value on fluid side (= Finite Volume scheme)
  size_t numtris = aerodata.size()/13;
  aerocoords.resize(numtris*3);
  aeroloads.resize(numtris*3);
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

        // pull out fluid forces
        tmp2(d) = aerodata[tri*13 + 9 + d] / 3.0 * forcescalingfac;
      }
      aerocoords[nodecounter] = tmp1;

      // pull out fluid heat flux
      tmp2(3) = aerodata[tri*13 + 12] / 3.0 * heatscalingfac;

      aeroloads[nodecounter] = tmp2;
      flux_in += tmp2(3);

      nodecounter++;
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

#endif

  return;
}


/*----------------------------------------------------------------------*
 | write output of TSI system and additional restart info   ghamm 12/13 |
 *----------------------------------------------------------------------*/
void FS3I::AeroTFSI::Output()
{
  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn = DRT::Problem::Instance()->TSIDynamicParams();
  int uprestart = tsidyn.get<int>("RESTARTEVRY");
  if(stopflag_ != 2)
  {
    // do either normal output or forced output when stopflag == 1
    tsi_->Output(stopflag_);
  }
  else
  {
    if((uprestart != 0) and (tsi_->Step()%uprestart == 0))
    {
      // normal restart output
      tsi_->Output(false);
    }
    else
    {
      // forced output is written
      tsi_->Output(true);
    }
  }

  // write latest step to file for proper restarting
  if(((uprestart != 0) and (tsi_->Step()%uprestart == 0)) or stopflag_!=0)
  {
    if(lcomm_.MyPID() == 0)
    {
      FILE *outFile;
      outFile = fopen("result_status_baci.inp", "w");
      fprintf(outFile, "%s\t%d", DRT::Problem::Instance()->OutputControlFile()->FileName().c_str(), tsi_->Step());
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
  if(lcomm_.MyPID() == 0)
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
    default:
      dserror("Coupling strategy %d is not yet implemented.", tfsi_coupling_);
      break;
    }

    IO::cout << "\n using the following scaling parameters for length: " << aerocoupling_->LengthScaling()
        << " and time: " << aerocoupling_->TimeScaling() << "\n" << IO::endl;
  }

  return;
}
