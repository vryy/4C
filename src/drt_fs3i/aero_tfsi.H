/*----------------------------------------------------------------------*/
/*!
\file aero_tfsi.H

\brief The aero code INCA is coupled with the monolithic TSI system
       based on communication with MPI. A loosely coupled approach is used.
       The time step starts with applying the time step from INCA to BACI.
       Then, INCA is solved and the interface loads (forces and heat fluxes)
       are transferred to BACI. After applying these data, the
       TSI system is solved. The time step is concluded by sending the resulting
       interface displacements and temperatures to INCA.

\level 3

\maintainer Georg Hammerl

*/


/*----------------------------------------------------------------------*
 | definitions                                              ghamm 12/11 |
 *----------------------------------------------------------------------*/
#ifndef AERO_TFSI_H
#define AERO_TFSI_H


/*----------------------------------------------------------------------*
 | headers                                                  ghamm 12/11 |
 *----------------------------------------------------------------------*/
#include "fs3i.H"
#include "Epetra_MpiComm.h"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <vector>
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_tsi.H"


// forward declaration
namespace TSI
{
  class Monolithic;
  class Partitioned;
}

namespace ADAPTER
{
  class Structure;
  class Thermo;
}


/*----------------------------------------------------------------------*
 |                                                          ghamm 12/11 |
 *----------------------------------------------------------------------*/
/// FS3I: Thermo-Fluid-Structure Interaction with Aero-code INCA
namespace FS3I
{
  //forward declaration
  namespace UTILS
  {
    class AeroCouplingUtils;
  }

  class AeroTFSI : public FS3I_Base
  {
  public:

    AeroTFSI(
    const Epetra_Comm& lcomm  ///< local epetra communicator
    );

    virtual ~AeroTFSI(){};

    /// write output of TSI system and additional restart info
    template<class A>
    void Output(Teuchos::RCP<A> solver);

    /// read restart
    virtual void ReadRestart();

    /// redistribute the  FPSI interface, if running on parallel. Just needed in the case of FPS3I
    void RedistributeInterface() {;} ;

    /// setup of the tsi system, mainly the two underlying discretizations
    virtual void SetupSystem();

    /// outer level TFSI time loop
    virtual void Timeloop();

    /// outer level TFSI time loop
    template<class A>
    void TimeloopT(Teuchos::RCP<A> solver);

    /// test results (if necessary)
    virtual void TestResults(
    const Epetra_Comm& comm  ///< local epetra communicator
    );


  protected:

    /// get information whether INCA will stop in this time step
    template<class A>
    void INCAfinshed(
    Teuchos::RCP<A> solver
    );

    /// Receive time step from INCA via MPI and set it in Baci
    template<class A>
    void GetAndSetDt(Teuchos::RCP<A> solver);

    /// Correct initial time in single fields
    void ResetInitialTimen();

    /// write status file to disk
    template<class A>
    void WriteStatusFile(Teuchos::RCP<A> solver);

    /// communication with INCA in order to receive forces and temperatures in the interface
    int ReceiveAeroData(
    std::vector<double>& aerodata  ///< raw data from INCA containing interface loads
    );

    /// communication with INCA in order to send disp, vel and temperatures in the interface
    void SendAeroData(
    std::vector<double>& aerosenddata  ///< interface data is stored in here
    );

    /// modifiy the received data into a more useful format
    double SplitData(
    std::vector<double>& aerodata,  ///< raw data from INCA
    std::vector<LINALG::Matrix<3,1> >& aerocoords,  ///< coordinates where loads are applied
    std::vector<LINALG::Matrix<4,1> >& aeroloads  ///< loads to be applied
    );

    /// received and processed data is applied to the structural/thermal solver
    void ApplyInterfaceData(
    Teuchos::RCP<Epetra_Vector> strfifc,  ///< interface forces due to the fluid (on full vector)
    Teuchos::RCP<Epetra_Vector> thrfifc  ///< interface thermal loads due to the fluid (on full vector)
    );

    /// get full vectors of structure
    void GetFullVectors(
    Teuchos::RCP<Epetra_Vector>& strfifc,
    Teuchos::RCP<Epetra_Vector>& thrfifc
    );

    /// get vectors with interface data of the solid at t_{n+1}
    void GetInterfaceDatanp(
    const int interf,
    Teuchos::RCP<const Epetra_Vector>& idispnp,
    Teuchos::RCP<const Epetra_Vector>& ivelnp,
    Teuchos::RCP<const Epetra_Vector>& itempnp
    );

    /// get vectors with interface data of the solid at t_{n}
    void GetInterfaceDatan(
    const int interf,
    Teuchos::RCP<const Epetra_Vector>& idispn,
    Teuchos::RCP<const Epetra_Vector>& iveln,
    Teuchos::RCP<const Epetra_Vector>& itempn
    );

    /// print chosen coupling strategy for TFSI to screen
    void PrintCouplingStrategy();

    /// compute wall clock timings
    void ComputeTiming(
    double& t_start,
    double& summation
    );

    /// summarize time monitoring and print it to screen
    void FinishTimeMonitoring(
    double CommTime,
    double INCATime,
    double BACITime,
    double BACICouplingTime
    );

    int localBACIleader_;  ///< leading processor in BACI that communicates with INCA in local communicator

    int INCAleader_;  ///< leading processor in INCA that communicates with BACI

    MPI_Comm intercomm_;  ///< intercommunicator between INCA and BACI for point to point communication

    const Epetra_Comm& lcomm_;  ///< local communicator in BACI stored as Epetra_Comm

    const MPI_Comm mpilcomm_;   ///< local communicator in BACI stored as MPI_Comm

  private:

    Teuchos::RCP<TSI::Monolithic> tsi_;  ///< Monolithic TSI system that is calculated
    Teuchos::RCP<ADAPTER::Structure> structure_;  ///< Structural field that is calculated (FSI)
    Teuchos::RCP<ADAPTER::Thermo> thermo_;  ///< Thermal field that is calculated
    Teuchos::RCP<Teuchos::ParameterList> probdyn_;   ///< parameter list of the problem

    INPAR::TSI::BaciIncaCoupling tfsi_coupling_;   ///< Coupling strategy for TFSI problem

    Teuchos::RCP<FS3I::UTILS::AeroCouplingUtils> aerocoupling_;  ///< Aero coupling helper class

    int stopflag_;   ///< stop flag that is sent from INCA

    double time_;  ///< time t_{n+1} of simulation

    bool diskoutput_;   ///< bool indicating if files are written to disk

  }; // AeroTFSI

} // namespace FS3I


/*----------------------------------------------------------------------*/
#endif  // AERO_TFSI_H
