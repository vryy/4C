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

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
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
  class Algorithm;
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

    /// read restart
    virtual void ReadRestart();

    /// setup of the tsi system, mainly the two underlying discretizations
    virtual void SetupSystem();

    /// outer level TFSI time loop
    virtual void Timeloop();

    /// test results (if necessary)
    virtual void TestResults(
    const Epetra_Comm& comm  ///< local epetra communicator
    );


  protected:

    /// get information whether INCA will stop in this time step
    bool INCAfinshed();

    /// get time step from INCA
    void GetTimeStep(
    std::vector<double>& timestep  ///< received time step size
    );

    /// apply time step from INCA to BACI
    void SetDt(
    const double dt  ///< time step size to be applied to BACI
    );

    /// apply time step from INCA to BACI
    void SetInitialDtAndTime(
    const double timestepsize  ///< time step size to be applied to BACI
    );

    /// communication with INCA in order to receive forces and temperatures in the interface
    void ReceiveAeroData(
    std::vector<double>& aerodata  ///< raw data from INCA containing interface loads
    );

    /// communication with INCA in order to send disp, vel and temperatures in the interface
    void SendAeroData(
    std::vector<double>& aerosenddata  ///< interface data is stored in here
    );

    /// modifiy the reveived data into a more useful format
    double SplitData(
    std::vector<double>& aerodata,  ///< raw data from INCA
    std::map<int, LINALG::Matrix<3,1> >& aerocoords,  ///< coordinates where loads are applied
    std::map<int, LINALG::Matrix<4,1> >& aeroforces,  ///< loads to be applied
    int startingvalue   ///< entry number in aerodata to start with
    );

    /// received and processed data is applied to the structural/thermal solver
    void ApplyInterfaceData(
    Teuchos::RCP<Epetra_Vector> strfifc,  ///< interface forces due to the fluid (on full vector)
    Teuchos::RCP<Epetra_Vector> thrfifc  ///< interface thermal loads due to the fluid (on full vector)
    );

    void PrintCouplingStrategy();   ///< print chosen coupling strategy for TFSI to screen

    /// call the solver for the current TSI time step
    void SolveTSIstep();


    int localBACIleader_;  ///< leading processor in BACI that communicates with INCA in local communicator

    int INCAleader_;  ///< leading processor in INCA that communicates with BACI

    MPI_Comm intercomm_;  ///< intercommunicator between INCA and BACI for point to point communication

    const Epetra_Comm& lcomm_;  ///< local communicator in BACI


  private:

    Teuchos::RCP<TSI::Algorithm> tsi_;  ///< Monolithic TSI system that is calculated

    INPAR::TSI::SolutionSchemeOverFields coupling_;  ///< Coupling type of the TSI system

    INPAR::TSI::BaciIncaCoupling tfsi_coupling_;   ///< Coupling strategy for TFSI problem

    Teuchos::RCP<FS3I::UTILS::AeroCouplingUtils> aerocoupling_;  ///< Aero coupling helper class

    bool additional_boundary_layer_;   ///< bool whether an additional boundary layer is used for mortar coupling

    bool fluid_forces_transfer_;   ///< bool whether fluid forces are transferred from INCA

    int stopflag_;   ///< stop flag that is sent from INCA

  }; // AeroTFSI

} // namespace FS3I


/*----------------------------------------------------------------------*/
#endif  // AERO_TFSI_H
