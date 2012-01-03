/*----------------------------------------------------------------------*/
/*!
\file aero_tfsi.H

\brief The aero code INCA is coupled with the monolithic TSI system
       based on communication with MPI. A loosely coupled approach is used.
       The time step starts with applying the time step from INCA to BACI.
       Then, INCA is solved and the interface loads (forces and heat fluxes)
       are transferred to BACI. After applying these data, the monolithic
       TSI system is solved. The time step is concluded by sending the resulting
       interface displacements and temperatures are send to INCA.

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-152537
</pre>
*/


/*----------------------------------------------------------------------*
 | definitions                                              ghamm 12/11 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef AERO_TFSI_H
#define AERO_TFSI_H


/*----------------------------------------------------------------------*
 | headers                                                  ghamm 12/11 |
 *----------------------------------------------------------------------*/
#include "fs3i.H"

#include "aero_tfsi_serv.H"
#include "../drt_tsi/tsi_monolithic.H"


/*----------------------------------------------------------------------*
 |                                                          ghamm 12/11 |
 *----------------------------------------------------------------------*/
/// FS3I: here --> Thermo-Fluid-Structure Interaction
namespace FS3I
{
  class AeroTFSI : public FS3I_Base
  {
  public:
    AeroTFSI(
    const Epetra_Comm& lcomm    ///< local epetra communicator
    );

    virtual ~AeroTFSI(){};

    /// outer level TFSI time loop
    virtual void Timeloop();

    /// test results (if necessary)
    virtual void TestResults(
    const Epetra_Comm& comm
    );

    /// read restart
    virtual void ReadRestart();

    /// setup of the tsi system, mainly the two underlying discretizations
    virtual void SetupSystem();

    /// get time step from INCA
    void GetTimeStep(
    vector<double>& timestep    ///< received time step size
    );

    /// apply time step from INCA to BACI
    void SetTimeStep(
    vector<double>& timestepsize    ///< time step size to be applied to BACI
    );

    /// apply time step from INCA to BACI
    void SetInitialTimeStepAndTime(
    vector<double>& timestepsize    ///< time step size to be applied to BACI
    );

    /// communication with INCA in order to receive forces and temperatures in the interface
    void ReceiveAeroData(
    vector<double>& aerodata    ///< raw data from INCA containing interface loads
    );

    /// communication with INCA in order to send disp, vel and temperatures in the interface
    void SendAeroData(
    vector<double>& aerosenddata    ///< interface data is stored in here
    );

    /// modifiy the reveived data into a more useful format
    void SplitData(
    vector<double> aerodata,    ///< raw data from INCA
    map<int, map<int, LINALG::Matrix<3,1> > >& aerocoords,   ///< coordinates where loads are applied
    map<int, map<int, LINALG::Matrix<4,1> > >& aeroforces    ///< loads to be applied
    );

    /// received and processed data is applied to the structural/thermal solver
    void ApplyInterfaceData(
    Teuchos::RCP<Epetra_Vector> iforce,   ///< interface forces due to the fluid
    Teuchos::RCP<Epetra_Vector> ithermoload   ///< interface thermal loads due to the fluid
    );


  protected:

  /// leading processor in BACI that communicates with INCA in local communicator
  int localBACIleader_;

  /// leading processor in INCA that communicates with BACI
  int INCAleader_;

  /// intercommunicator between INCA and BACI for point to point communication
  MPI_Comm intercomm_;

  /// local communicator in BACI
  const Epetra_Comm& lcomm_;


  private:

  Teuchos::RCP<TSI::Monolithic> tsi_;   ///< Monolithic TSI system that is calculated

  const Teuchos::ParameterList& sdynparams_;    ///< structural parameter list

  Teuchos::RCP<FS3I::UTILS::AeroCouplingUtils> aerocoupling_;  ///< Aero coupling helper class

  }; // AeroTFSI

} // namespace FS3I


/*----------------------------------------------------------------------*/
#endif  // AERO_TFSI_H
#endif  // CCADISCRET
