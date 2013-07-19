/*!-----------------------------------------------------------------------------------------------*
\file combust_algorithm.H

  \brief combustion algorithm

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef COMBUST_ALGORITHM_H
#define COMBUST_ALGORITHM_H

#include <Epetra_Comm.h>

#include "../drt_inpar/inpar_combust.H"
#include "../drt_adapter/adapter_scatra_fluid_coupling_algorithm.H"

/// namespace of the combustion module in BACI
namespace COMBUST
{

  class ReinitializationPDE;

  /*!

    \brief combustion algorithm class

    \detailed

    This combustion algorithm derives from the Scatra-Fluid Coupling Algorithm. Therefore it has
    both a Fluid Base Algorithm and a Scatra Base Algorithm.

    This algorithm handles the interaction between the different fields involved in a combustion
    problem, namely a fluid field and and a scalar transport field (level set/G-function). So far
    there is only one coupling algorithm.

    Until now there is only the implicit One-step-theta time integration scheme available for both
    fields.

    \author henke
    \date 06/08
   */
  class Algorithm : public ADAPTER::ScaTraFluidCouplingAlgorithm
  {
  public:

    /// constructor using a Epetra_Comm
    explicit Algorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& combustdyn, const Teuchos::ParameterList& solverparams);

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// time loop of algorithm for dynamic combustion problem
    void TimeLoop();

    /// algorithm for static combustion problem
    void SolveStationaryProblem();

    /// get initial field by solving stationary problem first
    void SolveInitialStationaryProblem();

    /// restart combustion problem
    void Restart(int step,
                    const bool restartscatrainput = false,
                    const bool restartturbinflow  = false);

  protected:

    /// reinitialize G-function
    void ReinitializeGfunc();

    /// call reinitialization
    void DoReinitialization();

    /// call reinitialization (old signed distance method)
    void ReinitializeGfuncSignedDistance();

    /// overwrite Navier-Stokes velocity field (used for test cases)
    const Teuchos::RCP<Epetra_Vector> OverwriteFluidVel();

    /// compute flame velocity
    const Teuchos::RCP<Epetra_Vector> ComputeFlameVel(const Teuchos::RCP<const Epetra_Vector>& velnp,
                                                      const Teuchos::RCP<const DRT::DofSet>& dofset
                                                      //const Teuchos::RCP<const Epetra_Map >& dbcmap
    );

    /// compute reinit velocity for advective scheme
    const Teuchos::RCP<Epetra_Vector> ComputeReinitVel(/*const Teuchos::RCP<Epetra_Vector>& velnp,*/
                                                      const Teuchos::RCP<const DRT::DofSet>& dofset
    );


    /// tests if the FGI iteration has converged
    bool NotConvergedFGI();

    /// prepare a new time step
    void PrepareTimeStep();

    /// start a new Fluid-G-function iteration
    void PrepareFGIteration();

    /// solve Navier-Stokes equations for current time step
    void DoFluidField();

    /// solve G-function equation for current time step
    void DoGfuncField();

    /// update fields in FGI iteration;
    void UpdateInterface();

    /// take current results for converged and save for next time step
    void UpdateTimeStep();

    /// write output
    void Output();

    /// compute volume on all processors
    double ComputeVolume();

    /// print mass conservation check on screen
    void printMassConservationCheck(const double volume_start, const double volume_end);

    /// manipulate the fluid field passed to the scatra
    const Teuchos::RCP<const Epetra_Vector> ManipulateFluidFieldForGfunc(const Teuchos::RCP<const Epetra_Vector>& convel,
                                                                         const Teuchos::RCP<const DRT::DofSet>& dofset);

    /// correct the volume of the minus(-) domain after reinitialization
    void CorrectVolume(const double targetvol, const double currentvol);

    /// Redistribute the scatra and fluid discretizations
    void Redistribute();

  private:

    //! @name Fluid-G-function iteration variables
    // Fluid-G-function iteration counter
    int fgiter_;
    // maximum number of Fluid-G-function iterations
    int fgitermax_;
    //
    double convtol_;

    // boolean indicating time step before reinitialization
    bool stepbeforereinit_;
    // boolean indicating time step of reinitialization
    bool stepreinit_;

    double volume_start_;


    // velocity vector of last FG iteration step
    Teuchos::RCP<Epetra_Vector> velnpi_;
    // g-function vector of last FG iteration step
    Teuchos::RCP<Epetra_Vector> phinpi_;
    //@}

    /// switch for type of combustion problem
    INPAR::COMBUST::CombustionType combusttype_;
    /// options for reinitialization of G-function;
    INPAR::COMBUST::ReInitialActionGfunc reinitaction_;
    /// initial reinitialization step
    bool reinitinitial_;
    /// interval for reinitialization (every 'reinitinterval_' time steps)
    int reinitinterval_;
    /// switch for reinitialization only within a band around the interface
    const bool reinitband_;
    /// band width for reinitialization (maximum G-function value)
    const double reinitbandwidth_;
    // boolian for reinitialization output
    const bool reinit_output_;
    /// reinitialization acception status
    bool reinitialization_accepted_;
    /// volume correction after reinitialization
    bool volcorrection_;
    /// the evaluation time ratio that triggers a redistribution
    double evaltimeratio_;

    /// switch for replacing the velocity of nodes at a certain distance of the interface by an approximate interface velocity
    const bool extract_interface_vel_;
    /// number of layers around the interface with unmodified velocity
    const int convel_layers_;

    /// boolean turning output to Gmsh on/off
    const bool gmshoutput_;

    /// combustion input parameter list
    const Teuchos::ParameterList& combustdyn_;
    //const Teuchos::ParameterList& xparams_;

    /// pointer to flame front holding information about interface (2D surface)
    Teuchos::RCP<COMBUST::FlameFront> flamefront_;

    /// pointer to PDE based reinitializer object
    Teuchos::RCP<COMBUST::ReinitializationPDE> reinit_pde_;

    /// type how the interface is transported
    INPAR::COMBUST::TransportVel transport_vel_;

    /// func no for interface transport
    int transport_vel_no_;

    /// indicate restarted problem
    bool restart_;

  };

} // namespace COMBUST

#endif // #ifndef COMBUST_ALGORITHM_H
