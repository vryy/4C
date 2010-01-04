/*!-----------------------------------------------------------------------------------------------*
\file combust_algorithm.H

  \brief combustion algorithm

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef COMBUST_ALGORITHM_H
#define COMBUST_ALGORITHM_H

#include <Epetra_Comm.h>

#include "../drt_inpar/inpar_combust.H"
#include "../drt_adapter/adapter_scatra_fluid_coupling_algorithm.H"
#include "combust_interface.H"
#include "combust_flamefront.H"


/// namespace of the combustion module in BACI
namespace COMBUST
{

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
    explicit Algorithm(Epetra_Comm& comm, const Teuchos::ParameterList& combustdyn);

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// time loop of algorithm for dynamic combustion problem
    void TimeLoop();

    /// algorithm for static combustion problem
    void SolveStationaryProblem();

  protected:

    /// FG iteration number
    int FGiter() const { return fgiter_; }

    /// reinitialize G-function
    void ReinitializeGfunc();

    /// overwrite Navier-Stokes velocity field (used for test cases)
    const Teuchos::RCP<Epetra_Vector> OverwriteFluidVel();

    /// compute flame velocity
    const Teuchos::RCP<Epetra_Vector> ComputeFlameVel(const Teuchos::RCP<Epetra_Vector>& velnp,
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
    void UpdateFGIteration();

    /// take current results for converged and save for next time step
    void UpdateTimeStep();

    /// write output
    void Output();

    /// print mass conservation check on screen
    void printMassConservationCheck(const double volume_start, const double volume_end);

  private:

    //! @name Fluid-G-function iteration variables
    // Fluid-G-function iteration counter
    int fgiter_;
    // maximum number of Fluid-G-function iterations
    int fgitermax_;
    // L2-norm of velocity in Fluid-G-function iteration
    double fgvelnormL2_;
    // L2-norm of g-function in Fluid-G-function iteration
    double fggfuncnormL2_;
    //
    double convtol_;
    // velocity vector of last iteration step
    Teuchos::RCP<Epetra_Vector> velnpi_;
    // velocity vector of new iteration step
    Teuchos::RCP<Epetra_Vector> velnpip_;
    // g-function vector of last iteration step
    Teuchos::RCP<Epetra_Vector> phinpi_;
    // g-function vector of new iteration step
    Teuchos::RCP<Epetra_Vector> phinpip_;
    //@}

    /// switch for type of combustion problem
    INPAR::COMBUST::CombustionType combusttype_;

    /// options for reinitialization of G-function;
    INPAR::COMBUST::ReInitialActionGfunc reinitializationaction_;

    /// interval for reinitialization (every "reinitinterval_" time steps)
    int reinitinterval_;

    /// combustion input parameter list
    const Teuchos::ParameterList& combustdyn_;
//    const Teuchos::ParameterList& xparams_;

    /// pointer to interfacehandle holding information about interface (3D integration cells)
    Teuchos::RCP<COMBUST::InterfaceHandleCombust> interfacehandle_;

    /// pointer to flame front holding information about interface (2D surface)
    Teuchos::RCP<COMBUST::FlameFront> flamefront_;

  };

} // namespace COMBUST

#endif // #ifndef COMBUST_ALGORITHM_H
#endif // #ifdef CCADISCRET
