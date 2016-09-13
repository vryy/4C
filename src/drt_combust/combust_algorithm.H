/*!-----------------------------------------------------------------------------------------------*
\file combust_algorithm.H

\brief combustion algorithm

\level 2

\maintainer Magnus Winter

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

    ///  Initialize this class
    virtual void Init(
        const Teuchos::ParameterList&   prbdyn,              ///< parameter list for global problem
        const Teuchos::ParameterList&   scatradyn,           ///< parameter list for scalar transport subproblem
        const Teuchos::ParameterList&   solverparams,        ///< parameter list for scalar transport solver
        const std::string&              disname = "scatra",  ///< name of scalar transport discretization
        const bool                      isale = false        ///< ALE flag
    );

    /// Setup this class
    virtual void Setup();

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

    /// perform result test
    void TestResults();

  protected:

    /// output of initial field
    void OutputInitialField();

    /// compute flame velocity
    const Teuchos::RCP<Epetra_Vector> ComputeFlameVel(const Teuchos::RCP<const Epetra_Vector>& velnp,
                                                      const Teuchos::RCP<const DRT::DofSet>& dofset
                                                      //const Teuchos::RCP<const Epetra_Map >& dbcmap
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

    /// Redistribute the scatra and fluid discretizations
    void Redistribute();

  private:

    /// compute center of mass of bubble etc
    void CenterOfMass();

    /// set velocity field in level-set algorithm
    void SetVelocityLevelSet(bool init=false);

    //! @name Fluid-G-function iteration variables
    // Fluid-G-function iteration counter
    int fgiter_;
    // maximum number of Fluid-G-function iterations
    int fgitermax_;
    //
    double convtol_;

    // velocity vector of last FG iteration step
    Teuchos::RCP<Epetra_Vector> velnpi_;
    // g-function vector of last FG iteration step
    Teuchos::RCP<Epetra_Vector> phinpi_;
    //@}

    /// switch for type of combustion problem
    INPAR::COMBUST::CombustionType combusttype_;

    /// the evaluation time ratio that triggers a redistribution
    double evaltimeratio_;

    /// boolean turning output to Gmsh on/off
    const bool gmshoutput_;

    /// combustion input parameter list
    const Teuchos::ParameterList& combustdyn_;
    //const Teuchos::ParameterList& xparams_;

    /// pointer to flame front holding information about interface (2D surface)
    Teuchos::RCP<COMBUST::FlameFront> flamefront_;

    /// indicate restarted problem
    bool restart_;

    /// merely perform flow generation and skip the soltion of the g-function field
    bool gen_flow_;

  };

} // namespace COMBUST

#endif // #ifndef COMBUST_ALGORITHM_H
