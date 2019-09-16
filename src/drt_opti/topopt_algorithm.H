/*---------------------------------------------------------------------*/
/*! \file

\brief base algorithm for topology optimization of fluid domains

\maintainer Martin Kronbichler

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef TOPOPT_ALGORITHM_H_
#define TOPOPT_ALGORITHM_H_


#include "../drt_adapter/ad_opt_fluid_algo.H"
#include "../drt_inpar/inpar_topopt.H"


namespace TOPOPT
{
  /*!

    \brief topology optimization algorithm class

    \detailed

    This topopt algorithm derives from the TopOpt Fluid Coupling Algorithm.
    Therefore it has a Fluid Base Algorithm and a TopOptForFluid Base Algorithm

    This algorithm handles the various interactions between the required steps
    during the topology optimization process of a fluid domain.

    Until now there is only the implicit One-step-theta time integration scheme available for both
    fields.

    \author winklmaier
    \date 12/11
   */
  class Algorithm : public ADAPTER::FluidTopOptCouplingAlgorithm
  {
   public:
    /// constructor using a Epetra_Comm
    explicit Algorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& topoptdyn);

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm() { return; };

    /// unused time loop
    void TimeLoop();

    /// perform the optimization
    void OptimizationLoop();

    /// restart topology optimization problem
    void Restart(const int step, const INPAR::TOPOPT::Restart type);

   protected:
    /// prepare the optimization process
    void PrepareOptimization();

    /// check if the optimization process has finished
    bool OptimizationFinished();

    /// prepare the computation of a new fluid field
    void PrepareFluidField();

    /// solve primary equations (stationary or instationary, for example NS-equations)
    void DoFluidField();

    /// prepare the computation of a new adjoint field
    void PrepareAdjointField();

    /// solve adjoint equations
    void DoAdjointField();

    /// compute gradient by finite differences
    void FDGradient(const int numFDPoints);

    /// evaluate the gradient of the objective
    void PrepareOptimizationStep();

    /// perform an optimization step
    void DoOptimizationStep();

    /// finish one optimization step
    void FinishOptimizationStep();

    /// write output
    void Output();

    /// update optimization procedure
    void Update();

   private:
    /// topology optimization input parameter list
    const Teuchos::ParameterList& topopt_;

    /// optimization algorithm
    Teuchos::RCP<TOPOPT::Optimizer>& optimizer_;

    bool doGradient_;  /// compute gradient(s) in current optimization step
    const INPAR::TOPOPT::GradientType gradienttype_;  /// type of gradient computation
    INPAR::TOPOPT::Restart restarttype_;              /// type of restart
  };

}  // namespace TOPOPT

#endif /* TOPOPT_ALGORITHM_H_ */
