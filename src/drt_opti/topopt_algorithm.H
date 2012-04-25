/*!------------------------------------------------------------------------------------------------*
\file topopt_algorithm.H

\brief base algorithm for topology optimization of fluid domains

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef TOPOPT_ALGORITHM_H_
#define TOPOPT_ALGORITHM_H_


#include "../drt_adapter/ad_opt_fluid_algo.H"
#include "../drt_inpar/inpar_topopt.H"


namespace TOPOPT
{

  /*!

    \brief topology optimization algorithm class

    \detailed

    This combustion algorithm derives from the TopOpt Fluid Coupling Algorithm.
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
    explicit Algorithm(
        const Epetra_Comm& comm,
        const Teuchos::ParameterList& topoptdyn
    );

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// unused time loop
    void TimeLoop();

    /// perform the optimization
    void OptimizationLoop();

    /// restart topology optimization problem
    void Restart(int step, const int type);

  protected:

    /// prepare the optimization process
    void PrepareOptimization();

    /// check if the optimization process has finished
    bool OptimizationNotFinished();

    /// prepare the computation of a new fluid field
    void PrepareFluidField();

    /// solve primary equations (stationary or instationary, for example NS-equations)
    void DoFluidField();

    /// prepare the computation of a new adjoint field
    void PrepareAdjointField();

    /// solve adjoint equations
    void DoAdjointField();

    /// evaluate the gradient of the objective
    void GetGradient();

    /// perform an optimization step
    void DoOptimizationStep();

    /// update the porosity field
    void UpdatePorosity();

    /// write output
    void Output();

    /// compute volume on all processors
    double ComputeVolume();

  private:

    /// topology optimization input parameter list
    const Teuchos::ParameterList& topopt_;

    /// parallel stuff
    const int myrank_;

    /// convergence check
    int iter_;
    const int max_iter_;
    const double res_tol_;
    const double inc_tol_;
    INPAR::TOPOPT::ConvCheck conv_check_type_;

    /// pseudo-porosity function for topology optimization
    Teuchos::RCP<Epetra_Vector> poro_;

    /// objective fields
    double objective_ip_; // objective at new iteration in row map
    double objective_i_; // objective at old iteration in row map
  };

} // namespace TOPOPT

#endif /* TOPOPT_ALGORITHM_H_ */
