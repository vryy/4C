/*!------------------------------------------------------------------------------------------------*
\file adapter_fluid_topopt_algo.H

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ADAPTER_FLUID_TOPOPT_ALGO_H_
#define ADAPTER_FLUID_TOPOPT_ALGO_H_

#include "adapter_fluid_base_algorithm.H"
#include "adapter_algorithmbase.H"
#include "adapter_topopt.H"

namespace ADAPTER
{

  /// basis coupling algorithm for the fields required for topology optimization
  /*!

    Base class coupling the fluid flow based on the Navier-Stokes equation with
    the optimization part containing the computation of the gradient of the
    objective with help of the adjoints and an optimization routine

    Various methods for the computation of the gradient via the adjoints
    (f.ex. discrete, continuous) and for the optimization routine are possible

    \author winklmaier
    \date 12/11
   */
  class FluidTopOptCouplingAlgorithm :
    public AlgorithmBase,
    public FluidBaseAlgorithm,
    public TopOptBaseAlgorithm
  {
  public:

    /// constructor using a Epetra_Comm
    explicit FluidTopOptCouplingAlgorithm(
        Epetra_Comm& comm,              ///< communicator
        const Teuchos::ParameterList& prbdyn ///< problem-specific parameters
        );

    /// virtual destructor to support polymorph destruction
    virtual ~FluidTopOptCouplingAlgorithm();

    /// outer level time loop (to be implemented by deriving classes)
    virtual void TimeLoop() = 0;

    /// read restart data
    virtual void ReadRestart(
        int step ///< step number where the calculation is continued
        );

    /// provide access to algorithm parameters
    virtual const Teuchos::ParameterList& AlgoParameters(){
        return params_;
    }

  private:

    /// problem-specific parameter list
    const Teuchos::ParameterList& params_;
  };

} // namespace ADAPTER

#endif /* ADAPTER_FLUID_TOPOPT_ALGO_H_ */
#endif  // #ifdef CCADISCRET
