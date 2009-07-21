/*----------------------------------------------------------------------*/
/*!
\file adapter_scatra_fluid_coupling_algorithm.H

\brief Basis of all algorithms that perform a coupling between Navier-Stokes
       and (active or passive) scalar transport equations

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ADAPTER_SCATRA_FLUID_COUPLING_ALGORITHM_H
#define ADAPTER_SCATRA_FLUID_COUPLING_ALGORITHM_H

#include "../drt_adapter/adapter_fluid_base_algorithm.H"
#include "../drt_adapter/adapter_ale.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "adapter_coupling.H"

#include "../drt_fsi/fsi_algorithmbase.H"

namespace ADAPTER
{

  /// basis coupling algorithm for scalar transport with fluid velocity field
  /*!

    Base class for scalar transport problems coupled to Navier-Stokes velocity field.
    Derives from FluidBaseAlgorithm and ScaTraBaseAlgorithm.
    There are different subclasses that implement different coupling schemes
    (one-way coupling for the transport of passive scalars and
     fully coupled schemes for other multiphysics applications like combustion or electrochemistry).

    \author gjb
    \date 07/08
   */
  class ScaTraFluidCouplingAlgorithm : public FSI::AlgorithmBase,
        public ADAPTER::FluidBaseAlgorithm,
        public ADAPTER::ScaTraBaseAlgorithm
  {
  public:

    /// constructor using a Epetra_Comm
    explicit ScaTraFluidCouplingAlgorithm(
        Epetra_Comm& comm,              ///< communicator
        const Teuchos::ParameterList& prbdyn, ///< problem-specific parameters
        bool isale  ///< do we need an ALE formulation of the fields?
        );

    /// virtual destructor to support polymorph destruction
    virtual ~ScaTraFluidCouplingAlgorithm();

    /// outer level time loop (to be implemented by deriving classes)
    virtual void TimeLoop() = 0;

    /// read restart data
    virtual void ReadRestart(
        int step ///< step number where the calculation is continued
        );

  private:

    /// problem-specific parameter list
    const Teuchos::ParameterList& params_;
  };

} // namespace ADAPTER

#endif // ADAPTER_SCATRA_FLUID_COUPLING_ALGORITHM_H
#endif // CCADISCRET
