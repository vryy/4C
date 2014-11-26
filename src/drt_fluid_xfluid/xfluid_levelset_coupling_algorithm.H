/*----------------------------------------------------------------------*/
/*!
\file xfluid_levelset_coupling_algorithm.H

\brief Base algorithm of Xfluid utilizing level set for cut

<pre>
Maintainer: Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef XFLUID_LEVELSET_ALGORITHM_H
#define XFLUID_LEVELSET_ALGORITHM_H

#include "../drt_adapter/adapter_scatra_fluid_coupling_algorithm.H"

namespace XFLUIDLEVELSET
{

  /// Xfluid level set coupling
  /*!

    Base class of xfluid and level set coupling. Derives from FluidBaseAlgorithm
    and ScatraBaseAlgorithm.

    \author winter
    \date 09/14
   */
  class Algorithm : public ADAPTER::ScaTraFluidCouplingAlgorithm
  {
  public:

    /// constructor
    explicit Algorithm(
        const Epetra_Comm& comm,             ///< communicator
        const Teuchos::ParameterList& prbdyn,///< problem-specific parameter list
        const Teuchos::ParameterList& solverparams ///<solver-specific parameter list
        );

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// Two phase flow time loop
    void TimeLoop(){dserror("Timeloop() not implemented yet.");};

    /// Two phase flow stationary problem
    void SolveStationaryProblem();

    /// Test the results aquired
    void TestResults();

  protected:

    /// set scatra values required in fluid
    void SetScaTraValuesInFluid();

    /// write output
    void Output();



  private:

  };

} // namespace XFLUIDLEVELSET

#endif // XFLUID_LEVELSET_ALGORITHM_H
