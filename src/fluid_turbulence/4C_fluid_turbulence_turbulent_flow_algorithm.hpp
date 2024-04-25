/*----------------------------------------------------------------------*/
/*! \file

\brief algorithm for turbulent flows with separate inflow section


\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TURBULENCE_TURBULENT_FLOW_ALGORITHM_HPP
#define FOUR_C_FLUID_TURBULENCE_TURBULENT_FLOW_ALGORITHM_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_base_algorithm.hpp"
#include "4C_fluid_result_test.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class FluidDiscretExtractor;

  class TurbulentFlowAlgorithm
  {
   public:
    /// constructor
    explicit TurbulentFlowAlgorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& fdyn);

    /// virtual destructor
    virtual ~TurbulentFlowAlgorithm() = default;

    /// do time loop
    void TimeLoop();

    /// read restart
    /// only during inflow generation
    void ReadRestart(const int restart);

    /// do result check
    Teuchos::RCP<DRT::ResultTest> DoResultCheck()
    {
      return fluidalgo_->FluidField()->CreateFieldTest();
    };

   private:
    /// method to transfer inflow velocity from inflow discret to compete discret
    void TransferInflowVelocity();

    /// discretization of the compete domain
    Teuchos::RCP<DRT::Discretization> fluiddis_;
    /// discretization of the separate part
    Teuchos::RCP<DRT::Discretization> inflowdis_;
    /// object for a redistributed evaluation of of the separated part
    Teuchos::RCP<FluidDiscretExtractor> inflowgenerator_;
    /// instance of fluid algorithm
    Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluidalgo_;
    /// instance of fluid inflow algorithm
    Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> inflowfluidalgo_;
    /// number of time steps
    int step_;
    /// number of development steps
    int numtimesteps_;
    /// velocity/pressure at time n+1 to be transferred to the complete fluid field
    Teuchos::RCP<Epetra_Vector> velnp_;
  };

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
