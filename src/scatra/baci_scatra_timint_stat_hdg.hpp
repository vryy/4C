/*----------------------------------------------------------------------*/
/*! \file
\brief solution algorithm for stationary problems

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_TIMINT_STAT_HDG_HPP
#define BACI_SCATRA_TIMINT_STAT_HDG_HPP

#include "baci_config.hpp"

#include "baci_linalg_utils_sparse_algebra_create.hpp"
#include "baci_scatra_timint_hdg.hpp"

BACI_NAMESPACE_OPEN

namespace SCATRA
{
  class TimIntStationaryHDG : public TimIntHDG
  {
   public:
    /// Standard Constructor
    TimIntStationaryHDG(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    /// initialize time integration scheme
    void Init() override;

    /// compute values at intermediate time steps (required for generalized-alpha)
    void ComputeIntermediateValues() override { return; };

    /// routine to return time integration specific parameters
    Teuchos::RCP<Teuchos::ParameterList> ScatraTimeParameterList() override
    {
      dserror("Not yet implemented!");
      return Teuchos::null;
    }

   protected:
    /// don't want = operator and cctor
    TimIntStationaryHDG operator=(const TimIntStationaryHDG& old);

    /// copy constructor
    TimIntStationaryHDG(const TimIntStationaryHDG& old);

    /// set time parameter for element evaluation
    void SetElementTimeParameter(bool forcedincrementalsolver = false) const override;

    /// set time for evaluation of Neumann boundary conditions
    void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params) override;

    /// calculate consistent initial conditions in compliance with initial scalar field
    /// this is not necessary for stationary calculations
    void CalcInitialTimeDerivative() override { return; };

    /// do explicit predictor step (nothing to predict for stationary problems!)
    void ExplicitPredictor() const override { return; };

    /// dynamic Smagorinsky model
    void DynamicComputationOfCs() override
    {
      dserror("no turbulence in stationary flows!");
      return;
    };

    /// dynamic Vreman model
    void DynamicComputationOfCv() override
    {
      dserror("no turbulence in stationary flows!");
      return;
    };

  };  // class TimIntStationaryHDG

}  // namespace SCATRA

BACI_NAMESPACE_CLOSE

#endif  // SCATRA_TIMINT_STAT_HDG_H