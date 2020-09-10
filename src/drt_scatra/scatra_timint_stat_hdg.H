/*----------------------------------------------------------------------*/
/*! \file
\brief solution algorithm for stationary problems

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_STAT_HDG_H
#define SCATRA_TIMINT_STAT_HDG_H

#include "scatra_timint_hdg.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"

namespace SCATRA
{
  class TimIntStationaryHDG : public TimIntHDG
  {
   public:
    /// Standard Constructor
    TimIntStationaryHDG(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    /// initialize time integration scheme
    void Init() override;

    /// compute values at intermediate time steps (required for generalized-alpha)
    void ComputeIntermediateValues() { return; };

    /// routine to return time integration specific parameters
    Teuchos::RCP<Teuchos::ParameterList> ScatraTimeParameterList()
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
    void SetElementTimeParameter(bool forcedincrementalsolver = false) const;

    /// set time for evaluation of Neumann boundary conditions
    void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params);

    /// calculate consistent initial conditions in compliance with initial scalar field
    /// this is not necessary for stationary calculations
    void CalcInitialTimeDerivative() { return; };

    /// do explicit predictor step (nothing to predict for stationary problems!)
    void ExplicitPredictor() const { return; };

    /// dynamic Smagorinsky model
    void DynamicComputationOfCs()
    {
      dserror("no turbulence in stationary flows!");
      return;
    };

    /// dynamic Vreman model
    void DynamicComputationOfCv()
    {
      dserror("no turbulence in stationary flows!");
      return;
    };

  };  // class TimIntStationaryHDG

}  // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_STAT_HDG_H