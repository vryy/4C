/*!
\file fluid_timint_ada.H

\brief Fluid wrapper for adaptive time stepping

<pre>
Maintainer: Raffaela Kruse
            kruse@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_TIMINT_ADA_H_
#define FLUID_TIMINT_ADA_H_

#include "../drt_adapter/ad_fld_wrapper.H"

#include <Teuchos_ParameterList.hpp>

namespace FLD
{

class FluidImplicitTimeInt;


class TimIntAda : public ADAPTER::FluidWrapper
{
public:

  explicit TimIntAda(
    Teuchos::ParameterList& params,
    Teuchos::RCP<Fluid> fluid) :
    FluidWrapper(fluid),
    dt_estimation_(DRT::INPUT::IntegralValue<INPAR::FLUID::AdaptiveTimeStepEstimator>(params,"ADAPTIVE_TIME_STEP_ESTIMATOR"))
  {
    if (dt_estimation_ == INPAR::FLUID::cfl_number)
      estimator_ = Teuchos::rcp(new EstimatorCFL(fluid));
    else
      dserror("No time step estimation active!");
  }

  ~TimIntAda(){}

  /// determine new time step and check, whether it is acceptable
  void Indicate
  (
    bool& accepted,
    double& dt_new
  );

  /*!
   * \brief Base class for time step estimation
   */
  class Estimator
  {
  public:

    Estimator(Teuchos::RCP<Fluid> fluid) :
      fluid_(fluid)
    {}

    virtual ~Estimator(){}

    /// check whether time step size is acceptable
    /// (error calculations - if necessary - in derived classes
    virtual bool TimeStepAccepted() = 0;

    /// calculate time step estimate
    virtual double CalculateDt(
      const double norm = 0.0 ///< current norm of local discretization error
    ) = 0;

    virtual void ResetStep(){}

    virtual void CalculateLocalError(){}

  protected:

    Teuchos::RCP<Fluid> fluid_;
  };

  /*!
   * \brief Time for time step estimation via CFL number
   */
  class EstimatorCFL : public Estimator
  {

  public:
    EstimatorCFL(Teuchos::RCP<Fluid> fluid) :
      Estimator(fluid),
      fluidimplicitintegration_(Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid))
    {
      if(fluidimplicitintegration_ == Teuchos::null)
        dserror("no fluidimplicitintegration available");
    }

    virtual ~EstimatorCFL(){}


    /// always accept time step size
    virtual bool TimeStepAccepted(){ return true; }

    /// calculate time step estimate
    virtual double CalculateDt(
      const double norm = 0.0 ///< current norm of local discretization error
    );

    virtual void ResetStep(){}

    virtual void CalculateLocalError(){}

    Teuchos::RCP<FluidImplicitTimeInt> fluidimplicitintegration_;
  };



private:

  INPAR::FLUID::AdaptiveTimeStepEstimator dt_estimation_;
  Teuchos::RCP<Estimator> estimator_;



};
}

#endif /* FLUID_TIMINT_ADA_H_ */
