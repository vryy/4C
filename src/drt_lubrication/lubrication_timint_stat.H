/*--------------------------------------------------------------------------*/
/*!
\file lubrication_timint_stat.H

\brief solution algorithm for stationary problems

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_LUBRICATION_LUBRICATION_TIMINT_STAT_H_
#define SRC_DRT_LUBRICATION_LUBRICATION_TIMINT_STAT_H_


#include "lubrication_timint_implicit.H"

namespace LUBRICATION
{
  class TimIntStationary : public TimIntImpl
  {
   public:
    /// Standard Constructor
    TimIntStationary(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    /// Destructor
    virtual ~TimIntStationary();

    /// initialize time integration scheme
    virtual void Init();

    /// Update the solution after convergence of the nonlinear iteration.
    /// Current solution becomes old solution of next timestep.
    virtual void Update(const int num = 0);

    /// read restart data
    void ReadRestart(int step);

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #raten_
    //! from scratch by using the newly updated #prenp_. The method
    //! respects the Dirichlet DOFs which are not touched.
    //! This method is necessary for certain predictors
    //! (like #PredictConstTempConsistRate)
    void UpdateIterIncrementally();

   protected:
    /// don't want = operator and cctor
    TimIntStationary operator=(const TimIntStationary& old);

    /// copy constructor
    TimIntStationary(const TimIntStationary& old);

    /// set time parameter for element evaluation
    void SetElementTimeParameter() const;

    //! set time for evaluation of Neumann boundary conditions
    void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params);

    /// add actual Neumann loads with time factor
    void AddNeumannToResidual();

    /// add parameters specific for time-integration scheme
    void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false);

    /// return the right time-scaling-factor for the true residual
    double ResidualScaling() const { return 1.0; }

  };  // class TimIntStationary

}  // namespace LUBRICATION

#endif /* SRC_DRT_LUBRICATION_LUBRICATION_TIMINT_STAT_H_ */
