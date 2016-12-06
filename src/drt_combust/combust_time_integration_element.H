/*!----------------------------------------------------------------------
\file combust_time_integration_element.H
\brief routines for fluid (in)stationary time-integration,

     including instationary formulations

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

*----------------------------------------------------------------------*/

#ifndef COMBUST_TIME_INTEGRATION_SCHEME_H
#define COMBUST_TIME_INTEGRATION_SCHEME_H

#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_fluid.H"

namespace COMBUST
{

//! namespace for all time discretization routines of the fluid field
namespace TIMEINT
{

  /*!
   * \brief Get the part of the righthandside belonging to the last timestep.
   *
   * set part of the rhs vector belonging to the old timestep
   *
   *   Stationary:
   *   Generalized Alpha:
   *
   *     hist = 0.0
   *
   *   One-step-Theta:
   *
   *     hist = veln + dta*(1-Theta)*accn
   *
   *   BDF2: for constant time step:
   *
   *     hist = 4/3 veln - 1/3 velnm
   *
   * \return hist          linear combination of old values
   *
   */
  LINALG::Matrix<3,1> GetOldPartOfRighthandside(
      const LINALG::Matrix<3,1>&                veln,     ///< velocity      n
      const LINALG::Matrix<3,1>&                velnm,    ///< velocity      n-1
      const LINALG::Matrix<3,1>&                accn,     ///< acceleration  n
      const INPAR::FLUID::TimeIntegrationScheme timealgo, ///< time discretization type
      const double                              dta,      ///< delta t (time step size)
      const double                              theta     ///< factor for one step theta scheme
  );

  /*!
   * Compute constant used in the weak form integration
   * \return timefac time integration factor
   */
  double ComputeTimeFac(
      const INPAR::FLUID::TimeIntegrationScheme timealgo, ///< time discretization type
      const double                              dt,       ///< delta t (time step size)
      const double                              theta,    ///< factor for one step theta scheme
      const double                              ga_alphaF,
      const double                              ga_alphaM,
      const double                              ga_gamma
      );
}

//! namespace for all time discretization routines of the fluid field for One-Step-Theta and BDF2 schemes
namespace UTILS
{

  /*!
   * \brief Get the part of the righthandside belonging to the last timestep.
   *
   * set part of the rhs vector belonging to the old timestep
   *
   *   Stationary:
   *
   *     hist = 0.0
   *
   *   One-step-Theta:
   *
   *     hist = veln + dta*(1-Theta)*accn
   *
   *   BDF2: for constant time step:
   *
   *     hist = 4/3 veln - 1/3 velnm
   *
   * \return hist          linear combination of old values
   *
   */
  LINALG::Matrix<3,1> GetOldPartOfRighthandside(
      const LINALG::Matrix<3,1>&        veln,       ///< velocity      n
      const LINALG::Matrix<3,1>&        velnm,      ///< velocity      n-1
      const LINALG::Matrix<3,1>&        accn,       ///< acceleration  n
      const INPAR::FLUID::TimeIntegrationScheme    timealgo,   ///< time discretization type
      const double               dta,        ///< delta t (time step size)
      const double               theta       ///< factor for one step theta scheme
  );

  /*!
   * Compute constant used in the weak form integration
   * \return timefac time integration factor
   */
  double ComputeTimeFac(
      const INPAR::FLUID::TimeIntegrationScheme    timealgo,   ///< time discretization type
      const double               dt,         ///< delta t (time step size)
      const double               theta       ///< factor for one step theta scheme
      );
}
}

#endif  // #ifndef TIME_INTEGRATION_SCHEME_H
