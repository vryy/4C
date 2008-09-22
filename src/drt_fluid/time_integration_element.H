/*!----------------------------------------------------------------------
\file time_integration_element.H
\brief routines for fluid (in)stationary time-integration,

     including instationary formulations

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef TIME_INTEGRATION_SCHEME_H
#define TIME_INTEGRATION_SCHEME_H

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_geometry/vector_definitions.H"


namespace FLD
{

//! namespace for all time discretization routines of the fluid field
namespace TIMEINT_THETA_BDF2
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
  BlitzVec3 GetOldPartOfRighthandside(
      const BlitzVec3&           veln,       ///< velocity      n
      const BlitzVec3&           velnm,      ///< velocity      n-1
      const BlitzVec3&           accn,       ///< acceleration  n
      const FLUID_TIMEINTTYPE    timealgo,   ///< time discretization type
      const double               dta,        ///< delta t (time step size)
      const double               theta       ///< factor for one step theta scheme
  );
  
  /*!
   * Compute constant used in the weak form integration 
   * \return timefac time integration factor
   */
  double ComputeTimeFac(
      const FLUID_TIMEINTTYPE    timealgo,   ///< time discretization type
      const double               dt,         ///< delta t (time step size)
      const double               theta       ///< factor for one step theta scheme
      );
}
}

#endif  // #ifndef TIME_INTEGRATION_SCHEME_H
#endif  // #ifdef CCADISCRET
