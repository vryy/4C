/*----------------------------------------------------------------------*/
/*!
\file inpar_cavitation.H

\brief Input parameters for cavitation

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*----------------------------------------------------------------------*/
#ifndef INPAR_CAVITATION_H
#define INPAR_CAVITATION_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/

namespace INPAR
{
  /// Input parameters for particle/cavitation problems
	namespace PARTICLE
	{
	  /// Contact type for particles
		enum ContactStrategy
		{
			None,
			Normal,
			NormalAndTang
		};

		/// Normal contact law
		enum NormalContact
		{
			LinSpring,
			Hertz,
			LinSpringDamp,
			LeeHerrmann,
			KuwabaraKono,
			Tsuji
		};
	}

  /// Input parameters for cavitation
  namespace CAVITATION
  {

    /// Type of coupling strategy for cavitation problems
    enum CouplingStrategyOverFields
    {
      OneWay,
      TwoWayMomentum,
      TwoWayFull
    };

    /// Strategy for calculation of void fraction of fluid elements
    enum VoidFractionCalculation
    {
      analytical_constpoly,
      analytical_quadraticpoly,
      gaussian_integration
    };


  }  // namespace CAVITATION

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_CAVITATION_H
