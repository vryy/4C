/*----------------------------------------------------------------------*/
/*!
\file inpar_potential.H

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_POTENTIAL_H
#define INPAR_POTENTIAL_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT { namespace INPUT { class ConditionDefinition; } }

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for interaction potentials
  namespace POTENTIAL
  {
    /// Type of effects which should be considered
    /// (this enum represents the input file parameter POTENTIAL_TYPE)
    enum PotentialType
    {
      potential_surface,          ///< only surface effects are considered
      potential_volume,           ///< only volume effects are considered
      potential_surfacevolume,    ///< surface and volume effects are considered
      potential_surface_fsi,      ///< only surface and fluid effects are considered
      potential_volume_fsi,       ///< only volume and fluid effects are considered
      potential_surfacevolume_fsi ///< surface, volume and fluid effects are considered
    };

    /// Type of approximation
    /// (this enum represents the input file parameter APPROXIMATION_TYPE)
    enum ApproximationType
    {
      approximation_none,       ///< no approximation
      approximation_surface,    ///< surface approximation
      approximation_point       ///< point approximation
    };

    /// Type of analytical solution
    /// (this enum represents the input file parameter ANALYTICALSOLUTION)
    enum SolutionType
    {
      solution_none,          ///< no analytical solution
      solution_sphere,        ///< sphere analytical solution
      solution_membrane       ///< membrane analytical solution
    };

    /// set the potential parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific potential conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

  }  // namespace POTENTIAL

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_POTENTIAL_H
