/*----------------------------------------------------------------------*/
/*!
\file inpar_ale.H

<pre>
 Maintainer: Matthias Mayr
             mayr@mhpc.mw.tum.de
             089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_ALE_H
#define INPAR_ALE_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

/*----------------------------------------------------------------------*/
namespace INPAR
{
  /// Input parameters for ALE
  namespace ALE
  {
    /// possible types of moving boundary simulation
    enum AleDynamic
    {
      none,              ///< default
      laplace_material,  ///< Laplacian smoothing based on material configuration
      laplace_spatial,   ///< Laplacian smoothing based on spatial configuration
      springs_material,  ///< use a spring analogy based on material configuration
      springs_spatial,   ///< use a spring analogy based on spatial configuration
      solid,             ///< nonlinear pseudo-structure approach
      solid_linear       ///< linear pseudo-structure approach
    };
    /// Handling of non-converged nonlinear solver
    enum DivContAct
    {
      divcont_stop,     ///< abort simulation
      divcont_continue  ///< continue nevertheless
    };

    /// mesh tying and mesh sliding algorithm
    enum MeshTying
    {
      no_meshtying,
      meshtying,
      meshsliding
    };

    /// Defines all valid parameters for acou problem
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// Defines ale specific conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace ALE

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_ALE_H
