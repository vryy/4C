/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for immersed

\level 1


*/

/*----------------------------------------------------------------------*/

#ifndef INPAR_IMMERSED_H_
#define INPAR_IMMERSED_H_

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

namespace INPAR
{
  namespace IMMERSED
  {
    /*----------------------------------------------------------------------*
     | Coupling Methods                                                     |
     *----------------------------------------------------------------------*/
    enum _IMMERSED_COUPLING
    {
      partitioned,
      monolithic
    };

    typedef enum _PARITIONED_SCHEME
    {
      cell_coupling_undefined = 0,
      cell_basic_sequ_stagg = 1,
      cell_iter_stagg_fixed_rel_param = 2,
      cell_iter_stagg_AITKEN_rel_param = 3,
    } PARITIONED_SCHEME;

    enum _IMMERSED_COUPLING_SCHEME
    {
      neumannneumann,
      dirichletneumann
    };

    enum _IMMERSED_PROJECTION
    {
      shapefunctions,
      mortar
    };

    enum _IMMERSED_RELAXATION
    {
      globally,
      selectively
    };

    enum _IMMERSED_NLNSOLVER
    {
      nlnsolver_stop,
      nlnsolver_continue
    };

    enum _IMMERSED_RELAXATIONPARAM
    {
      fixed,
      aitken
    };


    /// set the immersed parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific immersed conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace IMMERSED
}  // namespace INPAR
/*----------------------------------------------------------------------*/
#endif /* INPAR_IMMERSED_H_ */
