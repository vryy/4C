/*----------------------------------------------------------------------*/
/*!
\file inpar_fpsi.H
\brief fpsi parameters
\level 1
\maintainer Andreas Rauch
*/
/*----------------------------------------------------------------------*/

#ifndef INPAR_FPSI_H_
#define INPAR_FPSI_H_

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

/*----------------------------------------------------------------------*
 | Coupling Methods                                                                  |
 *----------------------------------------------------------------------*/
enum _FPSI_COUPLING
{
  fpsi_monolithic_plain,
  partitioned
};

namespace INPAR
{
  /// Input parameters for fluid porous structure interaction
  namespace FPSI
  {
    // type of partitioned coupling for FPSI problems
    enum PartitionedCouplingMethod
    {
      RobinNeumann,
      monolithic,
      nocoupling
    };

    // type of norm to check for convergence
    enum ConvergenceNorm
    {
      absoluteconvergencenorm,            // compare absolute value with single tolerance
      absoluteconvergencenorm_sys_split,  // compare absolute value with correction of systemsize
                                          // with different tolerances for each field
      relativconvergencenorm_sys  // compare relative value with correction of systemsize with
                                  // single tolerance
    };

    // type of norm to check for convergence
    enum BinaryOp
    {
      bop_and,
      bop_or
    };

    enum FluidFieldHierachy
    {
      fluid,
      porofluid
    };

    /// set the fpsi parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific fpsi conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace FPSI

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif /* INPAR_FPSI_H_ */
