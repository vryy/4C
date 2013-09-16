/*----------------------------------------------------------------------*/
/*!
\file inpar_fpsi.H

<pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef INPAR_FPSI_H_
#define INPAR_FPSI_H_

#include "inpar_parameterlist_utils.H"

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
   DirichletNeumann,
   RobinNeumann,
   monolithic,
   nocoupling
  };

  // type of norm to check for convergence
  enum ConvergenceNorm
  {
    absoluteconvergencenorm
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

}  // namespace FPSI

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif /* INPAR_FPSI_H_ */
