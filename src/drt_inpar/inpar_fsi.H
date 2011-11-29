/*----------------------------------------------------------------------*/
/*!
\file inpar_fsi.H

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef INPAR_FSI_H
#define INPAR_FSI_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for FSI
  namespace FSI
  {

    /// Type of partitioned coupling for FSI problems
    enum PartitionedCouplingMethod
    {
      DirichletNeumannSlideale,
      DirichletNeumann,
      RobinNeumann,
      DirichletRobin,
      RobinRobin
    };

    /// Linear solver algorithm for monolithic block system in monolithic FSI
    enum LinearBlockSolver
    {
      PreconditionedKrylov,
      FSIAMG,
      PartitionedAitken,
      PartitionedVectorExtrapolation,
      PartitionedJacobianFreeNewtonKrylov,
      BGSAitken,
      BGSVectorExtrapolation,
      BGSJacobianFreeNewtonKrylov
    };

    /// Projection methods for sliding fluid-structure interface
    enum SlideALEProj
    {
      ALEprojection_none,
      ALEprojection_curr,
      ALEprojection_ref,
      ALEprojection_rot_z,
      ALEprojection_rot_zsphere
    };

    /// Preconditioner for constraint fsi problem
    enum PrecConstr
    {
      Simple,
      Simplec
    };

    /// @name Solution technique and related

    /// type of norm to check for convergence of newton loop
    enum ConvNorm
    {
      convnorm_abs,  ///< absolute norm
      convnorm_rel,  ///< relative norm
      convnorm_mix   ///< mixed absolute-relative norm
    };

    /// type of norm to check for convergence
    enum BinaryOp
    {
      bop_and  ///<  and
    };

    //@}

  }  // namespace FSI

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_FSI_H
#endif  // #ifdef CCADISCRET
