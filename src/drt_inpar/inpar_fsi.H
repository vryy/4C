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
#ifndef INPAR_FSI_H
#define INPAR_FSI_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* The coupling methods for FSI. */
/*----------------------------------------------------------------------*/
// ToDo: put into the namespace INPAR::FSI ! No typedef?
typedef enum _FSI_COUPLING
{
  fsi_coupling_freesurface=-1,
  fsi_coupling_undefined=0,
  fsi_basic_sequ_stagg=1,
  fsi_iter_stagg_fixed_rel_param=4,
  fsi_iter_stagg_AITKEN_rel_param=5,
  fsi_iter_stagg_steep_desc=6,
  fsi_iter_stagg_CHEB_rel_param=7,
  fsi_iter_stagg_AITKEN_rel_force=8,
  fsi_iter_stagg_steep_desc_force=9,
  fsi_iter_stagg_Newton_FD=10,
  fsi_iter_stagg_Newton_I=11,
  fsi_iter_monolithicfluidsplit=13,
  fsi_iter_monolithicstructuresplit,
  fsi_iter_lung_monolithicstructuresplit,
  fsi_iter_lung_monolithicfluidsplit,
  fsi_iter_mortar_monolithicstructuresplit,
  fsi_iter_mortar_monolithicfluidsplit,
  fsi_iter_constr_monolithicstructuresplit,
  fsi_iter_constr_monolithicfluidsplit,
  fsi_iter_xfem_monolithic,
  fsi_iter_stagg_NLCG,
  fsi_iter_stagg_MFNK_FD,
  fsi_iter_stagg_MFNK_FSI,
  fsi_iter_stagg_MPE,
  fsi_iter_stagg_RRE,
  fsi_pseudo_structureale,
  fsi_iter_fluidfluid_monolithicstructuresplit,
  fsi_iter_fluidfluid_monolithicfluidsplit,
  fsi_iter_fluidfluid_monolithicstructuresplit_nox
} FSI_COUPLING;

namespace INPAR
{

  /// Input parameters for FSI
  namespace FSI
  {

    /// Type of partitioned coupling for FSI problems
    enum PartitionedCouplingMethod
    {
      DirichletNeumannSlideale,
      DirichletNeumann
    };

    /// Linear solver algorithm for monolithic block system in monolithic FSI
    enum LinearBlockSolver
    {
      PreconditionedKrylov,
      FSIAMG
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

    /// type of fluid field auxiliary time integrator for time adaptivity
    enum FluidMethod
    {
      timada_fld_none,            ///< no time adaptivity based on fluid field
      timada_fld_expleuler,       ///< use explicit Euler as auxiliary time integrator
      timada_fld_adamsbashforth2  ///< use Adams-Bashforth-2 as auxiliary time integrator
    };

    /// type of structure field auxiliary time integrator for time adaptivity
    enum SructureMethod
    {
      timada_str_none,            ///< no time adaptivity based on structure field
      timada_str_expleuler,       ///< use explicit Euler as auxiliary time integrator
      timada_str_adamsbashforth2  ///< use Adams-Bashforth-2 as auxiliary time integrator
    };

    //@}

  }  // namespace FSI

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_FSI_H
