/*----------------------------------------------------------------------*/
/*!
\file inpar_fsi.H

\level 1

\maintainer Matthias Mayr

\brief Input parameters for fluid structure interaction
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_FSI_H
#define INPAR_FSI_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT { namespace INPUT { class ConditionDefinition; } }

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
  fsi_iter_fluidfluid_monolithicstructuresplit_nonox,
  fsi_iter_fluidfluid_monolithicfluidsplit_nonox,
  fsi_iter_sliding_monolithicfluidsplit,
  fsi_iter_sliding_monolithicstructuresplit
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
      DirichletNeumann,
      DirichletNeumannVolCoupl
    };

    /// Linear preconditioning algorithm for monolithic block system
    enum LinearBlockSolver
    {
      PreconditionedKrylov, ///< BGS(AMG)
      FSIAMG, ///< AMG(BGS)
      AMGnxn, ///< general n-by-n block preconditioning framework
      HybridSchwarz ///< hybrid additive/multiplicative Schwarz
    };

    enum HybridASType
    {
      hybrid_as_type_ILU, ///< do processor-local ILU
      hybrid_as_type_Amesos_LU ///< do processor-local direct solve
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

    /// Handling of non-converged nonlinear solver
    enum DivContAct
    {
      divcont_stop,               ///< abort simulation
      divcont_continue,           ///< continue nevertheless
      divcont_halve_step,         ///< halve time step and carry on with simulation
      divcont_revert_dt           ///< revert time step size to previous one
    };

    /// Verbosity of FSI algorithm
    enum Verbosity
    {
      verbosity_subproblem, ///< output for FSI as a subproblem
      verbosity_low, ///< write only nonlinear solver status
      verbosity_medium, ///< write nonlinear and linear solver status
      verbosity_full, ///< write everything
    };

    /// options for parallel domain redistribution
    enum Redistribute
    {
      Redistribute_off, ///< no redistribution
      Redistribute_structure, ///< redistribute structure field only
      Redistribute_fluid, ///< redistribute fluid field only
      Redistribute_both, ///< redistribute both structure and fluid field
      Redistribute_monolithic ///< redistribute monolithic graph (fluid + structure)
    };

    //@}

    /// set the fsi parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific fsi conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

  }  // namespace FSI

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_FSI_H
