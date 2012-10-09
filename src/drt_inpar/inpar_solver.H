/*----------------------------------------------------------------------*/
/*!
\file inpar_solver.H

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_SOLVER_H
#define INPAR_SOLVER_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for solver
  namespace SOLVER
  {

    /// Available solvers (a few of these are likely out-dated)
    ///
    /// \author m.gee
    /// \date 7/01,11/08
    enum SolverType
    {
      aztec_msr,         ///< solver package aztec 2.1, matrix in dmsr format
      lapack_sym,        ///< symmetric lapack LU decomposition, matrix is dense
      lapack_nonsym,     ///< unsymmetric lapack LU decomposition, matrix is dense
      umfpack,           ///< solver package umfpack, matrix in compressed column format
      amesos_klu_sym,    ///< Trilinos' solver interface Amesos using KLU
      amesos_klu_nonsym, ///< Trilinos' solver interface Amesos using KLU
      superlu,           ///< Trilinos' solver interface Amesos using SuperLU_Dist
      vm3,               ///< VM3 solver
      belos,             ///< Belos solver
      stratimikos_amesos,///< Stratimikos solver framework (Amesos KLU)
      stratimikos_aztec, ///< Stratimikos solver framework (Aztec)
      stratimikos_belos, ///< Stratimikos solver framework (Belos)
      undefined          ///< undefined solver
    };

    /// Different solvers within the Aztec2.1 library
    ///
    /// \author m.gee
    /// \date 9/01,11/08
    enum AzSolverType
    {
      azsolv_CG,            ///< cg-solver
      azsolv_GMRES,         ///< gmres-solver
      azsolv_GMRESR,         ///< restarted gmres-solver
      azsolv_CGS,           ///< cg-squared-solver (can handle unsymmetric problems)
      azsolv_BiCGSTAB,      ///< bicgstab-solver
      azsolv_LU,            ///< LU-solver (I think this doesn't work for some reasons)
      azsolv_TFQMR,         ///< quasi-minimum residual-solver (never used it)
      belos_FGMRES          ///< flexible GMRES in Belos
    };

    /// Different preconditioners within the Aztec package
    ///
    /// \author m.gee
    /// \date 9/01,11/08
    enum AzPrecType
    {
      azprec_none,               ///< no preconditioning
      azprec_ILUT,               ///< incomplete LU fact. with numerical drop tolerance
      azprec_ILU,                ///< incomplete LU fact. with fill in levels
      azprec_Jacobi,             ///< Jacobi preconditioning
      azprec_Neumann,            ///< neumann polynomials
      azprec_Least_Squares,      ///< least squares something
      azprec_SymmGaussSeidel,    ///< symmetric n-step gauss-seidel preconditioner
      azprec_GaussSeidel,        ///< symmetric n-step gauss-seidel preconditioner
      azprec_DownwindGaussSeidel,///< symmetric n-step gauss-seidel preconditioner
      azprec_LU,                 ///< ?
      azprec_RILU,               ///< relaxed incomplete LU
      azprec_BILU,               ///< block incomplete LU (only with matrix in DVBR format
      azprec_ICC,                ///< incomplete cholesky
      azprec_MLfluid,            ///< ML for Fluids
      azprec_MLfluid2,           ///< energy optimal unsymmetric ML for Fluids
      azprec_MLAPI,              ///< experimental mlapi preconditioner
      azprec_ML,                 ///< standard ML for structures
      azprec_AMGBS,              ///< AMG(Braess-Sarazin) preconditioner for Fluids
      azprec_AMG,                ///< AMG preconditioner for Fluids (without ML)
      azprec_BGS2x2,             ///< Block Gauss-Seidel for 2x2 system
      azprec_TekoSIMPLE,         ///< SIMPLE preconditioner for 2x2 block system (Teko)
      azprec_BGSnxn,             ///< Block Gauss-Seidel for nxn system (Teko)
      azprec_CheapSIMPLE,        ///< CheapSIMPLE (BACI implementation for contact/meshtying/constraints)
      azprec_MueLuAMG_sym,	     ///< MueLu AMG preconditioner (Smoothed aggregation for symmetric problems)
      azprec_MueLuAMG_nonsym,	   ///< MueLu AMG preconditioner (PG-AMG for nonsymmetric problems ~ MLFLUID2)
      azprec_MueLuAMG_contact,   ///< MueLu AMG preconditioner (special handling for condensed contact formulation)
      azprec_MueLuAMG_contact2,  ///< MueLu AMG preconditioner (special enhanced handling for condensed contact formulation)
      azprec_MueLuAMG_contactSP  ///< MueLu AMG preconditioner for saddlepoint formulation
    };

    /// Different prolongation smoothing strategies for AMGBS preconditioner
    ///
    /// \author tw
    /// \date 03/10
    enum ProlongatorSmoothingStrategy
    {
        PA_AMG,         ///< plain aggregation: no smoothing
        SA_AMG,         ///< SA-AMG (smoothed aggregation)
        PG_AMG,         ///< PG-AMG (Petrov Galerkin energy minimization)
        PG2_AMG         ///< PG-AMG (Petrov Galerkin energy minimization), same as MATLAB variant
    };
  }  // namespace SOLVER

}  // namespace INPAR


/*----------------------------------------------------------------------*/

#endif  // #ifndef INPAR_SOLVER_H
