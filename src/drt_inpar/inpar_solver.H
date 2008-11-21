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
#ifdef CCADISCRET
#ifndef INPAR_SOLVER_H
#define INPAR_SOLVER_H

#ifdef __cplusplus /* stuff which is c++ specific */

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
      mlib_d_sp,        ///< solver package, hp's mlib, ect-sparse-symmetric-unsymmetric
      aztec_msr,        ///< solver package aztec 2.1, matrix in dmsr format 
      hypre_amg,        ///< solver package hypre, amg-solver, matrix in parcsr format 
      hypre_pcg,        ///< solver package hypre, cg, matrix in parcsr format 
      hypre_gmres,      ///< solver package hypre, gmres, matrix in parcsr format 
      hypre_bicgstab,   ///< solver package hypre, bicgstab, matrix in parcsr format 
      parsuperlu,       ///< solver package superlu, direkt parallel LU, matrix in redundant harwell-boeing 
      lapack_sym,       ///< symmetric lapack LU decomposition, matrix is dense 
      lapack_nonsym,    ///< unsymmetric lapack LU decomposition, matrix is dense 
      mumps_sym,        ///< solver package mumps, multifrontal parallel LU, matrixin row/column pointer format 
      mumps_nonsym,     ///< same but unsymmetric 
      colsol_solver,    ///< colsol 
      SPOOLES_sym,      ///< spooles parallel direct solver 
      SPOOLES_nonsym,   ///< spooles parallel direct solver 
      umfpack,          ///< solver package umfpack, matrix in compressed column format 
      mlpcg,            ///< internal solver package mlpcg 
      amesos_klu_sym,   ///< Trilinos' solver interface Amesos using KLU 
      amesos_klu_nonsym,///< Trilinos' solver interface Amesos using KLU 
      superlu,          ///< Trilinos' solver interface Amesos using SuperLU_Dist 
      vm3               ///< VM3 solver 
    };

    /// Different solvers within the Aztec2.1 library
    ///
    /// \author m.gee
    /// \date 9/01,11/08
    enum AzSolverType
    {
      azsolv_CG,            ///< cg-solver 
      azsolv_GMRES,         ///< gmres-solver 
      azsolv_CGS,           ///< cg-squared-solver (can handle unsymmetric problems) 
      azsolv_BiCGSTAB,      ///< bicgstab-solver 
      azsolv_LU,            ///< LU-solver (I think this doesn't work for some reasons) 
      azsolv_TFQMR          ///< quasi-minimum residual-solver (never used it) 
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
      azprec_ML                  ///< standard ML for structures 
    };

  }  // namespace SOLVER

}  // namespace INPAR


/*----------------------------------------------------------------------*/
#endif  // #ifdef __cplusplus

#endif  // #ifndef INPAR_SOLVER_H
#endif  // #ifdef CCADISCRET
