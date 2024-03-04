/*----------------------------------------------------------------------*/
/*! \file

\brief Smoothers for FSI AMG preconditioners

\level 1

*/
/*----------------------------------------------------------------------*/

#include "baci_fsi_overlapprec_fsiamg.hpp"
#include "baci_linear_solver_method_linalg.hpp"
#include "baci_linear_solver_preconditioner_linalg.hpp"

#include <EpetraExt_SolverMap_CrsMatrix.h>
#include <MLAPI_Expressions.h>
#include <MLAPI_LoadBalanceInverseOperator.h>
#include <MLAPI_LoadBalanceOperator.h>
#include <MLAPI_MultiVector.h>
#include <MLAPI_Workspace.h>
#include <Teuchos_Time.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double FSI::OverlappingBlockMatrixFSIAMG::RichardsonBGS_Mixed(const int myrank, const int sweeps,
    const double damp, std::vector<int>& blocksweeps, std::vector<double>& blockdamps,
    const bool sisamg, const bool fisamg, const bool aisamg, AnalyzeBest& sbest, AnalyzeBest& fbest,
    AnalyzeBest& abest, MLAPI::MultiVector& sy, MLAPI::MultiVector& fy, MLAPI::MultiVector& ay,
    const MLAPI::MultiVector& sf, const MLAPI::MultiVector& ff, const MLAPI::MultiVector& af,
    std::vector<MLAPI::Operator>& Ass, std::vector<MLAPI::Operator>& Pss,
    std::vector<MLAPI::Operator>& Rss, std::vector<MLAPI::Operator>& Aff,
    std::vector<MLAPI::Operator>& Pff, std::vector<MLAPI::Operator>& Rff,
    std::vector<MLAPI::Operator>& Aaa, std::vector<MLAPI::Operator>& Paa,
    std::vector<MLAPI::Operator>& Raa, std::vector<MLAPI::Operator>& Asf,
    std::vector<MLAPI::Operator>& Afs, std::vector<MLAPI::Operator>& Afa,
    std::vector<MLAPI::Operator>& Aaf, bool initiguesszero, bool analysis, bool silent) const
{
#if (FSIAMG_ANALYSIS >= 1)
  analysis = true;
  silent = false;
#endif

  MLAPI::MultiVector stmpf;
  MLAPI::MultiVector ftmpf;
  MLAPI::MultiVector atmpf;
  MLAPI::MultiVector sz(sy.GetVectorSpace(), 1, false);
  MLAPI::MultiVector fz(fy.GetVectorSpace(), 1, false);
  MLAPI::MultiVector az(ay.GetVectorSpace(), 1, false);

  double initsrl2 = 0.0;
  double initsrinf = 0.0;
  double initarl2 = 0.0;
  double initarinf = 0.0;
  double initfrl2 = 0.0;
  double initfrinf = 0.0;
  double initrl2 = 0.0;
  double initrinf = 0.0;
  double rl2 = 0.0;
  double rinf = 0.0;

  if (analysis)
  {
    stmpf = sf - Ass[0] * sy;
    stmpf = stmpf - Asf[0] * fy;
    atmpf = af - Aaa[0] * ay;
    if (structuresplit_)
      atmpf = atmpf - Aaf[0] * fy;
    else
      atmpf = atmpf - Aaf[0] * sy;
    ftmpf = ff - Aff[0] * fy;
    ftmpf = ftmpf - Afs[0] * sy;
    ftmpf = ftmpf - Afa[0] * ay;
    initsrl2 = stmpf.DotProduct(stmpf);
    initsrinf = stmpf.NormInf();
    initarl2 = atmpf.DotProduct(atmpf);
    initarinf = atmpf.NormInf();
    initfrl2 = ftmpf.DotProduct(ftmpf);
    initfrinf = ftmpf.NormInf();
    initrl2 = sqrt(initsrl2 + initarl2 + initfrl2);
    initrinf = std::max(initsrinf, initarinf);
    initrinf = std::max(initrinf, initfrinf);
  }
  Teuchos::Time timer("", true);
  double t1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;

  for (int i = 1; i <= sweeps; ++i)
  {
    //--------------------- structure block
    {
      if (analysis) timer.reset();
      // compute ( f - A x ) for structure row
      stmpf = sf - Ass[0] * sy;
      stmpf = stmpf - Asf[0] * fy;
      // zero initial guess
      sz = 0.0;

      if (sisamg)
        RichardsonV("(s)", myrank, blocksweeps[0], blockdamps[0], sbest.Sweeps(), sbest.Damp(), Ass,
            sbest.S(), Pss, Rss, 0, sbest.Nlevel(), sz, stmpf, true, false, true);
      else
        RichardsonMixed("(s)", myrank, 0, blocksweeps[0], blockdamps[0], Ass[0], Matrix(0, 0),
            structuresolver_, sz, stmpf, const_cast<int&>(srun_), true, false, true);

      sy.Update(damp, sz, 1.0);
      if (analysis) t1 += timer.totalElapsedTime(true);
    }
    //---------------------- ale block
    {
      if (analysis) timer.reset();
      // compute ( f - A x ) for ale row
      atmpf = af - Aaa[0] * ay;
      if (structuresplit_)
        atmpf = atmpf - Aaf[0] * fy;
      else
        atmpf = atmpf - Aaf[0] * sy;
      // zero initial guess
      az = 0.0;

      if (aisamg)
        RichardsonV("(a)", myrank, blocksweeps[2], blockdamps[2], abest.Sweeps(), abest.Damp(), Aaa,
            abest.S(), Paa, Raa, 0, abest.Nlevel(), az, atmpf, true, false, true);
      else
        RichardsonMixed("(a)", myrank, 0, blocksweeps[2], blockdamps[2], Aaa[0], Matrix(2, 2),
            alesolver_, az, atmpf, const_cast<int&>(arun_), true, false, true);

      ay.Update(damp, az, 1.0);
      if (analysis) t2 += timer.totalElapsedTime(true);
    }
    //------------------------ fluid block
    {
      if (analysis) timer.reset();
      // compute ( f - A x ) for fluid row
      ftmpf = ff - Aff[0] * fy;
      ftmpf = ftmpf - Afs[0] * sy;
      ftmpf = ftmpf - Afa[0] * ay;
      // zero initial guess
      fz = 0.0;

      if (fisamg)
        RichardsonV("(f)", myrank, blocksweeps[1], blockdamps[1], fbest.Sweeps(), fbest.Damp(), Aff,
            fbest.S(), Pff, Rff, 0, fbest.Nlevel(), fz, ftmpf, true, false, true);
      else
        RichardsonMixed("(f)", myrank, 0, blocksweeps[1], blockdamps[1], Aff[0], Matrix(1, 1),
            fluidsolver_, fz, ftmpf, const_cast<int&>(frun_), true, false, true);

      fy.Update(damp, fz, 1.0);
      if (analysis) t3 += timer.totalElapsedTime(true);
    }
  }  // iterations

  if (analysis)
  {
    double t = timer.totalElapsedTime(true);

    // final residuals
    stmpf = sf - Ass[0] * sy;
    stmpf = stmpf - Asf[0] * fy;
    double srl2 = stmpf.DotProduct(stmpf);
    double srinf = stmpf.NormInf();

    atmpf = af - Aaa[0] * ay;
    if (structuresplit_)
      atmpf = atmpf - Aaf[0] * fy;
    else
      atmpf = atmpf - Aaf[0] * sy;
    double arl2 = atmpf.DotProduct(atmpf);
    double arinf = atmpf.NormInf();

    ftmpf = ff - Aff[0] * fy;
    ftmpf = ftmpf - Afs[0] * sy;
    ftmpf = ftmpf - Afa[0] * ay;
    double frl2 = ftmpf.DotProduct(ftmpf);
    double frinf = ftmpf.NormInf();

    rl2 = sqrt(srl2 + arl2 + frl2);
    rinf = std::max(srinf, arinf);
    rinf = std::max(rinf, frinf);

    double rl2rate = Rate(myrank, t, rl2, initrl2,
        stmpf.GetGlobalLength() + atmpf.GetGlobalLength() + ftmpf.GetGlobalLength());
    double rinfrate = Rate(myrank, t, rinf, initrinf,
        stmpf.GetGlobalLength() + atmpf.GetGlobalLength() + ftmpf.GetGlobalLength());

    if (!myrank && !silent)
      printf(
          "RichardsonBGS_Mixed  (level %2d) r0 %10.5e rl_2 %10.5e rinf0 %10.5e rinf %10.5e t "
          "%10.5e damp %10.5e sweeps %2d l2rate %10.5e rinfrate %10.5e\n",
          0, initrl2, rl2, initrinf, rinf, t, damp, sweeps, rl2rate, rinfrate);
    rl2rate = std::max(rl2rate, rinfrate);

    MLAPI::GetEpetra_Comm().Broadcast(&rl2rate, 1, 0);
    return rl2rate;
  }
  else
    return 0.0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double FSI::OverlappingBlockMatrixFSIAMG::RichardsonMixed(const std::string field, const int myrank,
    const int level, const int sweeps, const double damp, const MLAPI::Operator& A,
    const CORE::LINALG::SparseMatrix& matrix,
    const Teuchos::RCP<CORE::LINALG::Preconditioner>& solver, MLAPI::MultiVector& x,
    const MLAPI::MultiVector& f, int& run, bool initiguesszero, bool analysis, bool silent) const
{
#if (FSIAMG_ANALYSIS >= 3)
  analysis = true;
  silent = false;
#endif

  MLAPI::MultiVector r(A.GetRangeSpace(), 1, true);
  MLAPI::MultiVector tmpx(A.GetDomainSpace(), 1, true);
  // MLAPI::MultiVector r(f.GetVectorSpace(),1,true);
  // MLAPI::MultiVector tmpx(x.GetVectorSpace(),1,true);

  if (initiguesszero)
    r = f;
  else
    r = f - A * x;

  // a view to these vectors (MUST be AFTER the above statement!)
  Teuchos::RCP<Epetra_Vector> er =
      Teuchos::rcp(new Epetra_Vector(View, matrix.RangeMap(), r.GetValues(0)));
  Teuchos::RCP<Epetra_Vector> etmpx =
      Teuchos::rcp(new Epetra_Vector(View, matrix.DomainMap(), tmpx.GetValues(0)));

  double initrinf = 0.0;
  double initrl2 = 0.0;
  double rinf = 0.0;
  double rl2 = 0.0;
  if (analysis)
  {
    initrl2 = r.Norm2();
    initrinf = r.NormInf();
  }

  Teuchos::Time timer("", true);

  for (int i = 1; i <= sweeps; ++i)
  {
    bool refactor = false;
    if (!run) refactor = true;
    tmpx = 0.0;
    solver->Solve(matrix.EpetraMatrix(), etmpx, er, refactor, false);
    x = x + damp * tmpx;
    if (i < sweeps || analysis)
    {
      r = f - A * x;
      // r is recreated here, so we need a fresh view (took me a while to find this one :-( )
      er = Teuchos::rcp(new Epetra_Vector(View, matrix.RangeMap(), r.GetValues(0)));
    }
    run++;
  }

  if (analysis)
  {
    rl2 = r.Norm2();
    rinf = r.NormInf();
    double t = timer.totalElapsedTime(true);
    double rl2rate = Rate(myrank, t, rl2, initrl2, r.GetGlobalLength());
    double rinfrate = Rate(myrank, t, rinf, initrinf, r.GetGlobalLength());

    if (!myrank && !silent)
      printf(
          "RichardsonMixed %s  (level %2d) r0 %10.5e rl_2 %10.5e t %10.5e damp %10.5e sweeps %2d "
          "l2rate %10.5e\n",
          field.c_str(), level, initrl2, rl2, t, damp, sweeps, rl2rate);
    rl2rate = std::max(rl2rate, rinfrate);

    MLAPI::GetEpetra_Comm().Broadcast(&rl2rate, 1, 0);
    return rl2rate;
  }
  else
    return 0.0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double FSI::OverlappingBlockMatrixFSIAMG::RichardsonS(const std::string field, const int myrank,
    const int level, const int sweeps, const double damp, const MLAPI::Operator& A,
    const MLAPI::InverseOperator& S, MLAPI::MultiVector& x, const MLAPI::MultiVector& f,
    bool initiguesszero, bool analysis, bool silent) const
{
#if (FSIAMG_ANALYSIS >= 3)
  analysis = true;
  silent = false;
#endif

  MLAPI::MultiVector r(f.GetVectorSpace(), 1, false);
  if (initiguesszero)
    r = f;
  else
    r = f - A * x;

  double initrinf = 0.0;
  double initrl2 = 0.0;
  double rinf = 0.0;
  double rl2 = 0.0;
  if (analysis)
  {
    initrl2 = r.Norm2();
    initrinf = r.NormInf();
  }

  MLAPI::MultiVector tmpx(x.GetVectorSpace(), 1, false);

  Teuchos::Time timer("", true);

  for (int i = 1; i <= sweeps; ++i)
  {
    tmpx = 0.0;
    S.Apply(r, tmpx);
    // x = x + damp * tmpx;
    x.Update(damp, tmpx, 1.0);
    if (i < sweeps || analysis) r = f - A * x;
  }

  if (analysis)
  {
    rl2 = r.Norm2();
    rinf = r.NormInf();
    double t = timer.totalElapsedTime(true);
    double rl2rate = Rate(myrank, t, rl2, initrl2, r.GetGlobalLength());
    double rinfrate = Rate(myrank, t, rinf, initrinf, r.GetGlobalLength());
    if (!myrank && !silent)
      printf(
          "RichardsonS %s      (level %2d) r0 %10.5e rl_2 %10.5e t %10.5e damp %10.5e sweeps %2d "
          "l2rate %10.5e\n",
          field.c_str(), level, initrl2, rl2, t, damp, sweeps, rl2rate);
    rl2rate = std::max(rl2rate, rinfrate);

    MLAPI::GetEpetra_Comm().Broadcast(&rl2rate, 1, 0);
    return rl2rate;
  }
  else
    return 0.0;
}

BACI_NAMESPACE_CLOSE
