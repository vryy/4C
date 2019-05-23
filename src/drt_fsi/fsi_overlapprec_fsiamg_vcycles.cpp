/*----------------------------------------------------------------------*/
/*!

\brief V-Cycles for FSI AMG preconditioners

\level 1

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/

#include "fsi_overlapprec_fsiamg.H"
#include "../linalg/linalg_solver.H"

#include <Epetra_Time.h>
#include <ml_MultiLevelPreconditioner.h>
#include "MLAPI_LoadBalanceOperator.h"
#include "MLAPI_LoadBalanceInverseOperator.h"
#include "MLAPI_Operator_Utils.h"
#include "MLAPI_CompObject.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_Workspace.h"

#include "EpetraExt_SolverMap_CrsMatrix.h"

/*----------------------------------------------------------------------*
A Richardson iteration wrapper of a single field V-cycle
 *----------------------------------------------------------------------*/
double FSI::OverlappingBlockMatrixFSIAMG::RichardsonV(const std::string field, const int myrank,
    int sweeps, const double damp, std::vector<int>& levelsweeps, std::vector<double>& leveldamps,
    std::vector<MLAPI::Operator>& A, std::vector<Teuchos::RCP<MLAPI::InverseOperator>>& S,
    std::vector<MLAPI::Operator>& P, std::vector<MLAPI::Operator>& R, const int level,
    const int nlevel, MLAPI::MultiVector& x, const MLAPI::MultiVector& f, bool initiguesszero,
    bool analysis, bool silent) const
{
#if (FSIAMG_ANALYSIS >= 3)
  analysis = true;
  silent = false;
#endif

  MLAPI::MultiVector r(f.GetVectorSpace(), 1, false);
  if (initiguesszero)
    r = f;
  else
    r = f - A[level] * x;

  double initrl2 = 0.0;
  double initrinf = 0.0;
  if (analysis)
  {
    initrl2 = r.Norm2();
    initrinf = r.NormInf();
  }

  MLAPI::MultiVector tmpx(x.GetVectorSpace(), 1, false);

  Epetra_Time timer(MLAPI::GetEpetra_Comm());

  for (int i = 1; i <= sweeps; ++i)
  {
    tmpx = 0.0;
    Vcycle(field, myrank, levelsweeps, leveldamps, level, nlevel, tmpx, r, A, S, P, R);
    x.Update(damp, tmpx, 1.0);
    // x = x + damp * tmpx;
    if (i < sweeps || analysis) r = f - A[level] * x;
  }

  if (analysis)
  {
    double t = timer.ElapsedTime();
    double rl2 = r.Norm2();
    double rinf = r.NormInf();
    double rl2rate = Rate(myrank, t, rl2, initrl2, r.GetGlobalLength());
    double rinfrate = Rate(myrank, t, rinf, initrinf, r.GetGlobalLength());
    if (!myrank && !silent)
      printf(
          "RichardsonV %s      (level %2d) r0 %10.5e rl_2 %10.5e t %10.5e damp %10.5e sweeps %2d "
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
A single field V-cycle
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::Vcycle(const std::string field, const int myrank,
    std::vector<int>& sweeps, std::vector<double>& damps, const int level, const int nlevel,
    MLAPI::MultiVector& z, const MLAPI::MultiVector& b, const std::vector<MLAPI::Operator>& A,
    const std::vector<Teuchos::RCP<MLAPI::InverseOperator>>& S,
    const std::vector<MLAPI::Operator>& P, const std::vector<MLAPI::Operator>& R) const
{
  // coarse solve
  if (level == nlevel - 1)
  {
    z = 0.0;
    RichardsonS(field, myrank, level, sweeps[level], damps[level], A[level], *(S[level]), z, b,
        true, false, true);
    return;
  }

  // presmoothing (initial guess = 0)
  z = 0.0;
  RichardsonS(field, myrank, level, sweeps[level], damps[level], A[level], *(S[level]), z, b, true,
      false, true);

  // coarse level residual and correction
  MLAPI::MultiVector bc;
  MLAPI::MultiVector zc(P[level].GetDomainSpace(), 1, false);
  bc = R[level] * (b - A[level] * z);

  // solve coarse problem
  Vcycle(field, myrank, sweeps, damps, level + 1, nlevel, zc, bc, A, S, P, R);

  // prolongate correction
  z = z + P[level] * zc;

  // postsmoothing (initial guess != 0 !!)
  MLAPI::MultiVector r(b.GetVectorSpace(), 1, false);
  MLAPI::MultiVector dz(b.GetVectorSpace(), 1, true);
  r = A[level] * z;
  r.Update(1.0, b, -1.0);
  dz = 0.0;
  RichardsonS(field, myrank, level, sweeps[level], damps[level], A[level], *(S[level]), dz, r, true,
      false, true);
  z.Update(1.0, dz, 1.0);

  return;
}


/*----------------------------------------------------------------------*
A Richardson iteration of a FSI block Gauss Seidel
with V-cycles for individual fields
 *----------------------------------------------------------------------*/
double FSI::OverlappingBlockMatrixFSIAMG::RichardsonBGS_V(const int myrank, const int sweeps,
    const double damp, std::vector<int>& blocksweeps, std::vector<double>& blockdamps,
    AnalyzeBest& sbest, AnalyzeBest& fbest, AnalyzeBest& abest, MLAPI::MultiVector& sy,
    MLAPI::MultiVector& fy, MLAPI::MultiVector& ay, const MLAPI::MultiVector& sf,
    const MLAPI::MultiVector& ff, const MLAPI::MultiVector& af, std::vector<MLAPI::Operator>& Ass,
    std::vector<MLAPI::Operator>& Pss, std::vector<MLAPI::Operator>& Rss,
    std::vector<MLAPI::Operator>& Aff, std::vector<MLAPI::Operator>& Pff,
    std::vector<MLAPI::Operator>& Rff, std::vector<MLAPI::Operator>& Aaa,
    std::vector<MLAPI::Operator>& Paa, std::vector<MLAPI::Operator>& Raa,
    std::vector<MLAPI::Operator>& Asf, std::vector<MLAPI::Operator>& Afs,
    std::vector<MLAPI::Operator>& Afa, std::vector<MLAPI::Operator>& Aaf, bool initiguesszero,
    bool analysis, bool silent) const
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
  Epetra_Time timer(MLAPI::GetEpetra_Comm());
  double t1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;

  for (int i = 1; i <= sweeps; ++i)
  {
    //--------------------- structure block
    {
      if (analysis) timer.ResetStartTime();
      // compute ( f - A x ) for structure row
      stmpf = sf - Ass[0] * sy;
      stmpf = stmpf - Asf[0] * fy;
      // zero initial guess
      sz = 0.0;
      RichardsonV("(s)", myrank, blocksweeps[0], blockdamps[0], sbest.Sweeps(), sbest.Damp(), Ass,
          sbest.S(), Pss, Rss, 0, sbest.Nlevel(), sz, stmpf, true, false, true);
      sy.Update(damp, sz, 1.0);
      if (analysis) t1 += timer.ElapsedTime();
    }
    //---------------------- ale block
    {
      if (analysis) timer.ResetStartTime();
      // compute ( f - A x ) for ale row
      atmpf = af - Aaa[0] * ay;
      if (structuresplit_)
        atmpf = atmpf - Aaf[0] * fy;
      else
        atmpf = atmpf - Aaf[0] * sy;
      // zero initial guess
      az = 0.0;
      RichardsonV("(a)", myrank, blocksweeps[2], blockdamps[2], abest.Sweeps(), abest.Damp(), Aaa,
          abest.S(), Paa, Raa, 0, abest.Nlevel(), az, atmpf, true, false, true);
      ay.Update(damp, az, 1.0);
      if (analysis) t2 += timer.ElapsedTime();
    }
    //------------------------ fluid block
    {
      if (analysis) timer.ResetStartTime();
      // compute ( f - A x ) for fluid row
      ftmpf = ff - Aff[0] * fy;
      ftmpf = ftmpf - Afs[0] * sy;
      ftmpf = ftmpf - Afa[0] * ay;
      // zero initial guess
      fz = 0.0;
      RichardsonV("(f)", myrank, blocksweeps[1], blockdamps[1], fbest.Sweeps(), fbest.Damp(), Aff,
          fbest.S(), Pff, Rff, 0, fbest.Nlevel(), fz, ftmpf, true, false, true);
      fy.Update(damp, fz, 1.0);
      if (analysis) t3 += timer.ElapsedTime();
    }
  }  // iterations

  if (analysis)
  {
    double t = timer.ElapsedTime();

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
          "RichardsonBGS_V      (level %2d) r0 %10.5e rl_2 %10.5e rinf0 %10.5e rinf %10.5e t "
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
A Richardson iteration wrapper of a block FSI Vcycle
 *----------------------------------------------------------------------*/
double FSI::OverlappingBlockMatrixFSIAMG::Richardson_BlockV(const int myrank, const int sweeps,
    const double damp, std::vector<int>& Vsweeps, std::vector<double>& Vdamps,
    std::vector<int>* blocksweeps, std::vector<double>* blockdamps,
    const std::vector<std::string>& blocksmoother, AnalyzeBest& sbest, AnalyzeBest& fbest,
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

  // initial residuals
  if (initiguesszero)
  {
    stmpf = sf;
    atmpf = af;
    ftmpf = ff;
  }
  else
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
  }
  double initsrl2 = 0.0;
  double initsrinf = 0.0;
  double initarl2 = 0.0;
  double initarinf = 0.0;
  double initfrl2 = 0.0;
  double initfrinf = 0.0;
  double initrl2 = 0.0;
  double initrinf = 0.0;

  if (analysis)
  {
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

  Epetra_Time timer(MLAPI::GetEpetra_Comm());

  for (int i = 1; i <= sweeps; ++i)
  {
    sz = 0.0;
    az = 0.0;
    fz = 0.0;

    BlockVcycle(myrank, Vsweeps, Vdamps, blocksweeps, blockdamps, blocksmoother, 0, minnlevel_,
        sbest, fbest, abest, sz, fz, az, stmpf, ftmpf, atmpf, Ass, Pss, Rss, Aff, Pff, Rff, Aaa,
        Paa, Raa, Asf, Afs, Afa, Aaf);

    sy.Update(damp, sz, 1.0);
    ay.Update(damp, az, 1.0);
    fy.Update(damp, fz, 1.0);

    if (i < sweeps || analysis)
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
    }
  }  // iterations

  if (analysis)
  {
    double t = timer.ElapsedTime();

    double srl2 = stmpf.DotProduct(stmpf);
    double srinf = stmpf.NormInf();
    double arl2 = atmpf.DotProduct(atmpf);
    double arinf = atmpf.NormInf();
    double frl2 = ftmpf.DotProduct(ftmpf);
    double frinf = ftmpf.NormInf();
    double rl2 = sqrt(srl2 + arl2 + frl2);
    double rinf = std::max(srinf, arinf);
    rinf = std::max(rinf, frinf);

    double rl2rate = Rate(myrank, t, rl2, initrl2,
        stmpf.GetGlobalLength() + atmpf.GetGlobalLength() + ftmpf.GetGlobalLength());
    double rinfrate = Rate(myrank, t, rinf, initrinf,
        stmpf.GetGlobalLength() + atmpf.GetGlobalLength() + ftmpf.GetGlobalLength());

    if (!myrank && !silent)
      printf(
          "Richardson_BlockV    (level %2d) r0 %10.5e rl_2 %10.5e rinf0 %10.5e rinf %10.5e t "
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
A block FSI Vcycle ( which is AMG(BGS/Schur) )
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::BlockVcycle(const int myrank, std::vector<int>& Vsweeps,
    std::vector<double>& Vdamps, std::vector<int>* blocksweeps, std::vector<double>* blockdamps,
    const std::vector<std::string>& blocksmoother, const int level, const int minnlevel,
    AnalyzeBest& sbest, AnalyzeBest& fbest, AnalyzeBest& abest, MLAPI::MultiVector& sy,
    MLAPI::MultiVector& fy, MLAPI::MultiVector& ay, const MLAPI::MultiVector& sf,
    const MLAPI::MultiVector& ff, const MLAPI::MultiVector& af, std::vector<MLAPI::Operator>& Ass,
    std::vector<MLAPI::Operator>& Pss, std::vector<MLAPI::Operator>& Rss,
    std::vector<MLAPI::Operator>& Aff, std::vector<MLAPI::Operator>& Pff,
    std::vector<MLAPI::Operator>& Rff, std::vector<MLAPI::Operator>& Aaa,
    std::vector<MLAPI::Operator>& Paa, std::vector<MLAPI::Operator>& Raa,
    std::vector<MLAPI::Operator>& Asf, std::vector<MLAPI::Operator>& Afs,
    std::vector<MLAPI::Operator>& Afa, std::vector<MLAPI::Operator>& Aaf) const
{
  sy = 0.0;
  ay = 0.0;
  fy = 0.0;

  // coarsest level Richardson
  if (level == minnlevel - 1)
  {
    if (blocksmoother[level] == "Schur")
    {
      if (structuresplit_)
        RichardsonSchurStructureSplit_S(myrank, Vsweeps[level], Vdamps[level], blocksweeps,
            blockdamps, level, sbest, fbest, abest, sy, fy, ay, sf, ff, af, Ass, Pss, Rss, Aff, Pff,
            Rff, Aaa, Paa, Raa, Asf[level], Afs[level], Afa[level], Aaf[level], true, false, true);
      else
        RichardsonSchurFluidSplit_S(myrank, Vsweeps[level], Vdamps[level], blocksweeps, blockdamps,
            level, sbest, fbest, abest, sy, fy, ay, sf, ff, af, Ass, Pss, Rss, Aff, Pff, Rff, Aaa,
            Paa, Raa, Asf[level], Afs[level], Afa[level], Aaf[level], true, false, true);
    }
    else  // block Gauss Seidel
      RichardsonBGS_SV(myrank, Vsweeps[level], Vdamps[level], blocksweeps, blockdamps, level, sbest,
          fbest, abest, sy, fy, ay, sf, ff, af, Ass, Pss, Rss, Aff, Pff, Rff, Aaa, Paa, Raa,
          Asf[level], Afs[level], Afa[level], Aaf[level], true, false, true);
    return;
  }

  // presmoothing
  if (blocksmoother[level] == "Schur")
  {
    if (structuresplit_)
      RichardsonSchurStructureSplit_S(myrank, Vsweeps[level], Vdamps[level], blocksweeps,
          blockdamps, level, sbest, fbest, abest, sy, fy, ay, sf, ff, af, Ass, Pss, Rss, Aff, Pff,
          Rff, Aaa, Paa, Raa, Asf[level], Afs[level], Afa[level], Aaf[level], true, false, true);
    else
      RichardsonSchurFluidSplit_S(myrank, Vsweeps[level], Vdamps[level], blocksweeps, blockdamps,
          level, sbest, fbest, abest, sy, fy, ay, sf, ff, af, Ass, Pss, Rss, Aff, Pff, Rff, Aaa,
          Paa, Raa, Asf[level], Afs[level], Afa[level], Aaf[level], true, false, true);
  }
  else  // block Gauss Seidel
    RichardsonBGS_SV(myrank, Vsweeps[level], Vdamps[level], blocksweeps, blockdamps, level, sbest,
        fbest, abest, sy, fy, ay, sf, ff, af, Ass, Pss, Rss, Aff, Pff, Rff, Aaa, Paa, Raa,
        Asf[level], Afs[level], Afa[level], Aaf[level], true, false, true);

  // create coarse level residual and correction vectors
  // structure:
  // sfc = Rss[level] * ( sf - Ass_[level] * sy - Asf[level] * fy)
  MLAPI::MultiVector sfc;
  {
    MLAPI::MultiVector tmp;
    tmp = sf - Ass[level] * sy;
    tmp = tmp - Asf[level] * fy;
    sfc.Reshape(Rss[level].GetOperatorRangeSpace());
    sfc = Rss_[level] * tmp;
  }
  // ale:
  // afc = Raa[level] * ( af - Aaa[level] * ay - Aaf[level] * fy)
  MLAPI::MultiVector afc;
  {
    MLAPI::MultiVector tmp;
    tmp = af - Aaa[level] * ay;
    if (structuresplit_)
      tmp = tmp - Aaf[level] * fy;
    else
      tmp = tmp - Aaf[level] * sy;
    afc.Reshape(Raa[level].GetOperatorRangeSpace());
    afc = Raa[level] * tmp;
  }
  // fluid:
  // ffc = Rff[level] * ( ff - Aff[level] * fy - Afs[level] * sy - Afa[level] * ay)
  MLAPI::MultiVector ffc;
  {
    MLAPI::MultiVector tmp;
    tmp = ff - Aff[level] * fy;
    tmp = tmp - Afs[level] * sy;
    tmp = tmp - Afa[level] * ay;
    ffc.Reshape(Rff[level].GetOperatorRangeSpace());
    ffc = Rff[level] * tmp;
  }

  // coarse level correction vectors (get set to zero inside V-cycle on next level)
  MLAPI::MultiVector syc(sfc.GetVectorSpace(), 1, false);
  MLAPI::MultiVector ayc(afc.GetVectorSpace(), 1, false);
  MLAPI::MultiVector fyc(ffc.GetVectorSpace(), 1, false);

  // solve coarser problem
  BlockVcycle(myrank, Vsweeps, Vdamps, blocksweeps, blockdamps, blocksmoother, level + 1, minnlevel,
      sbest, fbest, abest, syc, fyc, ayc, sfc, ffc, afc, Ass, Pss, Rss, Aff, Pff, Rff, Aaa, Paa,
      Raa, Asf, Afs, Afa, Aaf);

  // prolongate correction
  {
    MLAPI::MultiVector tmp;
    tmp = Pss[level] * syc;
    sy.Update(1.0, tmp, 1.0);
    tmp = Paa[level] * ayc;
    ay.Update(1.0, tmp, 1.0);
    tmp = Pff[level] * fyc;
    fy.Update(1.0, tmp, 1.0);
  }

  // postsmoothing (do NOT zero out initial guess!)
  if (blocksmoother[level] == "Schur")
  {
    if (structuresplit_)
      RichardsonSchurStructureSplit_S(myrank, Vsweeps[level], Vdamps[level], blocksweeps,
          blockdamps, level, sbest, fbest, abest, sy, fy, ay, sf, ff, af, Ass, Pss, Rss, Aff, Pff,
          Rff, Aaa, Paa, Raa, Asf[level], Afs[level], Afa[level], Aaf[level], false, false, true);
    else
      RichardsonSchurFluidSplit_S(myrank, Vsweeps[level], Vdamps[level], blocksweeps, blockdamps,
          level, sbest, fbest, abest, sy, fy, ay, sf, ff, af, Ass, Pss, Rss, Aff, Pff, Rff, Aaa,
          Paa, Raa, Asf[level], Afs[level], Afa[level], Aaf[level], false, false, true);
  }
  else  // block Gauss Seidel
    RichardsonBGS_SV(myrank, Vsweeps[level], Vdamps[level], blocksweeps, blockdamps, level, sbest,
        fbest, abest, sy, fy, ay, sf, ff, af, Ass, Pss, Rss, Aff, Pff, Rff, Aaa, Paa, Raa,
        Asf[level], Afs[level], Afa[level], Aaf[level], false, false, true);

  return;
}
