/*----------------------------------------------------------------------*/
/*!

\brief Analysis tool for FSI preconditioners

\level 3

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/

#include "fsi_overlapprec_fsiamg.H"
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
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::AnalyzeFSIAMG(const int myrank, const int snlevel,
    const int fnlevel, const int anlevel, Teuchos::ParameterList& sparams,
    Teuchos::ParameterList& fparams, Teuchos::ParameterList& aparams,
    std::vector<MLAPI::Operator>& Ass, std::vector<MLAPI::Operator>& Aff,
    std::vector<MLAPI::Operator>& Aaa, std::vector<MLAPI::Operator>& Pss,
    std::vector<MLAPI::Operator>& Rss, std::vector<MLAPI::Operator>& Pff,
    std::vector<MLAPI::Operator>& Rff, std::vector<MLAPI::Operator>& Paa,
    std::vector<MLAPI::Operator>& Raa, std::vector<MLAPI::Operator>& Asf,
    std::vector<MLAPI::Operator>& Afs, std::vector<MLAPI::Operator>& Afa,
    std::vector<MLAPI::Operator>& Aaf, ML* sml, ML* fml, ML* aml)
{
  if (!myrank)
  {
    printf("VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV\n");
    printf("running FSIAMG preconditioner analysis\n");
    printf("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
    fflush(stdout);
  }

  //-----------------------------------------analysis of structure field
  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("STRUCTURE field analysis\n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  AnalyzeBest sbest(snlevel);
  Analyse_SingleField("(s)", myrank, "STRUCTURE", snlevel, sparams, Ass, Pss, Rss, sbest);

  //-----------------------------------------analysis of fluid field
  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("FLUID field analysis\n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  AnalyzeBest fbest(fnlevel);
  Analyse_SingleField("(f)", myrank, "FLUID", fnlevel, fparams, Aff, Pff, Rff, fbest);

  //-----------------------------------------analysis of ale field
  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("ALE field analysis\n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  AnalyzeBest abest(anlevel);
  Analyse_SingleField("(a)", myrank, "ALE", anlevel, aparams, Aaa, Paa, Raa, abest);


  //---------------------------------------- run analysis of BGS(AMG)
  // which is called PreconditionedKrylov in our input file
  Analyse_BGSAMG(
      myrank, sbest, fbest, abest, Ass, Pss, Rss, Aff, Pff, Rff, Aaa, Paa, Raa, Asf, Afs, Afa, Aaf);

  //---------------------------------------- run analysis of AMG(BGS)
  // which is called FSIAMG in our input file
  Analyse_AMGBGS(
      myrank, sbest, fbest, abest, Ass, Pss, Rss, Aff, Pff, Rff, Aaa, Paa, Raa, Asf, Afs, Afa, Aaf);

  exit(EXIT_SUCCESS);


  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::Analyse_SingleField(const std::string fieldname,
    const int myrank, std::string field, const int nlevel, Teuchos::ParameterList& params,
    std::vector<MLAPI::Operator>& A, std::vector<MLAPI::Operator>& P,
    std::vector<MLAPI::Operator>& R, AnalyzeBest& best)
{
  Teuchos::RCP<MLAPI::InverseOperator> S;
  // Teuchos::RCP<MLAPI::LoadBalanceInverseOperator> lbS;

  // measure r_l2 decrease per time in bestrate
  std::vector<double> bestrate(nlevel - 1, 10.0);

  // loop levels in single field and test smoothers
  for (int level = 0; level < nlevel - 1; ++level)
  {
    if (!myrank)
    {
      printf("--------------------------------------------\n");
      printf("level %d\n", level);
      printf("--------------------------------------------\n");
      fflush(stdout);
    }
    Teuchos::ParameterList pushlist(params.sublist("smoother: ifpack list"));
    char levelstr[19];
    sprintf(levelstr, "(level %d)", level);
    Teuchos::ParameterList subp = params.sublist("smoother: list " + (std::string)levelstr);
    // std::cout << "pushlist\n" << pushlist << "subp\n" << subp << std::endl;

    MLAPI::Space rspace(A[level].GetOperatorRangeSpace());
    MLAPI::Space dspace(A[level].GetOperatorDomainSpace());
    MLAPI::MultiVector xref(dspace);
    xref.Random();
    xref.Scale(10000.0);
    MLAPI::MultiVector f(rspace, true);

    double localrate = 10.0;
    std::string localtype = "";
    double localdamp = 1.0;
    int localpoly = 1;
    Teuchos::RCP<MLAPI::InverseOperator> localS;

    //----------------------------------------------------------- test SGS
    if (!myrank)
    {
      printf("--------------------------------------------\n");
      printf("symmetric Gauss-Seidel\n");
      printf("--------------------------------------------\n");
      fflush(stdout);
    }
    localrate = 10.0;
    localtype = "";
    localdamp = 1.0;
    localpoly = 1;
    for (int i = 0; i < 99; ++i)
    {
      double damp = 1.0 - i * 0.01;  // reduce damping factor in steps of 0.01;
      std::string type = "";
      subp.set("smoother: type", "symmetric Gauss-Seidel");
      subp.set("smoother: sweeps", 1);
      subp.set("smoother: damping factor", damp);
      Teuchos::ParameterList p;
      SelectMLAPISmoother(type, level, subp, p, pushlist);
      S = Teuchos::rcp(new MLAPI::InverseOperator());
      S->Reshape(A[level], type, p, &pushlist);
      MLAPI::MultiVector x(dspace, 1, false);
      x.Update(1.0, xref, 0.0);
      double rate =
          RichardsonS(fieldname, myrank, level, 4, damp, A[level], *S, x, f, false, true, false);
      if (rate < localrate)
      {
        localrate = rate;
        localtype = "symmetric Gauss-Seidel";
        localdamp = damp;
        localS = S;
      }
    }
    if (localrate < bestrate[level])
    {
      bestrate[level] = localrate;
      best.Type()[level] = "symmetric Gauss-Seidel";
      best.Damp()[level] = localdamp;
      best.S()[level] = localS;
    }
    if (!myrank)
    {
      printf("--------------------------------------------\n");
      printf("BEST SGS: %s level %d %s damp %f poly %d\n", field.c_str(), level, localtype.c_str(),
          localdamp, localpoly);
      printf("--------------------------------------------\n");
    }



    //----------------------------------------------------------- test Chebychev
    if (!myrank)
    {
      printf("--------------------------------------------\n");
      printf("Chebychev\n");
      printf("--------------------------------------------\n");
      fflush(stdout);
    }
    localrate = 10.0;
    localtype = "";
    localdamp = 1.0;
    localpoly = 1;
    for (int poly = 1; poly <= 12; ++poly)
    {
      if (!myrank)
      {
        printf("--------------------------------------------\n");
        printf("polynomial order %d\n", poly);
        printf("--------------------------------------------\n");
        fflush(stdout);
      }
      double damp = 1.0;
      std::string type = "";
      subp.set("smoother: type", "MLS");
      subp.set("smoother: sweeps", 1);
      subp.set("smoother: MLS polynomial order", poly);
      Teuchos::ParameterList p;
      SelectMLAPISmoother(type, level, subp, p, pushlist);
      S = Teuchos::rcp(new MLAPI::InverseOperator());
      S->Reshape(A[level], type, p, &pushlist);
      MLAPI::MultiVector x(dspace, 1, false);
      x.Update(1.0, xref, 0.0);
      double rate =
          RichardsonS(fieldname, myrank, level, 4, damp, A[level], *S, x, f, false, true, false);
      if (rate < localrate)
      {
        localrate = rate;
        localtype = "Chebychev";
        localdamp = damp;
        localpoly = poly;
        localS = S;
      }
    }
    if (localrate < bestrate[level])
    {
      bestrate[level] = localrate;
      best.Type()[level] = "Chebychev";
      best.Damp()[level] = localdamp;
      best.Poly()[level] = localpoly;
      best.S()[level] = localS;
    }
    if (!myrank)
    {
      printf("--------------------------------------------\n");
      printf("BEST Chebychev: %s level %d %s damp %f poly %d\n", field.c_str(), level,
          localtype.c_str(), localdamp, localpoly);
      printf("--------------------------------------------\n");
    }

    //----------------------------------------------------------- test ILU
    if (!myrank)
    {
      printf("--------------------------------------------\n");
      printf("ILU\n");
      printf("--------------------------------------------\n");
      fflush(stdout);
    }
    localrate = 10.0;
    localtype = "";
    localdamp = 1.0;
    localpoly = 1;
    for (int fill = 0; fill <= 4; ++fill)
    {
      if (!myrank)
      {
        printf("--------------------------------------------\n");
        printf("fill-in level %d\n", fill);
        printf("--------------------------------------------\n");
        fflush(stdout);
      }
      std::string type = "";
      subp.set("smoother: type", "IFPACK");
      subp.set("smoother: ifpack type", "ILU");
      subp.set("smoother: ifpack level-of-fill", fill * 1.0);
      subp.set("smoother: sweeps", 1);
      subp.set("smoother: damping factor", 1.0);
      Teuchos::ParameterList p;
      SelectMLAPISmoother(type, level, subp, p, pushlist);
      S = Teuchos::rcp(new MLAPI::InverseOperator());
      S->Reshape(A[level], type, p, NULL);
      for (int i = 0; i < 99; ++i)
      {
        double damp = 1.0 - i * 0.01;  // reduce damping factor in steps of 0.01;
        MLAPI::MultiVector x(dspace, 1, false);
        x.Update(1.0, xref, 0.0);
        double rate =
            RichardsonS(fieldname, myrank, level, 4, damp, A[level], *S, x, f, false, true, false);
        if (rate < localrate)
        {
          localrate = rate;
          localtype = "ILU";
          localdamp = damp;
          localpoly = fill;
          localS = S;
        }
      }  // i
    }    // fill
    if (localrate < bestrate[level])
    {
      bestrate[level] = localrate;
      best.Type()[level] = "ILU";
      best.Damp()[level] = localdamp;
      best.Poly()[level] = localpoly;
      best.S()[level] = localS;
    }
    if (!myrank)
    {
      printf("--------------------------------------------\n");
      printf("BEST ILU: %s level %d %s damp %f poly %d\n", field.c_str(), level, localtype.c_str(),
          localdamp, localpoly);
      printf("--------------------------------------------\n");
    }

    //------------------------------------------------------------------
    if (!myrank)
    {
      printf("============================================\n");
      printf("BEST OVERALL: %s level %d %s damp %f poly %d\n", field.c_str(), level,
          best.Type()[level].c_str(), best.Damp()[level], best.Poly()[level]);
      printf("============================================\n");
    }
  }  // loop levels

  //==================================================== analyze V-cycle
  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("%s V cycle using best smoothers determining sweeps\n", field.c_str());
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  // create coarse level solve
  best.S()[nlevel - 1] = Teuchos::rcp(new MLAPI::InverseOperator());
  best.S()[nlevel - 1]->Reshape(A[nlevel - 1], "Amesos-KLU");

  MLAPI::Space rspace(A[0].GetOperatorRangeSpace());
  MLAPI::Space dspace(A[0].GetOperatorDomainSpace());
  MLAPI::MultiVector xref(dspace);
  xref.Random();
  xref.Scale(10000.0);
  MLAPI::MultiVector f(rspace, true);

  // number of sweeps per level
  std::vector<int> localVsweeps(6, 1);
  std::vector<int> loops(6, 1);
  for (int i = 0; i < nlevel - 1; ++i) loops[i] = 5;
  double bestVrate = 10.0;

  for (int i = 1; i <= loops[0]; ++i)
    for (int j = 1; j <= loops[1]; ++j)
      for (int k = 1; k <= loops[2]; ++k)
        for (int l = 1; l <= loops[3]; ++l)
          for (int m = 1; m <= loops[4]; ++m)
            for (int n = 1; n <= loops[5]; ++n)
            {
              localVsweeps[0] = i;
              localVsweeps[1] = j;
              localVsweeps[2] = k;
              localVsweeps[3] = l;
              localVsweeps[4] = m;
              localVsweeps[5] = n;
              if (!myrank)
              {
                printf("--------------------------------------------\n");
                printf("V cycle testing sweeps ");
                for (int o = 0; o < nlevel - 1; ++o) printf("%d ", localVsweeps[o]);
                printf("\n");
                fflush(stdout);
              }
              MLAPI::MultiVector x(dspace, 1, false);
              x.Update(1.0, xref, 0.0);
              double rate = RichardsonV(fieldname, myrank, 3, 1.0, localVsweeps, best.Damp(), A,
                  best.S(), P, R, 0, nlevel, x, f, false, true, false);
              if (rate < bestVrate)
              {
                if (!myrank)
                {
                  printf("Current best\n");
                  fflush(stdout);
                }
                bestVrate = rate;
                for (int o = 0; o < 6; ++o) best.Sweeps()[o] = localVsweeps[o];
              }
            }
  if (!myrank)
  {
    printf("============================================\n");
    printf("BEST V cycle sweeps: ");
    for (int i = 0; i < nlevel - 1; ++i) printf("%d ", best.Sweeps()[i]);
    printf("\n");
    printf("============================================\n");
    fflush(stdout);
  }

  if (!myrank)
  {
    printf("============================================\n");
    printf("============================================\n");
    printf("Recommendation Summary %s\n", field.c_str());
    for (int level = 0; level < nlevel - 1; ++level)
    {
      printf("Level %d type %s damping %8.4e polynomial order %d sweeps %d\n", level,
          best.Type()[level].c_str(), best.Damp()[level], best.Poly()[level], best.Sweeps()[level]);
    }
    printf("Level %d direct solve\n", nlevel - 1);
    printf("============================================\n");
    printf("============================================\n");
    fflush(stdout);
  }


  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::Analyse_BGSAMG(const int myrank, AnalyzeBest& sbest,
    AnalyzeBest& fbest, AnalyzeBest& abest, std::vector<MLAPI::Operator>& Ass,
    std::vector<MLAPI::Operator>& Pss, std::vector<MLAPI::Operator>& Rss,
    std::vector<MLAPI::Operator>& Aff, std::vector<MLAPI::Operator>& Pff,
    std::vector<MLAPI::Operator>& Rff, std::vector<MLAPI::Operator>& Aaa,
    std::vector<MLAPI::Operator>& Paa, std::vector<MLAPI::Operator>& Raa,
    std::vector<MLAPI::Operator>& Asf, std::vector<MLAPI::Operator>& Afs,
    std::vector<MLAPI::Operator>& Afa, std::vector<MLAPI::Operator>& Aaf)
{
  // determine: sweeps and damping of Vcycle of individual field
  //            sweeps of entire BlockGaussSeidel on finest level
  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("Analyze BGS(AMG) (called 'PreconditionedKrylov' in input file)\n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  MLAPI::Space srspace(Ass[0].GetOperatorRangeSpace());
  MLAPI::Space sdspace(Ass[0].GetOperatorDomainSpace());
  MLAPI::Space frspace(Aff[0].GetOperatorRangeSpace());
  MLAPI::Space fdspace(Aff[0].GetOperatorDomainSpace());
  MLAPI::Space arspace(Aaa[0].GetOperatorRangeSpace());
  MLAPI::Space adspace(Aaa[0].GetOperatorDomainSpace());
  MLAPI::MultiVector sxref(sdspace);
  MLAPI::MultiVector fxref(fdspace);
  MLAPI::MultiVector axref(adspace);
  sxref.Random();
  fxref.Random();
  axref.Random();
  sxref.Scale(10000.0);
  fxref.Scale(10000.0);
  axref.Scale(10000.0);
  MLAPI::MultiVector sf(srspace, true);
  MLAPI::MultiVector ff(frspace, true);
  MLAPI::MultiVector af(arspace, true);

  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("Analyze BGS(AMG) sweeps\n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  std::vector<int> bestsweeps(3, 1);
  std::vector<double> bestdamps(3, 1.0);
  std::vector<int> localsweeps(3, 1);
  std::vector<double> localdamps(3, 1.0);
  double bestrate = 100.0;
  for (int sweeps = 1; sweeps <= 5; ++sweeps)
    for (int fweeps = 1; fweeps <= 5; ++fweeps)
      for (int aweeps = 1; aweeps <= 5; ++aweeps)
      {
        localsweeps[0] = sweeps;
        localsweeps[1] = fweeps;
        localsweeps[2] = aweeps;
        if (!myrank)
        {
          printf("--------------------------------------------\n");
          printf(
              "BGS(AMG) sweeps S/F/A %d/%d/%d\n", localsweeps[0], localsweeps[1], localsweeps[2]);
          fflush(stdout);
        }
        MLAPI::MultiVector sx(sdspace, 1, false);
        MLAPI::MultiVector fx(fdspace, 1, false);
        MLAPI::MultiVector ax(adspace, 1, false);
        sx.Update(1.0, sxref, 0.0);
        fx.Update(1.0, fxref, 0.0);
        ax.Update(1.0, axref, 0.0);
        MLAPI::MultiVector sf(srspace, true);
        MLAPI::MultiVector ff(frspace, true);
        MLAPI::MultiVector af(arspace, true);
        double rate = RichardsonBGS_V(myrank, 3, 1.0, localsweeps, localdamps, sbest, fbest, abest,
            sx, fx, ax, sf, ff, af, Ass, Pss, Rss, Aff, Pff, Rff, Aaa, Paa, Raa, Asf, Afs, Afa, Aaf,
            false, true, false);
        if (rate < bestrate)
        {
          if (!myrank)
          {
            printf("** Current best **\n");
            fflush(stdout);
          }
          bestrate = rate;
          for (int k = 0; k < 3; ++k) bestsweeps[k] = localsweeps[k];
        }
      }  // sweeps, fsweeps, asweeps

  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("Analyze BGS(AMG) damps\n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  for (int j = 0; j < 3; ++j)
  {
    for (int k = 0; k < 3; ++k) localdamps[k] = bestdamps[k];
    for (int i = 0; i < 50; ++i)
    {
      double damp = 1.0 - i * 0.02;
      localdamps[j] = damp;
      if (!myrank)
      {
        printf("--------------------------------------------\n");
        printf("BGS(AMG) damps S/F/A %10.5e/%10.5e/%10.5e\n", localdamps[0], localdamps[1],
            localdamps[2]);
        fflush(stdout);
      }
      MLAPI::MultiVector sx(sdspace, 1, false);
      MLAPI::MultiVector fx(fdspace, 1, false);
      MLAPI::MultiVector ax(adspace, 1, false);
      sx.Update(1.0, sxref, 0.0);
      fx.Update(1.0, fxref, 0.0);
      ax.Update(1.0, axref, 0.0);
      MLAPI::MultiVector sf(srspace, true);
      MLAPI::MultiVector ff(frspace, true);
      MLAPI::MultiVector af(arspace, true);
      double rate = RichardsonBGS_V(myrank, 3, 1.0, bestsweeps, localdamps, sbest, fbest, abest, sx,
          fx, ax, sf, ff, af, Ass, Pss, Rss, Aff, Pff, Rff, Aaa, Paa, Raa, Asf, Afs, Afa, Aaf,
          false, true, false);
      if (rate < bestrate)
      {
        if (!myrank)
        {
          printf("** Current best **\n");
          fflush(stdout);
        }
        bestrate = rate;
        for (int k = 0; k < 3; ++k) bestdamps[k] = localdamps[k];
      }
    }  // i
  }    // j


  if (!myrank)
  {
    printf("============================================\n");
    printf("============================================\n");
    printf("BGS(AMG) Recommendation Summary\n");
    printf("--------------------------------------------\n");
    printf("Structure AMG hierarchy:\n");
    for (int level = 0; level < sbest.Nlevel() - 1; ++level)
    {
      printf("Level %d %s damping %8.4e polynomial order %d sweeps %d\n", level,
          sbest.Type()[level].c_str(), sbest.Damp()[level], sbest.Poly()[level],
          sbest.Sweeps()[level]);
    }
    printf("Level %d direct solve\n", sbest.Nlevel() - 1);
    printf("--------------------------------------------\n");
    printf("Fluid AMG hierarchy:\n");
    for (int level = 0; level < fbest.Nlevel() - 1; ++level)
    {
      printf("Level %d %s damping %8.4e polynomial order %d sweeps %d\n", level,
          fbest.Type()[level].c_str(), fbest.Damp()[level], fbest.Poly()[level],
          fbest.Sweeps()[level]);
    }
    printf("Level %d direct solve\n", fbest.Nlevel() - 1);
    printf("--------------------------------------------\n");
    printf("Ale AMG hierarchy:\n");
    for (int level = 0; level < abest.Nlevel() - 1; ++level)
    {
      printf("Level %d %s damping %8.4e polynomial order %d sweeps %d\n", level,
          abest.Type()[level].c_str(), abest.Damp()[level], abest.Poly()[level],
          abest.Sweeps()[level]);
    }
    printf("Level %d direct solve\n", abest.Nlevel() - 1);
    printf("********************************************\n");
    printf("BGS(AMG):\n");
    printf("Structure spciter %d spcomega %10.5e\n", bestsweeps[0], bestdamps[0]);
    printf("Fluid     fpciter %d fpcomega %10.5e\n", bestsweeps[1], bestdamps[1]);
    printf("Ale       apciter %d apcomega %10.5e\n", bestsweeps[2], bestdamps[2]);
    printf("============================================\n");
    printf("============================================\n");
    fflush(stdout);
  }

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::Analyse_AMGBGS(const int myrank, AnalyzeBest& sbest,
    AnalyzeBest& fbest, AnalyzeBest& abest, std::vector<MLAPI::Operator>& Ass,
    std::vector<MLAPI::Operator>& Pss, std::vector<MLAPI::Operator>& Rss,
    std::vector<MLAPI::Operator>& Aff, std::vector<MLAPI::Operator>& Pff,
    std::vector<MLAPI::Operator>& Rff, std::vector<MLAPI::Operator>& Aaa,
    std::vector<MLAPI::Operator>& Paa, std::vector<MLAPI::Operator>& Raa,
    std::vector<MLAPI::Operator>& Asf, std::vector<MLAPI::Operator>& Afs,
    std::vector<MLAPI::Operator>& Afa, std::vector<MLAPI::Operator>& Aaf)
{
  // determine: sweeps and damping of Vcycle of individual field
  //            sweeps of entire BlockGaussSeidel on finest level
  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("Analyze AMG(BGS) (called 'FSIAMG' in input file)\n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  std::vector<int> bestsweeps[3];
  for (int i = 0; i < 3; ++i) bestsweeps[i].resize(minnlevel_, 1);
  std::vector<double> bestdamps[3];
  for (int i = 0; i < 3; ++i) bestdamps[i].resize(minnlevel_, 1);
  std::vector<int> localsweeps[3];
  for (int i = 0; i < 3; ++i) localsweeps[i].resize(minnlevel_, 1);
  std::vector<double> localdamps[3];
  for (int i = 0; i < 3; ++i) localdamps[i].resize(minnlevel_, 1.0);
  std::vector<double> bestrate(minnlevel_, 100.0);
  std::vector<int> bestVsweeps(6, 1);
  std::vector<int> localVsweeps(6, 1);
  std::vector<double> bestVdamps(6, 1.0);
  std::vector<double> localVdamps(6, 1.0);
  std::vector<int> loops(6, 1);
  for (int i = 0; i < minnlevel_; ++i) loops[i] = 5;
  double bestVrate = 10.0;



  //------------------------------------------------------------------------------------
  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("Analyze AMG(BGS) sweeps for individual fields\n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  for (int level = 0; level < minnlevel_; ++level)
  {
    bestrate[level] = 100.0;
    MLAPI::Space srspace(Ass[level].GetOperatorRangeSpace());
    MLAPI::Space sdspace(Ass[level].GetOperatorDomainSpace());
    MLAPI::Space frspace(Aff[level].GetOperatorRangeSpace());
    MLAPI::Space fdspace(Aff[level].GetOperatorDomainSpace());
    MLAPI::Space arspace(Aaa[level].GetOperatorRangeSpace());
    MLAPI::Space adspace(Aaa[level].GetOperatorDomainSpace());
    MLAPI::MultiVector sxref(sdspace);
    MLAPI::MultiVector fxref(fdspace);
    MLAPI::MultiVector axref(adspace);
    sxref.Random();
    fxref.Random();
    axref.Random();
    sxref.Scale(10000.0);
    fxref.Scale(10000.0);
    axref.Scale(10000.0);
    MLAPI::MultiVector sf(srspace, true);
    MLAPI::MultiVector ff(frspace, true);
    MLAPI::MultiVector af(arspace, true);
    for (int sweeps = 1; sweeps <= 5; ++sweeps)
      for (int fweeps = 1; fweeps <= 5; ++fweeps)
        for (int aweeps = 1; aweeps <= 5; ++aweeps)
        {
          localsweeps[0][level] = sweeps;
          localsweeps[1][level] = fweeps;
          localsweeps[2][level] = aweeps;
          if (!myrank)
          {
            printf("--------------------------------------------\n");
            printf("AMG(BGS) level %d sweeps S/F/A %d/%d/%d\n", level, localsweeps[0][level],
                localsweeps[1][level], localsweeps[2][level]);
            fflush(stdout);
          }
          MLAPI::MultiVector sx(sdspace, 1, false);
          MLAPI::MultiVector fx(fdspace, 1, false);
          MLAPI::MultiVector ax(adspace, 1, false);
          sx.Update(1.0, sxref, 0.0);
          fx.Update(1.0, fxref, 0.0);
          ax.Update(1.0, axref, 0.0);
          MLAPI::MultiVector sf(srspace, true);
          MLAPI::MultiVector ff(frspace, true);
          MLAPI::MultiVector af(arspace, true);
          double rate = RichardsonBGS_SV(myrank, 3, 1.0, localsweeps, localdamps, level, sbest,
              fbest, abest, sx, fx, ax, sf, ff, af, Ass, Pss, Rss, Aff, Pff, Rff, Aaa, Paa, Raa,
              Asf[level], Afs[level], Afa[level], Aaf[level], false, true, false);
          if (rate < bestrate[level])
          {
            if (!myrank)
            {
              printf("** Current best **\n");
              fflush(stdout);
            }
            bestrate[level] = rate;
            for (int k = 0; k < 3; ++k) bestsweeps[k][level] = localsweeps[k][level];
          }
        }  // sweeps, fsweeps, asweeps
  }

  //------------------------------------------------------------------------------------
  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("Analyze AMG(BGS) damps for individual fields\n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  for (int level = 0; level < minnlevel_; ++level)
  {
    bestrate[level] = 1.0;
    MLAPI::Space srspace(Ass[level].GetOperatorRangeSpace());
    MLAPI::Space sdspace(Ass[level].GetOperatorDomainSpace());
    MLAPI::Space frspace(Aff[level].GetOperatorRangeSpace());
    MLAPI::Space fdspace(Aff[level].GetOperatorDomainSpace());
    MLAPI::Space arspace(Aaa[level].GetOperatorRangeSpace());
    MLAPI::Space adspace(Aaa[level].GetOperatorDomainSpace());
    MLAPI::MultiVector sxref(sdspace);
    MLAPI::MultiVector fxref(fdspace);
    MLAPI::MultiVector axref(adspace);
    sxref.Random();
    fxref.Random();
    axref.Random();
    sxref.Scale(10000.0);
    fxref.Scale(10000.0);
    axref.Scale(10000.0);
    MLAPI::MultiVector sf(srspace, true);
    MLAPI::MultiVector ff(frspace, true);
    MLAPI::MultiVector af(arspace, true);
    for (int s = 0; s < 12; ++s)
      for (int f = 0; f < 12; ++f)
        for (int a = 0; a < 12; ++a)
        {
          localdamps[0][level] = 1.0 - s * 0.08;
          localdamps[1][level] = 1.0 - f * 0.08;
          localdamps[2][level] = 1.0 - a * 0.08;
          if (!myrank)
          {
            printf("--------------------------------------------\n");
            printf("AMG(BGS) (level %d) damps S/F/A %10.5e/%10.5e/%10.5e\n", level,
                localdamps[0][level], localdamps[1][level], localdamps[2][level]);
            fflush(stdout);
          }
          MLAPI::MultiVector sx(sdspace, 1, false);
          MLAPI::MultiVector fx(fdspace, 1, false);
          MLAPI::MultiVector ax(adspace, 1, false);
          sx.Update(1.0, sxref, 0.0);
          fx.Update(1.0, fxref, 0.0);
          ax.Update(1.0, axref, 0.0);
          MLAPI::MultiVector sf(srspace, true);
          MLAPI::MultiVector ff(frspace, true);
          MLAPI::MultiVector af(arspace, true);
          double rate = RichardsonBGS_SV(myrank, 3, 1.0, bestsweeps, localdamps, level, sbest,
              fbest, abest, sx, fx, ax, sf, ff, af, Ass, Pss, Rss, Aff, Pff, Rff, Aaa, Paa, Raa,
              Asf[level], Afs[level], Afa[level], Aaf[level], false, true, false);
          if (rate < bestrate[level])
          {
            if (!myrank)
            {
              printf("** Current best **\n");
              fflush(stdout);
            }
            bestrate[level] = rate;
            for (int k = 0; k < 3; ++k) bestdamps[k][level] = localdamps[k][level];
          }
        }  // s,f,a
  }        // level


  //------------------------------------------------------------------------------------
  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("Analyze AMG(BGS) block-Vcycle sweeps on individual levels \n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  MLAPI::Space srspace(Ass[0].GetOperatorRangeSpace());
  MLAPI::Space sdspace(Ass[0].GetOperatorDomainSpace());
  MLAPI::Space frspace(Aff[0].GetOperatorRangeSpace());
  MLAPI::Space fdspace(Aff[0].GetOperatorDomainSpace());
  MLAPI::Space arspace(Aaa[0].GetOperatorRangeSpace());
  MLAPI::Space adspace(Aaa[0].GetOperatorDomainSpace());
  MLAPI::MultiVector sxref(sdspace);
  MLAPI::MultiVector fxref(fdspace);
  MLAPI::MultiVector axref(adspace);
  sxref.Random();
  fxref.Random();
  axref.Random();
  sxref.Scale(10000.0);
  fxref.Scale(10000.0);
  axref.Scale(10000.0);
  MLAPI::MultiVector sf(srspace, true);
  MLAPI::MultiVector ff(frspace, true);
  MLAPI::MultiVector af(arspace, true);
  bestVrate = 10.0;
  for (int i = 1; i <= loops[0]; ++i)
    for (int j = 1; j <= loops[1]; ++j)
      for (int k = 1; k <= loops[2]; ++k)
        for (int l = 1; l <= loops[3]; ++l)
          for (int m = 1; m <= loops[4]; ++m)
            for (int n = 1; n <= loops[5]; ++n)
            {
              localVsweeps[0] = i;
              localVsweeps[1] = j;
              localVsweeps[2] = k;
              localVsweeps[3] = l;
              localVsweeps[4] = m;
              localVsweeps[5] = n;
              if (!myrank)
              {
                printf("--------------------------------------------\n");
                printf("AMG(BGS) block-V-cycle testing sweeps ");
                for (int o = 0; o < minnlevel_; ++o) printf("%d ", localVsweeps[o]);
                printf("\n");
                fflush(stdout);
              }
              MLAPI::MultiVector sx(sdspace, 1, false);
              MLAPI::MultiVector fx(fdspace, 1, false);
              MLAPI::MultiVector ax(adspace, 1, false);
              sx.Update(1.0, sxref, 0.0);
              fx.Update(1.0, fxref, 0.0);
              ax.Update(1.0, axref, 0.0);
              MLAPI::MultiVector sf(srspace, true);
              MLAPI::MultiVector ff(frspace, true);
              MLAPI::MultiVector af(arspace, true);
              double rate = Richardson_BlockV(myrank, 3, 1.0, localVsweeps, localVdamps, bestsweeps,
                  bestdamps, blocksmoother_, sbest, fbest, abest, sx, fx, ax, sf, ff, af, Ass, Pss,
                  Rss, Aff, Pff, Rff, Aaa, Paa, Raa, Asf, Afs, Afa, Aaf, false, true, false);
              if (rate < bestVrate)
              {
                if (!myrank)
                {
                  printf("** Current best **\n");
                  fflush(stdout);
                }
                bestVrate = rate;
                for (int o = 0; o < minnlevel_; ++o) bestVsweeps[o] = localVsweeps[o];
              }


            }  // i,j,k,l,m,n


  //------------------------------------------------------------------------------------
  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("Analyze AMG(BGS) block-Vcycle damps on individual levels \n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  bestVrate = 10.0;
  for (int level = 0; level < minnlevel_; ++level)
  {
    for (int i = 0; i < minnlevel_; ++i) localVdamps[i] = bestVdamps[i];
    for (int j = 0; j < 50; ++j)
    {
      localVdamps[level] = 1.0 - j * 0.02;
      if (!myrank)
      {
        printf("--------------------------------------------\n");
        printf("AMG(BGS) block-V-cycle testing damps ");
        for (int i = 0; i < minnlevel_; ++i) printf("(level %d) %10.5e ", i, localVdamps[i]);
        printf("\n");
        fflush(stdout);
      }  // !myrank
      MLAPI::MultiVector sx(sdspace, 1, false);
      MLAPI::MultiVector fx(fdspace, 1, false);
      MLAPI::MultiVector ax(adspace, 1, false);
      sx.Update(1.0, sxref, 0.0);
      fx.Update(1.0, fxref, 0.0);
      ax.Update(1.0, axref, 0.0);
      MLAPI::MultiVector sf(srspace, true);
      MLAPI::MultiVector ff(frspace, true);
      MLAPI::MultiVector af(arspace, true);
      double rate = Richardson_BlockV(myrank, 3, 1.0, bestVsweeps, localVdamps, bestsweeps,
          bestdamps, blocksmoother_, sbest, fbest, abest, sx, fx, ax, sf, ff, af, Ass, Pss, Rss,
          Aff, Pff, Rff, Aaa, Paa, Raa, Asf, Afs, Afa, Aaf, false, true, false);
      if (rate < bestVrate)
      {
        if (!myrank)
        {
          printf("** Current best **\n");
          fflush(stdout);
        }
        bestVrate = rate;
        for (int i = 0; i < minnlevel_; ++i) bestVdamps[i] = localVdamps[i];
      }

    }  // j
  }    // level
  // test damping parameters once more backwards starting from the previous optimal parameters
  for (int level = minnlevel_ - 1; level >= 0; --level)
  {
    for (int i = 0; i < minnlevel_; ++i) localVdamps[i] = bestVdamps[i];
    for (int j = 0; j < 50; ++j)
    {
      localVdamps[level] = 1.0 - j * 0.02;
      if (!myrank)
      {
        printf("--------------------------------------------\n");
        printf("AMG(BGS) block-V-cycle testing damps ");
        for (int i = 0; i < minnlevel_; ++i) printf("(level %d) %10.5e ", i, localVdamps[i]);
        printf("\n");
        fflush(stdout);
      }  // !myrank
      MLAPI::MultiVector sx(sdspace, 1, false);
      MLAPI::MultiVector fx(fdspace, 1, false);
      MLAPI::MultiVector ax(adspace, 1, false);
      sx.Update(1.0, sxref, 0.0);
      fx.Update(1.0, fxref, 0.0);
      ax.Update(1.0, axref, 0.0);
      MLAPI::MultiVector sf(srspace, true);
      MLAPI::MultiVector ff(frspace, true);
      MLAPI::MultiVector af(arspace, true);
      double rate = Richardson_BlockV(myrank, 3, 1.0, bestVsweeps, localVdamps, bestsweeps,
          bestdamps, blocksmoother_, sbest, fbest, abest, sx, fx, ax, sf, ff, af, Ass, Pss, Rss,
          Aff, Pff, Rff, Aaa, Paa, Raa, Asf, Afs, Afa, Aaf, false, true, false);
      if (rate < bestVrate)
      {
        if (!myrank)
        {
          printf("** Current best **\n");
          fflush(stdout);
        }
        bestVrate = rate;
        for (int i = 0; i < minnlevel_; ++i) bestVdamps[i] = localVdamps[i];
      }

    }  // j
  }    // level


  //------------------------------------------------------------------------------------
  if (!myrank)
  {
    printf("============================================\n");
    printf("============================================\n");
    printf("AMG(BGS) Recommendation Summary\n");
    printf("--------------------------------------------\n");
    printf("Structure AMG hierarchy:\n");
    for (int level = 0; level < sbest.Nlevel() - 1; ++level)
    {
      printf("Level %d %s damping %8.4e polynomial order %d sweeps %d\n", level,
          sbest.Type()[level].c_str(), sbest.Damp()[level], sbest.Poly()[level],
          sbest.Sweeps()[level]);
    }
    printf("Level %d direct solve\n", sbest.Nlevel() - 1);
    printf("--------------------------------------------\n");
    printf("Fluid AMG hierarchy:\n");
    for (int level = 0; level < fbest.Nlevel() - 1; ++level)
    {
      printf("Level %d %s damping %8.4e polynomial order %d sweeps %d\n", level,
          fbest.Type()[level].c_str(), fbest.Damp()[level], fbest.Poly()[level],
          fbest.Sweeps()[level]);
    }
    printf("Level %d direct solve\n", fbest.Nlevel() - 1);
    printf("--------------------------------------------\n");
    printf("Ale AMG hierarchy:\n");
    for (int level = 0; level < abest.Nlevel() - 1; ++level)
    {
      printf("Level %d %s damping %8.4e polynomial order %d sweeps %d\n", level,
          abest.Type()[level].c_str(), abest.Damp()[level], abest.Poly()[level],
          abest.Sweeps()[level]);
    }
    printf("Level %d direct solve\n", abest.Nlevel() - 1);
    printf("********************************************\n");
    printf("AMG(BGS):\n");
    for (int level = 0; level < minnlevel_; ++level)
    {
      printf("Structure (level %d) spciter %d spcomega %10.5e\n", level, bestsweeps[0][level],
          bestdamps[0][level]);
      printf("Fluid     (level %d) fpciter %d fpcomega %10.5e\n", level, bestsweeps[1][level],
          bestdamps[1][level]);
      printf("Ale       (level %d) apciter %d apcomega %10.5e\n", level, bestsweeps[2][level],
          bestdamps[2][level]);
      printf("BGS       (level %d)  pciter %d  pcomega %10.5e\n", level, bestVsweeps[level],
          bestVdamps[level]);
    }
    printf("============================================\n");
    printf("============================================\n");
    fflush(stdout);
  }



  return;
}
