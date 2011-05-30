#ifdef CCADISCRET

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
void FSI::OverlappingBlockMatrixFSIAMG::AnalyzeFSIAMG(
                       const int myrank,
                       const int snlevel,
                       const int fnlevel,
                       const int anlevel,
                       Teuchos::ParameterList& sparams,
                       Teuchos::ParameterList& fparams,
                       Teuchos::ParameterList& aparams,
                       vector<MLAPI::Operator>& Ass,
                       vector<MLAPI::Operator>& Aff,
                       vector<MLAPI::Operator>& Aaa,
                       vector<MLAPI::Operator>& Pss, vector<MLAPI::Operator>& Rss,
                       vector<MLAPI::Operator>& Pff, vector<MLAPI::Operator>& Rff,
                       vector<MLAPI::Operator>& Paa, vector<MLAPI::Operator>& Raa,
                       vector<MLAPI::Operator>& Asf,
                       vector<MLAPI::Operator>& Afs,
                       vector<MLAPI::Operator>& Afa,
                       vector<MLAPI::Operator>& Aaf,
                       ML* sml,
                       ML* fml,
                       ML* aml)
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
  Analyse_SingleField("(s)",myrank,"STRUCTURE",snlevel,sparams,Ass,Pss,Rss,sbest);

  //-----------------------------------------analysis of fluid field
  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("FLUID field analysis\n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  AnalyzeBest fbest(fnlevel);
  Analyse_SingleField("(f)",myrank,"FLUID",fnlevel,fparams,Aff,Pff,Rff,fbest);

  //-----------------------------------------analysis of ale field
  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("ALE field analysis\n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  AnalyzeBest abest(anlevel);
  Analyse_SingleField("(a)",myrank,"ALE",anlevel,aparams,Aaa,Paa,Raa,abest);


#if 1
  //---------------------------------------- run analysis of BGS(AMG)
  // which is called PreconditionedKrylov in our input file
  Analyse_BGSAMG(myrank,sbest,fbest,abest,
                 Ass,Pss,Rss,
                 Aff,Pff,Rff,
                 Aaa,Paa,Raa,
                 Asf,Afs,Afa,Aaf);
#endif  
  


#if 1
  //---------------------------------------- run analysis of AMG(BGS)
  // which is called FSIAMG in our input file
  Analyse_AMGBGS(myrank,sbest,fbest,abest,
                 Ass,Pss,Rss,
                 Aff,Pff,Rff,
                 Aaa,Paa,Raa,
                 Asf,Afs,Afa,Aaf);
#endif  



  exit(EXIT_SUCCESS);
  
  
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::Analyse_SingleField(
                            const string fieldname,
                            const int myrank,
                            string field,
                            const int nlevel,
                            Teuchos::ParameterList& params,
                            vector<MLAPI::Operator>& A,
                            vector<MLAPI::Operator>& P, 
                            vector<MLAPI::Operator>& R,
                            AnalyzeBest& best)
{
  RCP<MLAPI::InverseOperator> S;
  //RCP<MLAPI::LoadBalanceInverseOperator> lbS;

  // measure r_l2 decrease per time in bestrate
  vector<double> bestrate(nlevel-1,10.0);

  // loop levels in single field and test smoothers
  for (int level=0; level<nlevel-1; ++level)
  {
    if (!myrank)
    {
      printf("--------------------------------------------\n");
      printf("level %d\n",level);
      printf("--------------------------------------------\n");
      fflush(stdout);
    }
    Teuchos::ParameterList pushlist(params.sublist("smoother: ifpack list"));
    char levelstr[11];
    sprintf(levelstr,"(level %d)",level);
    Teuchos::ParameterList subp =  params.sublist("smoother: list "+(string)levelstr);
    //cout << "pushlist\n" << pushlist << "subp\n" << subp << endl;
        
    MLAPI::Space rspace(A[level].GetOperatorRangeSpace());
    MLAPI::Space dspace(A[level].GetOperatorDomainSpace());
    MLAPI::MultiVector xref(dspace);
    xref.Random();
    xref.Scale(10000.0);
    MLAPI::MultiVector f(rspace,true);
    
    double localrate = 10.0;
    string localtype = "";
    double localdamp = 1.0;
    int    localpoly = 1;
    RCP<MLAPI::InverseOperator> localS;

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
    for (int i=0; i<99; ++i)
    {
      double damp = 1.0 - i * 0.01; // reduce damping factor in steps of 0.01;
      string type = "";
      subp.set("smoother: type","symmetric Gauss-Seidel");
      subp.set("smoother: sweeps",1);
      subp.set("smoother: damping factor",damp);
      Teuchos::ParameterList p;
      SelectMLAPISmoother(type,level,subp,p,pushlist);
      S = rcp(new MLAPI::InverseOperator());
      S->Reshape(A[level],type,p,&pushlist);
      MLAPI::MultiVector x(dspace,1,false);
      x.Update(1.0,xref,0.0);
      double rate = RichardsonS(fieldname,myrank,level,4,damp,A[level],*S,x,f,false,true,false);
      if (rate<localrate)
      {
        localrate = rate;
        localtype = "symmetric Gauss-Seidel";
        localdamp = damp;
        localS    = S;
      }
    }
    if (localrate<bestrate[level])
    {
      bestrate[level] = localrate;
      best.Type()[level] = "symmetric Gauss-Seidel";
      best.Damp()[level] = localdamp;
      best.S()[level]    = localS;
    }
    if (!myrank) 
    {
      printf("--------------------------------------------\n");
      printf("BEST SGS: %s level %d %s damp %f poly %d\n",
              field.c_str(),level,localtype.c_str(),localdamp,localpoly);
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
    for (int poly=1; poly<=12; ++poly)
    {
      if (!myrank)
      {
        printf("--------------------------------------------\n");
        printf("polynomial order %d\n",poly);
        printf("--------------------------------------------\n");
        fflush(stdout);
      }
      double damp = 1.0; 
      string type = "";
      subp.set("smoother: type","MLS");
      subp.set("smoother: sweeps",1);
      subp.set("smoother: MLS polynomial order",poly);
      Teuchos::ParameterList p;
      SelectMLAPISmoother(type,level,subp,p,pushlist);
      S = rcp(new MLAPI::InverseOperator());
      S->Reshape(A[level],type,p,&pushlist);
      MLAPI::MultiVector x(dspace,1,false);
      x.Update(1.0,xref,0.0);
      double rate = RichardsonS(fieldname,myrank,level,4,damp,A[level],*S,x,f,false,true,false);
      if (rate<localrate)
      {
        localrate = rate;
        localtype = "Chebychev";
        localdamp = damp;
        localpoly = poly;
        localS    = S;
      }
    }
    if (localrate<bestrate[level])
    {
      bestrate[level] = localrate;
      best.Type()[level] = "Chebychev";
      best.Damp()[level] = localdamp;
      best.Poly()[level] = localpoly;
      best.S()[level]    = localS;
    }
    if (!myrank) 
    {
      printf("--------------------------------------------\n");
      printf("BEST Chebychev: %s level %d %s damp %f poly %d\n",
             field.c_str(),level,localtype.c_str(),localdamp,localpoly);
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
    for (int fill=0; fill<=4; ++fill)
    {
      if (!myrank)
      {
        printf("--------------------------------------------\n");
        printf("fill-in level %d\n",fill);
        printf("--------------------------------------------\n");
        fflush(stdout);
      }
      string type = "";
      subp.set("smoother: type","IFPACK");
      subp.set("smoother: ifpack type","ILU");
      subp.set("smoother: ifpack level-of-fill",fill*1.0);
      subp.set("smoother: sweeps",1);
      subp.set("smoother: damping factor",1.0);
      Teuchos::ParameterList p; 
      SelectMLAPISmoother(type,level,subp,p,pushlist);
      S = rcp(new MLAPI::InverseOperator());
      S->Reshape(A[level],type,p,NULL);
      for (int i=0; i<99; ++i)
      {
        double damp = 1.0 - i * 0.01; // reduce damping factor in steps of 0.01;
        MLAPI::MultiVector x(dspace,1,false);
        x.Update(1.0,xref,0.0);
        double rate = RichardsonS(fieldname,myrank,level,4,damp,A[level],*S,x,f,false,true,false);
        if (rate<localrate)
        {
          localrate = rate;
          localtype = "ILU";
          localdamp = damp;
          localpoly = fill;
          localS    = S;
        }
      } // i
    } // fill
    if (localrate<bestrate[level])
    {
      bestrate[level] = localrate;
      best.Type()[level] = "ILU";
      best.Damp()[level] = localdamp;
      best.Poly()[level] = localpoly;
      best.S()[level]    = localS;
    }
    if (!myrank) 
    {
      printf("--------------------------------------------\n");
      printf("BEST ILU: %s level %d %s damp %f poly %d\n",
             field.c_str(),level,localtype.c_str(),localdamp,localpoly);
      printf("--------------------------------------------\n");
    }

    //------------------------------------------------------------------
    if (!myrank) 
    {
      printf("============================================\n");
      printf("BEST OVERALL: %s level %d %s damp %f poly %d\n",
             field.c_str(),level,best.Type()[level].c_str(),best.Damp()[level],best.Poly()[level]);
      printf("============================================\n");
    }
  } // loop levels
  
  //==================================================== analyze V-cycle
  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("%s V cycle using best smoothers determining sweeps\n",field.c_str());
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  // create coarse level solve
  best.S()[nlevel-1] = rcp(new MLAPI::InverseOperator());
  best.S()[nlevel-1]->Reshape(A[nlevel-1],"Amesos-KLU");

  MLAPI::Space rspace(A[0].GetOperatorRangeSpace());
  MLAPI::Space dspace(A[0].GetOperatorDomainSpace());
  MLAPI::MultiVector xref(dspace);
  xref.Random();
  xref.Scale(10000.0);
  MLAPI::MultiVector f(rspace,true);
  
  // number of sweeps per level
  vector<int> localVsweeps(6,1);
  vector<int> loops(6,1); for (int i=0; i<nlevel-1; ++i) loops[i] = 5;
  double bestVrate = 10.0;
  
  for (int i=1; i<=loops[0]; ++i)
    for (int j=1; j<=loops[1]; ++j)
      for (int k=1; k<=loops[2]; ++k)
        for (int l=1; l<=loops[3]; ++l)
          for (int m=1; m<=loops[4]; ++m)
            for (int n=1; n<=loops[5]; ++n)
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
                for (int o=0; o<nlevel-1; ++o) printf("%d ",localVsweeps[o]);
                printf("\n");
                fflush(stdout);
              }
              MLAPI::MultiVector x(dspace,1,false);
              x.Update(1.0,xref,0.0);
              double rate = RichardsonV(fieldname,myrank,3,1.0,localVsweeps,best.Damp(),A,best.S(),P,R,0,nlevel,x,f,false,true,false);
              if (rate<bestVrate)
              {
                if (!myrank) printf("Current best\n"); fflush(stdout);
                bestVrate = rate;
                for (int o=0; o<6; ++o) best.Sweeps()[o] = localVsweeps[o];
              }
            }
  if (!myrank)
  {
    printf("============================================\n");
    printf("BEST V cycle sweeps: ");
    for (int i=0; i<nlevel-1; ++i) printf("%d ",best.Sweeps()[i]);
    printf("\n");
    printf("============================================\n");
    fflush(stdout);
  }

  if (!myrank)
  {
    printf("============================================\n");
    printf("============================================\n");
    printf("Recommendation Summary %s\n",field.c_str());
    for (int level=0; level<nlevel-1; ++level)
    {
    printf("Level %d type %s damping %8.4e polynomial order %d sweeps %d\n",
           level,best.Type()[level].c_str(),best.Damp()[level],best.Poly()[level],best.Sweeps()[level]);
    }
    printf("Level %d direct solve\n",nlevel-1);
    printf("============================================\n");
    printf("============================================\n");
    fflush(stdout);
  }

  
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::Analyse_BGSAMG(
                       const int myrank,
                       AnalyzeBest& sbest,
                       AnalyzeBest& fbest,
                       AnalyzeBest& abest,  
                       vector<MLAPI::Operator>& Ass,
                       vector<MLAPI::Operator>& Pss, vector<MLAPI::Operator>& Rss,
                       vector<MLAPI::Operator>& Aff,
                       vector<MLAPI::Operator>& Pff, vector<MLAPI::Operator>& Rff,
                       vector<MLAPI::Operator>& Aaa,
                       vector<MLAPI::Operator>& Paa, vector<MLAPI::Operator>& Raa,
                       vector<MLAPI::Operator>& Asf,
                       vector<MLAPI::Operator>& Afs,
                       vector<MLAPI::Operator>& Afa,
                       vector<MLAPI::Operator>& Aaf
                       )
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
  MLAPI::MultiVector sf(srspace,true);
  MLAPI::MultiVector ff(frspace,true);
  MLAPI::MultiVector af(arspace,true);

  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("Analyze BGS(AMG) sweeps\n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  vector<int>    bestsweeps(3,1);
  vector<double> bestdamps(3,1.0);
  vector<int>    localsweeps(3,1);
  vector<double> localdamps(3,1.0);
  double bestrate = 100.0;
  for (int sweeps=1; sweeps<=5; ++sweeps)
    for (int fweeps=1; fweeps<=5; ++fweeps)
      for (int aweeps=1; aweeps<=5; ++aweeps)
      {
        localsweeps[0] = sweeps;
        localsweeps[1] = fweeps;
        localsweeps[2] = aweeps;
        if (!myrank)
        {
          printf("--------------------------------------------\n");
          printf("BGS(AMG) sweeps S/F/A %d/%d/%d\n",
                 localsweeps[0],localsweeps[1],localsweeps[2]);
          fflush(stdout);
        }
        MLAPI::MultiVector sx(sdspace,1,false);
        MLAPI::MultiVector fx(fdspace,1,false);
        MLAPI::MultiVector ax(adspace,1,false);
        sx.Update(1.0,sxref,0.0);
        fx.Update(1.0,fxref,0.0);
        ax.Update(1.0,axref,0.0);
        MLAPI::MultiVector sf(srspace,true);
        MLAPI::MultiVector ff(frspace,true);
        MLAPI::MultiVector af(arspace,true);
        double rate = RichardsonBGS_V(myrank,3,1.0,localsweeps,localdamps,
                                      sbest,fbest,abest,
                                      sx,fx,ax,sf,ff,af,
                                      Ass,Pss,Rss,Aff,Pff,Rff,Aaa,Paa,Raa,
                                      Asf,Afs,Afa,Aaf,
                                      false,true,false);
        if (rate<bestrate)
        {
          if (!myrank) printf("** Current best **\n"); fflush(stdout);
          bestrate = rate;
          for (int k=0; k<3; ++k)
            bestsweeps[k] = localsweeps[k];
        }
      } // sweeps, fsweeps, asweeps

  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("Analyze BGS(AMG) damps\n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  for (int j=0; j<3; ++j)
  {
    for (int k=0; k<3; ++k) localdamps[k] = bestdamps[k];
    for (int i=0; i<50; ++i)
    {
      double damp = 1.0-i*0.02;
      localdamps[j] = damp;
      if (!myrank)
      {
        printf("--------------------------------------------\n");
        printf("BGS(AMG) damps S/F/A %10.5e/%10.5e/%10.5e\n",
               localdamps[0],localdamps[1],localdamps[2]);
        fflush(stdout);
      }
      MLAPI::MultiVector sx(sdspace,1,false);
      MLAPI::MultiVector fx(fdspace,1,false);
      MLAPI::MultiVector ax(adspace,1,false);
      sx.Update(1.0,sxref,0.0);
      fx.Update(1.0,fxref,0.0);
      ax.Update(1.0,axref,0.0);
      MLAPI::MultiVector sf(srspace,true);
      MLAPI::MultiVector ff(frspace,true);
      MLAPI::MultiVector af(arspace,true);
      double rate = RichardsonBGS_V(myrank,3,1.0,bestsweeps,localdamps,
                                    sbest,fbest,abest,
                                    sx,fx,ax,sf,ff,af,
                                    Ass,Pss,Rss,Aff,Pff,Rff,Aaa,Paa,Raa,
                                    Asf,Afs,Afa,Aaf,
                                    false,true,false);
      if (rate<bestrate)
      {
        if (!myrank) printf("** Current best **\n"); fflush(stdout);
        bestrate = rate;
        for (int k=0; k<3; ++k)
          bestdamps[k] = localdamps[k];
      }
    } // i
  } // j


  if (!myrank)
  {
    printf("============================================\n");
    printf("============================================\n");
    printf("BGS(AMG) Recommendation Summary\n");
    printf("--------------------------------------------\n");
    printf("Structure AMG hierarchy:\n");
    for (int level=0; level<sbest.Nlevel()-1; ++level)
    {
    printf("Level %d %s damping %8.4e polynomial order %d sweeps %d\n",
           level,sbest.Type()[level].c_str(),sbest.Damp()[level],sbest.Poly()[level],sbest.Sweeps()[level]);
    }
    printf("Level %d direct solve\n",sbest.Nlevel()-1);
    printf("--------------------------------------------\n");
    printf("Fluid AMG hierarchy:\n");
    for (int level=0; level<fbest.Nlevel()-1; ++level)
    {
    printf("Level %d %s damping %8.4e polynomial order %d sweeps %d\n",
           level,fbest.Type()[level].c_str(),fbest.Damp()[level],fbest.Poly()[level],fbest.Sweeps()[level]);
    }
    printf("Level %d direct solve\n",fbest.Nlevel()-1);
    printf("--------------------------------------------\n");
    printf("Ale AMG hierarchy:\n");
    for (int level=0; level<abest.Nlevel()-1; ++level)
    {
    printf("Level %d %s damping %8.4e polynomial order %d sweeps %d\n",
           level,abest.Type()[level].c_str(),abest.Damp()[level],abest.Poly()[level],abest.Sweeps()[level]);
    }
    printf("Level %d direct solve\n",abest.Nlevel()-1);
    printf("********************************************\n");
    printf("BGS(AMG):\n");
    printf("Structure spciter %d spcomega %10.5e\n",bestsweeps[0],bestdamps[0]);
    printf("Fluid     fpciter %d fpcomega %10.5e\n",bestsweeps[1],bestdamps[1]);
    printf("Ale       apciter %d apcomega %10.5e\n",bestsweeps[2],bestdamps[2]);
    printf("============================================\n");
    printf("============================================\n");
    fflush(stdout);
  }

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::Analyse_AMGBGS(
                       const int myrank,
                       AnalyzeBest& sbest,
                       AnalyzeBest& fbest,
                       AnalyzeBest& abest,  
                       vector<MLAPI::Operator>& Ass,
                       vector<MLAPI::Operator>& Pss, vector<MLAPI::Operator>& Rss,
                       vector<MLAPI::Operator>& Aff,
                       vector<MLAPI::Operator>& Pff, vector<MLAPI::Operator>& Rff,
                       vector<MLAPI::Operator>& Aaa,
                       vector<MLAPI::Operator>& Paa, vector<MLAPI::Operator>& Raa,
                       vector<MLAPI::Operator>& Asf,
                       vector<MLAPI::Operator>& Afs,
                       vector<MLAPI::Operator>& Afa,
                       vector<MLAPI::Operator>& Aaf
                       )
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
  vector<int> bestsweeps[3];   for (int i=0; i<3; ++i) bestsweeps[i].resize(minnlevel_,1);
  vector<double> bestdamps[3]; for (int i=0; i<3; ++i) bestdamps[i].resize(minnlevel_,1);
  vector<int> localsweeps[3];   for (int i=0; i<3; ++i) localsweeps[i].resize(minnlevel_,1);
  vector<double> localdamps[3]; for (int i=0; i<3; ++i) localdamps[i].resize(minnlevel_,1.0);
  vector<double> bestrate(minnlevel_,100.0);
  vector<int> bestVsweeps(6,1);
  vector<int> localVsweeps(6,1);
  vector<double> bestVdamps(6,1.0);
  vector<double> localVdamps(6,1.0);
  vector<int> loops(6,1); for (int i=0; i<minnlevel_; ++i) loops[i] = 5;
  double bestVrate = 10.0;



  //------------------------------------------------------------------------------------
  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("Analyze AMG(BGS) sweeps for individual fields\n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  for (int level=0; level<minnlevel_; ++level)
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
    MLAPI::MultiVector sf(srspace,true);
    MLAPI::MultiVector ff(frspace,true);
    MLAPI::MultiVector af(arspace,true);
    for (int sweeps=1; sweeps<=5; ++sweeps)
      for (int fweeps=1; fweeps<=5; ++fweeps)
        for (int aweeps=1; aweeps<=5; ++aweeps)
        {
          localsweeps[0][level] = sweeps;
          localsweeps[1][level] = fweeps;
          localsweeps[2][level] = aweeps;
          if (!myrank)
          {
            printf("--------------------------------------------\n");
            printf("AMG(BGS) level %d sweeps S/F/A %d/%d/%d\n",
                   level,localsweeps[0][level],localsweeps[1][level],localsweeps[2][level]);
            fflush(stdout);
          }
          MLAPI::MultiVector sx(sdspace,1,false);
          MLAPI::MultiVector fx(fdspace,1,false);
          MLAPI::MultiVector ax(adspace,1,false);
          sx.Update(1.0,sxref,0.0);
          fx.Update(1.0,fxref,0.0);
          ax.Update(1.0,axref,0.0);
          MLAPI::MultiVector sf(srspace,true);
          MLAPI::MultiVector ff(frspace,true);
          MLAPI::MultiVector af(arspace,true);
          double rate = RichardsonBGS_SV(myrank,3,1.0,localsweeps,localdamps,level,
                                         sbest,fbest,abest,sx,fx,ax,sf,ff,af,
                                         Ass,Pss,Rss,Aff,Pff,Rff,Aaa,Paa,Raa,
                                         Asf[level],Afs[level],Afa[level],Aaf[level],
                                         false,true,false);
          if (rate<bestrate[level])
          {
            if (!myrank) printf("** Current best **\n"); fflush(stdout);
            bestrate[level] = rate;
            for (int k=0; k<3; ++k)
              bestsweeps[k][level] = localsweeps[k][level];
          }
        } // sweeps, fsweeps, asweeps
  }

  //------------------------------------------------------------------------------------
  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("Analyze AMG(BGS) damps for individual fields\n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  for (int level=0; level<minnlevel_; ++level)
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
    MLAPI::MultiVector sf(srspace,true);
    MLAPI::MultiVector ff(frspace,true);
    MLAPI::MultiVector af(arspace,true);
    for (int s=0; s<12; ++s)
      for (int f=0; f<12; ++f)  
        for (int a=0; a<12; ++a)
        {
          localdamps[0][level] = 1.0-s*0.08;
          localdamps[1][level] = 1.0-f*0.08;
          localdamps[2][level] = 1.0-a*0.08;
          if (!myrank)
          {
            printf("--------------------------------------------\n");
            printf("AMG(BGS) (level %d) damps S/F/A %10.5e/%10.5e/%10.5e\n",
                   level,localdamps[0][level],localdamps[1][level],localdamps[2][level]);
            fflush(stdout);
          }
          MLAPI::MultiVector sx(sdspace,1,false);
          MLAPI::MultiVector fx(fdspace,1,false);
          MLAPI::MultiVector ax(adspace,1,false);
          sx.Update(1.0,sxref,0.0);
          fx.Update(1.0,fxref,0.0);
          ax.Update(1.0,axref,0.0);
          MLAPI::MultiVector sf(srspace,true);
          MLAPI::MultiVector ff(frspace,true);
          MLAPI::MultiVector af(arspace,true);
            double rate = RichardsonBGS_SV(myrank,3,1.0,bestsweeps,localdamps,
                                           level,sbest,fbest,abest,sx,fx,ax,sf,ff,af,
                                           Ass,Pss,Rss,Aff,Pff,Rff,Aaa,Paa,Raa,
                                           Asf[level],Afs[level],Afa[level],Aaf[level],
                                           false,true,false);
          if (rate<bestrate[level])
          {
            if (!myrank) printf("** Current best **\n"); fflush(stdout);
            bestrate[level] = rate;
            for (int k=0; k<3; ++k)
              bestdamps[k][level] = localdamps[k][level];
          }
        } // s,f,a
  } // level


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
  MLAPI::MultiVector sf(srspace,true);
  MLAPI::MultiVector ff(frspace,true);
  MLAPI::MultiVector af(arspace,true);
  bestVrate = 10.0;
  for (int i=1; i<=loops[0]; ++i)
    for (int j=1; j<=loops[1]; ++j)
      for (int k=1; k<=loops[2]; ++k)
        for (int l=1; l<=loops[3]; ++l)
          for (int m=1; m<=loops[4]; ++m)
            for (int n=1; n<=loops[5]; ++n)
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
                for (int o=0; o<minnlevel_; ++o) printf("%d ",localVsweeps[o]);
                printf("\n");
                fflush(stdout);
              }
              MLAPI::MultiVector sx(sdspace,1,false);
              MLAPI::MultiVector fx(fdspace,1,false);
              MLAPI::MultiVector ax(adspace,1,false);
              sx.Update(1.0,sxref,0.0);
              fx.Update(1.0,fxref,0.0);
              ax.Update(1.0,axref,0.0);
              MLAPI::MultiVector sf(srspace,true);
              MLAPI::MultiVector ff(frspace,true);
              MLAPI::MultiVector af(arspace,true);
              double rate = Richardson_BlockV(myrank,3,1.0,localVsweeps,localVdamps,bestsweeps,bestdamps,
                                              sbest,fbest,abest,sx,fx,ax,sf,ff,af,
                                              Ass,Pss,Rss,Aff,Pff,Rff,Aaa,Paa,Raa,Asf,Afs,Afa,Aaf,
                                              false,true,false);
              if (rate<bestVrate)
              {
                if (!myrank) printf("** Current best **\n"); fflush(stdout);
                bestVrate = rate;
                for (int o=0; o<minnlevel_; ++o)
                  bestVsweeps[o] = localVsweeps[o];
              }
              
              
  } // i,j,k,l,m,n  


  //------------------------------------------------------------------------------------
  if (!myrank)
  {
    printf("--------------------------------------------\n");
    printf("Analyze AMG(BGS) block-Vcycle damps on individual levels \n");
    printf("--------------------------------------------\n");
    fflush(stdout);
  }
  bestVrate = 10.0;
  for (int level=0; level<minnlevel_; ++level)
  {
    for (int i=0; i<minnlevel_; ++i) localVdamps[i] = bestVdamps[i];
    for (int j=0; j<50; ++j)
    {
      localVdamps[level] = 1.0 - j * 0.02;
      if (!myrank)
      {
        printf("--------------------------------------------\n");
        printf("AMG(BGS) block-V-cycle testing damps ");
        for (int i=0; i<minnlevel_; ++i) printf("(level %d) %10.5e ",i,localVdamps[i]);
        printf("\n");
        fflush(stdout);
      } // !myrank
      MLAPI::MultiVector sx(sdspace,1,false);
      MLAPI::MultiVector fx(fdspace,1,false);
      MLAPI::MultiVector ax(adspace,1,false);
      sx.Update(1.0,sxref,0.0);
      fx.Update(1.0,fxref,0.0);
      ax.Update(1.0,axref,0.0);
      MLAPI::MultiVector sf(srspace,true);
      MLAPI::MultiVector ff(frspace,true);
      MLAPI::MultiVector af(arspace,true);
      double rate = Richardson_BlockV(myrank,3,1.0,bestVsweeps,localVdamps,bestsweeps,bestdamps,
                                      sbest,fbest,abest,sx,fx,ax,sf,ff,af,
                                      Ass,Pss,Rss,Aff,Pff,Rff,Aaa,Paa,Raa,Asf,Afs,Afa,Aaf,
                                      false,true,false);
      if (rate<bestVrate)
      {
        if (!myrank) printf("** Current best **\n"); fflush(stdout);
        bestVrate = rate;
        for (int i=0; i<minnlevel_; ++i) bestVdamps[i] = localVdamps[i];
      }
      
    } // j
  } // level
  // test damping parameters once more backwards starting from the previous optimal parameters
  for (int level=minnlevel_-1; level>=0; --level)
  {
    for (int i=0; i<minnlevel_; ++i) localVdamps[i] = bestVdamps[i];
    for (int j=0; j<50; ++j)
    {
      localVdamps[level] = 1.0 - j * 0.02;
      if (!myrank)
      {
        printf("--------------------------------------------\n");
        printf("AMG(BGS) block-V-cycle testing damps ");
        for (int i=0; i<minnlevel_; ++i) printf("(level %d) %10.5e ",i,localVdamps[i]);
        printf("\n");
        fflush(stdout);
      } // !myrank
      MLAPI::MultiVector sx(sdspace,1,false);
      MLAPI::MultiVector fx(fdspace,1,false);
      MLAPI::MultiVector ax(adspace,1,false);
      sx.Update(1.0,sxref,0.0);
      fx.Update(1.0,fxref,0.0);
      ax.Update(1.0,axref,0.0);
      MLAPI::MultiVector sf(srspace,true);
      MLAPI::MultiVector ff(frspace,true);
      MLAPI::MultiVector af(arspace,true);
      double rate = Richardson_BlockV(myrank,3,1.0,bestVsweeps,localVdamps,bestsweeps,bestdamps,
                                      sbest,fbest,abest,sx,fx,ax,sf,ff,af,
                                      Ass,Pss,Rss,Aff,Pff,Rff,Aaa,Paa,Raa,Asf,Afs,Afa,Aaf,
                                      false,true,false);
      if (rate<bestVrate)
      {
        if (!myrank) printf("** Current best **\n"); fflush(stdout);
        bestVrate = rate;
        for (int i=0; i<minnlevel_; ++i) bestVdamps[i] = localVdamps[i];
      }
      
    } // j
  } // level


  //------------------------------------------------------------------------------------
  if (!myrank)
  {
    printf("============================================\n");
    printf("============================================\n");
    printf("AMG(BGS) Recommendation Summary\n");
    printf("--------------------------------------------\n");
    printf("Structure AMG hierarchy:\n");
    for (int level=0; level<sbest.Nlevel()-1; ++level)
    {
    printf("Level %d %s damping %8.4e polynomial order %d sweeps %d\n",
           level,sbest.Type()[level].c_str(),sbest.Damp()[level],sbest.Poly()[level],sbest.Sweeps()[level]);
    }
    printf("Level %d direct solve\n",sbest.Nlevel()-1);
    printf("--------------------------------------------\n");
    printf("Fluid AMG hierarchy:\n");
    for (int level=0; level<fbest.Nlevel()-1; ++level)
    {
    printf("Level %d %s damping %8.4e polynomial order %d sweeps %d\n",
           level,fbest.Type()[level].c_str(),fbest.Damp()[level],fbest.Poly()[level],fbest.Sweeps()[level]);
    }
    printf("Level %d direct solve\n",fbest.Nlevel()-1);
    printf("--------------------------------------------\n");
    printf("Ale AMG hierarchy:\n");
    for (int level=0; level<abest.Nlevel()-1; ++level)
    {
    printf("Level %d %s damping %8.4e polynomial order %d sweeps %d\n",
           level,abest.Type()[level].c_str(),abest.Damp()[level],abest.Poly()[level],abest.Sweeps()[level]);
    }
    printf("Level %d direct solve\n",abest.Nlevel()-1);
    printf("********************************************\n");
    printf("AMG(BGS):\n");
    for (int level=0; level<minnlevel_; ++level)
    {
    printf("Structure (level %d) spciter %d spcomega %10.5e\n",level,bestsweeps[0][level],bestdamps[0][level]);
    printf("Fluid     (level %d) fpciter %d fpcomega %10.5e\n",level,bestsweeps[1][level],bestdamps[1][level]);
    printf("Ale       (level %d) apciter %d apcomega %10.5e\n",level,bestsweeps[2][level],bestdamps[2][level]);
    printf("BGS       (level %d)  pciter %d  pcomega %10.5e\n",level,bestVsweeps[level],bestVdamps[level]);
    }
    printf("============================================\n");
    printf("============================================\n");
    fflush(stdout);
  }



  
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double FSI::OverlappingBlockMatrixFSIAMG::RichardsonS(
                    const string field,
                    const int myrank,
                    const int level,
                    const int sweeps,
                    const double damp,
                    const MLAPI::Operator& A,
                    const MLAPI::InverseOperator& S,
                    MLAPI::MultiVector& x,
                    const MLAPI::MultiVector& f,
                    bool initiguesszero,
                    bool analysis,
                    bool silent) const
{
#if (FSIAMG_ANALYSIS>=3)
  analysis = true;
  silent = false;
#endif

  MLAPI::MultiVector r(f.GetVectorSpace(),1,false);
  if (initiguesszero) r = f;
  else                r = f - A * x;
  
  double initrinf = 0.0;
  double initrl2 = 0.0;
  double rinf=0.0;
  double rl2=0.0;
  if (analysis) 
  {
    initrl2 = r.Norm2();
    initrinf = r.NormInf();
  }
  
  MLAPI::MultiVector tmpx(x.GetVectorSpace(),1,false);
  
  Epetra_Time timer(MLAPI::GetEpetra_Comm());
  
  for (int i=1; i<=sweeps; ++i)
  {
    tmpx = 0.0;
    S.Apply(r,tmpx);
    //x = x + damp * tmpx;
    x.Update(damp,tmpx,1.0);
    if (i<sweeps || analysis) r = f - A * x;
  }

  if (analysis)
  {
    rl2 = r.Norm2();
    rinf = r.NormInf();
    double t = timer.ElapsedTime();
    double rl2rate = Rate(myrank,t,rl2,initrl2,r.GetGlobalLength());
    double rinfrate = Rate(myrank,t,rinf,initrinf,r.GetGlobalLength());
    if (!myrank && !silent) printf("RichardsonS %s      (level %2d) r0 %10.5e rl_2 %10.5e t %10.5e damp %10.5e sweeps %2d l2rate %10.5e\n",
                                   field.c_str(),level,initrl2,rl2,t,damp,sweeps,rl2rate);
    rl2rate = max(rl2rate,rinfrate);
  
  MLAPI::GetEpetra_Comm().Broadcast(&rl2rate,1,0);
  return rl2rate;
  }
  else return 0.0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double FSI::OverlappingBlockMatrixFSIAMG::RichardsonMixed(
                    const string field,
                    const int myrank,
                    const int level,
                    const int sweeps,
                    const double damp,
                    const MLAPI::Operator& A,
                    const LINALG::SparseMatrix& matrix,
                    const Teuchos::RCP<LINALG::Preconditioner>& solver,
                    MLAPI::MultiVector& x,
                    const MLAPI::MultiVector& f,
                    int& run,
                    bool initiguesszero,
                    bool analysis,
                    bool silent) const
{
#if (FSIAMG_ANALYSIS>=3)
  analysis = true;
  silent = false;
#endif
  
  MLAPI::MultiVector r(A.GetRangeSpace(),1,true);
  MLAPI::MultiVector tmpx(A.GetDomainSpace(),1,true);
  //MLAPI::MultiVector r(f.GetVectorSpace(),1,true);
  //MLAPI::MultiVector tmpx(x.GetVectorSpace(),1,true);
  
  if (initiguesszero) r = f;
  else                r = f - A * x;
  
  // a view to these vectors (MUST be AFTER the above statement!)
  RCP<Epetra_Vector> er = rcp(new Epetra_Vector(View,matrix.RangeMap(),r.GetValues(0)));
  RCP<Epetra_Vector> etmpx = rcp(new Epetra_Vector(View,matrix.DomainMap(),tmpx.GetValues(0)));

  double initrinf = 0.0;
  double initrl2 = 0.0;
  double rinf=0.0;
  double rl2=0.0;
  if (analysis) 
  {
    initrl2 = r.Norm2();
    initrinf = r.NormInf();
  }
  
  Epetra_Time timer(MLAPI::GetEpetra_Comm());
  
  for (int i=1; i<=sweeps; ++i)
  {
    bool refactor = false;
    if (!run) refactor = true;
    tmpx = 0.0;
    solver->Solve(matrix.EpetraMatrix(),etmpx,er,refactor,false);
    x = x + damp * tmpx;
    if (i<sweeps || analysis) 
    {
      r = f - A * x;
      // r is recreated here, so we need a fresh view (took me a while to find this one :-( )
      er = rcp(new Epetra_Vector(View,matrix.RangeMap(),r.GetValues(0)));
    }
    run++;
  }

  if (analysis)
  {
    rl2 = r.Norm2();
    rinf = r.NormInf();
    double t = timer.ElapsedTime();
    double rl2rate = Rate(myrank,t,rl2,initrl2,r.GetGlobalLength());
    double rinfrate = Rate(myrank,t,rinf,initrinf,r.GetGlobalLength());

    if (!myrank && !silent) printf("RichardsonMixed %s  (level %2d) r0 %10.5e rl_2 %10.5e t %10.5e damp %10.5e sweeps %2d l2rate %10.5e\n",
                                   field.c_str(),level,initrl2,rl2,t,damp,sweeps,rl2rate);
    rl2rate = max(rl2rate,rinfrate);
  
  MLAPI::GetEpetra_Comm().Broadcast(&rl2rate,1,0);
  return rl2rate;
  }
  else return 0.0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double FSI::OverlappingBlockMatrixFSIAMG::RichardsonV(
                    const string field,
                    const int myrank,
                    int sweeps,
                    const double damp,
                    vector<int>& levelsweeps,
                    vector<double>& leveldamps,
                    vector<MLAPI::Operator>& A,
                    vector<RCP<MLAPI::InverseOperator> >& S,
                    vector<MLAPI::Operator>& P, 
                    vector<MLAPI::Operator>& R,
                    const int level,
                    const int nlevel,
                    MLAPI::MultiVector& x,
                    const MLAPI::MultiVector& f,
                    bool initiguesszero,
                    bool analysis,
                    bool silent) const
{
#if (FSIAMG_ANALYSIS>=3)
  analysis = true;
  silent = false;
#endif

  MLAPI::MultiVector r(f.GetVectorSpace(),1,false);
  if (initiguesszero) r = f;
  else                r = f - A[level] * x;

  double initrl2 = 0.0;
  double initrinf = 0.0;
  if (analysis) 
  {
    initrl2 = r.Norm2();
    initrinf = r.NormInf();
  }
  
  MLAPI::MultiVector tmpx(x.GetVectorSpace(),1,false);
  
  Epetra_Time timer(MLAPI::GetEpetra_Comm());
  
  for (int i=1; i<=sweeps; ++i)
  {
    tmpx = 0.0;
    Vcycle(field,myrank,levelsweeps,leveldamps,level,nlevel,tmpx,r,A,S,P,R);
    x.Update(damp,tmpx,1.0);
    //x = x + damp * tmpx;
    if (i<sweeps || analysis) r = f - A[level] * x;
  }

  if (analysis)
  {
    double t = timer.ElapsedTime();
    double rl2 = r.Norm2();
    double rinf = r.NormInf();
    double rl2rate = Rate(myrank,t,rl2,initrl2,r.GetGlobalLength());
    double rinfrate = Rate(myrank,t,rinf,initrinf,r.GetGlobalLength());
    if (!myrank && !silent) printf("RichardsonV %s      (level %2d) r0 %10.5e rl_2 %10.5e t %10.5e damp %10.5e sweeps %2d l2rate %10.5e\n",
                                   field.c_str(),level,initrl2,rl2,t,damp,sweeps,rl2rate);
    rl2rate = max(rl2rate,rinfrate);
    MLAPI::GetEpetra_Comm().Broadcast(&rl2rate,1,0);
    return rl2rate;
  }
  else return 0.0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::Vcycle(
                const string field,
                const int myrank,
                vector<int>& sweeps,
                vector<double>& damps,
                const int level,
                const int nlevel,
                MLAPI::MultiVector& z,
                const MLAPI::MultiVector& b,
                const vector<MLAPI::Operator>& A,
                const vector<RCP<MLAPI::InverseOperator> >& S,
                const vector<MLAPI::Operator>& P,
                const vector<MLAPI::Operator>& R) const
{
  // coarse solve
  if (level==nlevel-1)
  {
    z = 0.0;
    RichardsonS(field,myrank,level,sweeps[level],damps[level],A[level],*(S[level]),z,b,true,false,true);
    return;
  }

  // presmoothing (initial guess = 0)
  z = 0.0;
  RichardsonS(field,myrank,level,sweeps[level],damps[level],A[level],*(S[level]),z,b,true,false,true);

  // coarse level residual and correction
  MLAPI::MultiVector bc;
  MLAPI::MultiVector zc(P[level].GetDomainSpace(),1,false);
  bc = R[level] * ( b - A[level] * z );

  // solve coarse problem
  Vcycle(field,myrank,sweeps,damps,level+1,nlevel,zc,bc,A,S,P,R);
      
  // prolongate correction
  z = z + P[level] * zc;
  
  // postsmoothing (initial guess != 0 !!)
  MLAPI::MultiVector r(b.GetVectorSpace(),1,false);
  MLAPI::MultiVector dz(b.GetVectorSpace(),1,true);
  r = A[level] * z;
  r.Update(1.0,b,-1.0);
  dz = 0.0;
  RichardsonS(field,myrank,level,sweeps[level],damps[level],A[level],*(S[level]),dz,r,true,false,true);
  z.Update(1.0,dz,1.0);
  
  return;
}





/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double FSI::OverlappingBlockMatrixFSIAMG::RichardsonBGS_V(
                          const int myrank,
                          const int sweeps,
                          const double damp,
                          vector<int>& blocksweeps,
                          vector<double>& blockdamps,
                          AnalyzeBest& sbest,
                          AnalyzeBest& fbest,
                          AnalyzeBest& abest,
                          MLAPI::MultiVector& sy,
                          MLAPI::MultiVector& fy,
                          MLAPI::MultiVector& ay,
                          const MLAPI::MultiVector& sf,
                          const MLAPI::MultiVector& ff,
                          const MLAPI::MultiVector& af,  
                          vector<MLAPI::Operator>& Ass,
                          vector<MLAPI::Operator>& Pss, vector<MLAPI::Operator>& Rss,
                          vector<MLAPI::Operator>& Aff,
                          vector<MLAPI::Operator>& Pff, vector<MLAPI::Operator>& Rff,
                          vector<MLAPI::Operator>& Aaa,
                          vector<MLAPI::Operator>& Paa, vector<MLAPI::Operator>& Raa,
                          vector<MLAPI::Operator>& Asf,
                          vector<MLAPI::Operator>& Afs,
                          vector<MLAPI::Operator>& Afa,
                          vector<MLAPI::Operator>& Aaf,
                          bool initiguesszero,
                          bool analysis,
                          bool silent) const
{
#if (FSIAMG_ANALYSIS>=1)
  analysis = true;
  silent = false;
#endif

  MLAPI::MultiVector stmpf;
  MLAPI::MultiVector ftmpf;
  MLAPI::MultiVector atmpf;
  MLAPI::MultiVector sz(sy.GetVectorSpace(),1,false);
  MLAPI::MultiVector fz(fy.GetVectorSpace(),1,false);
  MLAPI::MultiVector az(ay.GetVectorSpace(),1,false);

  double initsrl2  = 0.0;
  double initsrinf = 0.0;
  double initarl2  = 0.0;
  double initarinf = 0.0;
  double initfrl2  = 0.0;
  double initfrinf = 0.0;
  double initrl2 = 0.0;
  double initrinf = 0.0;
  double rl2  = 0.0;
  double rinf = 0.0;

  if (analysis)
  {
    stmpf = sf - Ass[0] * sy;
    stmpf = stmpf - Asf[0] * fy;
    atmpf = af - Aaa[0] * ay;
    if (structuresplit_) atmpf = atmpf - Aaf[0] * fy;
    else                 atmpf = atmpf - Aaf[0] * sy;
    ftmpf = ff - Aff[0] * fy;
    ftmpf = ftmpf - Afs[0] * sy;
    ftmpf = ftmpf - Afa[0] * ay;
    initsrl2  = stmpf.DotProduct(stmpf);
    initsrinf = stmpf.NormInf();
    initarl2  = atmpf.DotProduct(atmpf);
    initarinf = atmpf.NormInf();
    initfrl2  = ftmpf.DotProduct(ftmpf);
    initfrinf = ftmpf.NormInf();
    initrl2 = sqrt(initsrl2 + initarl2 + initfrl2);
    initrinf = max(initsrinf,initarinf);
    initrinf = max(initrinf,initfrinf);
  }
  Epetra_Time timer(MLAPI::GetEpetra_Comm());
  double t1=0.0;
  double t2=0.0;
  double t3=0.0;
  
  for (int i=1; i<=sweeps; ++i)
  {
    //--------------------- structure block
    {
      if (analysis) timer.ResetStartTime();
      // compute ( f - A x ) for structure row
      stmpf = sf - Ass[0] * sy;
      stmpf = stmpf - Asf[0] * fy;
      // zero initial guess
      sz = 0.0;
      RichardsonV("(s)",myrank,blocksweeps[0],blockdamps[0],sbest.Sweeps(),sbest.Damp(),
                  Ass,sbest.S(),Pss,Rss,0,sbest.Nlevel(),sz,stmpf,true,false,true);
      sy.Update(damp,sz,1.0); 
      if (analysis) t1 += timer.ElapsedTime();
    }
    //---------------------- ale block
    {
      if (analysis) timer.ResetStartTime();
      // compute ( f - A x ) for ale row
      atmpf = af - Aaa[0] * ay;
      if (structuresplit_) atmpf = atmpf - Aaf[0] * fy;
      else                 atmpf = atmpf - Aaf[0] * sy;
      // zero initial guess
      az = 0.0;
      RichardsonV("(a)",myrank,blocksweeps[2],blockdamps[2],abest.Sweeps(),abest.Damp(),
                  Aaa,abest.S(),Paa,Raa,0,abest.Nlevel(),az,atmpf,true,false,true);
      ay.Update(damp,az,1.0); 
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
      RichardsonV("(f)",myrank,blocksweeps[1],blockdamps[1],fbest.Sweeps(),fbest.Damp(),
                  Aff,fbest.S(),Pff,Rff,0,fbest.Nlevel(),fz,ftmpf,true,false,true);
      fy.Update(damp,fz,1.0); 
      if (analysis) t3 += timer.ElapsedTime();
    }
  } // iterations

  if (analysis)
  {
    double t = timer.ElapsedTime();

    // final residuals
    stmpf = sf - Ass[0] * sy;
    stmpf = stmpf - Asf[0] * fy;
    double srl2  = stmpf.DotProduct(stmpf);
    double srinf = stmpf.NormInf();
    
    atmpf = af - Aaa[0] * ay;
    if (structuresplit_) atmpf = atmpf - Aaf[0] * fy;
    else                 atmpf = atmpf - Aaf[0] * sy;
    double arl2  = atmpf.DotProduct(atmpf);
    double arinf = atmpf.NormInf();
    
    ftmpf = ff - Aff[0] * fy;
    ftmpf = ftmpf - Afs[0] * sy;
    ftmpf = ftmpf - Afa[0] * ay;
    double frl2  = ftmpf.DotProduct(ftmpf);
    double frinf = ftmpf.NormInf();

    rl2 = sqrt(srl2 + arl2 + frl2);
    rinf = max(srinf,arinf);
    rinf = max(rinf,frinf);

    double rl2rate =  Rate(myrank,t,rl2,initrl2,stmpf.GetGlobalLength()+atmpf.GetGlobalLength()+ftmpf.GetGlobalLength());
    double rinfrate = Rate(myrank,t,rinf,initrinf,stmpf.GetGlobalLength()+atmpf.GetGlobalLength()+ftmpf.GetGlobalLength());

    if (!myrank && !silent) printf("RichardsonBGS_V      (level %2d) r0 %10.5e rl_2 %10.5e rinf0 %10.5e rinf %10.5e t %10.5e damp %10.5e sweeps %2d l2rate %10.5e rinfrate %10.5e\n",
                                   0,initrl2,rl2,initrinf,rinf,t,damp,sweeps,rl2rate,rinfrate);
    rl2rate = max(rl2rate,rinfrate);

    MLAPI::GetEpetra_Comm().Broadcast(&rl2rate,1,0);
    return rl2rate;
  }
  else return 0.0;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double FSI::OverlappingBlockMatrixFSIAMG::RichardsonBGS_Mixed(
                          const int myrank,
                          const int sweeps,
                          const double damp,
                          vector<int>& blocksweeps,
                          vector<double>& blockdamps,
                          const bool sisamg,
                          const bool fisamg,
                          const bool aisamg,
                          AnalyzeBest& sbest,
                          AnalyzeBest& fbest,
                          AnalyzeBest& abest,
                          MLAPI::MultiVector& sy,
                          MLAPI::MultiVector& fy,
                          MLAPI::MultiVector& ay,
                          const MLAPI::MultiVector& sf,
                          const MLAPI::MultiVector& ff,
                          const MLAPI::MultiVector& af,  
                          vector<MLAPI::Operator>& Ass,
                          vector<MLAPI::Operator>& Pss, vector<MLAPI::Operator>& Rss,
                          vector<MLAPI::Operator>& Aff,
                          vector<MLAPI::Operator>& Pff, vector<MLAPI::Operator>& Rff,
                          vector<MLAPI::Operator>& Aaa,
                          vector<MLAPI::Operator>& Paa, vector<MLAPI::Operator>& Raa,
                          vector<MLAPI::Operator>& Asf,
                          vector<MLAPI::Operator>& Afs,
                          vector<MLAPI::Operator>& Afa,
                          vector<MLAPI::Operator>& Aaf,
                          bool initiguesszero,
                          bool analysis,
                          bool silent) const
{
#if (FSIAMG_ANALYSIS>=1)
  analysis = true;
  silent = false;
#endif

  MLAPI::MultiVector stmpf;
  MLAPI::MultiVector ftmpf;
  MLAPI::MultiVector atmpf;
  MLAPI::MultiVector sz(sy.GetVectorSpace(),1,false);
  MLAPI::MultiVector fz(fy.GetVectorSpace(),1,false);
  MLAPI::MultiVector az(ay.GetVectorSpace(),1,false);

  double initsrl2  = 0.0;
  double initsrinf = 0.0;
  double initarl2  = 0.0;
  double initarinf = 0.0;
  double initfrl2  = 0.0;
  double initfrinf = 0.0;
  double initrl2 = 0.0;
  double initrinf = 0.0;
  double rl2  = 0.0;
  double rinf = 0.0;

  if (analysis)
  {
    stmpf = sf - Ass[0] * sy;
    stmpf = stmpf - Asf[0] * fy;
    atmpf = af - Aaa[0] * ay;
    if (structuresplit_) atmpf = atmpf - Aaf[0] * fy;
    else                 atmpf = atmpf - Aaf[0] * sy;
    ftmpf = ff - Aff[0] * fy;
    ftmpf = ftmpf - Afs[0] * sy;
    ftmpf = ftmpf - Afa[0] * ay;
    initsrl2  = stmpf.DotProduct(stmpf);
    initsrinf = stmpf.NormInf();
    initarl2  = atmpf.DotProduct(atmpf);
    initarinf = atmpf.NormInf();
    initfrl2  = ftmpf.DotProduct(ftmpf);
    initfrinf = ftmpf.NormInf();
    initrl2 = sqrt(initsrl2 + initarl2 + initfrl2);
    initrinf = max(initsrinf,initarinf);
    initrinf = max(initrinf,initfrinf);
  }
  Epetra_Time timer(MLAPI::GetEpetra_Comm());
  double t1=0.0;
  double t2=0.0;
  double t3=0.0;
  
  for (int i=1; i<=sweeps; ++i)
  {
    //--------------------- structure block
    {
      if (analysis) timer.ResetStartTime();
      // compute ( f - A x ) for structure row
      stmpf = sf - Ass[0] * sy;
      stmpf = stmpf - Asf[0] * fy;
      // zero initial guess
      sz = 0.0;

      if (sisamg)
        RichardsonV("(s)",myrank,blocksweeps[0],blockdamps[0],sbest.Sweeps(),sbest.Damp(),
                    Ass,sbest.S(),Pss,Rss,0,sbest.Nlevel(),sz,stmpf,true,false,true);
      else
        RichardsonMixed("(s)",myrank,0,blocksweeps[0],blockdamps[0],Ass[0],
                        Matrix(0,0),structuresolver_,sz,stmpf,const_cast<int&>(srun_),true,false,true);

      sy.Update(damp,sz,1.0); 
      if (analysis) t1 += timer.ElapsedTime();
    }
    //---------------------- ale block
    {
      if (analysis) timer.ResetStartTime();
      // compute ( f - A x ) for ale row
      atmpf = af - Aaa[0] * ay;
      if (structuresplit_) atmpf = atmpf - Aaf[0] * fy;
      else                 atmpf = atmpf - Aaf[0] * sy;
      // zero initial guess
      az = 0.0;

      if (aisamg)
        RichardsonV("(a)",myrank,blocksweeps[2],blockdamps[2],abest.Sweeps(),abest.Damp(),
                    Aaa,abest.S(),Paa,Raa,0,abest.Nlevel(),az,atmpf,true,false,true);
      else
        RichardsonMixed("(a)",myrank,0,blocksweeps[2],blockdamps[2],Aaa[0],
                        Matrix(2,2),alesolver_,az,atmpf,const_cast<int&>(arun_),true,false,true);
      
      ay.Update(damp,az,1.0); 
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
      
      if (fisamg)
        RichardsonV("(f)",myrank,blocksweeps[1],blockdamps[1],fbest.Sweeps(),fbest.Damp(),
                    Aff,fbest.S(),Pff,Rff,0,fbest.Nlevel(),fz,ftmpf,true,false,true);
      else
        RichardsonMixed("(f)",myrank,0,blocksweeps[1],blockdamps[1],Aff[0],
                        Matrix(1,1),fluidsolver_,fz,ftmpf,const_cast<int&>(frun_),true,false,true);

      fy.Update(damp,fz,1.0); 
      if (analysis) t3 += timer.ElapsedTime();
    }
  } // iterations

  if (analysis)
  {
    double t = timer.ElapsedTime();

    // final residuals
    stmpf = sf - Ass[0] * sy;
    stmpf = stmpf - Asf[0] * fy;
    double srl2  = stmpf.DotProduct(stmpf);
    double srinf = stmpf.NormInf();
    
    atmpf = af - Aaa[0] * ay;
    if (structuresplit_) atmpf = atmpf - Aaf[0] * fy;
    else                 atmpf = atmpf - Aaf[0] * sy;
    double arl2  = atmpf.DotProduct(atmpf);
    double arinf = atmpf.NormInf();
    
    ftmpf = ff - Aff[0] * fy;
    ftmpf = ftmpf - Afs[0] * sy;
    ftmpf = ftmpf - Afa[0] * ay;
    double frl2  = ftmpf.DotProduct(ftmpf);
    double frinf = ftmpf.NormInf();

    rl2 = sqrt(srl2 + arl2 + frl2);
    rinf = max(srinf,arinf);
    rinf = max(rinf,frinf);

    double rl2rate =  Rate(myrank,t,rl2,initrl2,stmpf.GetGlobalLength()+atmpf.GetGlobalLength()+ftmpf.GetGlobalLength());
    double rinfrate = Rate(myrank,t,rinf,initrinf,stmpf.GetGlobalLength()+atmpf.GetGlobalLength()+ftmpf.GetGlobalLength());

    if (!myrank && !silent) printf("RichardsonBGS_Mixed  (level %2d) r0 %10.5e rl_2 %10.5e rinf0 %10.5e rinf %10.5e t %10.5e damp %10.5e sweeps %2d l2rate %10.5e rinfrate %10.5e\n",
                                   0,initrl2,rl2,initrinf,rinf,t,damp,sweeps,rl2rate,rinfrate);
    rl2rate = max(rl2rate,rinfrate);

    MLAPI::GetEpetra_Comm().Broadcast(&rl2rate,1,0);
    return rl2rate;
  }
  else return 0.0;
}






/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double FSI::OverlappingBlockMatrixFSIAMG::RichardsonBGS_SV(
                            const int myrank,
                            const int sweeps,
                            const double damp,
                            vector<int>* blocksweeps,
                            vector<double>* blockdamps,
                            const int level,
                            AnalyzeBest& sbest,
                            AnalyzeBest& fbest,
                            AnalyzeBest& abest,
                            MLAPI::MultiVector& sy,
                            MLAPI::MultiVector& fy,
                            MLAPI::MultiVector& ay,
                            const MLAPI::MultiVector& sf,
                            const MLAPI::MultiVector& ff,
                            const MLAPI::MultiVector& af,  
                            vector<MLAPI::Operator>& Ass,
                            vector<MLAPI::Operator>& Pss, vector<MLAPI::Operator>& Rss,
                            vector<MLAPI::Operator>& Aff,
                            vector<MLAPI::Operator>& Pff, vector<MLAPI::Operator>& Rff,
                            vector<MLAPI::Operator>& Aaa,
                            vector<MLAPI::Operator>& Paa, vector<MLAPI::Operator>& Raa,
                            MLAPI::Operator& Asf,
                            MLAPI::Operator& Afs,
                            MLAPI::Operator& Afa,
                            MLAPI::Operator& Aaf,
                            bool initiguesszero,
                            bool analysis,
                            bool silent) const
{
#if (FSIAMG_ANALYSIS>=2)
  analysis = true;
  silent = false;
#endif

  MLAPI::MultiVector stmpf;
  MLAPI::MultiVector ftmpf;
  MLAPI::MultiVector atmpf;
  MLAPI::MultiVector sz(sy.GetVectorSpace(),1,false);
  MLAPI::MultiVector fz(fy.GetVectorSpace(),1,false);
  MLAPI::MultiVector az(ay.GetVectorSpace(),1,false);

  double initsrl2  = 0.0;
  double initsrinf = 0.0;
  double initarl2  = 0.0;
  double initarinf = 0.0;
  double initfrl2  = 0.0;
  double initfrinf = 0.0;
  double initrl2 = 0.0;
  double initrinf = 0.0;
  double rl2  = 0.0;
  double rinf = 0.0;

  if (analysis)
  {
    stmpf = sf - Ass[level] * sy;
    stmpf = stmpf - Asf * fy;
    atmpf = af - Aaa[level] * ay;
    if (structuresplit_) atmpf = atmpf - Aaf * fy;
    else                 atmpf = atmpf - Aaf * sy;
    ftmpf = ff - Aff[level] * fy;
    ftmpf = ftmpf - Afs * sy;
    ftmpf = ftmpf - Afa * ay;
    initsrl2  = stmpf.DotProduct(stmpf);
    initsrinf = stmpf.NormInf();
    initarl2  = atmpf.DotProduct(atmpf);
    initarinf = atmpf.NormInf();
    initfrl2  = ftmpf.DotProduct(ftmpf);
    initfrinf = ftmpf.NormInf();
    initrl2 = sqrt(initsrl2 + initarl2 + initfrl2);
    initrinf = max(initsrinf,initarinf);
    initrinf = max(initrinf,initfrinf);
  }
  Epetra_Time timer(MLAPI::GetEpetra_Comm());
  double t1=0.0;
  double t2=0.0;
  double t3=0.0;

  for (int i=1; i<=sweeps; ++i)
  {
    //--------------------- structure block
    {
      if (analysis) timer.ResetStartTime();
      // compute ( f - A x ) for structure row
      stmpf = sf - Ass[level] * sy;
      stmpf = stmpf - Asf * fy;
      // zero initial guess
      sz = 0.0;

      if (level<minnlevel_-1)
        RichardsonS("(s)",myrank,level,blocksweeps[0][level],blockdamps[0][level],
                    Ass[level],*(sbest.S()[level]),sz,stmpf,true,false,true);
      else
        RichardsonV("(s)",myrank,blocksweeps[0][level],blockdamps[0][level],
                    sbest.Sweeps(),sbest.Damp(),Ass,sbest.S(),Pss,Rss,level,sbest.Nlevel(),
                    sz,stmpf,true,false,true);

      sy.Update(damp,sz,1.0); 
      if (analysis) t1 += timer.ElapsedTime();
    }
    //---------------------- ale block
    {
      if (analysis) timer.ResetStartTime();
      // compute ( f - A x ) for ale row
      atmpf = af - Aaa[level] * ay;
      if (structuresplit_) atmpf = atmpf - Aaf * fy;
      else                 atmpf = atmpf - Aaf * sy;
      // zero initial guess
      az = 0.0;

      if (level<minnlevel_-1)
        RichardsonS("(a)",myrank,level,blocksweeps[2][level],blockdamps[2][level],
                    Aaa[level],*(abest.S()[level]),az,atmpf,true,false,true);
      else
        RichardsonV("(a)",myrank,blocksweeps[2][level],blockdamps[2][level],
                    abest.Sweeps(),abest.Damp(),Aaa,abest.S(),Paa,Raa,level,abest.Nlevel(),
                    az,atmpf,true,false,true);

      ay.Update(damp,az,1.0); 
      if (analysis) t2 += timer.ElapsedTime();
    }
    //------------------------ fluid block
    {
      if (analysis) timer.ResetStartTime();
      // compute ( f - A x ) for fluid row
      ftmpf = ff - Aff[level] * fy;
      ftmpf = ftmpf - Afs * sy;
      ftmpf = ftmpf - Afa * ay;
      // zero initial guess
      fz = 0.0;

      if (level<minnlevel_-1)
        RichardsonS("(f)",myrank,level,blocksweeps[1][level],blockdamps[1][level],
                    Aff[level],*(fbest.S()[level]),fz,ftmpf,true,false,true);
      else
        RichardsonV("(f)",myrank,blocksweeps[1][level],blockdamps[1][level],
                    fbest.Sweeps(),fbest.Damp(),Aff,fbest.S(),Pff,Rff,level,fbest.Nlevel(),
                    fz,ftmpf,true,false,true);

      fy.Update(damp,fz,1.0); 
      if (analysis) t3 += timer.ElapsedTime();
    }
  } // iterations

  if (analysis)
  {
    double t = timer.ElapsedTime();

    // final residuals
    stmpf = sf - Ass[level] * sy;
    stmpf = stmpf - Asf * fy;
    double srl2  = stmpf.DotProduct(stmpf);
    double srinf = stmpf.NormInf();
    
    atmpf = af - Aaa[level] * ay;
    if (structuresplit_) atmpf = atmpf - Aaf * fy;
    else                 atmpf = atmpf - Aaf * sy;
    double arl2  = atmpf.DotProduct(atmpf);
    double arinf = atmpf.NormInf();
    
    ftmpf = ff - Aff[level] * fy;
    ftmpf = ftmpf - Afs * sy;
    ftmpf = ftmpf - Afa * ay;
    double frl2  = ftmpf.DotProduct(ftmpf);
    double frinf = ftmpf.NormInf();

    rl2 = sqrt(srl2 + arl2 + frl2);
    rinf = max(srinf,arinf);
    rinf = max(rinf,frinf);

    double rl2rate =  Rate(myrank,t,rl2,initrl2,stmpf.GetGlobalLength()+atmpf.GetGlobalLength()+ftmpf.GetGlobalLength());
    double rinfrate = Rate(myrank,t,rinf,initrinf,stmpf.GetGlobalLength()+atmpf.GetGlobalLength()+ftmpf.GetGlobalLength());

    if (!myrank && !silent) printf("RichardsonBGS_SV     (level %2d) r0 %10.5e rl_2 %10.5e rinf0 %10.5e rinf %10.5e t %10.5e damp %10.5e sweeps %2d l2rate %10.5e rinfrate %10.5e\n",
                                   level,initrl2,rl2,initrinf,rinf,t,damp,sweeps,rl2rate,rinfrate);
    rl2rate = max(rl2rate,rinfrate);

    MLAPI::GetEpetra_Comm().Broadcast(&rl2rate,1,0);
    return rl2rate;
  }
  else return 0.0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double FSI::OverlappingBlockMatrixFSIAMG::Richardson_BlockV(
                             const int myrank,
                             const int sweeps,
                             const double damp,
                             vector<int>& Vsweeps,
                             vector<double>& Vdamps,
                             vector<int>* blocksweeps,
                             vector<double>* blockdamps,
                             AnalyzeBest& sbest,
                             AnalyzeBest& fbest,
                             AnalyzeBest& abest,
                             MLAPI::MultiVector& sy,
                             MLAPI::MultiVector& fy,
                             MLAPI::MultiVector& ay,
                             const MLAPI::MultiVector& sf,
                             const MLAPI::MultiVector& ff,
                             const MLAPI::MultiVector& af,  
                             vector<MLAPI::Operator>& Ass,
                             vector<MLAPI::Operator>& Pss, vector<MLAPI::Operator>& Rss,
                             vector<MLAPI::Operator>& Aff,
                             vector<MLAPI::Operator>& Pff, vector<MLAPI::Operator>& Rff,
                             vector<MLAPI::Operator>& Aaa,
                             vector<MLAPI::Operator>& Paa, vector<MLAPI::Operator>& Raa,
                             vector<MLAPI::Operator>& Asf,
                             vector<MLAPI::Operator>& Afs,
                             vector<MLAPI::Operator>& Afa,
                             vector<MLAPI::Operator>& Aaf,
                             bool initiguesszero,
                             bool analysis,
                             bool silent) const
{
#if (FSIAMG_ANALYSIS>=1)
  analysis = true;
  silent = false;
#endif

  MLAPI::MultiVector stmpf;
  MLAPI::MultiVector ftmpf;
  MLAPI::MultiVector atmpf;
  MLAPI::MultiVector sz(sy.GetVectorSpace(),1,false);
  MLAPI::MultiVector fz(fy.GetVectorSpace(),1,false);
  MLAPI::MultiVector az(ay.GetVectorSpace(),1,false);

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
    if (structuresplit_) atmpf = atmpf - Aaf[0] * fy;
    else                 atmpf = atmpf - Aaf[0] * sy;
    ftmpf = ff - Aff[0] * fy;
    ftmpf = ftmpf - Afs[0] * sy;
    ftmpf = ftmpf - Afa[0] * ay;
  }
  double initsrl2  = 0.0;
  double initsrinf = 0.0;
  double initarl2  = 0.0;
  double initarinf = 0.0;
  double initfrl2  = 0.0;
  double initfrinf = 0.0;
  double initrl2 = 0.0;
  double initrinf = 0.0;
  
  if (analysis)
  {
    initsrl2  = stmpf.DotProduct(stmpf);
    initsrinf = stmpf.NormInf();
    initarl2  = atmpf.DotProduct(atmpf);
    initarinf = atmpf.NormInf();
    initfrl2  = ftmpf.DotProduct(ftmpf);
    initfrinf = ftmpf.NormInf();
    initrl2 = sqrt(initsrl2 + initarl2 + initfrl2); 
    initrinf = max(initsrinf,initarinf);
    initrinf = max(initrinf,initfrinf);
  }
  
  Epetra_Time timer(MLAPI::GetEpetra_Comm());
  
  for (int i=1; i<=sweeps; ++i)
  {
    sz = 0.0;
    az = 0.0;
    fz = 0.0;
    
    BlockVcycle(myrank,Vsweeps,Vdamps,blocksweeps,blockdamps,0,minnlevel_,
                sbest,fbest,abest,sz,fz,az,stmpf,ftmpf,atmpf,
                Ass,Pss,Rss,Aff,Pff,Rff,Aaa,Paa,Raa,Asf,Afs,Afa,Aaf);
    
    sy.Update(damp,sz,1.0);
    ay.Update(damp,az,1.0);
    fy.Update(damp,fz,1.0);
    
    if (i<sweeps || analysis)
    {
      stmpf = sf - Ass[0] * sy; 
      stmpf = stmpf - Asf[0] * fy;
      atmpf = af - Aaa[0] * ay;
      if (structuresplit_) atmpf = atmpf - Aaf[0] * fy;
      else                 atmpf = atmpf - Aaf[0] * sy;
      ftmpf = ff - Aff[0] * fy;
      ftmpf = ftmpf - Afs[0] * sy;
      ftmpf = ftmpf - Afa[0] * ay;
    }
  } // iterations    
    
  if (analysis)
  {
    double t = timer.ElapsedTime();

    double srl2  = stmpf.DotProduct(stmpf);
    double srinf = stmpf.NormInf();
    double arl2  = atmpf.DotProduct(atmpf);
    double arinf = atmpf.NormInf();
    double frl2  = ftmpf.DotProduct(ftmpf);
    double frinf = ftmpf.NormInf();
    double rl2 = sqrt(srl2 + arl2 + frl2);
    double rinf = max(srinf,arinf);
    rinf = max(rinf,frinf);

    double rl2rate =  Rate(myrank,t,rl2,initrl2,stmpf.GetGlobalLength()+atmpf.GetGlobalLength()+ftmpf.GetGlobalLength());
    double rinfrate = Rate(myrank,t,rinf,initrinf,stmpf.GetGlobalLength()+atmpf.GetGlobalLength()+ftmpf.GetGlobalLength());

    if (!myrank && !silent) printf("Richardson_BlockV    (level %2d) r0 %10.5e rl_2 %10.5e rinf0 %10.5e rinf %10.5e t %10.5e damp %10.5e sweeps %2d l2rate %10.5e rinfrate %10.5e\n",
                                   0,initrl2,rl2,initrinf,rinf,t,damp,sweeps,rl2rate,rinfrate);
    rl2rate = max(rl2rate,rinfrate);

    MLAPI::GetEpetra_Comm().Broadcast(&rl2rate,1,0);
    return rl2rate;
  }
  else return 0.0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::BlockVcycle(
                     const int myrank,
                     vector<int>& Vsweeps,
                     vector<double>& Vdamps,
                     vector<int>* blocksweeps,
                     vector<double>* blockdamps,
                     const int level,
                     const int minnlevel,
                     AnalyzeBest& sbest,
                     AnalyzeBest& fbest,
                     AnalyzeBest& abest,
                     MLAPI::MultiVector& sy,
                     MLAPI::MultiVector& fy,
                     MLAPI::MultiVector& ay,
                     const MLAPI::MultiVector& sf,
                     const MLAPI::MultiVector& ff,
                     const MLAPI::MultiVector& af,  
                     vector<MLAPI::Operator>& Ass,
                     vector<MLAPI::Operator>& Pss, vector<MLAPI::Operator>& Rss,
                     vector<MLAPI::Operator>& Aff,
                     vector<MLAPI::Operator>& Pff, vector<MLAPI::Operator>& Rff,
                     vector<MLAPI::Operator>& Aaa,
                     vector<MLAPI::Operator>& Paa, vector<MLAPI::Operator>& Raa,
                     vector<MLAPI::Operator>& Asf,
                     vector<MLAPI::Operator>& Afs,
                     vector<MLAPI::Operator>& Afa,
                     vector<MLAPI::Operator>& Aaf) const
{
  sy = 0.0;
  ay = 0.0;
  fy = 0.0;

  // coarsest level Richardson
  if (level==minnlevel-1)
  {
    RichardsonBGS_SV(myrank,Vsweeps[level],Vdamps[level],blocksweeps,blockdamps,
                     level,sbest,fbest,abest,sy,fy,ay,sf,ff,af,
                     Ass,Pss,Rss,Aff,Pff,Rff,Aaa,Paa,Raa,Asf[level],Afs[level],Afa[level],Aaf[level],
                     true,false,true);
    return;
  }
  
  // presmoothing
  RichardsonBGS_SV(myrank,Vsweeps[level],Vdamps[level],blocksweeps,blockdamps,
                   level,sbest,fbest,abest,sy,fy,ay,sf,ff,af,
                   Ass,Pss,Rss,Aff,Pff,Rff,Aaa,Paa,Raa,Asf[level],Afs[level],Afa[level],Aaf[level],
                   true,false,true);
  

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
    if (structuresplit_) tmp = tmp - Aaf[level] * fy;
    else                 tmp = tmp - Aaf[level] * sy;
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
  MLAPI::MultiVector syc(sfc.GetVectorSpace(),1,false);
  MLAPI::MultiVector ayc(afc.GetVectorSpace(),1,false);
  MLAPI::MultiVector fyc(ffc.GetVectorSpace(),1,false);
  
  // solve coarser problem
  BlockVcycle(myrank,Vsweeps,Vdamps,blocksweeps,blockdamps,level+1,minnlevel,sbest,fbest,abest,
              syc,fyc,ayc,sfc,ffc,afc,Ass,Pss,Rss,Aff,Pff,Rff,Aaa,Paa,Raa,Asf,Afs,Afa,Aaf);
              
  // prolongate correction
  {
    MLAPI::MultiVector tmp;
    tmp = Pss[level] * syc;
    sy.Update(1.0,tmp,1.0);
    tmp = Paa[level] * ayc;
    ay.Update(1.0,tmp,1.0);
    tmp = Pff[level] * fyc;
    fy.Update(1.0,tmp,1.0);
  }
  
  // postsmoothing (do NOT zero out initial guess!)
  RichardsonBGS_SV(myrank,Vsweeps[level],Vdamps[level],blocksweeps,blockdamps,
                   level,sbest,fbest,abest,sy,fy,ay,sf,ff,af,
                   Ass,Pss,Rss,Aff,Pff,Rff,Aaa,Paa,Raa,Asf[level],Afs[level],Afa[level],Aaf[level],
                   false,false,true);
  
  return;
}







#endif
