/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 | input of solver control variables  structure           m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inpctrsol(SOLVAR *solv)
{
INT           ierr;
char          buffer[50];

#ifdef MLIB_PACKAGE
MLVAR        *mlvar = NULL;
#endif

#ifdef AZTEC_PACKAGE
AZVAR        *azvar = NULL;
#endif

#ifdef HYPRE_PACKAGE
HYPREVARS    *hyprevars = NULL;
#endif

#ifdef PARSUPERLU_PACKAGE
PSUPERLUVARS *psuperluvars = NULL;
#endif

LAPACKVARS   *lapackvars = NULL;

#ifdef MUMPS_PACKAGE
MUMPSVARS    *mumpsvars = NULL;
#endif
COLSOLVARS   *colsolvars = NULL;

#ifdef MLPCG
MLPCGVARS    *mlpcgvars = NULL;
#endif

#ifdef DEBUG
dstrc_enter("inpctrsol");
#endif

solv->matrixtyp = matrix_none;
switch(solv->fieldtyp)
{
case structure:
if (frfind("-STRUCT SOLVER")==0) goto end;
break;
case fluid:
if (frfind("-FLUID SOLVER")==0) goto end;
break;
case ale:
if (frfind("-ALE SOLVER")==0) goto end;
break;
#ifdef D_TSI
case thermal:
if (frfind("-THERMAL SOLVER")==0) goto end;
break;
#endif /* D_TSI */
case pressure:
if (frfind("-PRESSURE SOLVER")==0) goto end;
break;
default:
  dserror("Unknown field typ in reading solver %d",solv->fieldtyp);
break;
}
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*------------------------------------ check typ of solver and allocate */
   frchar("SOLVER"    ,buffer,&ierr);
   if (ierr==1)
   {

      /*--------------------------------------------------- Amesos solver KLU */
      /* note that all Amesos solvers work only with Trilinos enabled */
      if (strncmp("Amesos_KLU_sym",buffer,14)==0)
      {
#ifndef TRILINOS_PACKAGE
         dserror("TRILINOS_PACKAGE package is not compiled in");
#else
         if (!genprob.usetrilinosalgebra)
           dserror("You have to use ALGEBRA Trilinos to use Amesos_KLU");
         solv->solvertyp = amesos_klu_sym;
#endif
      }
      if (strncmp("Amesos_KLU_nonsym",buffer,17)==0)
      {
#ifndef TRILINOS_PACKAGE
         dserror("TRILINOS_PACKAGE package is not compiled in");
#else
         if (!genprob.usetrilinosalgebra)
           dserror("You have to use ALGEBRA Trilinos to use Amesos_KLU");
         solv->solvertyp = amesos_klu_nonsym;
#endif
      }
      if (strncmp("Superlu",buffer,7)==0)
      {
#ifndef TRILINOS_PACKAGE
         dserror("TRILINOS_PACKAGE package is not compiled in");
#else
#ifndef PARALLEL
         dserror("Superlu can only be used with -DPARALLEL");
#endif
         if (!genprob.usetrilinosalgebra)
           dserror("You have to use ALGEBRA Trilinos to use Superlu");
         solv->solvertyp = superlu;
#endif
      }
      /*--------------------------------------------------- MLIB solver */
      if (strncmp("MLIB_D_SP",buffer,9)==0)
      {
#ifndef MLIB_PACKAGE
         dserror("MLIB package is not compiled in");
#else
         solv->solvertyp = mlib_d_sp;
         solv->mlvar = (MLVAR*)CCACALLOC(1,sizeof(MLVAR));
         mlvar = solv->mlvar;
#endif
      }
      /*---------------------------------------------------- VM3 solver */
      if (strncmp("vm3",buffer,3)==0)
      {
        solv->solvertyp = vm3;
        solv->azvar = (AZVAR*)CCACALLOC(1,sizeof(AZVAR));
        azvar = solv->azvar;
        azvar->azconv = AZ_noscaled;
      }
      /*-------------------------------------------------- Aztec solver */
      if (strncmp("Aztec_MSR",buffer,9)==0)
      {
#ifndef AZTEC_PACKAGE
         dserror("Aztec package is not compiled in");
#else
         solv->solvertyp = aztec_msr;
         solv->azvar = (AZVAR*)CCACALLOC(1,sizeof(AZVAR));
         azvar = solv->azvar;
	 azvar->azconv = AZ_noscaled; /* default value */
#endif
      }
      /*---------------------------------------- HYPRE solver BoomerAMG */
      if (strncmp("HYPRE_BoomerAMG",buffer,15)==0)
      {
#ifndef HYPRE_PACKAGE
         dserror("Hypre package is not compiled in");
#else
         solv->solvertyp = hypre_amg;
         solv->hyprevar = (HYPREVARS*)CCACALLOC(1,sizeof(HYPREVARS));
         hyprevars = solv->hyprevar;
#endif
      }
      /*---------------------------------------------- HYPRE solver PCG */
      if (strncmp("HYPRE_PCG",buffer,9)==0)
      {
#ifndef HYPRE_PACKAGE
         dserror("Hypre package is not compiled in");
#else
         solv->solvertyp = hypre_pcg;
         solv->hyprevar = (HYPREVARS*)CCACALLOC(1,sizeof(HYPREVARS));
         hyprevars = solv->hyprevar;
#endif
      }
      /*-------------------------------------------- HYPRE solver Gmres */
      if (strncmp("HYPRE_GMRES",buffer,11)==0)
      {
#ifndef HYPRE_PACKAGE
         dserror("Hypre package is not compiled in");
#else
         solv->solvertyp = hypre_gmres;
         solv->hyprevar = (HYPREVARS*)CCACALLOC(1,sizeof(HYPREVARS));
         hyprevars = solv->hyprevar;
#endif
      }
      /*----------------------------------------- HYPRE solver BiCGstab */
      if (strncmp("HYPRE_BiCGStab",buffer,14)==0)
      {
#ifndef HYPRE_PACKAGE
         dserror("Hypre package is not compiled in");
#else
         solv->solvertyp = hypre_bicgstab;
         solv->hyprevar = (HYPREVARS*)CCACALLOC(1,sizeof(HYPREVARS));
         hyprevars = solv->hyprevar;
#endif
      }
      /*------------------------------------------------ solver SuperLU */
      if (strncmp("ParSuperLU",buffer,10)==0)
      {
#ifndef PARSUPERLU_PACKAGE
         dserror("ParSuperLU package is not compiled in");
#else
         solv->solvertyp = parsuperlu;
         solv->psuperluvars = (PSUPERLUVARS*)CCACALLOC(1,sizeof(PSUPERLUVARS));
         psuperluvars = solv->psuperluvars;
#endif
      }
      /*------------------------------------------------- solver Lapack */
      if (strncmp("LAPACK_sym",buffer,10)==0)
      {
         solv->solvertyp = lapack_sym;
         solv->lapackvars = (LAPACKVARS*)CCACALLOC(1,sizeof(LAPACKVARS));
         lapackvars = solv->lapackvars;
      }
      /*------------------------------------------------- solver Lapack */
      if (strncmp("LAPACK_nonsym",buffer,13)==0)
      {
         solv->solvertyp = lapack_nonsym;
         solv->lapackvars = (LAPACKVARS*)CCACALLOC(1,sizeof(LAPACKVARS));
         lapackvars = solv->lapackvars;
      }
      /*-------------------------------------------------- solver Mumps */
      if (strncmp("MUMPS_sym",buffer,9)==0)
      {
#ifndef MUMPS_PACKAGE
         dserror("MUMPS package is not compiled in");
#else
         solv->solvertyp = mumps_sym;
         solv->mumpsvars = (MUMPSVARS*)CCACALLOC(1,sizeof(MUMPSVARS));
         mumpsvars = solv->mumpsvars;
#endif
      }
      /*-------------------------------------------------- solver Mumps */
      if (strncmp("MUMPS_nonsym",buffer,12)==0)
      {
#ifndef MUMPS_PACKAGE
         dserror("MUMPS package is not compiled in");
#else
         solv->solvertyp = mumps_nonsym;
         solv->mumpsvars = (MUMPSVARS*)CCACALLOC(1,sizeof(MUMPSVARS));
         mumpsvars = solv->mumpsvars;
#endif
      }
      /*------------------------------------------------ solver Umfpack */
      if (strncmp("UMFPACK",buffer,7)==0)
      {
#ifndef UMFPACK
         dserror("UMFPACK is not compiled in");
#else
         solv->solvertyp = umfpack;
#endif
      }
      /*------------------------------------------------- solver Colsol */
      if (strncmp("Colsol",buffer,6)==0)
      {
         solv->solvertyp = colsol_solver;
         solv->colsolvars = (COLSOLVARS*)CCACALLOC(1,sizeof(COLSOLVARS));
         colsolvars = solv->colsolvars;
      }
      /*------------------------------------------------- solver Spooles */
      if (strncmp("SPOOLES_nonsym",buffer,11)==0)
      {
#ifndef SPOOLES_PACKAGE
         dserror("SPOOLES package is not compiled in");
#else
         solv->solvertyp = SPOOLES_nonsym;
#endif
      }
      /*------------------------------------------------- solver Spooles */
      if (strncmp("SPOOLES_sym",buffer,11)==0)
      {
#ifndef SPOOLES_PACKAGE
         dserror("SPOOLES package is not compiled in");
#else
         solv->solvertyp = SPOOLES_sym;
#endif
      }
      /*--------------------------------------------------- solver MLPCG */
      if (strncmp("MLPCG",buffer,5)==0)
      {
#ifndef MLPCG
        dserror("MLPCG package is not compiled in");
#else
         solv->solvertyp = mlpcg;
         solv->mlpcgvars = (MLPCGVARS*)CCACALLOC(1,sizeof(MLPCGVARS));
         mlpcgvars = solv->mlpcgvars;
#endif
      }
   }/* end of (ierr==1) */
/*------------------------------------------- read solver specific data */
   switch (solv->solvertyp)
   {
#ifdef TRILINOS_PACKAGE
   case amesos_klu_sym:/*----------------------- read solver Amesos_KLU */
   case amesos_klu_nonsym:/*-------------------- read solver Amesos_KLU */
   case superlu:/*-------------------------- read solver Amesos_SuperLU */
   break;
#endif
#ifdef AZTEC_PACKAGE
   case vm3:/*----------------------------------------- read solver vm3 */
   case aztec_msr:/*--------------------------------- read solver aztec */
      frchar("AZSOLVE"   ,buffer,&ierr);
      if (ierr==1)
      {
         if (strncmp("CG",buffer,2)==0)       azvar->azsolvertyp = azsolv_CG;
         if (strncmp("GMRES",buffer,5)==0)    azvar->azsolvertyp = azsolv_GMRES;
         if (strncmp("CGS",buffer,3)==0)      azvar->azsolvertyp = azsolv_CGS;
         if (strncmp("BiCGSTAB",buffer,8)==0) azvar->azsolvertyp = azsolv_BiCGSTAB;
         if (strncmp("LU",buffer,2)==0)       azvar->azsolvertyp = azsolv_LU;
         if (strncmp("TFQMR",buffer,5)==0)    azvar->azsolvertyp = azsolv_TFQMR;
      }
      frchar("AZPREC"    ,buffer,&ierr);
      if (ierr==1)
      {
         if (strncmp("none",buffer,4)==0)             azvar->azprectyp = azprec_none;
         if (strncmp("LU",buffer,2)==0)               azvar->azprectyp = azprec_LU;
         if (strncmp("ILU",buffer,3)==0)              azvar->azprectyp = azprec_ILU;
         if (strncmp("ILUT",buffer,4)==0)             azvar->azprectyp = azprec_ILUT;
         if (strncmp("RILU",buffer,4)==0)             azvar->azprectyp = azprec_RILU;
         if (strncmp("BILU",buffer,4)==0)             azvar->azprectyp = azprec_BILU;
         if (strncmp("Jacobi",buffer,6)==0)           azvar->azprectyp = azprec_Jacobi;
         if (strncmp("SymmGaussSeidel",buffer,15)==0) azvar->azprectyp = azprec_SymmGaussSeidel;
         if (strncmp("Least_Squares",buffer,13)==0)   azvar->azprectyp = azprec_Least_Squares;
         if (strncmp("Neumann",buffer,7)==0)          azvar->azprectyp = azprec_Neumann;
         if (strncmp("ICC",buffer,3)==0)              azvar->azprectyp = azprec_ICC;
#ifdef TRILINOS_PACKAGE
         if (strncmp("ML",buffer,2)==0)               azvar->azprectyp = azprec_ML;
         if (strncmp("MLFLUID",buffer,7)==0)          azvar->azprectyp = azprec_MLfluid;
         if (strncmp("MLFLUID2",buffer,8)==0)         azvar->azprectyp = azprec_MLfluid2;
#endif
      }
      /* the corresponding numbers are defined in az_aztec_defs.h ... */
      frint("AZOUTPUT"  ,&(azvar->azoutput) ,&ierr);
      frint("AZREUSE"   ,&(azvar->azreuse)  ,&ierr);
      frint("AZGFILL"   ,&(azvar->azgfill)  ,&ierr);
      frint("AZITER"    ,&(azvar->aziter)   ,&ierr);
      frint("AZSUB"     ,&(azvar->azsub)    ,&ierr);
      frint("AZGRAPH"   ,&(azvar->azgraph)  ,&ierr);
      frint("AZPOLY"    ,&(azvar->azpoly)   ,&ierr);
      frint("AZBDIAG"   ,&(azvar->blockdiag),&ierr);
      frint("AZOVERLAP" ,&(azvar->azoverlap),&ierr);
      frdouble("AZDROP" ,&(azvar->azdrop)   ,&ierr);
      frdouble("AZFILL" ,&(azvar->azfill)   ,&ierr);
      frdouble("AZTOL"  ,&(azvar->aztol)    ,&ierr);
      frchar("AZCONV"    ,buffer,&ierr);
      if(ierr==1)
      {
             if (strncmp("AZ_r0"             ,buffer, 5)==0) azvar->azconv = AZ_r0;
        else if (strncmp("AZ_rhs"            ,buffer, 6)==0) azvar->azconv = AZ_rhs;
        else if (strncmp("AZ_Anorm"          ,buffer, 8)==0) azvar->azconv = AZ_Anorm;
	else if (strncmp("AZ_sol"            ,buffer, 6)==0) azvar->azconv = AZ_sol;
        else if (strncmp("AZ_weighted"       ,buffer,11)==0) azvar->azconv = AZ_weighted;
        else if (strncmp("AZ_expected_values",buffer,18)==0) azvar->azconv = AZ_expected_values;
        else if (strncmp("AZ_noscaled"       ,buffer,11)==0) azvar->azconv = AZ_noscaled;
#ifdef TRILINOS_PACKAGE
        else if (strncmp("AZTECOO_conv_test" ,buffer,17)==0) azvar->azconv = AZTECOO_conv_test;
        else if (strncmp("AZ_inf_noscaled"   ,buffer,15)==0) azvar->azconv = AZ_inf_noscaled;
#endif
        else
          dserror("unsupported AZCONV value");
      }
      frdouble("AZOMEGA",&(azvar->azomega)  ,&ierr);
#ifdef TRILINOS_PACKAGE
      frchar("AZSCAL"    ,buffer,&ierr);
      if (ierr==1)
      {
        if (strncmp("none"   ,buffer,4)==0) azvar->azscal = 0;
        if (strncmp("sym"    ,buffer,3)==0) azvar->azscal = 1;
        if (strncmp("infnorm",buffer,7)==0) azvar->azscal = 2;
      }
      /* parameters of ML preconditioner */
      frint("ML_PRINT"           ,&(azvar->mlprint),&ierr);
      frint("ML_MAXCOARSESIZE"   ,&(azvar->mlcsize),&ierr);
      frint("ML_MAXLEVEL"        ,&(azvar->mlmaxlevel),&ierr);
      frint("ML_AGG_SIZE"        ,&(azvar->mlaggsize),&ierr);
      frdouble("ML_DAMPFINE"     ,&(azvar->mldamp_fine),&ierr);
      frdouble("ML_DAMPMED"      ,&(azvar->mldamp_med),&ierr);
      frdouble("ML_DAMPCOARSE"   ,&(azvar->mldamp_coarse),&ierr);
      frdouble("ML_PROLONG_SMO"  ,&(azvar->mldamp_prolong),&ierr);
      frdouble("ML_PROLONG_THRES",&(azvar->ml_threshold),&ierr);
      if (azvar->mlmaxlevel)
      frint_n("ML_SMOTIMES",azvar->mlsmotimes,azvar->mlmaxlevel,&ierr);
      frchar("ML_COARSEN",buffer,&ierr);
      if (ierr==1)
      {
        if (strncmp("UC",buffer,2)==0)      azvar->mlcoarsentype = 0;
        if (strncmp("METIS",buffer,5)==0)   azvar->mlcoarsentype = 1;
        if (strncmp("VBMETIS",buffer,7)==0) azvar->mlcoarsentype = 2;
        if (strncmp("MIS",buffer,3)==0)     azvar->mlcoarsentype = 3;
      }
      frchar("ML_SMOOTHERFINE",buffer,&ierr);
      if (ierr==1)
      {
        if (strncmp("SGS",buffer,3)==0)       azvar->mlsmotype_fine = 0;
        if (strncmp("Jacobi",buffer,6)==0)    azvar->mlsmotype_fine = 1;
        if (strncmp("Chebychev",buffer,9)==0) azvar->mlsmotype_fine = 2;
        if (strncmp("MLS",buffer,3)==0)       azvar->mlsmotype_fine = 3;
        if (strncmp("ILU",buffer,3)==0)       azvar->mlsmotype_fine = 4;
        if (strncmp("KLU",buffer,3)==0)       azvar->mlsmotype_fine = 5;
        if (strncmp("Superlu",buffer,7)==0)   azvar->mlsmotype_fine = 6;
      }
      frchar("ML_SMOOTHERMED",buffer,&ierr);
      if (ierr==1)
      {
        if (strncmp("SGS",buffer,3)==0)       azvar->mlsmotype_med = 0;
        if (strncmp("Jacobi",buffer,6)==0)    azvar->mlsmotype_med = 1;
        if (strncmp("Chebychev",buffer,9)==0) azvar->mlsmotype_med = 2;
        if (strncmp("MLS",buffer,3)==0)       azvar->mlsmotype_med = 3;
        if (strncmp("ILU",buffer,3)==0)       azvar->mlsmotype_med = 4;
        if (strncmp("KLU",buffer,3)==0)       azvar->mlsmotype_med = 5;
        if (strncmp("Superlu",buffer,7)==0)   azvar->mlsmotype_med = 6;
      }
      frchar("ML_SMOOTHERCOARSE",buffer,&ierr);
      if (ierr==1)
      {
        if (strncmp("SGS",buffer,3)==0)       azvar->mlsmotype_coarse = 0;
        if (strncmp("Jacobi",buffer,6)==0)    azvar->mlsmotype_coarse = 1;
        if (strncmp("Chebychev",buffer,9)==0) azvar->mlsmotype_coarse = 2;
        if (strncmp("MLS",buffer,3)==0)       azvar->mlsmotype_coarse = 3;
        if (strncmp("ILU",buffer,3)==0)       azvar->mlsmotype_coarse = 4;
        if (strncmp("KLU",buffer,3)==0)       azvar->mlsmotype_coarse = 5;
        if (strncmp("Superlu",buffer,7)==0)   azvar->mlsmotype_coarse = 6;
      }
#endif
   break;
#endif

#ifdef HYPRE_PACKAGE
   case hypre_amg:/*----------------------------- read solver boomeramg */
      hyprevars->hypre_prectyp = hypreprec_none;

      frchar("HYPRE_IO",buffer,&ierr);
      if (ierr)
      {
       if (strncmp("full",buffer,4)==0)
          hyprevars->io=1;
       if (strncmp("none",buffer,4)==0)
          hyprevars->io=0;
      }

      frint("HYPRE_ITER"  ,&(hyprevars->maxiter),&ierr);
      frint("HYPRE_SFINE" ,&(hyprevars->sweep[0]),&ierr);
      frint("HYPRE_SDOWN" ,&(hyprevars->sweep[1]),&ierr);
      frint("HYPRE_SUP"   ,&(hyprevars->sweep[2]),&ierr);
      frint("HYPRE_SCOARS",&(hyprevars->sweep[3]),&ierr);

      frdouble("HYPRE_TOL"   ,&(hyprevars->tol),&ierr);
      frdouble("HYPRE_THREAS",&(hyprevars->threshold),&ierr);
   break;
   case hypre_pcg:/*----------------------------------- read solver CG */
      frchar("HYPRE_PREC",buffer,&ierr);
      if (ierr)
      {
         if (strncmp("none",buffer,4)==0)      hyprevars->hypre_prectyp = hypreprec_none;
         if (strncmp("Euclid",buffer,6)==0)    hyprevars->hypre_prectyp = hypreprec_euclid;
         if (strncmp("Parasails",buffer,9)==0) hyprevars->hypre_prectyp = hypreprec_parasails;
         if (strncmp("BoomerAMG",buffer,9)==0) hyprevars->hypre_prectyp = hypreprec_amg;
      }

      frchar("HYPRE_IO",buffer,&ierr);
      if (ierr)
      {
       if (strncmp("full",buffer,4)==0)
          hyprevars->io=1;
       if (strncmp("none",buffer,4)==0)
          hyprevars->io=0;
      }

      frchar("HYPRE_PARASY",buffer,&ierr);
      if (ierr)
      {
       if (strncmp("symm_pos_def",buffer,12)==0) hyprevars->parasymm = 1;
       if (strncmp("nonsym_def",buffer,10)==0)   hyprevars->parasymm = 2;
       if (strncmp("non_sym_indef",buffer,13)==0)hyprevars->parasymm = 0;
      }

      frint("HYPRE_ITER"   ,&(hyprevars->maxiter)  ,&ierr);
      frint("HYPRE_IFILL"  ,&(hyprevars->ifill)    ,&ierr);
      frint("HYPRE_PARALEV",&(hyprevars->paralevel),&ierr);
      frint("HYPRE_SFINE" ,&(hyprevars->sweep[0]),&ierr);
      frint("HYPRE_SDOWN" ,&(hyprevars->sweep[1]),&ierr);
      frint("HYPRE_SUP"   ,&(hyprevars->sweep[2]),&ierr);
      frint("HYPRE_SCOARS",&(hyprevars->sweep[3]),&ierr);

      frdouble("HYPRE_TOL"     ,&(hyprevars->tol)       ,&ierr);
      frdouble("HYPRE_DFILL"   ,&(hyprevars->dfill)     ,&ierr);
      frdouble("HYPRE_THREAS",&(hyprevars->threshold),&ierr);
      frdouble("HYPRE_PARATHR" ,&(hyprevars->parathresh),&ierr);
      frdouble("HYPRE_PARAFILT",&(hyprevars->parafilter),&ierr);
   break;
   case hypre_gmres:/*--------------------------------- read solver GMRES */

      frchar("HYPRE_PREC",buffer,&ierr);
      if (ierr)
      {
         if (strncmp("none",buffer,4)==0)      hyprevars->hypre_prectyp = hypreprec_none;
         if (strncmp("Euclid",buffer,6)==0)    hyprevars->hypre_prectyp = hypreprec_euclid;
         if (strncmp("Parasails",buffer,9)==0) hyprevars->hypre_prectyp = hypreprec_parasails;
         if (strncmp("BoomerAMG",buffer,9)==0) hyprevars->hypre_prectyp = hypreprec_amg;
      }

      frchar("HYPRE_IO",buffer,&ierr);
      if (ierr)
      {
       if (strncmp("full",buffer,4)==0)
          hyprevars->io=1;
       if (strncmp("none",buffer,4)==0)
          hyprevars->io=0;
      }

      frchar("HYPRE_PARASY",buffer,&ierr);
      if (ierr)
      {
       if (strncmp("symm_pos_def",buffer,12)==0) hyprevars->parasymm = 1;
       if (strncmp("nonsym_def",buffer,10)==0)   hyprevars->parasymm = 2;
       if (strncmp("non_sym_indef",buffer,13)==0)hyprevars->parasymm = 0;
      }

      frint("HYPRE_ITER"   ,&(hyprevars->maxiter)  ,&ierr);
      frint("HYPRE_IFILL"  ,&(hyprevars->ifill)    ,&ierr);
      frint("HYPRE_PARALEV",&(hyprevars->paralevel),&ierr);
      frint("HYPRE_SFINE"  ,&(hyprevars->sweep[0]) ,&ierr);
      frint("HYPRE_SDOWN"  ,&(hyprevars->sweep[1]) ,&ierr);
      frint("HYPRE_SUP"    ,&(hyprevars->sweep[2]) ,&ierr);
      frint("HYPRE_SCOARS" ,&(hyprevars->sweep[3]) ,&ierr);
      frint("HYPRE_KRYDIM" ,&(hyprevars->kryldim)  ,&ierr);

      frdouble("HYPRE_TOL"     ,&(hyprevars->tol)       ,&ierr);
      frdouble("HYPRE_DFILL"   ,&(hyprevars->dfill)     ,&ierr);
      frdouble("HYPRE_THREAS",&(hyprevars->threshold),&ierr);
      frdouble("HYPRE_PARATHR" ,&(hyprevars->parathresh),&ierr);
      frdouble("HYPRE_PARAFILT",&(hyprevars->parafilter),&ierr);
   break;
   case hypre_bicgstab:/*--------------------------- read solver BiCGStab */
      frchar("HYPRE_PREC",buffer,&ierr);
      if (ierr)
      {
         if (strncmp("none",buffer,4)==0)      hyprevars->hypre_prectyp = hypreprec_none;
         if (strncmp("Euclid",buffer,6)==0)    hyprevars->hypre_prectyp = hypreprec_euclid;
         if (strncmp("Parasails",buffer,9)==0) hyprevars->hypre_prectyp = hypreprec_parasails;
         if (strncmp("BoomerAMG",buffer,9)==0) hyprevars->hypre_prectyp = hypreprec_amg;
      }

      frchar("HYPRE_IO",buffer,&ierr);
      if (ierr)
      {
       if (strncmp("full",buffer,4)==0)
          hyprevars->io=1;
       if (strncmp("none",buffer,4)==0)
          hyprevars->io=0;
      }

      frchar("HYPRE_PARASY",buffer,&ierr);
      if (ierr)
      {
       if (strncmp("symm_pos_def",buffer,12)==0) hyprevars->parasymm = 1;
       if (strncmp("nonsym_def",buffer,10)==0)   hyprevars->parasymm = 2;
       if (strncmp("non_sym_indef",buffer,13)==0)hyprevars->parasymm = 0;
      }

      frint("HYPRE_ITER"   ,&(hyprevars->maxiter)  ,&ierr);
      frint("HYPRE_IFILL"  ,&(hyprevars->ifill)    ,&ierr);
      frint("HYPRE_PARALEV",&(hyprevars->paralevel),&ierr);
      frint("HYPRE_SFINE" ,&(hyprevars->sweep[0]),&ierr);
      frint("HYPRE_SDOWN" ,&(hyprevars->sweep[1]),&ierr);
      frint("HYPRE_SUP"   ,&(hyprevars->sweep[2]),&ierr);
      frint("HYPRE_SCOARS",&(hyprevars->sweep[3]),&ierr);

      frdouble("HYPRE_TOL"     ,&(hyprevars->tol)       ,&ierr);
      frdouble("HYPRE_DFILL"   ,&(hyprevars->dfill)     ,&ierr);
      frdouble("HYPRE_THREAS",&(hyprevars->threshold),&ierr);
      frdouble("HYPRE_PARATHR" ,&(hyprevars->parathresh),&ierr);
      frdouble("HYPRE_PARAFILT",&(hyprevars->parafilter),&ierr);
   break;
#endif

#ifdef MLIB_PACKAGE
   case mlib_d_sp:/*---------------------- read hp's direct mlib solver */
      frint("SYMM"   ,&(mlvar->symm  ),&ierr);
      frint("MSGLVL" ,&(mlvar->msglvl),&ierr);
      frint("MAXZER" ,&(mlvar->maxzer),&ierr);

      frdouble("PVTTOL",&(mlvar->pvttol) ,&ierr);

      frchar("ORDER" ,buffer,&ierr);
      if (ierr==1)
      {
         if (strncmp("NOO",buffer,5)==0) mlvar->order = 0;
         if (strncmp("MMD",buffer,5)==0) mlvar->order = 1;
         if (strncmp("NAT",buffer,5)==0) mlvar->order = 2;
         if (strncmp("CMD",buffer,5)==0) mlvar->order = 3;
         if (strncmp("RCM",buffer,5)==0) mlvar->order = 4;
         if (strncmp("1WD",buffer,5)==0) mlvar->order = 5;
         if (strncmp("MET",buffer,5)==0) mlvar->order = 6;
      }
   break;
#endif
   case parsuperlu:/*--------------------- read solver parallel SuperLU */
   break;
   case lapack_sym:/*------------------------------- read solver LAPACK */
   break;
   case lapack_nonsym:/*---------------------------- read solver LAPACK */
   break;
   case mumps_sym:/*--------------------------------- read solver MUMPS */
   break;
   case mumps_nonsym:/*------------------------------ read solver MUMPS */
   break;
   case umfpack:/*--------------------------------- read solver UMFPACK */
   break;
   case colsol_solver:/*---------------------------- read solver Colsol */
   break;
   case SPOOLES_nonsym:/*-------------------------- read solver Spooles */
   break;
   case SPOOLES_sym:/*----------------------------- read solver Spooles */
   break;

#ifdef MLPCG
   case mlpcg:/*------------------------------------- read solver MLPCG */
      frint("NUMLEV"      ,&(mlpcgvars->numlev)   ,&ierr);
      /*------------------------ there is a minimum of a 2-level method */
      if (ierr==1 && mlpcgvars->numlev < 2)
      dserror("Minimum number of levels to MLPCG is two");
      frint("ILU_N"       ,&(mlpcgvars->ilu_n)    ,&ierr);
      frint("COARSEILULEV",&(mlpcgvars->co_ilu_n) ,&ierr);
      frint("PRESWEEP"    ,&(mlpcgvars->presweep) ,&ierr);
      frint("POSTSWEEP"   ,&(mlpcgvars->postsweep),&ierr);
      frint("MAXITER"     ,&(mlpcgvars->maxiter)  ,&ierr);
      frint("COARSENUMDF" ,&(mlpcgvars->numdf)    ,&ierr);
      frint("OVERLAP"     ,&(mlpcgvars->overlap)  ,&ierr);
      frint("REUSE"       ,&(mlpcgvars->reuse)    ,&ierr);

      frdouble("PROLONGSMODAMP",&(mlpcgvars->p_omega),&ierr);
      frdouble("TOLERANCE"     ,&(mlpcgvars->tol)    ,&ierr);
      frdouble("GAMMA"         ,&(mlpcgvars->gamma)  ,&ierr);

      frchar("COARSESOLV" ,buffer,&ierr);
      if (ierr==1)
      {
         if (strncmp("SPOOLES",buffer,7)==0)
         {
            mlpcgvars->coarsesolv = co_spooles;
#ifndef SPOOLES_PACKAGE
            dserror("Spooles chosen as coarse solver but not compiled in");
#endif
         }
         if (strncmp("LAPACK",buffer,6)==0)  mlpcgvars->coarsesolv = co_lapack;
         if (strncmp("ILU",buffer,3)==0)     mlpcgvars->coarsesolv = co_ilu;
      }
      frchar("PRESMO" ,buffer,&ierr);
      if (ierr==1)
      {
         if (strncmp("fwdGS",buffer,5)==0)  mlpcgvars->presmoother = pre_fwdGS;
         if (strncmp("Jacobi",buffer,6)==0) mlpcgvars->presmoother = pre_Jacobi;
         if (strncmp("ILU",buffer,3)==0)    mlpcgvars->presmoother = pre_ilu;
      }
      frchar("POSTSMO" ,buffer,&ierr);
      if (ierr==1)
      {
         if (strncmp("bckGS",buffer,5)==0)  mlpcgvars->postsmoother = post_bckGS;
         if (strncmp("Jacobi",buffer,6)==0) mlpcgvars->postsmoother = post_Jacobi;
         if (strncmp("ILU",buffer,3)==0)    mlpcgvars->postsmoother = post_ilu;
      }
      frchar("TYPE" ,buffer,&ierr);
      if (ierr==1)
      {
         if (strncmp("Vanek",buffer,5)==0)  mlpcgvars->typ=2;
         if (strncmp("Fish",buffer,4)==0)   mlpcgvars->typ=1;
      }
   break;
#endif

   default:/*-------------------------------------------- reading error */
      dserror("Unknown solvertyp - error in reading");
   break;
   }/*------------------------------------ end of switch ofer solvertyp */
   /*------------------------------------------- read for typ of matrix */
   frchar("MATRIXTYP" ,buffer,&ierr);
   if (ierr==1)
   {
      if (strncmp("OLL",buffer,3)==0) solv->matrixtyp = oll_matrix;
   }

   /*------------------------------------- read for typ of partitioning */
   frchar("PARTITION" ,buffer,&ierr);
   if (ierr==1)
   {
      if (strncmp("Cut_Elements",buffer,12)==0) solv->parttyp = cut_elements;
      if (strncmp("Cut_Nodes"   ,buffer, 9)==0) solv->parttyp = cut_nodes;
   }

   frread();
}
frrewind();
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inpctrsol */



