#include "../headers/standardtypes.h"
#include "../headers/solution.h"

/*----------------------------------------------------------------------*
 | input of solver control variables  structure           m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inpctrsol(SOLVAR *solv)
{
int           ierr;
char          buffer[50];
AZVAR        *azvar;
HYPREVARS    *hyprevars;
PSUPERLUVARS *psuperluvars;
LAPACKVARS   *lapackvars;
MUMPSVARS    *mumpsvars;
#ifdef DEBUG 
dstrc_enter("inpctrsol");
#endif
/*----------------------------------------------------------------------*/
switch(solv->fieldtyp)
{
case structure:
frfind("-STRUCT SOLVER");
break;
case fluid:
frfind("-FLUID SOLVER");
break;
case ale:
frfind("-ALE SOLVER");
break;
default:
   dserror("Unknown field typ in reading solver");
break;
}
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*------------------------------------ check typ of solver and allocate */
   frchar("SOLVER"    ,buffer,&ierr);
   if (ierr==1)
   {
      if (strncmp("Aztec_MSR",buffer,9)==0) 
      {
#ifndef AZTEC_PACKAGE
         dserror("Aztec package is not compiled in");
#endif
         solv->solvertyp = aztec_msr;
         solv->azvar = (AZVAR*)calloc(1,sizeof(AZVAR));
         if (!(solv->azvar)) dserror("Allocation of AZVAR failed");
         azvar = solv->azvar;
      }
      if (strncmp("HYPRE_BoomerAMG",buffer,15)==0) 
      {
#ifndef HYPRE_PACKAGE
         dserror("Hypre package is not compiled in");
#endif
         solv->solvertyp = hypre_amg;
         solv->hyprevar = (HYPREVARS*)calloc(1,sizeof(HYPREVARS));
         if (!(solv->hyprevar)) dserror("Allocation of HYPREVARS failed");
         hyprevars = solv->hyprevar;
      }
      if (strncmp("HYPRE_PCG",buffer,9)==0) 
      {
#ifndef HYPRE_PACKAGE
         dserror("Hypre package is not compiled in");
#endif
         solv->solvertyp = hypre_pcg;
         solv->hyprevar = (HYPREVARS*)calloc(1,sizeof(HYPREVARS));
         if (!(solv->hyprevar)) dserror("Allocation of HYPREVARS failed");
         hyprevars = solv->hyprevar;
      }
      if (strncmp("HYPRE_GMRES",buffer,11)==0) 
      {
#ifndef HYPRE_PACKAGE
         dserror("Hypre package is not compiled in");
#endif
         solv->solvertyp = hypre_gmres;
         solv->hyprevar = (HYPREVARS*)calloc(1,sizeof(HYPREVARS));
         if (!(solv->hyprevar)) dserror("Allocation of HYPREVARS failed");
         hyprevars = solv->hyprevar;
      }
      if (strncmp("HYPRE_BiCGStab",buffer,14)==0) 
      {
#ifndef HYPRE_PACKAGE
         dserror("Hypre package is not compiled in");
#endif
         solv->solvertyp = hypre_bicgstab;
         solv->hyprevar = (HYPREVARS*)calloc(1,sizeof(HYPREVARS));
         if (!(solv->hyprevar)) dserror("Allocation of HYPREVARS failed");
         hyprevars = solv->hyprevar;
      }
      if (strncmp("ParSuperLU",buffer,10)==0) 
      {
#ifndef PARSUPERLU_PACKAGE
         dserror("ParSuperLU package is not compiled in");
#endif
         solv->solvertyp = parsuperlu;
         solv->psuperluvars = (PSUPERLUVARS*)calloc(1,sizeof(PSUPERLUVARS));
         if (!(solv->psuperluvars)) dserror("Allocation of PSUPERLUVARS failed");
         psuperluvars = solv->psuperluvars;
      }
      if (strncmp("LAPACK_sym",buffer,10)==0) 
      {
         solv->solvertyp = lapack_sym;
         solv->lapackvars = (LAPACKVARS*)calloc(1,sizeof(LAPACKVARS));
         if (!(solv->lapackvars)) dserror("Allocation of LAPACKVARS failed");
         lapackvars = solv->lapackvars;
      }
      if (strncmp("LAPACK_nonsym",buffer,13)==0) 
      {
         solv->solvertyp = lapack_nonsym;
         solv->lapackvars = (LAPACKVARS*)calloc(1,sizeof(LAPACKVARS));
         if (!(solv->lapackvars)) dserror("Allocation of LAPACKVARS failed");
         lapackvars = solv->lapackvars;
      }
      if (strncmp("MUMPS_sym",buffer,9)==0) 
      {
#ifndef MUMPS_PACKAGE
         dserror("MUMPS package is not compiled in");
#endif
         solv->solvertyp = mumps_sym;
         solv->mumpsvars = (MUMPSVARS*)calloc(1,sizeof(MUMPSVARS));
         if (!(solv->mumpsvars)) dserror("Allocation of MUMPSVARS failed");
         mumpsvars = solv->mumpsvars;
      }
      if (strncmp("MUMPS_nonsym",buffer,12)==0) 
      {
#ifndef MUMPS_PACKAGE
         dserror("MUMPS package is not compiled in");
#endif
         solv->solvertyp = mumps_nonsym;
         solv->mumpsvars = (MUMPSVARS*)calloc(1,sizeof(MUMPSVARS));
         if (!(solv->mumpsvars)) dserror("Allocation of MUMPSVARS failed");
         mumpsvars = solv->mumpsvars;
      }
   }/* end of (ierr==1) */
/*------------------------------------------- read solver specific data */
   switch (solv->solvertyp)
   {
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
      }
      frint("AZREUSE",&(azvar->azreuse),&ierr);
      frint("AZGFILL",&(azvar->azgfill),&ierr);
      frint("AZITER" ,&(azvar->aziter),&ierr);
      frint("AZSUB"  ,&(azvar->azsub),&ierr);
      frint("AZGRAPH",&(azvar->azgraph),&ierr);
      frint("AZPOLY" ,&(azvar->azpoly),&ierr);
      frdouble("AZDROP" ,&(azvar->azdrop)  ,&ierr);
      frdouble("AZFILL" ,&(azvar->azfill)  ,&ierr);
      frdouble("AZTOL"  ,&(azvar->aztol)   ,&ierr);
      frdouble("AZOMEGA",&(azvar->azomega) ,&ierr);
   break;
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
   default:/*-------------------------------------------- reading error */
      dserror("Unknown solvertyp - error in reading");
   break;
   }/*------------------------------------ end of switch ofer solvertyp */

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
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpctrsol */



