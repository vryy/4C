#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 | input of materials                                     m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_material()
{
int  ierr;
int  i;
char *colpointer;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("inp_material");
#endif
/*----------------------------------------------------------------------*/
mat = (MATERIAL*)calloc(genprob.nmat,sizeof(MATERIAL));
if (mat==NULL) dserror("Allocation of MATERIAL failed");
/*----------------------------------------------------------------------*/
frfind("--MATERIALS");
frread();
i=0;
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   if (i==genprob.nmat) dserror("number of materials incorrect");

   frint("MAT",&(mat[i].Id),&ierr);

   frchk("MAT_fluid",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_fluid;
      mat[i].m.fluid = (FLUID*)calloc(1,sizeof(FLUID));
      if (mat[i].m.fluid==NULL) dserror("Alloocation of FLUID material failed");
      frdouble("VISCOSITY",&(mat[i].m.fluid->viscosity),&ierr);
      frdouble("DENS"  ,&(mat[i].m.fluid->density)  ,&ierr);
   }
   frchk("MAT_Struct_StVenantKirchhoff",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_stvenant;
      mat[i].m.stvenant = (STVENANT*)calloc(1,sizeof(STVENANT));
      if (mat[i].m.stvenant==NULL) dserror("Alloocation of STVENANT material failed");
      frdouble("YOUNG"  ,&(mat[i].m.stvenant->youngs)      ,&ierr);
      frdouble("NUE"    ,&(mat[i].m.stvenant->possionratio),&ierr);
      frdouble("DENS",&(mat[i].m.stvenant->density)     ,&ierr);
   }
   frchk("MAT_Struct_NeoHooke",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_neohooke;
      mat[i].m.neohooke = (NEO_HOOKE*)calloc(1,sizeof(NEO_HOOKE));
      if (mat[i].m.neohooke==NULL) dserror("Alloocation of NEO_HOOKE material failed");
      frdouble("YOUNG",&(mat[i].m.neohooke->youngs)        ,&ierr);
      frdouble("NUE"  ,&(mat[i].m.neohooke->possionratio)  ,&ierr);
      frdouble("DENSITY",&(mat[i].m.neohooke->density)     ,&ierr);
   }
   frchk("MAT_MisesPlastic",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_pl_mises;
      mat[i].m.pl_mises = (PL_MISES*)calloc(1,sizeof(PL_MISES));
      if (mat[i].m.pl_mises==NULL) dserror("Alloocation of MISES material failed");
      frdouble("YOUNG",&(mat[i].m.pl_mises->youngs)        ,&ierr);
      frdouble("NUE"  ,&(mat[i].m.pl_mises->possionratio)  ,&ierr);
      frdouble("ALFAT",&(mat[i].m.pl_mises->ALFAT)         ,&ierr);
      frdouble("Sigy" ,&(mat[i].m.pl_mises->Sigy)          ,&ierr);
      mat[i].m.pl_mises->Hard = 0.; 
      mat[i].m.pl_mises->GF   = 0.; 
      frdouble("Hard" ,&(mat[i].m.pl_mises->Hard)          ,&ierr);
      frdouble("GF"   ,&(mat[i].m.pl_mises->GF)            ,&ierr);
   }
   frchk("MAT_DP_Plastic",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_pl_dp;
      mat[i].m.pl_dp = (PL_DP*)calloc(1,sizeof(PL_DP));
      if (mat[i].m.pl_dp==NULL) dserror("Alloocation of Drucker Prager material failed");
      frdouble("YOUNG",&(mat[i].m.pl_dp->youngs)        ,&ierr);
      frdouble("NUE"  ,&(mat[i].m.pl_dp->possionratio)  ,&ierr);
      frdouble("ALFAT",&(mat[i].m.pl_dp->ALFAT)         ,&ierr);
      frdouble("Sigy" ,&(mat[i].m.pl_dp->Sigy)          ,&ierr);
      frdouble("Hard" ,&(mat[i].m.pl_dp->Hard)          ,&ierr);
      frdouble("PHI"  ,&(mat[i].m.pl_dp->PHI)           ,&ierr);
   }
   frchk("MAT_Porous_MisesPlastic",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_pl_por_mises;
      mat[i].m.pl_por_mises = (PL_POR_MISES*)calloc(1,sizeof(PL_POR_MISES));
      if (mat[i].m.pl_por_mises==NULL) dserror("Alloocation of MISES material failed");
      frdouble("YOUNG"   ,&(mat[i].m.pl_por_mises->youngs)        ,&ierr);
      frdouble("DP_YM"   ,&(mat[i].m.pl_por_mises->DP_YM )        ,&ierr);
      frdouble("NUE"     ,&(mat[i].m.pl_por_mises->possionratio)  ,&ierr);
      frdouble("ALFAT"   ,&(mat[i].m.pl_por_mises->ALFAT)         ,&ierr);
      frdouble("Sigy"    ,&(mat[i].m.pl_por_mises->Sigy)          ,&ierr);
      frdouble("DP_Sigy" ,&(mat[i].m.pl_por_mises->DP_Sigy)       ,&ierr);
      frdouble("Hard"    ,&(mat[i].m.pl_por_mises->Hard)          ,&ierr);
      frdouble("DP_Hard" ,&(mat[i].m.pl_por_mises->DP_Hard)       ,&ierr);
   }
   i++;
/*----------------------------------------------------------------------*/
   frread();
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_material */
