#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
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
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
extern struct _MULTIMAT  *multimat;
/*----------------------------------------------------------------------*
 | input of materials                                     m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_material()
{
int  ierr, ierralloc, int_dummy;
int  i, j, k, ncm, num_klay, num_mlay;
struct    _KINLAY *actlay;           /*actual kinematic layer -> shell9 */ 
double    klay_sum;                  /*total hight or shell9*/
double    mlay_sum;                  /*hight of a kinematic layer*/

char *colpointer;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("inp_material");
#endif
/*----------------------------------------------------------------------*/
mat = (MATERIAL*)CCACALLOC(genprob.nmat,sizeof(MATERIAL));
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
      mat[i].m.fluid = (FLUID*)CCACALLOC(1,sizeof(FLUID));
      if (mat[i].m.fluid==NULL) dserror("Alloocation of FLUID material failed");
      frdouble("VISCOSITY",&(mat[i].m.fluid->viscosity),&ierr);
      frdouble("DENS"  ,&(mat[i].m.fluid->density)  ,&ierr);
   }
   frchk("MAT_Struct_StVenantKirchhoff",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_stvenant;
      mat[i].m.stvenant = (STVENANT*)CCACALLOC(1,sizeof(STVENANT));
      if (mat[i].m.stvenant==NULL) dserror("Alloocation of STVENANT material failed");
      frdouble("YOUNG"  ,&(mat[i].m.stvenant->youngs)      ,&ierr);
      frdouble("NUE"    ,&(mat[i].m.stvenant->possionratio),&ierr);
      frdouble("DENS",&(mat[i].m.stvenant->density)     ,&ierr);
   }
   frchk("MAT_Struct_STVENPOR",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_stvenpor;
      mat[i].m.stvenpor = (STVENPOR*)CCACALLOC(1,sizeof(STVENPOR));
      if (mat[i].m.stvenpor==NULL) dserror("Alloocation of Porous StVen. material failed");
      frdouble("YOUNG"  ,&(mat[i].m.stvenpor->youngs)        ,&ierr);
      frdouble("NUE"    ,&(mat[i].m.stvenpor->possionratio)  ,&ierr);
      frdouble("DENS"   ,&(mat[i].m.stvenpor->density)       ,&ierr);
      frdouble("REFDENS",&(mat[i].m.stvenpor->refdens)       ,&ierr);
      frdouble("EXPO"   ,&(mat[i].m.stvenpor->exponent)      ,&ierr);
   }
   frchk("MAT_Struct_NeoHooke",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_neohooke;
      mat[i].m.neohooke = (NEO_HOOKE*)CCACALLOC(1,sizeof(NEO_HOOKE));
      if (mat[i].m.neohooke==NULL) dserror("Alloocation of NEO_HOOKE material failed");
      frdouble("YOUNG",&(mat[i].m.neohooke->youngs)        ,&ierr);
      frdouble("NUE"  ,&(mat[i].m.neohooke->possionratio)  ,&ierr);
      frdouble("DENS",&(mat[i].m.neohooke->density)     ,&ierr);
   }
   frchk("MAT_3DMisesPlastic",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_pl_mises_3D;
      mat[i].m.pl_mises = (PL_MISES*)CCACALLOC(1,sizeof(PL_MISES));
      if (mat[i].m.pl_mises==NULL) dserror("Allocation of MISES_3D material failed");
      frdouble("YOUNG",&(mat[i].m.pl_mises->youngs)        ,&ierr);
      frdouble("NUE"  ,&(mat[i].m.pl_mises->possionratio)  ,&ierr);
      frdouble("ALFAT",&(mat[i].m.pl_mises->ALFAT)         ,&ierr);
      frdouble("Sigy" ,&(mat[i].m.pl_mises->Sigy)          ,&ierr);
      mat[i].m.pl_mises->Hard = 0.; 
      mat[i].m.pl_mises->GF   = 0.; 
      mat[i].m.pl_mises->betah= 1.; 
      frdouble("Hard" ,&(mat[i].m.pl_mises->Hard)          ,&ierr);
      frdouble("GF"   ,&(mat[i].m.pl_mises->GF)            ,&ierr);
      frdouble("BETAH",&(mat[i].m.pl_mises->betah)         ,&ierr);
   }
   frchk("MAT_MisesPlastic",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_pl_mises;
      mat[i].m.pl_mises = (PL_MISES*)CCACALLOC(1,sizeof(PL_MISES));
      if (mat[i].m.pl_mises==NULL) dserror("Alloocation of MISES material failed");
      frdouble("YOUNG",&(mat[i].m.pl_mises->youngs)        ,&ierr);
      frdouble("NUE"  ,&(mat[i].m.pl_mises->possionratio)  ,&ierr);
      frdouble("ALFAT",&(mat[i].m.pl_mises->ALFAT)         ,&ierr);
      frdouble("Sigy" ,&(mat[i].m.pl_mises->Sigy)          ,&ierr);
      mat[i].m.pl_mises->Hard = 0.; 
      mat[i].m.pl_mises->GF   = 0.; 
      mat[i].m.pl_mises->betah= 1.; 
      frdouble("Hard" ,&(mat[i].m.pl_mises->Hard)          ,&ierr);
      frdouble("GF"   ,&(mat[i].m.pl_mises->GF)            ,&ierr);
      frdouble("BETAH",&(mat[i].m.pl_mises->betah)         ,&ierr);
   }
   frchk("MAT_FoamPlastic",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_pl_foam;
      mat[i].m.pl_foam = (PL_FOAM*)CCACALLOC(1,sizeof(PL_FOAM));
      if (mat[i].m.pl_foam==NULL) dserror("Alloocation of FOAM material failed");
      frdouble("YOUNG",&(mat[i].m.pl_foam->youngs)        ,&ierr);
      frdouble("NUE"  ,&(mat[i].m.pl_foam->possionratio)  ,&ierr);
      frdouble("ALFAT",&(mat[i].m.pl_foam->ALFAT)         ,&ierr);
      frdouble("Sigy" ,&(mat[i].m.pl_foam->Sigy)          ,&ierr);
      mat[i].m.pl_foam->Hard = 0.; 
      mat[i].m.pl_foam->GF   = 0.; 
      frdouble("Hard" ,&(mat[i].m.pl_foam->Hard)          ,&ierr);
      frdouble("GF"   ,&(mat[i].m.pl_foam->GF)            ,&ierr);
   }
   frchk("MAT_DP_Plastic",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_pl_dp;
      mat[i].m.pl_dp = (PL_DP*)CCACALLOC(1,sizeof(PL_DP));
      if (mat[i].m.pl_dp==NULL) dserror("Alloocation of Drucker Prager material failed");
      frdouble("YOUNG",&(mat[i].m.pl_dp->youngs)        ,&ierr);
      frdouble("NUE"  ,&(mat[i].m.pl_dp->possionratio)  ,&ierr);
      frdouble("ALFAT",&(mat[i].m.pl_dp->ALFAT)         ,&ierr);
      frdouble("Sigy" ,&(mat[i].m.pl_dp->Sigy)          ,&ierr);
      frdouble("Hard" ,&(mat[i].m.pl_dp->Hard)          ,&ierr);
      frdouble("PHI"  ,&(mat[i].m.pl_dp->PHI)           ,&ierr);
   }
   frchk("MAT_ConcretePlastic",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_pl_epc;
      mat[i].m.pl_epc = (PL_EPC*)CCACALLOC(1,sizeof(PL_EPC));
      if (mat[i].m.pl_epc==NULL) dserror("Allocation of elpl-concrete material failed");
      /* initialize */
      mat[i].m.pl_epc->gamma1 = 3.;
      mat[i].m.pl_epc->gamma2 = 6./5.;
      
      
      frdouble("DENS"    ,&(mat[i].m.pl_epc->dens        )        ,&ierr);
      /* concrete */
      frdouble("YOUNG"   ,&(mat[i].m.pl_epc->youngs      )        ,&ierr);
      frdouble("NUE"     ,&(mat[i].m.pl_epc->possionratio)        ,&ierr);
      frdouble("ALFAT"   ,&(mat[i].m.pl_epc->alfat       )        ,&ierr);
      frdouble("XSI"     ,&(mat[i].m.pl_epc->xsi         )        ,&ierr);
      frdouble("Sigy"    ,&(mat[i].m.pl_epc->sigy        )        ,&ierr);
      frread();
      frdouble("FTM"     ,&(mat[i].m.pl_epc->ftm         )        ,&ierr);
      frdouble("FCM"     ,&(mat[i].m.pl_epc->fcm         )        ,&ierr);
      frdouble("GT"      ,&(mat[i].m.pl_epc->gt          )        ,&ierr);
      frdouble("GC"      ,&(mat[i].m.pl_epc->gc          )        ,&ierr);
      frdouble("GAMMA1"  ,&(mat[i].m.pl_epc->gamma1      )        ,&ierr);
      if(mat[i].m.pl_epc->gamma1<1.)mat[i].m.pl_epc->gamma1=3.; 
      frdouble("GAMMA2"  ,&(mat[i].m.pl_epc->gamma2      )        ,&ierr);
     /* tension stiffening - next line in input file!*/
      frread();
      frint(   "NSTIFF"  ,&(mat[i].m.pl_epc->nstiff      )        ,&ierr);
      /* number of rebars - next line in input file! */
      frread();
      frint(   "MAXREB"   ,&(mat[i].m.pl_epc->maxreb     )        ,&ierr);
      /* allocate memory */
      ncm       = mat[i].m.pl_epc->maxreb;
      ierralloc = 0;
      if ((mat[i].m.pl_epc->rebar=(int*)CCACALLOC(ncm,sizeof(int)))==NULL) ierralloc=1;
      if ((mat[i].m.pl_epc->reb_area  =(double*)CCACALLOC(ncm,sizeof(double)))==NULL) ierralloc=1;
      if ((mat[i].m.pl_epc->reb_ang   =(double*)CCACALLOC(ncm,sizeof(double)))==NULL) ierralloc=1;
      if ((mat[i].m.pl_epc->reb_so    =(double*)CCACALLOC(ncm,sizeof(double)))==NULL) ierralloc=1;
      if ((mat[i].m.pl_epc->reb_ds    =(double*)CCACALLOC(ncm,sizeof(double)))==NULL) ierralloc=1;
      if ((mat[i].m.pl_epc->reb_rgamma=(double*)CCACALLOC(ncm,sizeof(double)))==NULL) ierralloc=1;
      if ((mat[i].m.pl_epc->reb_dens  =(double*)CCACALLOC(ncm,sizeof(double)))==NULL) ierralloc=1;
      if ((mat[i].m.pl_epc->reb_alfat =(double*)CCACALLOC(ncm,sizeof(double)))==NULL) ierralloc=1;
      if ((mat[i].m.pl_epc->reb_emod  =(double*)CCACALLOC(ncm,sizeof(double)))==NULL) ierralloc=1;
      if ((mat[i].m.pl_epc->reb_rebnue=(double*)CCACALLOC(ncm,sizeof(double)))==NULL) ierralloc=1;
      if ((mat[i].m.pl_epc->reb_sigy  =(double*)CCACALLOC(ncm,sizeof(double)))==NULL) ierralloc=1;
      if ((mat[i].m.pl_epc->reb_hard  =(double*)CCACALLOC(ncm,sizeof(double)))==NULL) ierralloc=1;
      
      if (ierralloc) dserror("Allocation of elpl-concrete material failed");
      /* rebar data - next line in input file! */
      if(ncm==0)
      {
        frread();
        frread();
        frread();
      }
      for(j=0;j<ncm;j++)
      {
        frread();
        frint(   "REBAR"   ,&(mat[i].m.pl_epc->rebar[j]     ),&ierr);
        frdouble("REBAREA" ,&(mat[i].m.pl_epc->reb_area[j]  ),&ierr);
        frdouble("REBANG"  ,&(mat[i].m.pl_epc->reb_ang[j]   ),&ierr);
        frdouble("REBSO"   ,&(mat[i].m.pl_epc->reb_so[j]    ),&ierr);
        frdouble("REBDS"   ,&(mat[i].m.pl_epc->reb_ds[j]    ),&ierr);
        frdouble("REBGAMMA",&(mat[i].m.pl_epc->reb_rgamma[j]),&ierr);
        frread();
        frdouble("REBDENS" ,&(mat[i].m.pl_epc->reb_dens[j]  ),&ierr);
        frdouble("REBALFAT",&(mat[i].m.pl_epc->reb_alfat[j] ),&ierr);
        frdouble("REBEMOD" ,&(mat[i].m.pl_epc->reb_emod[j]  ),&ierr);
        frdouble("REBNUE"  ,&(mat[i].m.pl_epc->reb_rebnue[j]),&ierr);
        frread();                               
        frdouble("REBSIGY" ,&(mat[i].m.pl_epc->reb_sigy[j]  ),&ierr);
        frdouble("REBHARD" ,&(mat[i].m.pl_epc->reb_hard[j]  ),&ierr);
      }
   }
   frchk("MAT_Porous_MisesPlastic",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_pl_por_mises;
      mat[i].m.pl_por_mises = (PL_POR_MISES*)CCACALLOC(1,sizeof(PL_POR_MISES));
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
   /*multi layer material */
   frchk("MAT_Multilayer",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_multi_layer;
      mat[i].m.multi_layer = (MULTI_LAYER*)CCACALLOC(1,sizeof(MULTI_LAYER));
      if (mat[i].m.multi_layer==NULL) dserror("Alloocation of MULTI_LAYER material failed");
      /*read number of kinamtic layer*/
      frint("NUM_KLAY"   ,&(mat[i].m.multi_layer->num_klay)        ,&ierr);

      /*allocate memory for the different kinematic layers*/
      num_klay = mat[i].m.multi_layer->num_klay;
      ierralloc = 0;
      if ((mat[i].m.multi_layer->klayhgt=(double*)CCACALLOC(num_klay,sizeof(double)))==NULL) ierralloc=1;
      if (ierralloc) dserror("Allocation of klayhgt in MULTI_LAYER material failed");
      ierralloc = 0;
      if ((mat[i].m.multi_layer->kinlay=(KINLAY*)CCACALLOC(num_klay,sizeof(KINLAY)))==NULL) ierralloc=1;
      if (ierralloc) dserror("Allocation of KINLAYs in MULTI_LAYER material failed");

      /*----- read section data  -> hgt of different kinematic layers ---*/
      frdouble_n("SEC_KLAY",mat[i].m.multi_layer->klayhgt,num_klay,&ierr);
      if (ierr!=1) dserror("Reading of klayhgt in MULTI_LAYER material failed");
      /*------------------------ check if sectian data adds up to 100 ---*/
      klay_sum = 0.0;
      for(j=0; j<num_klay; j++) klay_sum += mat[i].m.multi_layer->klayhgt[j];
      if (FABS(klay_sum-100.0) > EPS5) dserror("klay_sum != 100 -> change Data in SEC_KLAY in MULTI_LAYER material");

      /*read one kinematic layer per line*/
      for(j=0; j<num_klay; j++)
      {
        actlay = &(mat[i].m.multi_layer->kinlay[j]);
        frread();
        
        frint("KINLAY"   ,&(int_dummy)        ,&ierr);
        /*read number of material layers to this kinematic layer*/
        frint("NUM_MLAY" ,&(actlay->num_mlay) ,&ierr);
        num_mlay = actlay->num_mlay;
        /*allocate memory for different material layers*/
        ierralloc = 0;
        if ((actlay->mlayhgt=(double*)CCACALLOC(num_mlay,sizeof(double)))==NULL) ierralloc=1;
        if (ierralloc) dserror("Allocation of mlayhgt in MULTI_LAYER material failed");
        ierralloc = 0;
        if ((actlay->mmatID=(int*)CCACALLOC(num_mlay,sizeof(int)))==NULL) ierralloc=1;
        if (ierralloc) dserror("Allocation of mmatID in MULTI_LAYER material failed");
        ierralloc = 0;
        if ((actlay->phi=(double*)CCACALLOC(num_mlay,sizeof(double)))==NULL) ierralloc=1;
        if (ierralloc) dserror("Allocation of phi in MULTI_LAYER material failed");
        ierralloc = 0;
        if ((actlay->rot=(int*)CCACALLOC(num_mlay,sizeof(int)))==NULL) ierralloc=1;
        if (ierralloc) dserror("Allocation of rot in MULTI_LAYER material failed");

        /*----- read section data of this kin layer -> hgt of different mat layers ---*/
        frdouble_n("SEC_MLAY",actlay->mlayhgt,num_mlay,&ierr);
        if (ierr!=1) dserror("Reading of mlayhgt (SEC_MLAY) in MULTI_LAYER material failed");
        frint_n("SEC_MAT",actlay->mmatID,num_mlay,&ierr);
        if (ierr!=1) dserror("Reading of mmatID (SEC_MAT) in MULTI_LAYER material failed");
        frdouble_n("SEC_PHI",actlay->phi,num_mlay,&ierr);
        if (ierr!=1) dserror("Reading of phi (SEC_PHI) in MULTI_LAYER material failed");
        frint_n("SEC_ROT",actlay->rot,num_mlay,&ierr);
        if (ierr!=1) dserror("Reading of rot (SEC_ROT) in MULTI_LAYER material failed");
        /*------------------------ check if sectian data adds up to 100 ---*/
        mlay_sum = 0.0;
        for(k=0; k<num_mlay; k++) mlay_sum += actlay->mlayhgt[k];
        if (FABS(mlay_sum-100.0) > EPS5) dserror("mlay_sum != 100 -> change Data in SEC_MLAY in MULTI_LAYER material");
                
      }
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


/*----------------------------------------------------------------------*
 | input of multilayer materials                            sh 10/02    |
 *----------------------------------------------------------------------*/
void inp_multimat()
{
int  counter=0;
int  check = 0;
int  ierr, ierralloc;
int  i, j, ncm;
char *colpointer;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("inp_multimat");
#endif
/*---- read multilayer material only if a mattyp == m_multi_layer ------*/
for (i=0; i<genprob.nmat; i++)     
{ 
   if  (mat[i].mattyp == m_multi_layer) check = 1;
}
if (check == 0 ) goto end;
/*--------------------------------- count number of multilayer materials */
frrewind();
frfind("--MULTILAYER MATERIALS");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   counter++;
   frread();
}
if (counter == 0) goto end;       /* no multilayer material set */

/*--------------------------------------------------- allocate MULTIMAT */
multimat = (MULTIMAT*)CCACALLOC(counter,sizeof(MULTIMAT));
if (multimat==NULL) dserror("Allocation of MULTIMAT failed");
/*----------------------------------------------------------------------*/
frrewind();
frfind("--MULTILAYER MATERIALS");
frread();
i=0;
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frint("MULTIMAT",&(multimat[i].Id),&ierr);

   frchk("MAT_Struct_StVenantKirchhoff",&ierr);
   if (ierr==1)
   {
      multimat[i].mattyp = m_stvenant;
      multimat[i].m.stvenant = (STVENANT*)CCACALLOC(1,sizeof(STVENANT));
      if (multimat[i].m.stvenant==NULL) dserror("Allocation of STVENANT material failed");
      frdouble("YOUNG"  ,&(multimat[i].m.stvenant->youngs)      ,&ierr);
      frdouble("NUE"    ,&(multimat[i].m.stvenant->possionratio),&ierr);
      frdouble("DENS"   ,&(multimat[i].m.stvenant->density)     ,&ierr);
   }
   frchk("MAT_Struct_NeoHooke",&ierr);
   if (ierr==1)
   {
      multimat[i].mattyp = m_neohooke;
      multimat[i].m.neohooke = (NEO_HOOKE*)CCACALLOC(1,sizeof(NEO_HOOKE));
      if (multimat[i].m.neohooke==NULL) dserror("Allocation of NEO_HOOKE material failed");
      frdouble("YOUNG"   ,&(multimat[i].m.neohooke->youngs)        ,&ierr);
      frdouble("NUE"     ,&(multimat[i].m.neohooke->possionratio)  ,&ierr);
      frdouble("DENS"    ,&(multimat[i].m.neohooke->density)       ,&ierr);
   }
   frchk("MAT_Struct_Orthotropic",&ierr);
   if (ierr==1)
   {
      multimat[i].mattyp = m_orthotropic;
      multimat[i].m.orthotropic = (ORTHOTROPIC*)CCACALLOC(1,sizeof(ORTHOTROPIC));
      if (multimat[i].m.orthotropic==NULL) dserror("Allocation of ORTHOTROPIC material failed");
      frdouble("EMOD1"   ,&(multimat[i].m.orthotropic->emod1)        ,&ierr);
      frdouble("EMOD2"   ,&(multimat[i].m.orthotropic->emod2)        ,&ierr);
      frdouble("EMOD3"   ,&(multimat[i].m.orthotropic->emod3)        ,&ierr);
      frdouble("GMOD12"  ,&(multimat[i].m.orthotropic->gmod12)       ,&ierr);
      frdouble("GMOD13"  ,&(multimat[i].m.orthotropic->gmod13)       ,&ierr);
      frdouble("GMOD23"  ,&(multimat[i].m.orthotropic->gmod23)       ,&ierr);
      frdouble("XNUE12"  ,&(multimat[i].m.orthotropic->xnue12)       ,&ierr);
      frdouble("XNUE13"  ,&(multimat[i].m.orthotropic->xnue13)       ,&ierr);
      frdouble("XNUE23"  ,&(multimat[i].m.orthotropic->xnue23)       ,&ierr);
   }
   i++;
/*----------------------------------------------------------------------*/
   frread();
}
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_multimat */
