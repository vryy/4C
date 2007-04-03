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
INT  ierr, ierralloc, int_dummy;
INT  i, j, k, ncm, num_klay, num_mlay;
struct    _KINLAY *actlay;           /*actual kinematic layer -> shell9 */
DOUBLE    klay_sum;                  /*total hight or shell9*/
DOUBLE    mlay_sum;                  /*hight of a kinematic layer*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("inp_material");
#endif
/*----------------------------------------------------------------------*/
mat = (MATERIAL*)CCACALLOC(genprob.nmat,sizeof(MATERIAL));
/*----------------------------------------------------------------------*/
if (frfind("--MATERIALS")==0) dserror("frfind: MATERIALS is not in input file");
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
      frdouble("VISCOSITY",&(mat[i].m.fluid->viscosity),&ierr);
      frdouble("DENS"  ,&(mat[i].m.fluid->density)  ,&ierr);
      frdouble("GAMMA",&(mat[i].m.fluid->gamma)  ,&ierr);
   }
   frchk("MAT_Struct_StVenantKirchhoff",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_stvenant;
      mat[i].m.stvenant = (STVENANT*)CCACALLOC(1,sizeof(STVENANT));
      frdouble("YOUNG"  ,&(mat[i].m.stvenant->youngs)      ,&ierr);
      frdouble("NUE"    ,&(mat[i].m.stvenant->possionratio),&ierr);
      frdouble("DENS"   ,&(mat[i].m.stvenant->density)     ,&ierr);
      frdouble("THEXPANS",&(mat[i].m.stvenant->thermexpans) ,&ierr);
   }
   frchk("MAT_Struct_Orthotropic",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_el_orth;
      mat[i].m.el_orth = (EL_ORTH*)CCACALLOC(1,sizeof(EL_ORTH));
      frdouble("EMOD1"   ,&(mat[i].m.el_orth->emod1)        ,&ierr);
      frdouble("EMOD2"   ,&(mat[i].m.el_orth->emod2)        ,&ierr);
      frdouble("EMOD3"   ,&(mat[i].m.el_orth->emod3)        ,&ierr);
      frdouble("GMOD12"  ,&(mat[i].m.el_orth->gmod12)       ,&ierr);
      frdouble("GMOD13"  ,&(mat[i].m.el_orth->gmod13)       ,&ierr);
      frdouble("GMOD23"  ,&(mat[i].m.el_orth->gmod23)       ,&ierr);
      frdouble("XNUE12"  ,&(mat[i].m.el_orth->xnue12)       ,&ierr);
      frdouble("XNUE13"  ,&(mat[i].m.el_orth->xnue13)       ,&ierr);
      frdouble("XNUE23"  ,&(mat[i].m.el_orth->xnue23)       ,&ierr);
   }
   frchk("MAT_Struct_STVENPOR",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_stvenpor;
      mat[i].m.stvenpor = (STVENPOR*)CCACALLOC(1,sizeof(STVENPOR));
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
      frdouble("YOUNG",&(mat[i].m.neohooke->youngs)        ,&ierr);
      frdouble("NUE"  ,&(mat[i].m.neohooke->possionratio)  ,&ierr);
      frdouble("DENS",&(mat[i].m.neohooke->density)        ,&ierr);
   }
   frchk("MAT_MFOC",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_mfoc;
      mat[i].m.mfoc = (MFOC*)CCACALLOC(1,sizeof(MFOC));
      if (mat[i].m.mfoc==NULL) dserror("Alloocation of Open Cell foam material failed");
      frdouble("Es"     ,&(mat[i].m.mfoc->es)          ,&ierr);
      frdouble("pr"     ,&(mat[i].m.mfoc->pr)          ,&ierr);
      frdouble("dens"   ,&(mat[i].m.mfoc->dens)        ,&ierr);  /* por. density */
      frdouble("denss"  ,&(mat[i].m.mfoc->denss)       ,&ierr);  /* ref. density */
      frdouble("oce"    ,&(mat[i].m.mfoc->oce)         ,&ierr);
      frdouble("ocf"    ,&(mat[i].m.mfoc->ocf)         ,&ierr);
      frdouble("densmin",&(mat[i].m.mfoc->denmin)      ,&ierr);
      frdouble("densmax",&(mat[i].m.mfoc->denmax)      ,&ierr);
   }
   frchk("MAT_MFCC",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_mfcc;
      mat[i].m.mfcc = (MFCC*)CCACALLOC(1,sizeof(MFCC));
      if (mat[i].m.mfoc==NULL) dserror("Alloocation of Closed Cell foam material failed");
      frdouble("Es"     ,&(mat[i].m.mfcc->es)          ,&ierr);
      frdouble("pr"     ,&(mat[i].m.mfcc->pr)          ,&ierr);
      frdouble("dens"   ,&(mat[i].m.mfcc->dens)        ,&ierr);  /* por. density */
      frdouble("denss"  ,&(mat[i].m.mfcc->denss)       ,&ierr);  /* ref. density */
      frdouble("cce"    ,&(mat[i].m.mfcc->cce)         ,&ierr);
      frdouble("ccf"    ,&(mat[i].m.mfcc->ccf)         ,&ierr);
      frdouble("densmin",&(mat[i].m.mfcc->denmin)      ,&ierr);
      frdouble("densmax",&(mat[i].m.mfcc->denmax)      ,&ierr);
   }
   frchk("MAT_NeoHMFCC",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_nhmfcc;
      mat[i].m.nhmfcc = (NHMFCC*)CCACALLOC(1,sizeof(NHMFCC));
      if (mat[i].m.nhmfcc==NULL) dserror("Alloocation of foam material failed");
      frdouble("Es"     ,&(mat[i].m.nhmfcc->es)          ,&ierr);
      frdouble("pr"     ,&(mat[i].m.nhmfcc->pr)          ,&ierr);
      frdouble("dens"   ,&(mat[i].m.nhmfcc->dens)        ,&ierr);  /* por. density */
      frdouble("denss"  ,&(mat[i].m.nhmfcc->denss)       ,&ierr);  /* ref. density */
      frdouble("cce"    ,&(mat[i].m.nhmfcc->cce)         ,&ierr);
      frdouble("ccf"    ,&(mat[i].m.nhmfcc->ccf)         ,&ierr);
      frdouble("densmin",&(mat[i].m.nhmfcc->denmin)      ,&ierr);
      frdouble("densmax",&(mat[i].m.nhmfcc->denmax)      ,&ierr);
   }
   frchk("MAT_Struct_Ogden",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_compogden;
      mat[i].m.compogden = (COMPOGDEN*)CCACALLOC(1,sizeof(COMPOGDEN));
      frdouble("NUE"  ,&(mat[i].m.compogden->nue)     ,&ierr);
      frdouble("BETA" ,&(mat[i].m.compogden->beta)    ,&ierr);
      frdouble("ALFA1",&(mat[i].m.compogden->alfap[0]),&ierr);
      frdouble("ALFA2",&(mat[i].m.compogden->alfap[1]),&ierr);
      frdouble("ALFA3",&(mat[i].m.compogden->alfap[2]),&ierr);
      frdouble("NU1"  ,&(mat[i].m.compogden->mup[0])  ,&ierr);
      frdouble("NU2"  ,&(mat[i].m.compogden->mup[1])  ,&ierr);
      frdouble("NU3"  ,&(mat[i].m.compogden->mup[2])  ,&ierr);
      frdouble("DENS" ,&(mat[i].m.compogden->density) ,&ierr);
   }
   frchk("MAT_Struct_Viscohyper",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_viscohyper;
      mat[i].m.viscohyper = (VISCOHYPER*)CCACALLOC(1,sizeof(VISCOHYPER));
      frdouble("NUE"  ,&(mat[i].m.viscohyper->nue)     ,&ierr);
      frdouble("BETA" ,&(mat[i].m.viscohyper->beta)    ,&ierr);
      frdouble("ALFA1",&(mat[i].m.viscohyper->alfap[0]),&ierr);
      frdouble("ALFA2",&(mat[i].m.viscohyper->alfap[1]),&ierr);
      frdouble("ALFA3",&(mat[i].m.viscohyper->alfap[2]),&ierr);
      frdouble("NU1"  ,&(mat[i].m.viscohyper->mup[0])  ,&ierr);
      frdouble("NU2"  ,&(mat[i].m.viscohyper->mup[1])  ,&ierr);
      frdouble("NU3"  ,&(mat[i].m.viscohyper->mup[2])  ,&ierr);
      frdouble("DENS" ,&(mat[i].m.viscohyper->density) ,&ierr);
      frint   ("NMAXW",&(mat[i].m.viscohyper->nmaxw)   ,&ierr);
      frdouble("TAU1" ,&(mat[i].m.viscohyper->tau[0])  ,&ierr);
      frdouble("TAU2" ,&(mat[i].m.viscohyper->tau[1])  ,&ierr);
      frdouble("TAU3" ,&(mat[i].m.viscohyper->tau[2])  ,&ierr);
      frdouble("TAU4" ,&(mat[i].m.viscohyper->tau[3])  ,&ierr);
      frdouble("BETA1",&(mat[i].m.viscohyper->betas[0]),&ierr);
      frdouble("BETA2",&(mat[i].m.viscohyper->betas[1]),&ierr);
      frdouble("BETA3",&(mat[i].m.viscohyper->betas[2]),&ierr);
      frdouble("BETA4",&(mat[i].m.viscohyper->betas[3]),&ierr);
   }
   frchk("MAT_3DMisesPlastic",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_pl_mises_3D;
      mat[i].m.pl_mises = (PL_MISES*)CCACALLOC(1,sizeof(PL_MISES));
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
   frchk("MAT_Damage",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_damage;
      mat[i].m.damage = (DAMAGE*)CCACALLOC(1,sizeof(DAMAGE));
      if (mat[i].m.damage==NULL) dserror("Allocation of DAMAGE material failed");
      frdouble("YOUNG",&(mat[i].m.damage->youngs)        ,&ierr);
      frdouble("NUE"  ,&(mat[i].m.damage->possionratio)  ,&ierr);
      frint("Equival" ,&(mat[i].m.damage->Equival)       ,&ierr);
      frint("Damtyp"  ,&(mat[i].m.damage->Damtyp)        ,&ierr);
      frdouble("Kappa_0",&(mat[i].m.damage->Kappa_0)     ,&ierr);
      frdouble("Kappa_m",&(mat[i].m.damage->Kappa_m)     ,&ierr);
      frdouble("Alpha",  &(mat[i].m.damage->Alpha)       ,&ierr);
      frdouble("Beta" ,  &(mat[i].m.damage->Beta)        ,&ierr);
      frdouble("k_fac" ,  &(mat[i].m.damage->k_fac)      ,&ierr);
   }
   frchk("MAT_FoamPlastic",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_pl_foam;
      mat[i].m.pl_foam = (PL_FOAM*)CCACALLOC(1,sizeof(PL_FOAM));
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
      frdouble("YOUNG",&(mat[i].m.pl_dp->youngs)        ,&ierr);
      frdouble("NUE"  ,&(mat[i].m.pl_dp->possionratio)  ,&ierr);
      frdouble("ALFAT",&(mat[i].m.pl_dp->ALFAT)         ,&ierr);
      frdouble("Sigy" ,&(mat[i].m.pl_dp->Sigy)          ,&ierr);
      frdouble("PHI"  ,&(mat[i].m.pl_dp->PHI)           ,&ierr);
      mat[i].m.pl_dp->Hard = 0.;
      mat[i].m.pl_dp->GF   = 0.;
      mat[i].m.pl_dp->betah= 1.;
      frdouble("Hard" ,&(mat[i].m.pl_dp->Hard)          ,&ierr);
      frdouble("GF"   ,&(mat[i].m.pl_dp->GF)            ,&ierr);
      frdouble("BETAH",&(mat[i].m.pl_dp->betah)         ,&ierr);
   }
   frchk("MAT_ConcretePlastic",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_pl_epc;
      mat[i].m.pl_epc = (PL_EPC*)CCACALLOC(1,sizeof(PL_EPC));
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
      mat[i].m.pl_epc->maxreb = 0;
      frint(   "MAXREB"   ,&(mat[i].m.pl_epc->maxreb     )        ,&ierr);
      /* allocate memory */
      ncm       = mat[i].m.pl_epc->maxreb;
      ierralloc = 0;
      mat[i].m.pl_epc->rebar=(INT*)CCACALLOC(ncm,sizeof(INT));
      mat[i].m.pl_epc->reb_area  =(DOUBLE*)CCACALLOC(ncm,sizeof(DOUBLE));
      mat[i].m.pl_epc->reb_ang   =(DOUBLE*)CCACALLOC(ncm,sizeof(DOUBLE));
      mat[i].m.pl_epc->reb_so    =(DOUBLE*)CCACALLOC(ncm,sizeof(DOUBLE));
      mat[i].m.pl_epc->reb_ds    =(DOUBLE*)CCACALLOC(ncm,sizeof(DOUBLE));
      mat[i].m.pl_epc->reb_rgamma=(DOUBLE*)CCACALLOC(ncm,sizeof(DOUBLE));
      mat[i].m.pl_epc->reb_dens  =(DOUBLE*)CCACALLOC(ncm,sizeof(DOUBLE));
      mat[i].m.pl_epc->reb_alfat =(DOUBLE*)CCACALLOC(ncm,sizeof(DOUBLE));
      mat[i].m.pl_epc->reb_emod  =(DOUBLE*)CCACALLOC(ncm,sizeof(DOUBLE));
      mat[i].m.pl_epc->reb_rebnue=(DOUBLE*)CCACALLOC(ncm,sizeof(DOUBLE));
      mat[i].m.pl_epc->reb_sigy  =(DOUBLE*)CCACALLOC(ncm,sizeof(DOUBLE));
      mat[i].m.pl_epc->reb_hard  =(DOUBLE*)CCACALLOC(ncm,sizeof(DOUBLE));

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
   frchk("MAT_3DConcretePlastic",&ierr);  /* 3D formulation, the same as used for multilayer materials*/
   if (ierr==1)
   {
      /*write a warning to use an unsymmetric solver*/
      printf("|---------------------------------------------------------------------------------------| \n");
      printf("|    WARNING in input_material.c:    ================================================   | \n");
      printf("|    WARNING in input_material.c:    MAT_ConcretePlastic -> use an UNsymmetric solver   | \n");
      printf("|    WARNING in input_material.c:    if trial stresses could be in the apex region      | \n");
      printf("|    WARNING in input_material.c:    ================================================   | \n");
      printf("|---------------------------------------------------------------------------------------| \n");

      mat[i].mattyp = m_pl_epc3D;
      mat[i].m.pl_epc = (PL_EPC*)CCACALLOC(1,sizeof(PL_EPC));

      /* initialize */
      mat[i].m.pl_epc->dfac = 0.;
      mat[i].m.pl_epc->gamma1 = 3.;
      mat[i].m.pl_epc->gamma2 = 6./5.;
      mat[i].m.pl_epc->gamma3 = 1./3.;
      mat[i].m.pl_epc->gamma4 = 4./3.;
      mat[i].m.pl_epc->maxreb = 0;

      /* concrete */
      frdouble("YOUNG"   ,&(mat[i].m.pl_epc->youngs      )        ,&ierr);
      frdouble("NUE"     ,&(mat[i].m.pl_epc->possionratio)        ,&ierr);
      frdouble("FTM"     ,&(mat[i].m.pl_epc->ftm         )        ,&ierr);
      frdouble("FCM"     ,&(mat[i].m.pl_epc->fcm         )        ,&ierr);
      frdouble("GT"      ,&(mat[i].m.pl_epc->gt          )        ,&ierr);
      frdouble("GC"      ,&(mat[i].m.pl_epc->gc          )        ,&ierr);
      frdouble("DFAC"    ,&(mat[i].m.pl_epc->dfac        )        ,&ierr);
      frdouble("GAMMA1"  ,&(mat[i].m.pl_epc->gamma1      )        ,&ierr);
      frdouble("GAMMA2"  ,&(mat[i].m.pl_epc->gamma2      )        ,&ierr);
      frdouble("GAMMA3"  ,&(mat[i].m.pl_epc->gamma3      )        ,&ierr);
      frdouble("GAMMA4"  ,&(mat[i].m.pl_epc->gamma4      )        ,&ierr);
   }
   frchk("MAT_Porous_MisesPlastic",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_pl_por_mises;
      mat[i].m.pl_por_mises = (PL_POR_MISES*)CCACALLOC(1,sizeof(PL_POR_MISES));
      frdouble("YOUNG"   ,&(mat[i].m.pl_por_mises->youngs)        ,&ierr);
      frdouble("DP_YM"   ,&(mat[i].m.pl_por_mises->DP_YM )        ,&ierr);
      frdouble("NUE"     ,&(mat[i].m.pl_por_mises->possionratio)  ,&ierr);
      frdouble("ALFAT"   ,&(mat[i].m.pl_por_mises->ALFAT)         ,&ierr);
      frdouble("Sigy"    ,&(mat[i].m.pl_por_mises->Sigy)          ,&ierr);
      frdouble("DP_Sigy" ,&(mat[i].m.pl_por_mises->DP_Sigy)       ,&ierr);
      frdouble("Hard"    ,&(mat[i].m.pl_por_mises->Hard)          ,&ierr);
      frdouble("DP_Hard" ,&(mat[i].m.pl_por_mises->DP_Hard)       ,&ierr);
   }
   frchk("MAT_IF",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp  = m_ifmat;
      mat[i].m.ifmat = (IFMAT*)CCACALLOC(1,sizeof(IFMAT));
      frdouble("EMOD"   ,&(mat[i].m.ifmat->emod)     ,&ierr);
      frdouble("KMOD"   ,&(mat[i].m.ifmat->kmod )    ,&ierr);
      frdouble("GMOD"   ,&(mat[i].m.ifmat->gmod )    ,&ierr);
      frdouble("DICK"   ,&(mat[i].m.ifmat->dick)     ,&ierr);
      frdouble("QMOD"   ,&(mat[i].m.ifmat->qmod)     ,&ierr);
      frdouble("DN"     ,&(mat[i].m.ifmat->deltan)   ,&ierr);
      frdouble("DT"     ,&(mat[i].m.ifmat->deltat)   ,&ierr);
      frdouble("MU"     ,&(mat[i].m.ifmat->mu)       ,&ierr);
   }
   frchk("MAT_DAM_MP",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp   = m_dam_mp;
      mat[i].m.dam_mp = (DAM_MP*)CCACALLOC(1,sizeof(DAM_MP));
      frdouble("YOUNG"   ,&(mat[i].m.dam_mp->youngs)     ,&ierr);
      frdouble("NUE"     ,&(mat[i].m.dam_mp->nue)        ,&ierr);
      frdouble("KAPPA_0" ,&(mat[i].m.dam_mp->kappa_0)    ,&ierr);
      frdouble("ALPHA"   ,&(mat[i].m.dam_mp->alpha)      ,&ierr);
      frdouble("BETA"    ,&(mat[i].m.dam_mp->beta)       ,&ierr);
   }
   frchk("MAT_DAMAGE_GE",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp      = m_damage_ge;
      mat[i].m.damage_ge = (DAMAGE_GE*)CCACALLOC(1,sizeof(DAMAGE_GE));
      frint("EQU"        ,&(mat[i].m.damage_ge->equival)    ,&ierr);
      frint("DAMT"       ,&(mat[i].m.damage_ge->damtyp)     ,&ierr);
      frdouble("CRAD"    ,&(mat[i].m.damage_ge->crad)       ,&ierr);
      frdouble("YOUNG"   ,&(mat[i].m.damage_ge->youngs)     ,&ierr);
      frdouble("NUE"     ,&(mat[i].m.damage_ge->nue)        ,&ierr);
      frdouble("KAPPA_0" ,&(mat[i].m.damage_ge->kappa_0)    ,&ierr);
      frdouble("KAPPA_M" ,&(mat[i].m.damage_ge->kappa_m)    ,&ierr);
      frdouble("ALPHA"   ,&(mat[i].m.damage_ge->alpha)      ,&ierr);
      frdouble("BETA"    ,&(mat[i].m.damage_ge->beta)       ,&ierr);
      frdouble("K_FAC"   ,&(mat[i].m.damage_ge->k_fac)      ,&ierr);
   }
   frchk("MAT_HYPER_POLYCONVEX",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp      = m_hyper_polyconvex;
      mat[i].m.hyper_polyconvex = (HYPER_POLYCONVEX*)CCACALLOC(1,sizeof(HYPER_POLYCONVEX));
      frdouble("C"       ,&(mat[i].m.hyper_polyconvex->c)       ,&ierr);
      frdouble("K1"      ,&(mat[i].m.hyper_polyconvex->k1)     ,&ierr);
      frdouble("K2"      ,&(mat[i].m.hyper_polyconvex->k2)        ,&ierr);
      frdouble("EPSILON" ,&(mat[i].m.hyper_polyconvex->epsilon)    ,&ierr);
      frdouble("GAMMA"   ,&(mat[i].m.hyper_polyconvex->gamma)    ,&ierr);
      frdouble("DENS"    ,&(mat[i].m.hyper_polyconvex->density)     ,&ierr);
   }
   /* Fourier's law of isotropic heat conduction --> heat cond. coeff. */
   frchk("MAT_Therm_Fourier_iso",&ierr);
   if (ierr==1)
   {
     mat[i].mattyp = m_th_fourier_iso;
     mat[i].m.th_fourier_iso = (TH_FOURIER_ISO*)CCACALLOC(1,sizeof(TH_FOURIER_ISO));
     frdouble("CONDUCT"  ,&(mat[i].m.th_fourier_iso->conduct)   ,&ierr);
   }
   /* Fourier's law of general heat conduction --> heat cond. matrix */
   frchk("MAT_Therm_Fourier_gen",&ierr);
   if (ierr==1)
   {
     mat[i].mattyp = m_th_fourier_gen;
     mat[i].m.th_fourier_gen = (TH_FOURIER_GEN*)CCACALLOC(1,sizeof(TH_FOURIER_GEN));
     frdouble_n("CONDUCT"  ,&(mat[i].m.th_fourier_gen->conduct[0]), 9   ,&ierr);
     if (ierr == 1)
     {
       dserror("Insufficient conduction parameters!");
     }
   }
   /* Robinson's visco-plastic material for structures */
   frchk("MAT_Struct_Robinson",&ierr);
   if (ierr==1)
   {
     mat[i].mattyp = m_vp_robinson;
     mat[i].m.vp_robinson = (VP_ROBINSON*) CCACALLOC(1,sizeof(VP_ROBINSON));
     VP_ROBINSON* robin = mat[i].m.vp_robinson;
     frdouble("YOUNG", &(robin->youngs), &ierr);
     frdouble("NUE", &(robin->possionratio), &ierr);
     frdouble("DENS", &(robin->density), &ierr);
     frdouble("THEXPANS", &(robin->thermexpans), &ierr);
     frdouble("HRDN_FACT", &(robin->hrdn_fact), &ierr);
     frdouble("HRDN_EXPO", &(robin->hrdn_expo), &ierr);
     /* SHRTHRSHLD */
     robin->shrthrshld_ipl = vp_robinson_ipl_none;
     frchk("SHRTHRSHLD", &ierr);
     if (ierr == 1)
     {
       /* polynomial interpolation */
       frchk("SHRTHRSHLD POLY", &ierr);
       if ( (ierr == 1) && (robin->shrthrshld == NULL) )
       {
         robin->shrthrshld_ipl = vp_robinson_ipl_poly;
         frint("SHRTHRSHLD POLY", &(robin->shrthrshld_n), &ierr);
         if (robin->shrthrshld_n >= 1)
         {
           DOUBLE* robin_tmp 
             = (DOUBLE*) CCACALLOC(robin->shrthrshld_n+1, sizeof(DOUBLE));
           frdouble_n("SHRTHRSHLD POLY", &(robin_tmp[0]),
                      robin->shrthrshld_n+1, &ierr);
           robin->shrthrshld 
             = (DOUBLE*) CCACALLOC(robin->shrthrshld_n, sizeof(DOUBLE));
           INT i;
           for (i=0; i<robin->shrthrshld_n; ++i)
           {
             robin->shrthrshld[i] = robin_tmp[i+1];
           }
           CCAFREE(robin_tmp);
         }
       }
       /* piecewise linear interpolation */
       frchk("SHRTHRSHLD PCWS", &ierr);
       if ( (ierr == 1) && (robin->shrthrshld == NULL) )
       {
         robin->shrthrshld_ipl = vp_robinson_ipl_pcwslnr;
         frint("SHRTHRSHLD PCWS", &(robin->shrthrshld_n), &ierr);
         if (robin->shrthrshld_n >= 1)
         {
           DOUBLE* robin_tmp 
             = (DOUBLE*) CCACALLOC(robin->shrthrshld_n+1, sizeof(DOUBLE));
           frdouble_n("SHRTHRSHLD PCWS", &(robin_tmp[0]),
                      robin->shrthrshld_n+1, &ierr);
           robin->shrthrshld 
             = (DOUBLE*) CCACALLOC(robin->shrthrshld_n, sizeof(DOUBLE));
           INT i;
           for (i=0; i<robin->shrthrshld_n; ++i)
           {
             robin->shrthrshld[i] = robin_tmp[i+1];
           }
           CCAFREE(robin_tmp);
         }
       }
       /* constant */
       if (robin->shrthrshld == NULL)
       {
         robin->shrthrshld_ipl = vp_robinson_ipl_const;
         robin->shrthrshld_n = 1;
         robin->shrthrshld = (DOUBLE*) CCACALLOC(1, sizeof(DOUBLE));
         frdouble("SHRTHRSHLD", &(robin->shrthrshld[0]), &ierr);
       }
     }
     /* RCVRY */
     robin->rcvry_ipl = vp_robinson_ipl_none;
     frchk("RCVRY", &ierr);
     if (ierr == 1)
     {
       /* polynomial interpolation */
       frchk("RCVRY POLY", &ierr);
       if ( (ierr == 1) && (robin->rcvry == NULL) )
       {
         robin->rcvry_ipl = vp_robinson_ipl_poly;
         frint("RCVRY POLY", &(robin->rcvry_n), &ierr);
         if (robin->rcvry_n >= 1)
         {
           DOUBLE* robin_tmp 
             = (DOUBLE*) CCACALLOC(robin->rcvry_n+1, sizeof(DOUBLE));
           frdouble_n("RCVRY POLY", &(robin_tmp[0]), robin->rcvry_n+1, &ierr);
           robin->rcvry = (DOUBLE*) CCACALLOC(robin->rcvry_n, sizeof(DOUBLE));
           INT i;
           for (i=0; i<robin->rcvry_n; ++i)
           {
             robin->rcvry[i] = robin_tmp[i+1];
           }
           CCAFREE(robin_tmp);
         }
       }
       /* piecewise linear interpolation */
       frchk("RCVRY PCWS", &ierr);
       if ( (ierr == 1) && (robin->rcvry == NULL) )
       {
         robin->rcvry_ipl = vp_robinson_ipl_pcwslnr;
         frint("RCVRY PCWS", &(robin->rcvry_n), &ierr);
         if (robin->rcvry_n >= 1)
         {
           DOUBLE* robin_tmp 
             = (DOUBLE*) CCACALLOC(robin->rcvry_n+1, sizeof(DOUBLE));
           frdouble_n("RCVRY PCWS", &(robin_tmp[0]), robin->rcvry_n+1, &ierr);
           robin->rcvry 
             = (DOUBLE*) CCACALLOC(robin->rcvry_n, sizeof(DOUBLE));
           INT i;
           for (i=0; i<robin->rcvry_n; ++i)
           {
             robin->rcvry[i] = robin_tmp[i+1];
           }
           CCAFREE(robin_tmp);
         }
       }
       /* constant */
       if (robin->rcvry == NULL)
       {
         robin->rcvry_ipl = vp_robinson_ipl_const;
         robin->rcvry_n = 1;
         robin->rcvry = (DOUBLE*) CCACALLOC(1, sizeof(DOUBLE));
         frdouble("RCVRY", &(robin->rcvry[0]), &ierr);
       }
     }
     frdouble("ACTV_TMPR", &(robin->actv_tmpr), &ierr);
     frdouble("ACTV_ERGY", &(robin->actv_ergy), &ierr);
     frdouble("G0", &(robin->g0), &ierr);
     frdouble("M_EXPO", &(robin->m), &ierr);
     /* BETA */
     robin->beta_ipl = vp_robinson_ipl_none;
     frchk("BETA", &ierr);
     if (ierr == 1)
     {
       /* polynomial interpolation */
       frchk("BETA POLY", &ierr);
       if ( (ierr == 1) && (robin->beta == NULL) )
       {
         robin->beta_ipl = vp_robinson_ipl_poly;
         frint("BETA POLY", &(robin->beta_n), &ierr);
         if (robin->beta_n >= 1)
         {
           DOUBLE* robin_tmp 
             = (DOUBLE*) CCACALLOC(robin->beta_n+1, sizeof(DOUBLE));
           frdouble_n("BETA POLY", &(robin_tmp[0]), robin->beta_n+1, &ierr);
           robin->beta = (DOUBLE*) CCACALLOC(robin->beta_n, sizeof(DOUBLE));
           INT i;
           for (i=0; i<robin->beta_n; ++i)
           {
             robin->beta[i] = robin_tmp[i+1];
           }
           CCAFREE(robin_tmp);
         }
       }
       /* piecewise linear interpolation */
       frchk("BETA PCWS", &ierr);
       if ( (ierr == 1) && (robin->beta == NULL) )
       {
         robin->beta_ipl = vp_robinson_ipl_pcwslnr;
         frint("BETA PCWS", &(robin->beta_n), &ierr);
         if (robin->beta_n >= 1)
         {
           DOUBLE* robin_tmp 
             = (DOUBLE*) CCACALLOC(robin->beta_n+1, sizeof(DOUBLE));
           frdouble_n("BETA PCWS", &(robin_tmp[0]), robin->beta_n+1, &ierr);
           robin->beta = (DOUBLE*) CCACALLOC(robin->beta_n, sizeof(DOUBLE));
           INT i;
           for (i=0; i<robin->beta_n; ++i)
           {
             robin->beta[i] = robin_tmp[i+1];
           }
           CCAFREE(robin_tmp);
         }
       }
       /* constant */
       if (robin->beta == NULL)
       {
         robin->beta_ipl = vp_robinson_ipl_const;
         robin->beta_n = 1;
         robin->beta = (DOUBLE*) CCACALLOC(1, sizeof(DOUBLE));
         frdouble("BETA", &(robin->beta[0]), &ierr);
       }
     }
     /* H_FACT */
     robin->h_ipl = vp_robinson_ipl_none;
     frchk("H_FACT", &ierr);
     if (ierr == 1)
     {
       /* polynomial interpolation */
       frchk("H_FACT POLY", &ierr);
       if ( (ierr == 1) && (robin->h == NULL) )
       {
         robin->h_ipl = vp_robinson_ipl_poly;
         frint("H_FACT POLY", &(robin->h_n), &ierr);
         if (robin->h_n >= 1)
         {
           DOUBLE* robin_tmp 
             = (DOUBLE*) CCACALLOC(robin->h_n+1, sizeof(DOUBLE));
           frdouble_n("H_FACT POLY", &(robin_tmp[0]), robin->h_n+1, &ierr);
           robin->h = (DOUBLE*) CCACALLOC(robin->h_n, sizeof(DOUBLE));
           INT i;
           for (i=0; i<robin->h_n; ++i)
           {
             robin->h[i] = robin_tmp[i+1];
           }
           CCAFREE(robin_tmp);
         }
       }
       /* piecewise linear interpolation */
       frchk("H_FACT PCWS", &ierr);
       if ( (ierr == 1) && (robin->h == NULL) )
       {
         robin->h_ipl = vp_robinson_ipl_pcwslnr;
         frint("H_FACT PCWS", &(robin->h_n), &ierr);
         if (robin->h_n >= 1)
         {
           DOUBLE* robin_tmp 
             = (DOUBLE*) CCACALLOC(robin->h_n+1, sizeof(DOUBLE));
           frdouble_n("H_FACT PCWS", &(robin_tmp[0]), robin->h_n+1, &ierr);
           robin->h = (DOUBLE*) CCACALLOC(robin->h_n, sizeof(DOUBLE));
           INT i;
           for (i=0; i<robin->h_n; ++i)
           {
             robin->h[i] = robin_tmp[i+1];
           }
           CCAFREE(robin_tmp);
         }
       }
       /* constant */
       if (robin->h == NULL)
       {
         robin->h_ipl = vp_robinson_ipl_const;
         robin->h_n = 1;
         robin->h = (DOUBLE*) CCACALLOC(1, sizeof(DOUBLE));
         frdouble("H_FACT", &(robin->h[0]), &ierr);
       }
     }
     /* check if allocatables have indeed been allocated */
     if (robin->beta == NULL)
     {
       dserror("Beta was not found!");
     }
     if (robin->shrthrshld == NULL)
     {
       dserror("Shear threshold was not found!");
     }
     if (robin->rcvry == NULL)
     {
       dserror("Recovery factor was not found!");
     }
     if (robin->h == NULL)
     {
       dserror("H was not found!");
     }
   }
   /*multi layer material */
   frchk("MAT_Multilayer",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_multi_layer;
      mat[i].m.multi_layer = (MULTI_LAYER*)CCACALLOC(1,sizeof(MULTI_LAYER));
      /*read number of kinamtic layer*/
      frint("NUM_KLAY"   ,&(mat[i].m.multi_layer->num_klay)        ,&ierr);

      /*allocate memory for the different kinematic layers*/
      num_klay = mat[i].m.multi_layer->num_klay;
      mat[i].m.multi_layer->klayhgt=(DOUBLE*)CCACALLOC(num_klay,sizeof(DOUBLE));
      mat[i].m.multi_layer->kinlay=(KINLAY*)CCACALLOC(num_klay,sizeof(KINLAY));

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
        actlay->mlayhgt=(DOUBLE*)CCACALLOC(num_mlay,sizeof(DOUBLE));
        actlay->mmatID=(INT*)CCACALLOC(num_mlay,sizeof(INT));
        actlay->phi=(DOUBLE*)CCACALLOC(num_mlay,sizeof(DOUBLE));
        actlay->rot=(INT*)CCACALLOC(num_mlay,sizeof(INT));

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
char buffer[50];
INT  counter=0;
INT  check = 0;
INT  ierr;
INT  i;
/*----------------------------------------------------------------------*/
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
if (frfind("--MULTILAYER MATERIALS")==1)
{
  frread();
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    frchar("MULTIMAT",buffer,&ierr);
    if (ierr==1) counter++;
    frread();
  }
}
if (counter == 0) goto end;       /* no multilayer material set */

/*--------------------------------------------------- allocate MULTIMAT */
multimat = (MULTIMAT*)CCACALLOC(counter,sizeof(MULTIMAT));
/*----------------------------------------------------------------------*/
if (frfind("--MULTILAYER MATERIALS")==0) goto end;
frread();
i=0;
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frint("MULTIMAT",&(multimat[i].Id),&ierr);
   if (ierr==0)
   {
      frread();
      continue;
   }

   frchk("MAT_Struct_StVenantKirchhoff",&ierr);
   if (ierr==1)
   {
      multimat[i].mattyp = m_stvenant;
      multimat[i].m.stvenant = (STVENANT*)CCACALLOC(1,sizeof(STVENANT));
      frdouble("YOUNG"  ,&(multimat[i].m.stvenant->youngs)      ,&ierr);
      frdouble("NUE"    ,&(multimat[i].m.stvenant->possionratio),&ierr);
      frdouble("DENS"   ,&(multimat[i].m.stvenant->density)     ,&ierr);
   }
   frchk("MAT_Struct_NeoHooke",&ierr);
   if (ierr==1)
   {
      multimat[i].mattyp = m_neohooke;
      multimat[i].m.neohooke = (NEO_HOOKE*)CCACALLOC(1,sizeof(NEO_HOOKE));
      frdouble("YOUNG"   ,&(multimat[i].m.neohooke->youngs)        ,&ierr);
      frdouble("NUE"     ,&(multimat[i].m.neohooke->possionratio)  ,&ierr);
      frdouble("DENS"    ,&(multimat[i].m.neohooke->density)       ,&ierr);
   }
   frchk("MAT_Struct_Orthotropic",&ierr);
   if (ierr==1)
   {
      multimat[i].mattyp = m_el_orth;
      multimat[i].m.el_orth = (EL_ORTH*)CCACALLOC(1,sizeof(EL_ORTH));
      frdouble("EMOD1"   ,&(multimat[i].m.el_orth->emod1)        ,&ierr);
      frdouble("EMOD2"   ,&(multimat[i].m.el_orth->emod2)        ,&ierr);
      frdouble("EMOD3"   ,&(multimat[i].m.el_orth->emod3)        ,&ierr);
      frdouble("GMOD12"  ,&(multimat[i].m.el_orth->gmod12)       ,&ierr);
      frdouble("GMOD13"  ,&(multimat[i].m.el_orth->gmod13)       ,&ierr);
      frdouble("GMOD23"  ,&(multimat[i].m.el_orth->gmod23)       ,&ierr);
      frdouble("XNUE12"  ,&(multimat[i].m.el_orth->xnue12)       ,&ierr);
      frdouble("XNUE13"  ,&(multimat[i].m.el_orth->xnue13)       ,&ierr);
      frdouble("XNUE23"  ,&(multimat[i].m.el_orth->xnue23)       ,&ierr);
   }
   frchk("MAT_MisesPlastic",&ierr);
   if (ierr==1)
   {
      multimat[i].mattyp = m_pl_mises;
      multimat[i].m.pl_mises = (PL_MISES*)CCACALLOC(1,sizeof(PL_MISES));
      frdouble("YOUNG",&(multimat[i].m.pl_mises->youngs)        ,&ierr);
      frdouble("NUE"  ,&(multimat[i].m.pl_mises->possionratio)  ,&ierr);
      frdouble("ALFAT",&(multimat[i].m.pl_mises->ALFAT)         ,&ierr);
      frdouble("Sigy" ,&(multimat[i].m.pl_mises->Sigy)          ,&ierr);
      multimat[i].m.pl_mises->Hard = 0.;
      multimat[i].m.pl_mises->GF   = 0.;
      multimat[i].m.pl_mises->betah= 1.;
      frdouble("Hard" ,&(multimat[i].m.pl_mises->Hard)          ,&ierr);
      frdouble("GF"   ,&(multimat[i].m.pl_mises->GF)            ,&ierr);
      frdouble("BETAH",&(multimat[i].m.pl_mises->betah)         ,&ierr);
   }
   frchk("MAT_DP_Plastic",&ierr);
   if (ierr==1)
   {
      /*write a warning to use an unsymmetric solver*/
      printf("|---------------------------------------------------------------------------------------| \n");
      printf("|    WARNING in input_material.c:    ==============================================     | \n");
      printf("|    WARNING in input_material.c:    MAT_DP_Plastic -> use an UNsymmetric solver if     | \n");
      printf("|    WARNING in input_material.c:    trial stresses could be in the apex region and     | \n");
      printf("|    WARNING in input_material.c:    hardening law not fully kinematic (betah > 0)      | \n");
      printf("|    WARNING in input_material.c:    ==============================================     | \n");
      printf("|---------------------------------------------------------------------------------------| \n");

      multimat[i].mattyp = m_pl_dp;
      multimat[i].m.pl_dp = (PL_DP*)CCACALLOC(1,sizeof(PL_DP));
      frdouble("YOUNG",&(multimat[i].m.pl_dp->youngs)        ,&ierr);
      frdouble("NUE"  ,&(multimat[i].m.pl_dp->possionratio)  ,&ierr);
      frdouble("ALFAT",&(multimat[i].m.pl_dp->ALFAT)         ,&ierr);
      frdouble("Sigy" ,&(multimat[i].m.pl_dp->Sigy)          ,&ierr);
      frdouble("PHI"  ,&(multimat[i].m.pl_dp->PHI)           ,&ierr);
      multimat[i].m.pl_dp->Hard = 0.;
      multimat[i].m.pl_dp->GF   = 0.;
      multimat[i].m.pl_dp->betah= 1.;
      frdouble("Hard" ,&(multimat[i].m.pl_dp->Hard)          ,&ierr);
      frdouble("GF"   ,&(multimat[i].m.pl_dp->GF)            ,&ierr);
      frdouble("BETAH",&(multimat[i].m.pl_dp->betah)         ,&ierr);
   }
   frchk("MAT_HoffPlastic",&ierr);
   if (ierr==1)
   {
      /*write a warning to use an unsymmetric solver*/
      printf("|---------------------------------------------------------------------------------------| \n");
      printf("|    WARNING in input_material.c:    ============================================       | \n");
      printf("|    WARNING in input_material.c:    MAT_HoffPlastic -> use an UNsymmetric solver       | \n");
      printf("|    WARNING in input_material.c:    ============================================       | \n");
      printf("|---------------------------------------------------------------------------------------| \n");

      multimat[i].mattyp = m_pl_hoff;
      multimat[i].m.pl_hoff = (PL_HOFF*)CCACALLOC(1,sizeof(PL_HOFF));
      frdouble("EMOD1"   ,&(multimat[i].m.pl_hoff->emod1)        ,&ierr);
      frdouble("EMOD2"   ,&(multimat[i].m.pl_hoff->emod2)        ,&ierr);
      frdouble("EMOD3"   ,&(multimat[i].m.pl_hoff->emod3)        ,&ierr);
      frdouble("GMOD12"  ,&(multimat[i].m.pl_hoff->gmod12)       ,&ierr);
      frdouble("GMOD13"  ,&(multimat[i].m.pl_hoff->gmod13)       ,&ierr);
      frdouble("GMOD23"  ,&(multimat[i].m.pl_hoff->gmod23)       ,&ierr);
      frdouble("XNUE12"  ,&(multimat[i].m.pl_hoff->xnue12)       ,&ierr);
      frdouble("XNUE13"  ,&(multimat[i].m.pl_hoff->xnue13)       ,&ierr);
      frdouble("XNUE23"  ,&(multimat[i].m.pl_hoff->xnue23)       ,&ierr);
      frdouble("UNIAX"   ,&(multimat[i].m.pl_hoff->uniax)        ,&ierr);

      frread(); /*new line*/

      frdouble("S11T"    ,&(multimat[i].m.pl_hoff->s11T)         ,&ierr);
      frdouble("S11C"    ,&(multimat[i].m.pl_hoff->s11C)         ,&ierr);
      frdouble("S22T"    ,&(multimat[i].m.pl_hoff->s22T)         ,&ierr);
      frdouble("S22C"    ,&(multimat[i].m.pl_hoff->s22C)         ,&ierr);
      frdouble("S33T"    ,&(multimat[i].m.pl_hoff->s33T)         ,&ierr);
      frdouble("S33C"    ,&(multimat[i].m.pl_hoff->s33C)         ,&ierr);
      frdouble("S12"     ,&(multimat[i].m.pl_hoff->s12)          ,&ierr);
      frdouble("S23"     ,&(multimat[i].m.pl_hoff->s23)          ,&ierr);
      frdouble("S13"     ,&(multimat[i].m.pl_hoff->s13)          ,&ierr);

      frread(); /*new line*/

      frdouble("SH11T"   ,&(multimat[i].m.pl_hoff->sh11T)        ,&ierr);
      frdouble("SH11C"   ,&(multimat[i].m.pl_hoff->sh11C)        ,&ierr);
      frdouble("SH22T"   ,&(multimat[i].m.pl_hoff->sh22T)        ,&ierr);
      frdouble("SH22C"   ,&(multimat[i].m.pl_hoff->sh22C)        ,&ierr);
      frdouble("SH33T"   ,&(multimat[i].m.pl_hoff->sh33T)        ,&ierr);
      frdouble("SH33C"   ,&(multimat[i].m.pl_hoff->sh33C)        ,&ierr);
      frdouble("SH12"    ,&(multimat[i].m.pl_hoff->sh12)         ,&ierr);
      frdouble("SH23"    ,&(multimat[i].m.pl_hoff->sh23)         ,&ierr);
      frdouble("SH13"    ,&(multimat[i].m.pl_hoff->sh13)         ,&ierr);

      frread(); /*new line*/

      frdouble("HA11T"   ,&(multimat[i].m.pl_hoff->ha11T)        ,&ierr);
      frdouble("HA11C"   ,&(multimat[i].m.pl_hoff->ha11C)        ,&ierr);
      frdouble("HA22T"   ,&(multimat[i].m.pl_hoff->ha22T)        ,&ierr);
      frdouble("HA22C"   ,&(multimat[i].m.pl_hoff->ha22C)        ,&ierr);
      frdouble("HA33T"   ,&(multimat[i].m.pl_hoff->ha33T)        ,&ierr);
      frdouble("HA33C"   ,&(multimat[i].m.pl_hoff->ha33C)        ,&ierr);
      frdouble("HA12"    ,&(multimat[i].m.pl_hoff->ha12)         ,&ierr);
      frdouble("HA23"    ,&(multimat[i].m.pl_hoff->ha23)         ,&ierr);
      frdouble("HA13"    ,&(multimat[i].m.pl_hoff->ha13)         ,&ierr);
   }
   frchk("MAT_ConcretePlastic",&ierr);
   if (ierr==1)
   {
      /*write a warning to use an unsymmetric solver*/
      printf("|---------------------------------------------------------------------------------------| \n");
      printf("|    WARNING in input_material.c:    ================================================   | \n");
      printf("|    WARNING in input_material.c:    MAT_ConcretePlastic -> use an UNsymmetric solver   | \n");
      printf("|    WARNING in input_material.c:    if trial stresses could be in the apex region      | \n");
      printf("|    WARNING in input_material.c:    ================================================   | \n");
      printf("|---------------------------------------------------------------------------------------| \n");

      multimat[i].mattyp = m_pl_epc;
      multimat[i].m.pl_epc = (PL_EPC*)CCACALLOC(1,sizeof(PL_EPC));

      /* initialize */
      multimat[i].m.pl_epc->dfac = 0.;
      multimat[i].m.pl_epc->gamma1 = 3.;
      multimat[i].m.pl_epc->gamma2 = 6./5.;
      multimat[i].m.pl_epc->gamma3 = 1./3.;
      multimat[i].m.pl_epc->gamma4 = 4./3.;
      multimat[i].m.pl_epc->maxreb = 0;

      /* concrete */
      frdouble("YOUNG"   ,&(multimat[i].m.pl_epc->youngs      )        ,&ierr);
      frdouble("NUE"     ,&(multimat[i].m.pl_epc->possionratio)        ,&ierr);
      frdouble("FTM"     ,&(multimat[i].m.pl_epc->ftm         )        ,&ierr);
      frdouble("FCM"     ,&(multimat[i].m.pl_epc->fcm         )        ,&ierr);
      frdouble("GT"      ,&(multimat[i].m.pl_epc->gt          )        ,&ierr);
      frdouble("GC"      ,&(multimat[i].m.pl_epc->gc          )        ,&ierr);
      frdouble("GAMMA1"  ,&(multimat[i].m.pl_epc->gamma1      )        ,&ierr);
      frdouble("GAMMA2"  ,&(multimat[i].m.pl_epc->gamma2      )        ,&ierr);
      frdouble("GAMMA3"  ,&(multimat[i].m.pl_epc->gamma3      )        ,&ierr);
      frdouble("GAMMA4"  ,&(multimat[i].m.pl_epc->gamma4      )        ,&ierr);
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
