/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/wiechert
            089 - 28915303
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "global_inp_control2.H"

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}
/*----------------------------------------------------------------------*
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 | input of materials                                     m.gee 4/01    |
 *----------------------------------------------------------------------*/
void DRT::Problem::ReadMaterial()
{
_MATERIAL  localmat;
INT  ierr, ierralloc, int_dummy;
INT  i, j, k, ncm, num_klay, num_mlay;
struct    _KINLAY *actlay;           /*actual kinematic layer -> shell9 */
DOUBLE    klay_sum;                  /*total hight or shell9*/
DOUBLE    mlay_sum;                  /*hight of a kinematic layer*/

/*----------------------------------------------------------------------*/
if (frfind("--MATERIALS")==0) dserror("frfind: MATERIALS is not in input file");
frread();
i=0;
while(strncmp(allfiles.actplace,"------",6)!=0)
{
  localmat.m.fluid = NULL;

   frint("MAT",&(localmat.Id),&ierr);
   frchk("MAT_fluid",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_fluid;
      localmat.m.fluid = new _FLUID();
      frdouble("VISCOSITY",&(localmat.m.fluid->viscosity),&ierr);
      frdouble("DENS"  ,&(localmat.m.fluid->density)  ,&ierr);
      frdouble("GAMMA",&(localmat.m.fluid->gamma)  ,&ierr);
   }
   frchk("MAT_condif",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_condif;
      localmat.m.condif = new _CONDIF();
      frdouble("DIFFUSIVITY",&(localmat.m.condif->diffusivity),&ierr);
   }
   frchk("MAT_Struct_StVenantKirchhoff",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_stvenant;
      localmat.m.stvenant = new _STVENANT();
      frdouble("YOUNG"  ,&(localmat.m.stvenant->youngs)      ,&ierr);
      frdouble("NUE"    ,&(localmat.m.stvenant->possionratio),&ierr);
      frdouble("DENS"   ,&(localmat.m.stvenant->density)     ,&ierr);
      frdouble("THEXPANS",&(localmat.m.stvenant->thermexpans) ,&ierr);
   }
   frchk("MAT_Struct_Orthotropic",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_el_orth;
      localmat.m.el_orth = new EL_ORTH();
      frdouble("EMOD1"   ,&(localmat.m.el_orth->emod1)        ,&ierr);
      frdouble("EMOD2"   ,&(localmat.m.el_orth->emod2)        ,&ierr);
      frdouble("EMOD3"   ,&(localmat.m.el_orth->emod3)        ,&ierr);
      frdouble("GMOD12"  ,&(localmat.m.el_orth->gmod12)       ,&ierr);
      frdouble("GMOD13"  ,&(localmat.m.el_orth->gmod13)       ,&ierr);
      frdouble("GMOD23"  ,&(localmat.m.el_orth->gmod23)       ,&ierr);
      frdouble("XNUE12"  ,&(localmat.m.el_orth->xnue12)       ,&ierr);
      frdouble("XNUE13"  ,&(localmat.m.el_orth->xnue13)       ,&ierr);
      frdouble("XNUE23"  ,&(localmat.m.el_orth->xnue23)       ,&ierr);
   }
   frchk("MAT_Struct_STVENPOR",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_stvenpor;
      localmat.m.stvenpor = new STVENPOR();
      frdouble("YOUNG"  ,&(localmat.m.stvenpor->youngs)        ,&ierr);
      frdouble("NUE"    ,&(localmat.m.stvenpor->possionratio)  ,&ierr);
      frdouble("DENS"   ,&(localmat.m.stvenpor->density)       ,&ierr);
      frdouble("REFDENS",&(localmat.m.stvenpor->refdens)       ,&ierr);
      frdouble("EXPO"   ,&(localmat.m.stvenpor->exponent)      ,&ierr);
   }
   frchk("MAT_Struct_NeoHooke",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_neohooke;
      localmat.m.neohooke = new NEO_HOOKE();
      frdouble("YOUNG",&(localmat.m.neohooke->youngs)        ,&ierr);
      frdouble("NUE"  ,&(localmat.m.neohooke->possionratio)  ,&ierr);
      frdouble("DENS",&(localmat.m.neohooke->density)        ,&ierr);
   }
   frchk("MAT_MFOC",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_mfoc;
      localmat.m.mfoc = new MFOC();
      if (localmat.m.mfoc==NULL) dserror("Allocation of Open Cell foam material failed");
      frdouble("Es"     ,&(localmat.m.mfoc->es)          ,&ierr);
      frdouble("pr"     ,&(localmat.m.mfoc->pr)          ,&ierr);
      frdouble("dens"   ,&(localmat.m.mfoc->dens)        ,&ierr);  /* por. density */
      frdouble("denss"  ,&(localmat.m.mfoc->denss)       ,&ierr);  /* ref. density */
      frdouble("oce"    ,&(localmat.m.mfoc->oce)         ,&ierr);
      frdouble("ocf"    ,&(localmat.m.mfoc->ocf)         ,&ierr);
      frdouble("densmin",&(localmat.m.mfoc->denmin)      ,&ierr);
      frdouble("densmax",&(localmat.m.mfoc->denmax)      ,&ierr);
   }
   frchk("MAT_MFCC",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_mfcc;
      localmat.m.mfcc = new MFCC();
      if (localmat.m.mfoc==NULL) dserror("Alloocation of Closed Cell foam material failed");
      frdouble("Es"     ,&(localmat.m.mfcc->es)          ,&ierr);
      frdouble("pr"     ,&(localmat.m.mfcc->pr)          ,&ierr);
      frdouble("dens"   ,&(localmat.m.mfcc->dens)        ,&ierr);  /* por. density */
      frdouble("denss"  ,&(localmat.m.mfcc->denss)       ,&ierr);  /* ref. density */
      frdouble("cce"    ,&(localmat.m.mfcc->cce)         ,&ierr);
      frdouble("ccf"    ,&(localmat.m.mfcc->ccf)         ,&ierr);
      frdouble("densmin",&(localmat.m.mfcc->denmin)      ,&ierr);
      frdouble("densmax",&(localmat.m.mfcc->denmax)      ,&ierr);
   }
   frchk("MAT_NeoHMFCC",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_nhmfcc;
      localmat.m.nhmfcc = new NHMFCC();
      if (localmat.m.nhmfcc==NULL) dserror("Alloocation of foam material failed");
      frdouble("Es"     ,&(localmat.m.nhmfcc->es)          ,&ierr);
      frdouble("pr"     ,&(localmat.m.nhmfcc->pr)          ,&ierr);
      frdouble("dens"   ,&(localmat.m.nhmfcc->dens)        ,&ierr);  /* por. density */
      frdouble("denss"  ,&(localmat.m.nhmfcc->denss)       ,&ierr);  /* ref. density */
      frdouble("cce"    ,&(localmat.m.nhmfcc->cce)         ,&ierr);
      frdouble("ccf"    ,&(localmat.m.nhmfcc->ccf)         ,&ierr);
      frdouble("densmin",&(localmat.m.nhmfcc->denmin)      ,&ierr);
      frdouble("densmax",&(localmat.m.nhmfcc->denmax)      ,&ierr);
   }
   frchk("MAT_Struct_Ogden",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_compogden;
      localmat.m.compogden = new COMPOGDEN();
      frdouble("NUE"  ,&(localmat.m.compogden->nue)     ,&ierr);
      frdouble("BETA" ,&(localmat.m.compogden->beta)    ,&ierr);
      frdouble("ALFA1",&(localmat.m.compogden->alfap[0]),&ierr);
      frdouble("ALFA2",&(localmat.m.compogden->alfap[1]),&ierr);
      frdouble("ALFA3",&(localmat.m.compogden->alfap[2]),&ierr);
      frdouble("NU1"  ,&(localmat.m.compogden->mup[0])  ,&ierr);
      frdouble("NU2"  ,&(localmat.m.compogden->mup[1])  ,&ierr);
      frdouble("NU3"  ,&(localmat.m.compogden->mup[2])  ,&ierr);
      frdouble("DENS" ,&(localmat.m.compogden->density) ,&ierr);
   }
   frchk("MAT_Struct_Viscohyper",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_viscohyper;
      localmat.m.viscohyper = new VISCOHYPER();
      frdouble("NUE"  ,&(localmat.m.viscohyper->nue)     ,&ierr);
      frdouble("BETA" ,&(localmat.m.viscohyper->beta)    ,&ierr);
      frdouble("ALFA1",&(localmat.m.viscohyper->alfap[0]),&ierr);
      frdouble("ALFA2",&(localmat.m.viscohyper->alfap[1]),&ierr);
      frdouble("ALFA3",&(localmat.m.viscohyper->alfap[2]),&ierr);
      frdouble("NU1"  ,&(localmat.m.viscohyper->mup[0])  ,&ierr);
      frdouble("NU2"  ,&(localmat.m.viscohyper->mup[1])  ,&ierr);
      frdouble("NU3"  ,&(localmat.m.viscohyper->mup[2])  ,&ierr);
      frdouble("DENS" ,&(localmat.m.viscohyper->density) ,&ierr);
      frint   ("NMAXW",&(localmat.m.viscohyper->nmaxw)   ,&ierr);
      frdouble("TAU1" ,&(localmat.m.viscohyper->tau[0])  ,&ierr);
      frdouble("TAU2" ,&(localmat.m.viscohyper->tau[1])  ,&ierr);
      frdouble("TAU3" ,&(localmat.m.viscohyper->tau[2])  ,&ierr);
      frdouble("TAU4" ,&(localmat.m.viscohyper->tau[3])  ,&ierr);
      frdouble("BETA1",&(localmat.m.viscohyper->betas[0]),&ierr);
      frdouble("BETA2",&(localmat.m.viscohyper->betas[1]),&ierr);
      frdouble("BETA3",&(localmat.m.viscohyper->betas[2]),&ierr);
      frdouble("BETA4",&(localmat.m.viscohyper->betas[3]),&ierr);
   }
   frchk("MAT_3DMisesPlastic",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_pl_mises_3D;
      localmat.m.pl_mises = new PL_MISES();
      frdouble("YOUNG",&(localmat.m.pl_mises->youngs)        ,&ierr);
      frdouble("NUE"  ,&(localmat.m.pl_mises->possionratio)  ,&ierr);
      frdouble("ALFAT",&(localmat.m.pl_mises->ALFAT)         ,&ierr);
      frdouble("Sigy" ,&(localmat.m.pl_mises->Sigy)          ,&ierr);
      localmat.m.pl_mises->Hard = 0.;
      localmat.m.pl_mises->GF   = 0.;
      localmat.m.pl_mises->betah= 1.;
      frdouble("Hard" ,&(localmat.m.pl_mises->Hard)          ,&ierr);
      frdouble("GF"   ,&(localmat.m.pl_mises->GF)            ,&ierr);
      frdouble("BETAH",&(localmat.m.pl_mises->betah)         ,&ierr);
   }
   frchk("MAT_MisesPlastic",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_pl_mises;
      localmat.m.pl_mises = new PL_MISES();
      frdouble("YOUNG",&(localmat.m.pl_mises->youngs)        ,&ierr);
      frdouble("NUE"  ,&(localmat.m.pl_mises->possionratio)  ,&ierr);
      frdouble("ALFAT",&(localmat.m.pl_mises->ALFAT)         ,&ierr);
      frdouble("Sigy" ,&(localmat.m.pl_mises->Sigy)          ,&ierr);
      localmat.m.pl_mises->Hard = 0.;
      localmat.m.pl_mises->GF   = 0.;
      localmat.m.pl_mises->betah= 1.;
      frdouble("Hard" ,&(localmat.m.pl_mises->Hard)          ,&ierr);
      frdouble("GF"   ,&(localmat.m.pl_mises->GF)            ,&ierr);
      frdouble("BETAH",&(localmat.m.pl_mises->betah)         ,&ierr);
   }
   frchk("MAT_Damage",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_damage;
      localmat.m.damage = new DAMAGE();
      if (localmat.m.damage==NULL) dserror("Allocation of DAMAGE material failed");
      frdouble("YOUNG",&(localmat.m.damage->youngs)        ,&ierr);
      frdouble("NUE"  ,&(localmat.m.damage->possionratio)  ,&ierr);
      frint("Equival" ,&(localmat.m.damage->Equival)       ,&ierr);
      frint("Damtyp"  ,&(localmat.m.damage->Damtyp)        ,&ierr);
      frdouble("Kappa_0",&(localmat.m.damage->Kappa_0)     ,&ierr);
      frdouble("Kappa_m",&(localmat.m.damage->Kappa_m)     ,&ierr);
      frdouble("Alpha",  &(localmat.m.damage->Alpha)       ,&ierr);
      frdouble("Beta" ,  &(localmat.m.damage->Beta)        ,&ierr);
      frdouble("k_fac" ,  &(localmat.m.damage->k_fac)      ,&ierr);
   }
   frchk("MAT_FoamPlastic",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_pl_foam;
      localmat.m.pl_foam = new PL_FOAM();
      frdouble("YOUNG",&(localmat.m.pl_foam->youngs)        ,&ierr);
      frdouble("NUE"  ,&(localmat.m.pl_foam->possionratio)  ,&ierr);
      frdouble("ALFAT",&(localmat.m.pl_foam->ALFAT)         ,&ierr);
      frdouble("Sigy" ,&(localmat.m.pl_foam->Sigy)          ,&ierr);
      localmat.m.pl_foam->Hard = 0.;
      localmat.m.pl_foam->GF   = 0.;
      frdouble("Hard" ,&(localmat.m.pl_foam->Hard)          ,&ierr);
      frdouble("GF"   ,&(localmat.m.pl_foam->GF)            ,&ierr);
   }
   frchk("MAT_DP_Plastic",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_pl_dp;
      localmat.m.pl_dp = new PL_DP();
      frdouble("YOUNG",&(localmat.m.pl_dp->youngs)        ,&ierr);
      frdouble("NUE"  ,&(localmat.m.pl_dp->possionratio)  ,&ierr);
      frdouble("ALFAT",&(localmat.m.pl_dp->ALFAT)         ,&ierr);
      frdouble("Sigy" ,&(localmat.m.pl_dp->Sigy)          ,&ierr);
      frdouble("PHI"  ,&(localmat.m.pl_dp->PHI)           ,&ierr);
      localmat.m.pl_dp->Hard = 0.;
      localmat.m.pl_dp->GF   = 0.;
      localmat.m.pl_dp->betah= 1.;
      frdouble("Hard" ,&(localmat.m.pl_dp->Hard)          ,&ierr);
      frdouble("GF"   ,&(localmat.m.pl_dp->GF)            ,&ierr);
      frdouble("BETAH",&(localmat.m.pl_dp->betah)         ,&ierr);
   }
   frchk("MAT_ConcretePlastic",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_pl_epc;
      localmat.m.pl_epc = new PL_EPC();
      /* initialize */
      localmat.m.pl_epc->gamma1 = 3.;
      localmat.m.pl_epc->gamma2 = 6./5.;


      frdouble("DENS"    ,&(localmat.m.pl_epc->dens        )        ,&ierr);
      /* concrete */
      frdouble("YOUNG"   ,&(localmat.m.pl_epc->youngs      )        ,&ierr);
      frdouble("NUE"     ,&(localmat.m.pl_epc->possionratio)        ,&ierr);
      frdouble("ALFAT"   ,&(localmat.m.pl_epc->alfat       )        ,&ierr);
      frdouble("XSI"     ,&(localmat.m.pl_epc->xsi         )        ,&ierr);
      frdouble("Sigy"    ,&(localmat.m.pl_epc->sigy        )        ,&ierr);
      frread();
      frdouble("FTM"     ,&(localmat.m.pl_epc->ftm         )        ,&ierr);
      frdouble("FCM"     ,&(localmat.m.pl_epc->fcm         )        ,&ierr);
      frdouble("GT"      ,&(localmat.m.pl_epc->gt          )        ,&ierr);
      frdouble("GC"      ,&(localmat.m.pl_epc->gc          )        ,&ierr);
      frdouble("GAMMA1"  ,&(localmat.m.pl_epc->gamma1      )        ,&ierr);
      if(localmat.m.pl_epc->gamma1<1.)localmat.m.pl_epc->gamma1=3.;
      frdouble("GAMMA2"  ,&(localmat.m.pl_epc->gamma2      )        ,&ierr);
     /* tension stiffening - next line in input file!*/
      frread();
      frint(   "NSTIFF"  ,&(localmat.m.pl_epc->nstiff      )        ,&ierr);
      /* number of rebars - next line in input file! */
      frread();
      localmat.m.pl_epc->maxreb = 0;
      frint(   "MAXREB"   ,&(localmat.m.pl_epc->maxreb     )        ,&ierr);
      /* allocate memory */
      ncm       = localmat.m.pl_epc->maxreb;
      ierralloc = 0;
      localmat.m.pl_epc->rebar     = new INT[ncm];
      localmat.m.pl_epc->reb_area  = new DOUBLE[ncm];
      localmat.m.pl_epc->reb_ang   = new DOUBLE[ncm];
      localmat.m.pl_epc->reb_so    = new DOUBLE[ncm];
      localmat.m.pl_epc->reb_ds    = new DOUBLE[ncm];
      localmat.m.pl_epc->reb_rgamma= new DOUBLE[ncm];
      localmat.m.pl_epc->reb_dens  = new DOUBLE[ncm];
      localmat.m.pl_epc->reb_alfat = new DOUBLE[ncm];
      localmat.m.pl_epc->reb_emod  = new DOUBLE[ncm];
      localmat.m.pl_epc->reb_rebnue= new DOUBLE[ncm];
      localmat.m.pl_epc->reb_sigy  = new DOUBLE[ncm];
      localmat.m.pl_epc->reb_hard  = new DOUBLE[ncm];

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
        frint(   "REBAR"   ,&(localmat.m.pl_epc->rebar[j]     ),&ierr);
        frdouble("REBAREA" ,&(localmat.m.pl_epc->reb_area[j]  ),&ierr);
        frdouble("REBANG"  ,&(localmat.m.pl_epc->reb_ang[j]   ),&ierr);
        frdouble("REBSO"   ,&(localmat.m.pl_epc->reb_so[j]    ),&ierr);
        frdouble("REBDS"   ,&(localmat.m.pl_epc->reb_ds[j]    ),&ierr);
        frdouble("REBGAMMA",&(localmat.m.pl_epc->reb_rgamma[j]),&ierr);
        frread();
        frdouble("REBDENS" ,&(localmat.m.pl_epc->reb_dens[j]  ),&ierr);
        frdouble("REBALFAT",&(localmat.m.pl_epc->reb_alfat[j] ),&ierr);
        frdouble("REBEMOD" ,&(localmat.m.pl_epc->reb_emod[j]  ),&ierr);
        frdouble("REBNUE"  ,&(localmat.m.pl_epc->reb_rebnue[j]),&ierr);
        frread();
        frdouble("REBSIGY" ,&(localmat.m.pl_epc->reb_sigy[j]  ),&ierr);
        frdouble("REBHARD" ,&(localmat.m.pl_epc->reb_hard[j]  ),&ierr);
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

      localmat.mattyp = m_pl_epc3D;
      localmat.m.pl_epc = new PL_EPC();

      /* initialize */
      localmat.m.pl_epc->dfac = 0.;
      localmat.m.pl_epc->gamma1 = 3.;
      localmat.m.pl_epc->gamma2 = 6./5.;
      localmat.m.pl_epc->gamma3 = 1./3.;
      localmat.m.pl_epc->gamma4 = 4./3.;
      localmat.m.pl_epc->maxreb = 0;

      /* concrete */
      frdouble("YOUNG"   ,&(localmat.m.pl_epc->youngs      )        ,&ierr);
      frdouble("NUE"     ,&(localmat.m.pl_epc->possionratio)        ,&ierr);
      frdouble("FTM"     ,&(localmat.m.pl_epc->ftm         )        ,&ierr);
      frdouble("FCM"     ,&(localmat.m.pl_epc->fcm         )        ,&ierr);
      frdouble("GT"      ,&(localmat.m.pl_epc->gt          )        ,&ierr);
      frdouble("GC"      ,&(localmat.m.pl_epc->gc          )        ,&ierr);
      frdouble("DFAC"    ,&(localmat.m.pl_epc->dfac        )        ,&ierr);
      frdouble("GAMMA1"  ,&(localmat.m.pl_epc->gamma1      )        ,&ierr);
      frdouble("GAMMA2"  ,&(localmat.m.pl_epc->gamma2      )        ,&ierr);
      frdouble("GAMMA3"  ,&(localmat.m.pl_epc->gamma3      )        ,&ierr);
      frdouble("GAMMA4"  ,&(localmat.m.pl_epc->gamma4      )        ,&ierr);
   }
   frchk("MAT_Porous_MisesPlastic",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_pl_por_mises;
      localmat.m.pl_por_mises = new PL_POR_MISES();
      frdouble("YOUNG"   ,&(localmat.m.pl_por_mises->youngs)        ,&ierr);
      frdouble("DP_YM"   ,&(localmat.m.pl_por_mises->DP_YM )        ,&ierr);
      frdouble("NUE"     ,&(localmat.m.pl_por_mises->possionratio)  ,&ierr);
      frdouble("ALFAT"   ,&(localmat.m.pl_por_mises->ALFAT)         ,&ierr);
      frdouble("Sigy"    ,&(localmat.m.pl_por_mises->Sigy)          ,&ierr);
      frdouble("DP_Sigy" ,&(localmat.m.pl_por_mises->DP_Sigy)       ,&ierr);
      frdouble("Hard"    ,&(localmat.m.pl_por_mises->Hard)          ,&ierr);
      frdouble("DP_Hard" ,&(localmat.m.pl_por_mises->DP_Hard)       ,&ierr);
   }
   frchk("MAT_IF",&ierr);
   if (ierr==1)
   {
      localmat.mattyp  = m_ifmat;
      localmat.m.ifmat = new IFMAT();
      frdouble("EMOD"   ,&(localmat.m.ifmat->emod)     ,&ierr);
      frdouble("KMOD"   ,&(localmat.m.ifmat->kmod )    ,&ierr);
      frdouble("GMOD"   ,&(localmat.m.ifmat->gmod )    ,&ierr);
      frdouble("DICK"   ,&(localmat.m.ifmat->dick)     ,&ierr);
      frdouble("QMOD"   ,&(localmat.m.ifmat->qmod)     ,&ierr);
      frdouble("DN"     ,&(localmat.m.ifmat->deltan)   ,&ierr);
      frdouble("DT"     ,&(localmat.m.ifmat->deltat)   ,&ierr);
      frdouble("MU"     ,&(localmat.m.ifmat->mu)       ,&ierr);
   }
   frchk("MAT_DAM_MP",&ierr);
   if (ierr==1)
   {
      localmat.mattyp   = m_dam_mp;
      localmat.m.dam_mp = new DAM_MP();
      frdouble("YOUNG"   ,&(localmat.m.dam_mp->youngs)     ,&ierr);
      frdouble("NUE"     ,&(localmat.m.dam_mp->nue)        ,&ierr);
      frdouble("KAPPA_0" ,&(localmat.m.dam_mp->kappa_0)    ,&ierr);
      frdouble("ALPHA"   ,&(localmat.m.dam_mp->alpha)      ,&ierr);
      frdouble("BETA"    ,&(localmat.m.dam_mp->beta)       ,&ierr);
   }
   frchk("MAT_DAMAGE_GE",&ierr);
   if (ierr==1)
   {
      localmat.mattyp      = m_damage_ge;
      localmat.m.damage_ge = new DAMAGE_GE();
      frint("EQU"        ,&(localmat.m.damage_ge->equival)    ,&ierr);
      frint("DAMT"       ,&(localmat.m.damage_ge->damtyp)     ,&ierr);
      frdouble("CRAD"    ,&(localmat.m.damage_ge->crad)       ,&ierr);
      frdouble("YOUNG"   ,&(localmat.m.damage_ge->youngs)     ,&ierr);
      frdouble("NUE"     ,&(localmat.m.damage_ge->nue)        ,&ierr);
      frdouble("KAPPA_0" ,&(localmat.m.damage_ge->kappa_0)    ,&ierr);
      frdouble("KAPPA_M" ,&(localmat.m.damage_ge->kappa_m)    ,&ierr);
      frdouble("ALPHA"   ,&(localmat.m.damage_ge->alpha)      ,&ierr);
      frdouble("BETA"    ,&(localmat.m.damage_ge->beta)       ,&ierr);
      frdouble("K_FAC"   ,&(localmat.m.damage_ge->k_fac)      ,&ierr);
   }
   frchk("MAT_HYPER_POLYCONVEX",&ierr);
   if (ierr==1)
   {
      localmat.mattyp      = m_hyper_polyconvex;
      localmat.m.hyper_polyconvex = new HYPER_POLYCONVEX();
      frdouble("C"       ,&(localmat.m.hyper_polyconvex->c)       ,&ierr);
      frdouble("K1"      ,&(localmat.m.hyper_polyconvex->k1)     ,&ierr);
      frdouble("K2"      ,&(localmat.m.hyper_polyconvex->k2)        ,&ierr);
      frdouble("EPSILON" ,&(localmat.m.hyper_polyconvex->epsilon)    ,&ierr);
      frdouble("GAMMA"   ,&(localmat.m.hyper_polyconvex->gamma)    ,&ierr);
      frdouble("DENS"    ,&(localmat.m.hyper_polyconvex->density)     ,&ierr);
   }
   // Anisotropic Polyconvex Material Law based on Balzani et. al.
   frchk("MAT_ANISOTROPIC_BALZANI",&ierr);
   if (ierr==1)
   {
      localmat.mattyp      = m_anisotropic_balzani;
      localmat.m.anisotropic_balzani = new ANISOTROPIC_BALZANI();
      frdouble("C1"     ,&(localmat.m.anisotropic_balzani->c1)       ,&ierr);
      frdouble("EPS1"   ,&(localmat.m.anisotropic_balzani->eps1)     ,&ierr);
      frdouble("EPS2"   ,&(localmat.m.anisotropic_balzani->eps2)     ,&ierr);
      frdouble("ALPHA1" ,&(localmat.m.anisotropic_balzani->alpha1)   ,&ierr);
      frdouble("ALPHA2" ,&(localmat.m.anisotropic_balzani->alpha2)   ,&ierr);
      frdouble("DENS"   ,&(localmat.m.anisotropic_balzani->density)  ,&ierr);
      frint   ("ALOC"   ,&(localmat.m.anisotropic_balzani->aloc)     ,&ierr);
      frdouble("A1X"    ,&(localmat.m.anisotropic_balzani->a1[0])    ,&ierr);
      frdouble("A1Y"    ,&(localmat.m.anisotropic_balzani->a1[1])    ,&ierr);
      frdouble("A1Z"    ,&(localmat.m.anisotropic_balzani->a1[2])    ,&ierr);
      frdouble("ALPHA1_2" ,&(localmat.m.anisotropic_balzani->alpha1_2)   ,&ierr);
      frdouble("ALPHA2_2" ,&(localmat.m.anisotropic_balzani->alpha2_2)   ,&ierr);
      frdouble("A2X"    ,&(localmat.m.anisotropic_balzani->a2[0])    ,&ierr);
      frdouble("A2Y"    ,&(localmat.m.anisotropic_balzani->a2[1])    ,&ierr);
      frdouble("A2Z"    ,&(localmat.m.anisotropic_balzani->a2[2])    ,&ierr);
   }
   /* Fourier's law of isotropic heat conduction --> heat cond. coeff. */
   frchk("MAT_Therm_Fourier_iso",&ierr);
   if (ierr==1)
   {
     localmat.mattyp = m_th_fourier_iso;
     localmat.m.th_fourier_iso = new TH_FOURIER_ISO();
     frdouble("CONDUCT"  ,&(localmat.m.th_fourier_iso->conduct)   ,&ierr);
     frdouble("CAPACITY" ,&(localmat.m.th_fourier_iso->capacity)  ,&ierr);
   }
   /* Fourier's law of general heat conduction --> heat cond. matrix */
   frchk("MAT_Therm_Fourier_gen",&ierr);
   if (ierr==1)
   {
     localmat.mattyp = m_th_fourier_gen;
     localmat.m.th_fourier_gen = new TH_FOURIER_GEN();
     frdouble_n("CONDUCT"  ,&(localmat.m.th_fourier_gen->conduct[0]), 9   ,&ierr);
     if (ierr == 1)
     {
       dserror("Insufficient conduction parameters!");
     }
   }
   /* Robinson's visco-plastic material for structures */
   frchk("MAT_Struct_Robinson",&ierr);
   if (ierr==1)
   {
     dserror("Robinson's visco-plastic material is not yet implemented for new discretization!\n");
   }
   /*multi layer material */
   frchk("MAT_Multilayer",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_multi_layer;
      localmat.m.multi_layer = (MULTI_LAYER*)CCACALLOC(1,sizeof(MULTI_LAYER));
      /*read number of kinamtic layer*/
      frint("NUM_KLAY"   ,&(localmat.m.multi_layer->num_klay)        ,&ierr);

      /*allocate memory for the different kinematic layers*/
      num_klay = localmat.m.multi_layer->num_klay;
      localmat.m.multi_layer->klayhgt=(DOUBLE*)CCACALLOC(num_klay,sizeof(DOUBLE));
      localmat.m.multi_layer->kinlay=(KINLAY*)CCACALLOC(num_klay,sizeof(KINLAY));

      /*----- read section data  -> hgt of different kinematic layers ---*/
      frdouble_n("SEC_KLAY",localmat.m.multi_layer->klayhgt,num_klay,&ierr);
      if (ierr!=1) dserror("Reading of klayhgt in MULTI_LAYER material failed");
      /*------------------------ check if sectian data adds up to 100 ---*/
      klay_sum = 0.0;
      for(j=0; j<num_klay; j++) klay_sum += localmat.m.multi_layer->klayhgt[j];
      if (FABS(klay_sum-100.0) > EPS5) dserror("klay_sum != 100 -> change Data in SEC_KLAY in MULTI_LAYER material");

      /*read one kinematic layer per line*/
      for(j=0; j<num_klay; j++)
      {
        actlay = &(localmat.m.multi_layer->kinlay[j]);
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
   /* Structural micro-scale approach: material parameters are
    * calculated from microscale simulation */
   frchk("MAT_Struct_Multiscale",&ierr);
   if (ierr==1)
   {
      localmat.mattyp = m_struct_multiscale;
      localmat.m.struct_multiscale = new _STRUCT_MULTISCALE();
      localmat.m.struct_multiscale->microdis = 1; /* currently only one
                                                   * microscale discretization
                                                   * is used in all Gauss
                                                   * points */
      char buffer[500];
      frchar("MICROFILE", buffer, &ierr);
      if (ierr!=1) dserror("No inputfile for microstructure given!\n");
      int length = strlen(buffer);

      localmat.m.struct_multiscale->micro_inputfile_name =
        (char*)CCACALLOC(length+1, sizeof(char));
      strcpy(localmat.m.struct_multiscale->micro_inputfile_name, &buffer[0]);
   }
   i++;

   /*----------------------------------------------------------------------*/
   /* add local material vector to problem instance                        */
   /*----------------------------------------------------------------------*/
   AddMaterial(localmat);

   /*----------------------------------------------------------------------*/
   frread();
}

if (i==0)
  dserror("No material could be read from inputfile\n");
return;
} /* end of inp_material */

#endif  // #ifdef CCADISCRET
