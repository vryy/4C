#ifdef D_MLSTRUCT
/*!----------------------------------------------------------------------
\file
\brief input of submesh data

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../wall1/wall1.h"
#include "../struct2_ml/s2ml.h"
#include "../interf/interf.h"
/* #include "../headers/solution_mlpcg.h" */
#include "../solver/solver.h"
/*----------------------------------------------------------------------*
 |                                                       ah 03/04       |
 | vector of numfld submesh-FIELDs, defined in global_control.c         |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *sm_field;
/*- ist eine globale Variable und dafuer global in global_control.c  
    definiert, deshalb hier extern */
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
 | globale Variable, die hier definiert wird -> deshalb hier kein "extern"
 *----------------------------------------------------------------------*/
struct _MATERIAL  *sm_mat;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 |              Structure for solution process                          |
 | globale Variable, die hier definiert wird -> deshalb hier kein "extern"
 *----------------------------------------------------------------------*/
struct _SOLVAR  *sm_solv;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 |            one proc's info about his submesh partition               |
 | globale Variable, die hier definiert wird -> deshalb hier kein "extern"
 *----------------------------------------------------------------------*/
struct _PARTITION  *sm_part;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 |            ranks and communicators for submesh                       |
 | globale Variable, die hier definiert wird -> deshalb hier kein "extern"
 *----------------------------------------------------------------------*/
struct _PAR  *sm_par;
/*----------------------------------------------------------------------*
 | input of submesh                                          ah 4/04    |
 *----------------------------------------------------------------------*/
void inp_submesh()
{
INT  ierr=0;
INT  counter=0;
INT  i,j,k;
INT  nodeId;
INT  elenumber;
char *colpointer;

struct _ELEMENT    *actsmele;
struct _NODE       *actsmnode;

#ifdef DEBUG 
dstrc_enter("inp_submesh");
#endif


/*---------------------------------------- allocate the SUBMESH-FIELD */      
if(genprob.numfld != 1) dserror("NUMFIELD has to be 1 for structural multiscale problem");
sm_field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));

sm_field[0].fieldtyp = structure;
sm_field[0].ndis     = 1;
/*------------------------------------ allocate the SUBMESH-FIELD->DIS */      
sm_field[0].dis = (DISCRET*)CCACALLOC(sm_field[0].ndis,sizeof(DISCRET));

/*----------------------------------------------------------------------*
 |                   Einlesen der Submeshinformationen                  |
 *----------------------------------------------------------------------*/

if (frfind("--SMSIZE")==1)
{
  frread();
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    frint("SMELEMENTS", &(sm_field[0].dis[0].numele),&ierr);
    frint("SMNODES",    &(sm_field[0].dis[0].numnp),&ierr);
    frint("SMMATERIALS",&(genprob.nmat_sm),&ierr);

    frread();
  }
}

/*----------------------------------------------------------------------*
 |                   Einlesen der SubmeshELEMENTE                       |
 *----------------------------------------------------------------------*/

/*--------------------------- allocate the SUBMESH-FIELD->DIS->ELEMENTS */      
sm_field[0].dis[0].element=(ELEMENT*)CCACALLOC(sm_field[0].dis[0].numele,sizeof(ELEMENT));

if (frfind("--SMELE")==1)
{
  frread();
  counter=0;  
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    colpointer = allfiles.actplace;
    elenumber  = strtol(colpointer,&colpointer,10);
    sm_field[0].dis[0].element[counter].Id = --elenumber;  
   /*--------------------------------------------------- discretisation */
   frchk("QUAD4",&ierr);
   if (ierr==1)
   {
     sm_field[0].dis[0].element[counter].distyp = quad4;
     sm_field[0].dis[0].element[counter].numnp  = 4;
     /*------------------------------------------------ allocate the LM */      
     sm_field[0].dis[0].element[counter].lm  = (INT*)CCACALLOC(sm_field[0].dis[0].element[counter].numnp,sizeof(INT));
     /*-- da ich mehrere elemente in submech allociert habe, und mit [couter] die einzelnen
       dereferenziere, kommt vor lm ein . und kein ->, obwohl mit *lm in sm_element deklariert */
     if (sm_field[0].dis[0].element[counter].lm==NULL) dserror("Allocation of lm for submesh in ELEMENT failed");
     frint_n("QUAD4",&(sm_field[0].dis[0].element[counter].lm[0]),sm_field[0].dis[0].element[counter].numnp,&ierr);
     if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
   } 

   frchk("QUAD8",&ierr);
   if (ierr==1)
   {
     sm_field[0].dis[0].element[counter].distyp = quad8;
     sm_field[0].dis[0].element[counter].numnp  = 8;
     /*------------------------------------------------ allocate the LM */      
     sm_field[0].dis[0].element[counter].lm  = (INT*)CCACALLOC(sm_field[0].dis[0].element[counter].numnp,sizeof(INT));
     if (sm_field[0].dis[0].element[counter].lm==NULL) dserror("Allocation of lm for submesh in ELEMENT failed");
     frint_n("QUAD8",&(sm_field[0].dis[0].element[counter].lm[0]),sm_field[0].dis[0].element[counter].numnp,&ierr);
     if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
   } 
   /*-------------------- reduce node numbers in location matrix by one */
   for (i=0; i<sm_field[0].dis[0].element[counter].numnp; i++) (sm_field[0].dis[0].element[counter].lm[i])--;
   /*----------------------------------------- read the material number */
   frint("MAT",&(sm_field[0].dis[0].element[counter].mat),&ierr);
   if (ierr!=1) dserror("Reading of Submesh-material element failed");

   /*----------------------------------------------- elementtyp is WALL */
   frchk("WALL",&ierr);
   if (ierr==1)
   {
     sm_field[0].dis[0].element[counter].eltyp = el_wall1;
     /*--------------------------------------- allocate the WALLelement */      
     sm_field[0].dis[0].element[counter].e.w1 = (WALL1*)CCACALLOC(1,sizeof(WALL1));
     if (sm_field[0].dis[0].element[counter].e.w1==NULL) dserror("Allocation of submesh-WALLelement failed");
     /*------------------------------------ read WALL element thickness */
     frdouble("THICK",&(sm_field[0].dis[0].element[counter].e.w1->thick),&ierr);
     if (ierr!=1) dserror("Reading of Submesh WALL1-thickness element failed");
    /*-------------------------- read gaussian points for wall elements */
     frint_n("GP",&(sm_field[0].dis[0].element[counter].e.w1->nGP[0]),2,&ierr);
     if (ierr!=1) dserror("Reading of WALL1 Number of GP failed");
    /*-------------------------------------------- read 2D problem type */
     sm_field[0].dis[0].element[counter].e.w1->wtype = plane_stress;/* default */
     frchk("PLANE_STRESS",&ierr);
     if (ierr==1) sm_field[0].dis[0].element[counter].e.w1->wtype = plane_stress;
     frchk("PLANE_STRAIN",&ierr);
     if (ierr==1) sm_field[0].dis[0].element[counter].e.w1->wtype = plane_strain;
   }

   /*-------------------------------------------- elementtyp is Interf  */
   frchk("IF",&ierr);
   if (ierr==1)
   {
      sm_field[0].dis[0].element[counter].eltyp = el_interf;
      /*---------------------------------------- allocate the IFelement */      
      sm_field[0].dis[0].element[counter].e.interf = (INTERF*)CCACALLOC(1,sizeof(INTERF));
      if (sm_field[0].dis[0].element[counter].e.interf==NULL) dserror("Allocation of IFelement failed");
     /*-------------------------------------- read IFelement thickness */
     frdouble("THICK",&(sm_field[0].dis[0].element[counter].e.interf->thick),&ierr);
     if (ierr!=1) dserror("Reading of Submesh IF-thickness element failed");
    /*---------------------------- read gaussian points for IF elements */
     frint_n("GP",&(sm_field[0].dis[0].element[counter].e.interf->nGP),1,&ierr);
     if (ierr!=1) dserror("Reading of IF Number of GP failed");
   }

   counter++;
   frread();
  }/*----- end while (strncmp(allfiles.actplace,"------",6)!=0) --------*/
}
frrewind();

/*----------------------------------------------------------------------*
 |                   Einlesen der SubmeshNODES                          |
 *----------------------------------------------------------------------*/
/*------------------------------ allocate the SUBMESH-FIELD->DIS->NODES */      
sm_field[0].dis[0].node = (NODE*)CCACALLOC(sm_field[0].dis[0].numnp,sizeof(NODE));
/*-------------------------------------------------------------------- */      
if (frfind("--SMNODECOORD")==1)
{
  frread();
  counter=0;  
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /*------------------------------------------------------- NODE-Id */      
    frint("NODE",&(sm_field[0].dis[0].node[counter].Id),&ierr);
    if (ierr!=1) dserror("reading of submeshnodes failed");
    (sm_field[0].dis[0].node[counter].Id)--;
    /*---------------------------------------------- NODE-Coordinates */      
    frdouble_n("COORD",&(sm_field[0].dis[0].node[counter].x[0]),3,&ierr);
    if (ierr!=1) dserror("reading of submeshnodes failed");
    
    counter++;
    frread();
  }/*----- end while (strncmp(allfiles.actplace,"------",6)!=0) -------*/
}
frrewind();

/*----------------------------------------------------------------------*
 |                  Submesh-Topologie NODE-ELEMENT                      |
 *----------------------------------------------------------------------*/

/*-------------------------------------- pointer from ELEMENTS to NODES */      
for (i=0; i<sm_field[0].dis[0].numele; i++)
{
  actsmele = &(sm_field[0].dis[0].element[i]);
  /*------------------------------------------------ localID = globalID */      
  actsmele->Id_loc = actsmele->Id;
  actsmele->node = (NODE**)CCACALLOC(actsmele->numnp,sizeof(NODE*));

  for (j=0; j<actsmele->numnp; j++)
  {
    nodeId = actsmele->lm[j];
    actsmele->node[j] = &(sm_field[0].dis[0].node[nodeId]);
  } 
}

/*-------------------------------------- pointer from NODES to ELEMENTS */      
for (i=0; i<sm_field[0].dis[0].numnp; i++) 
{
   sm_field[0].dis[0].node[i].numele = 0;
   /*-------------------------------------------------- NODE-number DOF */      
   sm_field[0].dis[0].node[i].numdf = 2;
}

/*---------------------------- count the number of elements to one node */
for (i=0; i<sm_field[0].dis[0].numele; i++)
{
   actsmele = &(sm_field[0].dis[0].element[i]);
   for (j=0; j<actsmele->numnp; j++)
   {
      (actsmele->node[j]->numele)++; 
   }
}

/*------------------------- allocate space for element pointers in NODE */
for (i=0; i<sm_field[0].dis[0].numnp; i++)
{
   actsmnode = &(sm_field[0].dis[0].node[i]);
  /*------------------------------------------------ localID = globalID */      
   actsmnode->Id_loc = actsmnode->Id;
   actsmnode->element = (ELEMENT**)CCACALLOC(actsmnode->numele,sizeof(ELEMENT*));
   for (j=0; j<actsmnode->numele; j++) actsmnode->element[j] = NULL;
}

/*---------------- loop elements and point from their nodes to themself */
for (i=0; i<sm_field[0].dis[0].numele; i++)
{
   actsmele = &(sm_field[0].dis[0].element[i]);
   for (j=0; j<actsmele->numnp; j++)
   {
      actsmnode = actsmele->node[j];
      for (k=0; k<actsmnode->numele; k++)
      {
         if (actsmnode->element[k] == NULL) break;
      }
      actsmnode->element[k] = actsmele;
   }
}

/*----------------------------------------------------------------------*
 |          Erstellen der GNODES und NODE<->GNODE Topologie             |
 *----------------------------------------------------------------------*/
sm_field[0].dis[0].ngnode = sm_field[0].dis[0].numnp;

/*----------------------------- allocate the SUBMESH-FIELD->DIS->GNODES */      
sm_field[0].dis[0].gnode = (GNODE*)CCACALLOC(sm_field[0].dis[0].ngnode,sizeof(GNODE));
/*----------------------------------------- set pointers gnode <-> node */
for (i=0; i<sm_field[0].dis[0].numnp; i++)
{
  sm_field[0].dis[0].node[i].gnode = &(sm_field[0].dis[0].gnode[i]);
  sm_field[0].dis[0].gnode[i].node = &(sm_field[0].dis[0].node[i]);
}

/*----------------------------------------------------------------------*
 |                 Einlesen der Submesh-MATERIALS                      |
 *---------------------------------------------------------------------*/
/*------------------------------------------ allocate Submesh-Material */      
sm_mat = (MATERIAL*)CCACALLOC(genprob.nmat_sm,sizeof(MATERIAL));
/*---------------------------------------------------------------------*/
if (frfind("--SMMAT")==0) dserror("frfind: submesh-MATERIAL is not in input file");
frread();
i=0;
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frint("MAT",&(sm_mat[i].Id),&ierr);
  /*-------------------------------------------------------------------*/
   frchk("MAT_Struct_StVenantKirchhoff",&ierr);
   if (ierr==1)
   {
      sm_mat[i].mattyp = m_stvenant;
      sm_mat[i].m.stvenant = (STVENANT*)CCACALLOC(1,sizeof(STVENANT));
      frdouble("YOUNG"  ,&(sm_mat[i].m.stvenant->youngs)      ,&ierr);
      frdouble("NUE"    ,&(sm_mat[i].m.stvenant->possionratio),&ierr);
      frdouble("DENS"   ,&(sm_mat[i].m.stvenant->density)     ,&ierr);
   }
   frchk("MAT_MisesPlastic",&ierr);
   if (ierr==1)
   {
      sm_mat[i].mattyp = m_pl_mises;
      sm_mat[i].m.pl_mises = (PL_MISES*)CCACALLOC(1,sizeof(PL_MISES));
      frdouble("YOUNG",&(sm_mat[i].m.pl_mises->youngs)        ,&ierr);
      frdouble("NUE"  ,&(sm_mat[i].m.pl_mises->possionratio)  ,&ierr);
      frdouble("ALFAT",&(sm_mat[i].m.pl_mises->ALFAT)         ,&ierr);
      frdouble("Sigy" ,&(sm_mat[i].m.pl_mises->Sigy)          ,&ierr);
      sm_mat[i].m.pl_mises->Hard = 0.; 
      sm_mat[i].m.pl_mises->GF   = 0.; 
      sm_mat[i].m.pl_mises->betah= 1.; 
      frdouble("Hard" ,&(sm_mat[i].m.pl_mises->Hard)          ,&ierr);
      frdouble("GF"   ,&(sm_mat[i].m.pl_mises->GF)            ,&ierr);
      frdouble("BETAH",&(sm_mat[i].m.pl_mises->betah)         ,&ierr);
   }
   frchk("MAT_Damage",&ierr);
   if (ierr==1)
   {
      sm_mat[i].mattyp = m_damage;
      sm_mat[i].m.damage = (DAMAGE*)CCACALLOC(1,sizeof(DAMAGE));
      if (sm_mat[i].m.damage==NULL) dserror("Allocation of DAMAGE material failed");
      frdouble("YOUNG",&(sm_mat[i].m.damage->youngs)        ,&ierr);
      frdouble("NUE"  ,&(sm_mat[i].m.damage->possionratio)  ,&ierr);
      frint("Equival" ,&(sm_mat[i].m.damage->Equival)       ,&ierr);
      frint("Damtyp"  ,&(sm_mat[i].m.damage->Damtyp)        ,&ierr);
      frdouble("Kappa_0",&(sm_mat[i].m.damage->Kappa_0)     ,&ierr);
      frdouble("Kappa_m",&(sm_mat[i].m.damage->Kappa_m)     ,&ierr);
      frdouble("Alpha",  &(sm_mat[i].m.damage->Alpha)       ,&ierr);
      frdouble("Beta" ,  &(sm_mat[i].m.damage->Beta)        ,&ierr);
      frdouble("k_fac" ,  &(sm_mat[i].m.damage->k_fac)      ,&ierr);
   }
   frchk("MAT_IF",&ierr);
   if (ierr==1)
   {
      sm_mat[i].mattyp  = m_ifmat;
      sm_mat[i].m.ifmat = (IFMAT*)CCACALLOC(1,sizeof(IFMAT));
      frdouble("EMOD"   ,&(sm_mat[i].m.ifmat->emod)     ,&ierr);
      frdouble("KMOD"   ,&(sm_mat[i].m.ifmat->kmod )    ,&ierr);
      frdouble("GMOD"   ,&(sm_mat[i].m.ifmat->gmod )    ,&ierr);
      frdouble("DICK"   ,&(sm_mat[i].m.ifmat->dick)     ,&ierr);
      frdouble("QMOD"   ,&(sm_mat[i].m.ifmat->qmod)     ,&ierr);
      frdouble("DN"     ,&(sm_mat[i].m.ifmat->deltan)   ,&ierr);
      frdouble("DT"     ,&(sm_mat[i].m.ifmat->deltat)   ,&ierr);
      frdouble("MU"     ,&(sm_mat[i].m.ifmat->mu)       ,&ierr);
   }
   frchk("MAT_INTERF_THERM",&ierr);
   if (ierr==1)
   {
      sm_mat[i].mattyp  = m_interf_therm;
      sm_mat[i].m.interf_therm = (INTERF_THERM *)CCACALLOC(1,sizeof(INTERF_THERM ));
      frdouble("EMOD"     ,&(sm_mat[i].m.interf_therm->emod)     ,&ierr);
      frdouble("NU"       ,&(sm_mat[i].m.interf_therm->nu )      ,&ierr);
      frdouble("DICK"     ,&(sm_mat[i].m.interf_therm->dick )    ,&ierr);
      frint("EQUIVAL" ,&(sm_mat[i].m.interf_therm->equival)      ,&ierr);
      frint("DAMTYP"  ,&(sm_mat[i].m.interf_therm->damtyp)       ,&ierr);
      frdouble("KAPPA0_N" ,&(sm_mat[i].m.interf_therm->kappa0_n) ,&ierr);
      frdouble("ALPHA_N"  ,&(sm_mat[i].m.interf_therm->alpha_n)  ,&ierr);
      frdouble("BETA_N"   ,&(sm_mat[i].m.interf_therm->beta_n)   ,&ierr);
      frdouble("KAPPA0_T" ,&(sm_mat[i].m.interf_therm->kappa0_t) ,&ierr);
      frdouble("ALPHA_T"  ,&(sm_mat[i].m.interf_therm->alpha_t)  ,&ierr);
      frdouble("BETA_T"   ,&(sm_mat[i].m.interf_therm->beta_t)   ,&ierr);
   }
   frchk("MAT_DAM_MP",&ierr);
   if (ierr==1)
   {
      sm_mat[i].mattyp   = m_dam_mp;
      sm_mat[i].m.dam_mp = (DAM_MP*)CCACALLOC(1,sizeof(DAM_MP));
      frdouble("YOUNG"   ,&(sm_mat[i].m.dam_mp->youngs)     ,&ierr);
      frdouble("NUE"     ,&(sm_mat[i].m.dam_mp->nue)        ,&ierr);
      frdouble("KAPPA_0" ,&(sm_mat[i].m.dam_mp->kappa_0)    ,&ierr);
      frdouble("ALPHA"   ,&(sm_mat[i].m.dam_mp->alpha)      ,&ierr);
      frdouble("BETA"    ,&(sm_mat[i].m.dam_mp->beta)       ,&ierr);
   }
   frchk("MAT_DAMAGE_GE",&ierr);
   if (ierr==1)
   {
      sm_mat[i].mattyp      = m_damage_ge;
      sm_mat[i].m.damage_ge = (DAMAGE_GE*)CCACALLOC(1,sizeof(DAMAGE_GE));
      frint("EQU"        ,&(sm_mat[i].m.damage_ge->equival)    ,&ierr);
      frint("DAMT"       ,&(sm_mat[i].m.damage_ge->damtyp)     ,&ierr);
      frdouble("CRAD"    ,&(sm_mat[i].m.damage_ge->crad)       ,&ierr);
      frdouble("YOUNG"   ,&(sm_mat[i].m.damage_ge->youngs)     ,&ierr);
      frdouble("NUE"     ,&(sm_mat[i].m.damage_ge->nue)        ,&ierr);
      frdouble("KAPPA_0" ,&(sm_mat[i].m.damage_ge->kappa_0)    ,&ierr);
      frdouble("KAPPA_M" ,&(sm_mat[i].m.damage_ge->kappa_m)    ,&ierr);
      frdouble("ALPHA"   ,&(sm_mat[i].m.damage_ge->alpha)      ,&ierr);
      frdouble("BETA"    ,&(sm_mat[i].m.damage_ge->beta)       ,&ierr);
      frdouble("K_FAC"   ,&(sm_mat[i].m.damage_ge->k_fac)      ,&ierr);
   }
  /*-------------------------------------------------------------------*/
   i++;
  /*-------------------------------------------------------------------*/
   frread();
}/*------- end while (strncmp(allfiles.actplace,"------",6)!=0) -------*/
if (i != genprob.nmat_sm) dserror("number of submesh-materials incorrect");
frrewind();

/*----------------------------------------------------------------------*
 |  Dirichletrandbedingungen fuer G-Knoten am Rand des submeshes        |
 *----------------------------------------------------------------------*/

for (i=0; i<sm_field[0].dis[0].numnp; i++)
{
   actsmnode = &(sm_field[0].dis[0].node[i]);
    /*------------------------------------------------------ Randknoten */
   if(actsmnode->numele <= 2)
   {
    /*------------------------------------- Allokieren der Dirichlet-RB */
      actsmnode->gnode->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));

      amdef("onoff",&(actsmnode->gnode->dirich->dirich_onoff),6,1,"IV");
      amdef("val",  &(actsmnode->gnode->dirich->dirich_val),  6,1,"DV");
      amdef("curve",&(actsmnode->gnode->dirich->curve),       6,1,"IV");

    /*----------------------- Null-RB auf die beiden FHG'e vorschreiben */
      actsmnode->gnode->dirich->dirich_onoff.a.iv[0] = 1; 
      actsmnode->gnode->dirich->dirich_onoff.a.iv[1] = 1; 
      for (j=2; j<6; j++)
      actsmnode->gnode->dirich->dirich_onoff.a.iv[j] = 0;
   
      for (j=0; j<6; j++)
      {
        actsmnode->gnode->dirich->dirich_val.a.dv[j] = 0.0; 
        actsmnode->gnode->dirich->curve.a.iv[j] = 0; 
      }
   }
/*------------------------------- allocate macro info at submesh nodes */      
   actsmnode->sm_macroinfo = (SM_MACRO_INFO*)CCACALLOC(1,sizeof(SM_MACRO_INFO));
}

/*----------------------------------------------------------------------*
 |             Werte fuer SM_SOLVAR und SM_PAR,... festlegen            |
 |             (werden nicht im Eingeabefile eingelesen)                |
 *----------------------------------------------------------------------*/

/*--------------------------------------------- allocate Submesh-SOLVAR */      
sm_solv = (SOLVAR*)CCACALLOC(1,sizeof(SOLVAR));

sm_solv->matrixtyp = matrix_none;
sm_solv->solvertyp = umfpack;
sm_solv->parttyp   = cut_elements;
sm_solv->nsysarray = 1;

/*------------------------------------------ allocate Submesh-PARTITION */      
sm_part = (PARTITION*)CCACALLOC(1,sizeof(PARTITION));
sm_part->ndis = 1;
/*---------------------------- -> sm_part->pdis[0]->coupledofs.fdim==0; */
sm_part->pdis = (PARTDISCRET*)CCACALLOC(1,sizeof(PARTDISCRET));
/*- das was mit Makronetz in part_assignfield() getan wird fuer nprocs=1*/

for (i=0;i<sm_part->ndis;i++)
{
   sm_part->fieldtyp = sm_field[0].fieldtyp;
   sm_part->pdis[i].numnp       = sm_field[0].dis[i].numnp;
   sm_part->pdis[i].numele      = sm_field[0].dis[i].numele;
   sm_part->pdis[i].bou_numnp   = 0;
   sm_part->pdis[i].bou_numele  = 0;
   sm_part->pdis[i].bou_element = NULL;
   sm_part->pdis[i].bou_node    = NULL;
   sm_part->pdis[i].element = (ELEMENT**)CCACALLOC(sm_part->pdis[i].numele,sizeof(ELEMENT*));
   sm_part->pdis[i].node    = (NODE**)CCACALLOC(sm_part->pdis[i].numnp,sizeof(NODE*));
   if (sm_part->pdis[i].element==NULL) dserror("Allocation of element pointer in PARTITION failed");
   if (sm_part->pdis[i].node==NULL)    dserror("Allocation of node pointer in PARTITION failed");
   for (j=0; j<sm_field[0].dis[i].numele; j++) 
   sm_part->pdis[i].element[j] = &(sm_field[0].dis[i].element[j]);
   for (j=0; j<sm_field[0].dis[i].numnp; j++)  
   sm_part->pdis[i].node[j] = &(sm_field[0].dis[i].node[j]);      
}
/*------------------------------------------------ allocate Submesh-PAR */      
sm_par = (PAR*)CCACALLOC(1,sizeof(PAR));
sm_par->nprocs = 1;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_submesh */

#endif /* D_MLSTRUCT */



