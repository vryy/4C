#ifdef D_XFEM 
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid_full/fluid_prototypes.h"
#include "../ls/ls_prototypes.h"
#include "xfem_prototypes.h"



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA            *alldyn;   
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL   *mat;





static ARRAY     eveln_a;     /* element velocities at (n) */
static DOUBLE  **eveln;
static ARRAY     evelng_a;    /* element velocities at (n+gamma) */
static DOUBLE  **evelng;
static ARRAY     epren_a;     /* element pressures at (n) */
static DOUBLE   *epren;
static ARRAY     edeadn_a;    /* element dead load (selfweight) */
static DOUBLE   *edeadng;
static ARRAY     edeadng_a;   /* element dead load (selfweight) */
static DOUBLE   *edeadn;
static ARRAY     funct_a;     /* shape functions */
static DOUBLE   *funct;
static ARRAY     deriv_a;     /* first natural derivatives */
static DOUBLE  **deriv;
static ARRAY     deriv2_a;    /* second natural derivatives */
static DOUBLE  **deriv2;
static ARRAY     xyze_a;
static DOUBLE  **xyze;   
static ARRAY     xyzen_a;
static DOUBLE  **xyzen;   
static ARRAY     xjm_a;       /* Jocobian matrix */
static DOUBLE  **xjm;
static ARRAY     velint_a;    /* velocities at integration point */
static DOUBLE   *velint;
static ARRAY     vel2int_a;   /* velocities at integration point */
static DOUBLE   *vel2int;
static ARRAY	 covint_a;    /* convective velocities at integration point */
static DOUBLE   *covint;
static ARRAY     vderxy_a;    /* velocity derivatives */
static DOUBLE  **vderxy;
static ARRAY     pderxy_a;    /* pressure derivatives */
static DOUBLE   *pderxy;
static ARRAY     vderxy2_a;   /* velocity 2nd derivatives */
static DOUBLE  **vderxy2;
static ARRAY     derxy_a;     /* coordinate derivatives */
static DOUBLE  **derxy;
static ARRAY     derxy2_a;    /* 2nd coordinate derivatives */
static DOUBLE  **derxy2;
static ARRAY     ekappan_a;   /* surface curvature at (n) */
static DOUBLE   *ekappan;
static ARRAY     ekappang_a;  /* surface curvature at (n+1) */
static DOUBLE   *ekappang;
static ARRAY     w1_a;        /* working array of arbitrary chosen size */
static DOUBLE  **wa1;         /* used in different element routines */
static ARRAY     w2_a;        /* working array of arbitrary chosen size */
static DOUBLE  **wa2;         /* used in different element routines */
static DOUBLE  **estif;       /* pointer to global ele-stif */
static DOUBLE  **emass;       /* pointer to galerkin ele-stif */
static DOUBLE   *etforce;     /* pointer to Time RHS */
static DOUBLE   *eiforce;     /* pointer to Iteration RHS */
static DOUBLE   *edforce;     /* pointer to RHS due to dirichl. conditions */



static ARRAY     iarr_a;
static INT      *iarr;        /* local connectivity array */
static ARRAY     iand_a;
static INT      *iand;        /* index vector of active nodes */
static INT       nact;        /* nact => (iel     )  in case of standard formulation */
                              /*      => (iel+icnt)  in case of extended formulation */
static INT       ntotal;
static DOUBLE    thdt;        /* theta*dt  */
static INT       iel;         /* number of nodes per element  */
static ELEMENT  *myls2;       /* corresponding LS2 element */



static ARRAY     estif_temp_a;/* temporary arrays */
static DOUBLE  **estif_temp;
static ARRAY     emass_temp_a;
static DOUBLE  **emass_temp;
static ARRAY     eiforce_temp_a;
static DOUBLE   *eiforce_temp;
static ARRAY     etforce_temp_a;
static DOUBLE   *etforce_temp;

static FLUID_DYNAMIC   *fdyn;





/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_calele(
  FLUID_DATA     *data, 
  ELEMENT        *ele,             
  ARRAY          *estif_global,   
  ARRAY          *emass_global,   
  ARRAY          *etforce_global,       
  ARRAY          *eiforce_global, 
  ARRAY          *edforce_global,		
  INT            *hasdirich,      
  INT            *hasext,
  INT             imyrank,
  INT             is_relax,
  INT             init            
  )
{
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calele");
#endif
/*----------------------------------------------------------------------*/
  
  if (init==1) /* allocate working arrays and set pointers */
  {
    eveln     = amdef("eveln"    ,&eveln_a    ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
    evelng    = amdef("evelng"   ,&evelng_a   ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
    epren     = amdef("epren"    ,&epren_a    ,MAXNOD_F2,1,"DV");
    edeadn    = amdef("edeadn"   ,&edeadn_a   ,2,1,"DV");
    edeadng   = amdef("edeadng"  ,&edeadng_a  ,2,1,"DV");      
    funct     = amdef("funct"    ,&funct_a    ,MAXNOD_F2,1,"DV");
    deriv     = amdef("deriv"    ,&deriv_a    ,2,MAXNOD_F2,"DA");
    deriv2    = amdef("deriv2"   ,&deriv2_a   ,3,MAXNOD_F2,"DA");
    xjm       = amdef("xjm"      ,&xjm_a      ,2,2        ,"DA");
    xyze      = amdef("xyze"     ,&xyze_a     ,2,MAXNOD_F2,"DA");
    xyzen     = amdef("xyzen"    ,&xyzen_a    ,2,MAXNOD_F2,"DA");
    velint    = amdef("velint"   ,&velint_a   ,NUM_F2_VELDOF,1,"DV");
    vel2int   = amdef("vel2int"  ,&vel2int_a  ,NUM_F2_VELDOF,1,"DV");
    covint    = amdef("covint"   ,&covint_a   ,NUM_F2_VELDOF,1,"DV");
    vderxy    = amdef("vderxy"   ,&vderxy_a   ,2,2,"DA");
    pderxy    = amdef("pderxy"   ,&pderxy_a   ,2,1,"DV");
    vderxy2   = amdef("vderxy2"  ,&vderxy2_a  ,2,3,"DA");
    derxy     = amdef("derxy"    ,&derxy_a    ,2,MAXNOD_F2,"DA");
    derxy2    = amdef("derxy2"   ,&derxy2_a   ,3,MAXNOD_F2,"DA");
    ekappan   = amdef("ekappan"  ,&ekappan_a  ,MAXNOD_F2,1 ,"DV");
    ekappang  = amdef("ekappang" ,&ekappang_a ,MAXNOD_F2,1 ,"DV");
    wa1       = amdef("wa1"      ,&w1_a       ,50,50,"DA");
    wa2       = amdef("wa2"      ,&w2_a       ,50,50,"DA");

    /* allocate temp arrays and iarr and iand */
    iarr         = amdef("iarr"        , &iarr_a        , MAXNOD*MAXDOFPERNODE, ONE                 , "IV");
    iand         = amdef("iand"        , &iand_a        , MAXNOD              , ONE                 , "IV");
    estif_temp   = amdef("estif_temp"  , &estif_temp_a  , MAXNOD*MAXDOFPERNODE, MAXNOD*MAXDOFPERNODE, "DA");
    emass_temp   = amdef("emass_temp"  , &emass_temp_a  , MAXNOD*MAXDOFPERNODE, MAXNOD*MAXDOFPERNODE, "DA");
    eiforce_temp = amdef("eiforce_temp", &eiforce_temp_a, MAXNOD*MAXDOFPERNODE, ONE                 , "DV");
    etforce_temp = amdef("etforce_temp", &etforce_temp_a, MAXNOD*MAXDOFPERNODE, ONE                 , "DV");
    
    estif   = estif_global->a.da;
    emass   = emass_global->a.da;
    eiforce = eiforce_global->a.dv;
    etforce = etforce_global->a.dv;
    edforce = edforce_global->a.dv;

    fdyn = alldyn[genprob.numff].fdyn;

    goto end;
  }

  /* initialise with ZERO */
  amzero(estif_global);
  amzero(emass_global);
  amzero(eiforce_global);
  amzero(etforce_global);
  amzero(edforce_global);
  
  *hasdirich=0;
  *hasext=0;
  
  /* set some parameters */
  xfem_f2_init(ele);
  /* initialize iand and temp arrays */  
  xfem_f2_array_init();
  /* construct iand vector */
  xfem_f2_iand();
  /* construct local connectivity array */
  xfem_f2_loc_con();
  /* set element data */
  xfem_f2_calset(ele,xyze,eveln,evelng,epren,edeadn,edeadng,hasext);
  /* calculate element size and stability parameters */
  xfem_f2_calelesize(
    ele,data,xyze,funct,deriv,deriv2,xjm,derxy,
    vderxy,evelng,velint,wa1
    );
  /* calculate element tangent and mass matrices and "internal force" vectors */
  xfem_f2_calint(
    data,ele,hasext,estif_temp,emass_temp,etforce_temp,
    eiforce_temp,xyze,funct,deriv,deriv2,xjm,derxy,derxy2,
    eveln,evelng,epren,edeadn,edeadng,velint,vel2int,covint,
    vderxy,pderxy,vderxy2,wa1,wa2
    );
  /* perform local assembly for the tangent */
  xfem_f2_loc_ass_tangent();
  /* perform local assembly for the etforce (time RHS) */
  if (fdyn->nif!=0) xfem_f2_loc_ass_intforce(etforce,etforce_temp);
  /* perform local assembly for the eiforce (iteration RHS) */
  if (fdyn->nii+(*hasext)!=0) xfem_f2_loc_ass_intforce(eiforce,eiforce_temp);
  /* calculate element load vector edforce */
  fluid_caldirich(ele,edforce,estif,hasdirich,is_relax);
  
/*----------------------------------------------------------------------*/
 end:
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return; 
} /* end of xfem_f2_calele */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_loc_con()
{
  INT     i;
  INT     icnt;
  INT     piv,piv1;
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_loc_con");
#endif
/*----------------------------------------------------------------------*/

  /* set pivots */
  piv = THREE*iel;
  piv1 = TWO*iel;
  /*
   * => Part I "local connectivity for standard element formulation"
   */
  for (i=0; i<iel; i++)
  {
    iarr[2*(i+1)-2] = FIVE*(i+1)-4;
    iarr[2*(i+1)-1] = FIVE*(i+1)-3;
    iarr[piv1+(i+1)-1] = FIVE*(i+1)-2;          
  }
  /*
   * => Part II = Part I + "local connectivity for enriched formulation"
   */
  icnt = 0;
  for (i=0; i<iel; i++)
  {
    if (iand[i]==1)
    {
      iarr[piv+2*(i+1)-2] = FIVE*(i+1)-1;
      iarr[piv+2*(i+1)-1] = FIVE*(i+1);
      icnt++;
    }
  }
  nact = THREE*iel + TWO*icnt;
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_loc_con */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_loc_ass_tangent()
{
  INT        i,j;
  INT        i1,j1;
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_loc_ass_tangent");
#endif
/*----------------------------------------------------------------------*/

  /* perform local assembly */
  for (i=0; i<ntotal; i++)
  {
    i1 = iarr[i]-1;
    if (i1!=-1)
    {
      for (j=0; j<ntotal; j++)
      {
        j1 = iarr[j]-1;
        if (j1!=-1)
        {
          estif[i1][j1] = thdt*estif_temp[i][j] + emass_temp[i][j];
        }
      }
    }
  }
  
/*
 * => NOTE
 * to avoid singular system at least diagonal terms of the
 * tangent corresponding to dummy equations have to be set
 */
  /* finalize */
  for (i=0; i<ntotal; i++)
  {
    if (estif[i][i]==0.0) estif[i][i] = 1.0;
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_loc_ass_tangent */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_loc_ass_intforce(
  DOUBLE *intforce,
  DOUBLE *intforce_temp
  )
{
  INT        i;
  INT        i1;
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_loc_ass_intforce");
#endif
/*----------------------------------------------------------------------*/

  /* perform local assembly */
  for (i=0; i<ntotal; i++)
  {
    i1 = iarr[i]-1;
    if (i1!=-1)
    {
      intforce[i1] = intforce_temp[i];
    }
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_loc_ass_intforce */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_array_init()
{
  INT     i,j;
  
#ifdef DEBUG 
  dstrc_enter("xfem_array_init");
#endif
/*----------------------------------------------------------------------*/

  /* initialize iand */
  for (i=0; i<iel; i++) iand[i] = 0;
  /* initialize iarr */
  for (i=0; i<ntotal; i++) iarr[i] = 0;
  /* initialize temp arrays */
  for (i=0; i<ntotal; i++)
  {
    etforce_temp[i] = 0.0;
    eiforce_temp[i] = 0.0;    
  }
  
  for (i=0; i<ntotal; i++)
  {
    for (j=0; j<ntotal; j++)
    {
      estif_temp[i][j] = 0.0;
      emass_temp[i][j] = 0.0;    
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_array_init */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_init(
  ELEMENT *ele
  )
{
#ifdef DEBUG 
  dstrc_enter("xfem_f2_init");
#endif
/*----------------------------------------------------------------------*/

  /* set iel */
  iel = ele->numnp;
  /* set total number of equations */
  ntotal = FIVE*iel;
  /* compute theta*dt */
  thdt = fdyn->thsl;
  /* set ls2 element to me */
  myls2 = ele->e.f2->my_ls;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_init */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_iand()
{
  INT       i;
  INT       locid;
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_iand");
#endif
/*----------------------------------------------------------------------*/
  
/*
 * NOTE =>
 * be sure that the local node numbers are identical for both
 * fluid and levelset element
 */
  for (i=0; i<myls2->e.ls2->nenode; i++)
  {
    locid = myls2->e.ls2->enode[i]-1;
    iand[locid] = 1;
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_iand */
#endif
