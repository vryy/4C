/*-----------------------------------------------------------------------*/
/*!
\file
\brief Brief description.

  Very detailed description.

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

 */
/*-----------------------------------------------------------------------*/

/*!
\addtogroup Fluid3_fast
*//*! @{ (documentation module open)*/

#ifdef D_FLUID3_F

#include "../headers/standardtypes.h"
#include "../fluid_full/fluid_prototypes.h"
#include "../fluid3_fast/f3f_prototypes.h"
#include "../fluid3/fluid3.h"

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | pointer to allocate dynamic variables if needed                      |
  | dedfined in global_control.c                                         |
  | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


static DOUBLE   *evhist;
static DOUBLE   *evelng;
static DOUBLE   *ealecovn;
static DOUBLE   *ealecovng;
static DOUBLE   *egridv;
static DOUBLE   *epren;
static DOUBLE   *edeadng;
static DOUBLE   *edeadn;

static DOUBLE   *elecord;
static DOUBLE   *tau;

static DOUBLE   *funct;
static DOUBLE   *deriv;
static DOUBLE   *deriv2;
static DOUBLE   *xjm;
static DOUBLE   *velint;
static DOUBLE   *vel2int;
static DOUBLE   *covint;
static DOUBLE   *alecovint;
static DOUBLE   *gridvint;
static DOUBLE   *vderxy;
static DOUBLE   *pderxy;
static DOUBLE   *vderxy2;
static DOUBLE   *derxy;
static DOUBLE   *derxy2;
static DOUBLE   *wa1;
static DOUBLE   *wa2;

static DOUBLE   *sigint;
static DOUBLE   *nostr;

static DOUBLE   *estif;
static DOUBLE   *emass;
static DOUBLE   *eforce;
static DOUBLE   *edforce;

static FLUID_DYNAMIC   *fdyn;



/*-----------------------------------------------------------------------*/
/*!
  \brief control routine for element integration of fast fluid3 elements

  This routine controls the element evaluation:
  - actual vel. and pres. variables are set
  - stabilisation parameters are calculated
  - element integration is performed
    --> element stiffness matrix and
    --> element load vectors
  - stiffness matrix and load vectors are permuted for assembling
  - element load vector due to dirichlet conditions is calculated

  \param ele[LOOPL]      ELEMENT  (i) the set of elements
  \param estif          *DOUBLE   (i) element stiffness matrix
  \param emass          *DOUBLE   (i) element mass matrix
  \param eforce         *DOUBLE   (i) element force
  \param edforce        *DOUBLE   (i) element dirichlet force
  \param ipos                     (i) node array positions
  \param hasdirich      *INT      (i) flag if s.th. was written to edforce
  \param hasext         *INT      (i) flag if there are ext forces
  \param init            INT      (i) init flag
  \param aloopl          INT      (i) num of elements in ele[]

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3fcalele(
    ELEMENT        *ele[LOOPL],
    ARRAY          *estif_fast,
    ARRAY          *emass_fast,
    ARRAY          *eforce_fast,
    ARRAY          *edforce_fast,
    ARRAY_POSITION *ipos,
    INT            *hasdirich,
    INT            *hasext,
    INT             init,
    INT             aloopl
    )
{

  INT l;
  INT sizevec[6];

  DOUBLE thsl;
  INT    nis;

#ifdef DEBUG
  dstrc_enter("f3fcalele");
#endif

  if (init==1) /* allocate working arrays and set pointers */
  {
    estif   = estif_fast->a.dv;
    emass   = emass_fast->a.dv;
    eforce  = eforce_fast->a.dv;
    edforce = edforce_fast->a.dv;

    evhist    = (DOUBLE *)CCAMALLOC(LOOPL*NUM_F3_VELDOF*MAXNOD_F3*sizeof(DOUBLE));
    evelng    = (DOUBLE *)CCAMALLOC(LOOPL*NUM_F3_VELDOF*MAXNOD_F3*sizeof(DOUBLE));
    ealecovn  = (DOUBLE *)CCAMALLOC(LOOPL*NUM_F3_VELDOF*MAXNOD_F3*sizeof(DOUBLE));
    ealecovng = (DOUBLE *)CCAMALLOC(LOOPL*NUM_F3_VELDOF*MAXNOD_F3*sizeof(DOUBLE));
    egridv    = (DOUBLE *)CCAMALLOC(LOOPL*NUM_F3_VELDOF*MAXNOD_F3*sizeof(DOUBLE));
    epren     = (DOUBLE *)CCAMALLOC(LOOPL*MAXNOD_F3*sizeof(DOUBLE));
    edeadn    = (DOUBLE *)CCAMALLOC(LOOPL*3*sizeof(DOUBLE));
    edeadng   = (DOUBLE *)CCAMALLOC(LOOPL*3*sizeof(DOUBLE));

    elecord   = (DOUBLE *)CCAMALLOC(LOOPL*MAXNOD_F3*3*sizeof(DOUBLE));
    tau       = (DOUBLE *)CCAMALLOC(LOOPL*3*sizeof(DOUBLE));

    funct     = (DOUBLE *)CCAMALLOC(MAXNOD_F3*sizeof(DOUBLE));
    deriv     = (DOUBLE *)CCAMALLOC(3*MAXNOD_F3*sizeof(DOUBLE));
    deriv2    = (DOUBLE *)CCAMALLOC(6*MAXNOD_F3*sizeof(DOUBLE));

    xjm       = (DOUBLE *)CCAMALLOC(LOOPL*3*3*sizeof(DOUBLE));
    velint    = (DOUBLE *)CCAMALLOC(LOOPL*NUM_F3_VELDOF*sizeof(DOUBLE));
    vel2int   = (DOUBLE *)CCAMALLOC(LOOPL*NUM_F3_VELDOF*sizeof(DOUBLE));
    covint    = (DOUBLE *)CCAMALLOC(LOOPL*NUM_F3_VELDOF*sizeof(DOUBLE));
    alecovint = (DOUBLE *)CCAMALLOC(LOOPL*NUM_F3_VELDOF*sizeof(DOUBLE));
    gridvint  = (DOUBLE *)CCAMALLOC(LOOPL*NUM_F3_VELDOF*sizeof(DOUBLE));
    vderxy    = (DOUBLE *)CCAMALLOC(LOOPL*3*3*sizeof(DOUBLE));
    pderxy    = (DOUBLE *)CCAMALLOC(LOOPL*3*sizeof(DOUBLE));
    vderxy2   = (DOUBLE *)CCAMALLOC(LOOPL*3*6*sizeof(DOUBLE));
    derxy     = (DOUBLE *)CCAMALLOC(LOOPL*3*MAXNOD_F3*sizeof(DOUBLE));
    derxy2    = (DOUBLE *)CCAMALLOC(LOOPL*6*MAXNOD_F3*sizeof(DOUBLE));
    wa1       = (DOUBLE *)CCAMALLOC(LOOPL*32*32*sizeof(DOUBLE));
    wa2       = (DOUBLE *)CCAMALLOC(LOOPL*32*32*sizeof(DOUBLE));

    sigint    = (DOUBLE *)CCAMALLOC(LOOPL*MAXGAUSS*6*sizeof(DOUBLE));
    nostr     = (DOUBLE *)CCAMALLOC(LOOPL*MAXNOD_F3*6*sizeof(DOUBLE));

    fdyn      = alldyn[genprob.numff].fdyn;
    goto end;
  } /* endif (init==1) */


  /* initialise with ZERO */
  amzero(estif_fast);
  amzero(emass_fast);
  amzero(eforce_fast);
  amzero(edforce_fast);

  for(l=0;l<LOOPL;l++)
  {
    hasext[l]=0;
    hasdirich[l]=0;
  }


  /*Variables for fortran functions.*/
  sizevec[0] = MAXNOD_F3;
  sizevec[1] = ele[0]->numnp;
  sizevec[2] = MAXNOD*MAXDOFPERNODE;
  sizevec[3] = LOOPL;
  sizevec[4] = aloopl;
  sizevec[5] = MAXGAUSS;



  switch(ele[0]->e.f3->is_ale)
  {
    /* Euler */
    case 0:

      /* set element data */
      f3fcalset(ele,evhist,evelng,epren,edeadn,edeadng,ipos,hasext,sizevec);

      /*A "C-function", but elecord[3,8] used in fortran this can be out of
        the gauss point loop*/
      f3fcalelecord(ele,elecord,sizevec);


    switch (ele[0]->e.f3->stab_type)
    {
      case stab_gls:

      /* calculate element size and stab-parameter: */
      f3fcalelesize(ele,funct,deriv,deriv2,
          derxy,xjm,evelng,velint,wa1,elecord,tau,sizevec);


      /* calculate element stiffness matrices
         and element force vectors */
      f3fcalint(ele,elecord,tau,hasext,estif,
          emass,eforce,
          funct,deriv,deriv2,
          xjm,derxy,derxy2,evelng,
          epren,edeadn,edeadng,velint,vel2int,
          covint,vderxy,pderxy,vderxy2,wa1,
          wa2,sizevec);

      break; /* end WAW GLS */

      case stab_usfem: /* completely linearised version */
      /* stab-parameter */
      f3fcaltau(ele,elecord,funct,deriv,derxy,velint,xjm,evelng,
                tau,wa1,sizevec);

      /* perform element integration */
      f3fint_usfem(ele,elecord,tau,hasext,estif,eforce,
                   funct,deriv,deriv2,xjm,derxy,derxy2,evelng,
                   evhist,NULL,epren,edeadng,velint,vel2int,NULL,
                   vderxy,vderxy2,pderxy,wa1,wa2,sizevec);
      goto dirich;

      default: dserror("unknown stabilisation type");

    } /* end switch over stab_type */

    break;

#ifdef D_FSI

    /* ALE */
    case 1:

      /* set element data */
      f3fcalseta(ele,evhist,evelng,ealecovng,egridv,
          epren,edeadn,edeadng,ipos,hasext,sizevec);


      /*A "C-function", but elecord[3,8] used in fortran this can be out of
        the gauss point loop*/
      f3falecord(ele,elecord,sizevec);

    switch (ele[0]->e.f3->stab_type)
    {
      case stab_gls:
      /* calculate element size and stab-parameter: */
      f3fcalelesize(ele,funct,deriv,deriv2,
          derxy,xjm,evelng,velint,wa1,elecord,tau,sizevec);


      /* calculate element stiffness matrices
         and element force vectors */
      f3fcalinta(ele,elecord,tau,hasext,estif,
          emass,eforce,
          funct,deriv,deriv2,
          xjm,derxy,derxy2,evelng,
          ealecovn,ealecovng,egridv,
          epren,edeadn,edeadng,velint,vel2int,
          covint,alecovint,gridvint,vderxy,pderxy,vderxy2,wa1,
          wa2,sizevec);
        break;/* end WAW GLS */

        case stab_usfem: /* completely linearised version */
        /* stab-parameter */
        f3fcaltau(ele,elecord,funct,deriv,derxy,alecovint,xjm,ealecovng,
                  tau,wa1,sizevec);

        /* perform element integration */
        f3fint_usfem(ele,elecord,tau,hasext,estif,eforce,
                     funct,deriv,deriv2,xjm,derxy,derxy2,evelng,
                     evhist,egridv,epren,edeadng,velint,vel2int,gridvint,
                     vderxy,vderxy2,pderxy,wa1,wa2,sizevec);
        goto dirich;

        default: dserror("unknown stabilisation type");
      } /* end switch over stab_type */

      break;

#endif /* ifdef D_FSI */


    default:
      dserror("parameter is_ale not 0 or 1!\n");
  }



#if 0
  /* add emass and estif to estif */
  f3fmake_estif(estif,emass,sizevec);
#else
  nis  = fdyn->nis;
  thsl = fdyn->thsl;
  f3fmast(estif,emass,&thsl,&nis,sizevec);
#endif



#if 0
  /* local co-ordinate system */
  if(ele->locsys==locsys_yes)
    locsys_trans(ele,estif,NULL,NULL,eforce);
#endif


    /* calculate emass * hist */
    f3fmassrhs(emass,evhist,edeadng,eforce,hasext,&(fdyn->thsl),sizevec);


dirich:


  /* calculate element load vector edforce */
  f3fcaldirich(ele,edforce,estif,ipos,hasdirich,0,sizevec);



end:

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3_calele */




/*-----------------------------------------------------------------------*/
/*!
  \brief control routine for stabilisation calculation

  evaluation of stabilisation parameter at the end of a time step

  \param ele[LOOPL]      ELEMENT  (i) the set of elements
  \param ipos                     (i) node array positions
  \param aloopl          INT      (i) num of elements in ele[]

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3fcalstab(
    ELEMENT      *ele[LOOPL],
    ARRAY_POSITION* ipos,
    INT           aloopl
    )
{

  INT      i,j;
  NODE    *actnode;
  INT      sizevec[6];
  INT      velnp;


#ifdef DEBUG
  dstrc_enter("f3fcalstab");
#endif

  fdyn   = alldyn[genprob.numff].fdyn;
  velnp  = ipos->velnp;


  sizevec[0] = MAXNOD_F3;
  sizevec[1] = ele[0]->numnp;
  sizevec[2] = MAXNOD*MAXDOFPERNODE;
  sizevec[3] = LOOPL;
  sizevec[4] = aloopl;
  sizevec[5] = MAXGAUSS;


  /* get actual co-ordinates */
  if (ele[0]->e.f3->is_ale==0)
    f3fcalelecord(ele,elecord,sizevec);
#ifdef D_FSI
  else
    f3falecord(ele,elecord,sizevec);
#endif


  for(i=0; i<sizevec[1]; i++) /* loop nodes */
  {
    for(j=0;j<sizevec[4];j++)
    {
      actnode=ele[j]->node[i];
      /* get actual velocity */
      evelng[                  LOOPL*i+j] = actnode->sol_increment.a.da[velnp][0];
      evelng[  MAXNOD_F3*LOOPL+LOOPL*i+j] = actnode->sol_increment.a.da[velnp][1];
      evelng[2*MAXNOD_F3*LOOPL+LOOPL*i+j] = actnode->sol_increment.a.da[velnp][2];
    }/*loop*/
  } /* end of loop over nodes */


  /* calculate element size and stab-parameter: */
  f3fcalelesize(ele,funct,deriv,deriv2,derxy,xjm,evelng,
      velint,wa1,elecord,tau,sizevec);


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3fcalstab */



/*-----------------------------------------------------------------------*/
/*!
  \brief control routine for stress calculation


  \param str             FLUID_STRESS (i) flag for stress calculation
  \param viscstr         INT      (i) viscose stresses yes/no?
  \param ele[LOOPL]      ELEMENT  (i) the set of elements
  \param ipos                     (i) node array positions
  \param ele[LOOPL]      ELEMENT  (i) the set of elements

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3fstress(
    FLUID_STRESS  str,
    INT           viscstr,
    ELEMENT      *ele[LOOPL],
    ARRAY_POSITION *ipos,
    INT           aloopl
    )
{

#ifdef D_FSI
  INT       i,l;
  INT       coupled;      /* flag for fsi interface element */
  INT       iel;          /* number of nodes per element */
  GNODE    *actgnode;     /* actual gnode */
#endif


  INT sizevec[6];


#ifdef DEBUG
  dstrc_enter("f3fstress");
#endif

  /*Variables for fortran functions.*/
  sizevec[0] = MAXNOD_F3;
  sizevec[1] = ele[0]->numnp;
  sizevec[2] = MAXNOD*MAXDOFPERNODE;
  sizevec[3] = LOOPL;
  sizevec[4] = aloopl;
  sizevec[5] = MAXGAUSS;

  switch(str)
  {
    case str_none:
      break;

#ifdef D_FSI
    case str_fsicoupling:
      /* check if any fluid element is coupled to struct element */
      coupled=0;
      for (l=0;l<aloopl;l++)
      {
        iel=ele[l]->numnp;
        for (i=0;i<iel;i++)
        {
          actgnode = ele[l]->node[i]->gnode;
          /* check if there is a coupled struct node */
          if (actgnode->mfcpnode[genprob.numsf]==NULL)
            continue;

          coupled=1;
          break;
        }
      }

      if (coupled==1)
        f3fcalelestress(viscstr,ele,evelng,epren,funct,
            deriv,derxy,vderxy,xjm,wa1,elecord,sigint,nostr,ipos,sizevec);
      break;
#endif

    case str_liftdrag:
    case str_all:
      f3fcalelestress(viscstr,ele,evelng,epren,funct,
          deriv,derxy,vderxy,xjm,wa1,elecord,sigint,nostr,ipos,sizevec);
      break;

    default:
      dserror("stress calculation not possible!\n");
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3fstress */

/*!---------------------------------------------------------------------
\brief control routine for integration of element residual

<pre>                                                        chfoe 05/05

This routine controls the integration of the elemental residual for fast
elements which is required to compute consistent nodal forces. 
These are used to be FSI coupling forces

</pre>

\param  *ele	         ELEMENT	(i)   actual element
\param  *eforce_global   ARRAY	        (o)   ele iteration force
\param  *hasext          INT	        (o)   element flag
\param  *ipos            ARRAY_POSITION (i)   position flag for solution
\return void

------------------------------------------------------------------------*/
void f3fcaleleres(
    ELEMENT        *ele[LOOPL],
    ARRAY          *eforce_fast,
    INT            *hasext,
    ARRAY_POSITION *ipos,
    INT             aloopl
    )
{

  INT l;
  INT sizevec[6];

#ifdef DEBUG
  dstrc_enter("f3fcaleleres");
#endif

/* initialise with ZERO */
amzero(eforce_fast);

for(l=0;l<LOOPL;l++)
{
  hasext[l]=0;
}


/*Variables for fortran functions.*/
sizevec[0] = MAXNOD_F3;
sizevec[1] = ele[0]->numnp;
sizevec[2] = MAXNOD*MAXDOFPERNODE;
sizevec[3] = LOOPL;
sizevec[4] = aloopl;
sizevec[5] = MAXGAUSS;



switch(ele[0]->e.f3->is_ale)
{
  /* Euler */
  case 0:

  /* set element data */
  f3fcalset(ele,evhist,evelng,epren,edeadn,edeadng,ipos,hasext,sizevec);

  /*A "C-function", but elecord[3,8] used in fortran this can be out of
    the gauss point loop*/
  f3fcalelecord(ele,elecord,sizevec);

  /* stab-parameter */
  f3fcaltau(ele,elecord,funct,deriv,derxy,velint,xjm,evelng,
            tau,wa1,sizevec);

  /* perform element integration */
  f3fint_res(ele,elecord,tau,hasext,eforce,
             funct,deriv,deriv2,xjm,derxy,derxy2,evelng,
             evhist,NULL,epren,edeadng,velint,vel2int,NULL,
             vderxy,vderxy2,pderxy,wa1,wa2,sizevec);

  break;

#ifdef D_FSI

 /* ALE */
  case 1:

   /* set element data */
  f3fcalseta(ele,evhist,evelng,ealecovng,egridv,
      epren,edeadn,edeadng,ipos,hasext,sizevec);

  /*A "C-function", but elecord[3,8] used in fortran this can be out of
    the gauss point loop*/
  f3falecord(ele,elecord,sizevec);
  /* stab-parameter */
  f3fcaltau(ele,elecord,funct,deriv,derxy,alecovint,xjm,ealecovng,
            tau,wa1,sizevec);

    /* perform element integration */
  f3fint_res(ele,elecord,tau,hasext,eforce,
             funct,deriv,deriv2,xjm,derxy,derxy2,evelng,
             evhist,ealecovng,epren,edeadng,velint,vel2int,alecovint,
             vderxy,vderxy2,pderxy,wa1,wa2,sizevec);
  break;
#endif /* D_FSI */


  default:
    dserror("parameter is_ale not 0 or 1!\n");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3_caleleres */




/*-----------------------------------------------------------------------*/
/*!
  \brief control function for error calculation of f3f elements


  \param ele[LOOPL]  ELEMENT        (i) the set of elements
  \param hasext     *INT            (i) flag whether there are ext forces
  \param container  *CONTAINER      (i) contains variables defined in container.h
  \param ipos       *ARRAY_POSITION (i)
  \param aloopl      INT            (i) num of elements in ele[]

  \return void

  \author mn
  \date   08/05
 */
/*-----------------------------------------------------------------------*/
void f3f_calerror(
    ELEMENT          *ele[LOOPL],
    INT              *hasext,
    CONTAINER        *container,
    ARRAY_POSITION   *ipos,
    INT               aloopl
    )
{

  INT l;
  INT sizevec[6];

  DOUBLE thsl;
  INT    nis;

#ifdef DEBUG
  dstrc_enter("f3f_calerror");
#endif



  /*Variables for fortran functions.*/
  sizevec[0] = MAXNOD_F3;
  sizevec[1] = ele[0]->numnp;
  sizevec[2] = MAXNOD*MAXDOFPERNODE;
  sizevec[3] = LOOPL;
  sizevec[4] = aloopl;
  sizevec[5] = MAXGAUSS;



  switch(ele[0]->e.f3->is_ale)
  {
    /* Euler */
    case 0:

      /* set element data */
      f3fcalset(ele,evhist,evelng,epren,edeadn,edeadng,ipos,hasext,sizevec);


      /*A "C-function", but elecord[3,8] used in fortran this can be out of
        the gauss point loop*/
      f3fcalelecord(ele,elecord,sizevec);
      break;

#ifdef D_FSI

      /* ALE */
    case 1:

      /* set element data */
      f3fcalseta(ele,evhist,evelng,ealecovng,egridv,
          epren,edeadn,edeadng,ipos,hasext,sizevec);


      /*A "C-function", but elecord[3,8] used in fortran this can be out of
        the gauss point loop*/
      f3falecord(ele,elecord,sizevec);

      break;

#endif /* ifdef D_FSI */


    default:
      dserror("parameter is_ale not 0 or 1!\n");
  }


  /* perform element integration */
  f3f_int_error(ele, elecord, funct, deriv, xjm, evelng, epren, velint, vel2int,
      container, sizevec);

end:

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3f_calerror */





#endif


/*! @} (documentation module close)*/

