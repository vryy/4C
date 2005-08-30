/*!
\file
\brief Postprocessing utility that shows the results using visual2.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Filters like this one are special inhabitants of the ccarat
world. They are always single processor applications yet they share
some code with ccarat and are closely linked to ccarat internals.

For visual we indeed do load all the results at once.

A huge part of this file was taken from the file visual_vis2.c.

\author u.kue
\date 10/04

*/

#include <ctype.h>
#include "post_visual2.h"


/*----------------------------------------------------------------------*/
/* the result values */
/*----------------------------------------------------------------------*/
static POST_DISCRETIZATION* discret;
static FIELD_DATA* struct_field = NULL;
static FIELD_DATA* fluid_field = NULL;
static FIELD_DATA* ale_field = NULL;
static INT struct_idx;
static INT fluid_idx;
static INT ale_idx;

#ifdef D_FSI
static INT* fluid_struct_connect;
static INT* fluid_ale_connect;
#endif

/*----------------------------------------------------------------------*
 |                                                        genk 07/02    |
 | all variables needed by VISUAL2 are defined extern                   |
 *----------------------------------------------------------------------*/
static INT      MPTRI,MPPTRI,MFACE;   /* see VISUAL2 manual             */
static INT      MPFACE,MEDGE,MPEDGE;  /* see VISUAL2 manual             */
static INT      KCELL;                /* number of visualised cells -
                                         different from number of elements! */
static INT      IOPT;                 /* program mode                   */
static INT      IVORT;                /* flag for vorticity calculation */
static INT      ITURBU;               /* flag for turbulence calc.      */
static INT      SCAL;
static INT      DSTEP;                /* increment of visualised steps  */
static INT      INCRE;                /* increment of visualised steps  */
static INT      NUMA;                 /* number of ALE-FIELD            */
static INT      CMNCOL=300;           /* see VISUAL2 manual             */
static INT      CMUNIT=37;            /* see VISUAL2 manual             */
static INT      MNODE;                /* maximum number of nodes        */
static INT      NKEYS=14;             /* see VISUAL2 manual             */
static INT      XYPIX[2];             /* size of VISUAL2 window         */
static INT      IKEYS[]={112,115,118,102,120,121,97,116,102,92,109,107,100,110};
static INT      FKEYS[]={1,1,1,3,1,1,1,1,3,1,1,1,1,1};
static float    FLIMS[14][2];         /* data limits                    */
static float    XYMIN[2], XYMAX[2];   /* min. and max. coordinates      */
static struct  _ARRAY PCELL_A;
static struct  _ARRAY WCELL_A;
static struct  _ARRAY WNODE_A;
static struct  _ARRAY WFACE_A;
static struct  _ARRAY CEDGE_A;  /* arrays needed in qat2v2        */
static INT    **PCELL;
static INT     *WCELL;
static INT     *WNODE;
static INT     *WFACE;
static INT     *CEDGE;          /* pointers to arrays             */


static DOUBLE  minxy,maxxy;
static DOUBLE  centerx,centery;
static DOUBLE  minpre,maxpre;
static DOUBLE  minvx,maxvx;
static DOUBLE  minvy,maxvy;
static DOUBLE  minabsv,maxabsv;
static DOUBLE  minvort,maxvort;
static DOUBLE  minkappa,maxkappa;
static DOUBLE  mindissi,maxdissi;
static DOUBLE  minvisco,maxvisco;
static DOUBLE  mingv,maxgv;
static DOUBLE  minvort,maxvort;

static INT     isstr=0;         /* flag for streamline calculation	*/
static INT     icol=-1;         /* act. num. of sol step to be visual.  */
static INT     FIRSTSTEP;
static INT     LASTSTEP;
static INT     ACTSTEP;
static INT     hsize, vsize;    /* size for VISUAL2 window		*/
static INT     for_true=-1;     /* flag for v2_cursor			*/
static INT     for_false=0;     /* flag for v2_cursor			*/
static INT     STSTEP=-2;       /* stopping step                        */
static INT     IMOVIE=0;        /* counter to slow down for movie creat.*/
static INT     bgcolour;        /* background colour                    */
static float   sstrval=0.0;     /* stationary streamline value		*/
static float   dx,dy,dxy;       /* for determination of coord. limits*/
static float   ratio;	        /*					*/
static DOUBLE  yscale;	        /* scaling factors for y-direction	*/
static DOUBLE  velx;
static DOUBLE  vely;
static DOUBLE  pres;
static DOUBLE  absv;
static DOUBLE  vort;
static DOUBLE  kappa;
static DOUBLE  dissi;
static DOUBLE  visco;

static FIELDTYP actfieldtyp;
static INT    nsteps;
static ARRAY  time_a;           /* time array				*/
static ARRAY  step_a;           /* time array				*/

static CHAR* yesno_options[] = { "no", "yes", NULL };


/*----------------------------------------------------------------------*/
/*!
  \brief set VISUAL2's constants and structures

  use QAT2V2 to make VISUAL2 structures
  this routine is called by VISUAL2 during the visualisation

  </pre>
  \param  *KNODE     INT   (o)    Number of nodes
  \param	*PTRI	   INT	 (o)	PolyTri pointers to nodes
  \param	*CTRI	   INT	 (o)	PolyTri pointers to neighbouring cells
  \param	*ITRI	   INT	 (o)	PolyTri pointers to original cells
  \param	*KPTRI	   INT	 (o)	Number of PolyTri strips
  \param	*PPTRI	   INT	 (o)	Pointers to ends of strips
  \param	*PFACE	   INT	 (o)	Face PolyLine pointers to nodes
  \param	*KPFACE    INT	 (o)	Number of face PolyLines
  \param	*PPFACE    INT	 (o)	Pointers to ends of face PolyLines
  \param	*PEDGE	   INT	 (o)	Edge PolyLine pointers to nodes
  \param	*KPEDGE    INT	 (o)	Number of Edge PolyLines
  \param	*PPEDGE    INT	 (o)	Pointers to ends of edge PolyLines

  \author genk
  \date 07/02
*/
/*----------------------------------------------------------------------*/
void v2data(INT *KNODE, INT *PTRI, INT *CTRI, INT *ITRI, INT *KPPTRI,
            INT *PPTRI, INT *PFACE, INT *KPFACE, INT *PPFACE,
	    INT *PEDGE, INT *KPEDGE, INT *PPEDGE)
{
  INT KPTRI, KEDGE, KFACE;

#ifdef DEBUG
  dstrc_enter("V2DATA");
#endif

  *KNODE = fluid_field->numnp;

/*------------------------------------------------ generate polystrips */
  qat2v2(&PCELL[0][0],WCELL,&KCELL,WNODE,KNODE,WFACE,
         &MPTRI,&MPPTRI,&MEDGE,&MPEDGE,&MFACE,&MPFACE,
         PTRI,CTRI,ITRI,&KPTRI,PPTRI,KPPTRI,
         PEDGE,CEDGE,&KEDGE,PPEDGE,KPEDGE,
         PFACE,&KFACE,PPFACE,KPFACE);

/*-------------------------------------------------------------- output */
  printf("\n");
  printf("   Number of Nodes           = %d \n",MNODE);
  printf("   Number of PolyTri Strips  = %d \n",*KPPTRI);
  printf("   Length of PolyTri Strips  = %d \n",KPTRI);
  if (*KPPTRI!=0)
    printf("   Ave PolyTri Strip Length  = %d \n",KPTRI/(*KPPTRI));
  printf("   Number of PolyLine Strips = %d \n",*KPFACE);
  printf("   Length of PolyLine Strips = %d \n",KFACE);
  if (*KPFACE!=0)
    printf("   Ave PolyLine Strip Length = %d \n",KFACE/(*KPFACE));
  printf("   Number of edge PLn Strips = %d \n",*KPEDGE);
  printf("   Length of edge PLn Strips = %d \n",KEDGE);
  if (*KPEDGE!=0)
    printf("   Ave Edge PLn Strip Length = %d \n",KEDGE/(*KPEDGE));
  printf("\n");

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of V2DATA */


/*----------------------------------------------------------------------*/
/*!
  \brief set VISUAL2's grid coordinates

  this routine is called by VISUAL2 during the visualisation

  \param  *XY        float   (o)  xy coordinates of grid nodes
  \param	*XYDELT    float   (o)  Maximum width (dx) and heigth dy of any cell
  \param	*XYMINP    float   (o)  Minimum values of x,y in each PolyTri strip
  \param	*XYMAXP    float   (o)  Maximum values of x,y in each PolyTri strip
  \param	*XYMINF    float   (o)  Minimum values of x,y in each face PolyTri
  \param	*XYMAXF    float   (o)  Maximum values of x,y in each face PolyTri
  \param	*XYMINE    float   (o)  Minimum values of x,y in each edge PolyTri
  \param	*XYMAXE    float   (o)  Maximum values of x,y in each edge PolyTri

  \author genk
  \date 07/02
*/
/*----------------------------------------------------------------------*/
void v2grid(float XY[][2], float *XYDELT, float *XYMINP, float *XYMAXP,
            float *XYMINF, float *XYMAXF, float *XYMINE, float *XYMAXE)
{

  INT i;

#ifdef DEBUG
  dstrc_enter("V2GRID");
#endif

  switch (actfieldtyp) {
  case fluid:
    if (NUMA>=0 && IOPT==0)
    {
      for (;;) {
        v2_cursor(&for_true);
        printf("Give the column number (min. 0; max. %d): ?\n",nsteps-1);
        scanf("%d",&icol);
        v2_cursor(&for_false);
        if (icol<0 || icol>nsteps-1)
        {
          printf("Column number out of range. New input!\n");
          printf("\n");
        }
        else {
          break;
        }
      }
    }
    if (IOPT>=2)
    {
#ifdef D_FSI
      INT i;
      for (i=0;i<fluid_field->numnp;i++)
      {
        NODE* actnode = &(discret[fluid_idx].node[i]);
        if (fluid_ale_connect[i] == -1) {
          XY[i][0] = actnode->x[0];
          XY[i][1] = actnode->x[1];
        }
        else
        {
          NODE* actanode = &(discret[ale_idx].node[fluid_ale_connect[i]]);
          XY[i][0] = actanode->x[0] + actanode->sol.a.da[icol][0];
          XY[i][1] = actanode->x[1] + actanode->sol.a.da[icol][1];
        }
      }
#else
      dserror("FSI-functions not compiled in!\n");
#endif
    }
    else
    {
      for (i=0;i<fluid_field->numnp;i++)
      {
        NODE* actnode = &(discret[fluid_idx].node[i]);
        XY[i][0] = actnode->x[0];
        XY[i][1] = actnode->x[1];
      }
    }
    break;
  case structure:
    dserror("fieldtyp not implemented yet!");
  case ale:
    dserror("fieldtyp not implemented yet!");
  default:
    dserror("fieldtyp not valid!\n");
  } /* end switch(actfieldtyp) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of V2GRID */


/*----------------------------------------------------------------------*/
/*!
  \brief set VISUAL2's scalar function values

  This routine is called by VISUAL2 during the visualisation.
  It filters the scalar values out of the nodal sol-field at the actual
  position icol.
  icol is set every timestep in v2update

  \param  *JKEY      INT    (i)  Index of key
  \param	*S         float  (o)  Scalar function at nodes or cells
  \sa v2update

  \author genk
  \date 07/02
*/
/*----------------------------------------------------------------------*/
void v2scal(INT *JKEY, float *S)
{
  INT i;
  float vx,vy;

#ifdef DEBUG
  dstrc_enter("V2SCAL");
#endif

  switch (actfieldtyp)
  {
  case fluid:
    if (IOPT==0) /* read in the column number */
    {
      for (;;) {
        v2_cursor(&for_true);
        printf("Give the column number (min. 0; max. %d): ?\n",nsteps-1);
        scanf("%d",&icol);
        v2_cursor(&for_false);
        if (icol<0 || icol>nsteps-1)
        {
          printf("Column number out of range. New input!\n");
          printf("\n");
        }
        else {
          break;
        }
      }
    }
    /*-------------------------------------------------- get scalar data */
    switch(*JKEY)
    {
    case 1: /* pressure */
      for (i=0;i<fluid_field->numnp;i++)
      {
        NODE* actnode = &(discret[fluid_idx].node[i]);
        S[i] = actnode->sol.a.da[icol][2];
      }
      break;
      /*-------------------------------------------------------------------*/
    case 2: /* ???? */
      printf("NOT IMPLEMTED SO FAR!!!!!!!!!!!!!!!!!!!!\n");
      break;
      /*-------------------------------------------------------------------*/
    case 3: /* vorticity */
      if (IVORT==1)
      {
        for (i=0;i<fluid_field->numnp;i++)
        {
          NODE* actnode = &(discret[fluid_idx].node[i]);
          S[i]=actnode->sol.a.da[icol][3];
        }
      }
      else
      {
        printf("vorticity not calculated!!!! --> ZERO-field\n");
        for (i=0;i<fluid_field->numnp;i++)
          S[i]=ZERO;
      }
      break;
    case 4:
      printf("VECTOR DATA IN SCALAR ROUTINE ???????\n");
      break;
      /*-------------------------------------------------------------------*/
    case 5: /* horizontal velocity */
      for (i=0;i<fluid_field->numnp;i++)
      {
        NODE* actnode = &(discret[fluid_idx].node[i]);
        S[i]=actnode->sol.a.da[icol][0];
      }
      break;
      /*-------------------------------------------------------------------*/
    case 6: /* vertical velocity */
      for (i=0;i<fluid_field->numnp;i++)
      {
        NODE* actnode = &(discret[fluid_idx].node[i]);
        S[i]=actnode->sol.a.da[icol][1];
      }
      break;
      /*-------------------------------------------------------------------*/
    case 7: /* absolute velocity */
      for (i=0;i<fluid_field->numnp;i++)
      {
        NODE* actnode = &(discret[fluid_idx].node[i]);
        vx=actnode->sol.a.da[icol][0];
        vy=actnode->sol.a.da[icol][1];
        S[i] = sqrt(vx*vx+vy*vy);
      }
      break;
      /*-------------------------------------------------------------------*/
    case 8: /* switch streamline - stationary streamline */
      if (isstr==0)
      {
        isstr=1;
        v2_cursor(&for_true);
        printf("Subtraction value for stationary streamlines: ?\n");
        scanf("%f",&sstrval);
        v2_cursor(&for_false);
      }
      else
        isstr=0;
      for (i=0;i<fluid_field->numnp;i++) {
        NODE* actnode = &(discret[fluid_idx].node[i]);
        S[i] = actnode->sol.a.da[icol][2];
      }
      break;
      /*-------------------------------------------------------------------*/
    case 9:
      printf("VECTOR DATA IN SCALAR ROUTINE ???????\n");
      break;
      /*-------------------------------------------------------------------*/
    case 10: /* stopping at step */
      printf("   Which step do you want to stop?\n");
      scanf("%d",&STSTEP);
      break;
      /*-------------------------------------------------------------------*/
    case 11:
      if (IMOVIE==0)
      {
        printf("\n");
        printf("   Starting movie creation at time 0.0\n");
        printf("\n");
        IMOVIE=1;
      }
      break;
      /*-------------------------------------------------------------------*/
    case 12: /* kappa */
      if (ITURBU==1)
      {
#if 0
        for (i=0;i<fluid_field->numnp;i++)
        {
          NODE* actnode = &(discret[fluid_idx].node[i]);
          /*actnode=&(actfield->dis[1].node[i]);*/
	  S[i]=actnode->sol.a.da[icol][0];
        }
#endif
      }
      else
      {
        printf("kappa not calculated!!!! --> ZERO-field\n");
        for (i=0;i<fluid_field->numnp;i++)
          S[i]=ZERO;
      }
      break;
      /*-------------------------------------------------------------------*/
    case 13: /* dissi */
      if (ITURBU==1)
      {
#if 0
        for (i=0;i<fluid_field->numnp;i++)
        {
          actnode=&(actfield->dis[1].node[i]);
	  S[i] = actnode->sol.a.da[icol][2];
        }
#endif
      }
      else
      {
        printf("dissipation not calculated!!!! --> ZERO-field\n");
        for (i=0;i<fluid_field->numnp;i++)
          S[i]=ZERO;
      }
      break;
      /*-------------------------------------------------------------------*/
    case 14: /* visco */
      if (ITURBU==1)
      {
#if 0
        for (i=0;i<fluid_field->numnp;i++)
        {
          actnode=&(actfield->dis[1].node[i]);
	  S[i] = actnode->sol.a.da[icol][1];
        }
#endif
      }
      else
      {
        printf("eddy-viscosity with ke- or kw-model not calculated!!! --> ZERO-field\n");
        for (i=0;i<fluid_field->numnp;i++)
          S[i]=ZERO;
      }
      break;
      /*-------------------------------------------------------------------*/
    } /* end switch(*JKEY) */
    break;
  case structure:
    dserror("fieldtyp not implemented yet!\n");
  case ale:
    dserror("fieldtyp not implemented yet!\n");
  default:
    dserror("fieldtyp not valid!\n");
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of V2SCAL */


/*----------------------------------------------------------------------*/
/*!
  \brief set VISUAL2's vector function values

  This routine is called by VISUAL2 during the visualisation.
  It filters the vector values out of the nodal sol-field at the actual
  position icol.
  icol is set every timestep in v2update

  \param  *JKEY      INT    (i)  Index of key
  \param	*V         float  (o)  Vector function at nodes or cells
  \return void
  \sa v2update

  \author genk
  \date 07/02
*/
/*----------------------------------------------------------------------*/
void v2vect(INT *JKEY, float V[][2])
{
  INT i;

#ifdef DEBUG
  dstrc_enter("V2VECT");
#endif

  switch (actfieldtyp)
  {
  case fluid:
    /*-------------------------------------------------- get vector data */
    if (*JKEY==4)
    {
      switch(isstr)
      {
      case 0: /* streamlines */
        for (i=0;i<fluid_field->numnp;i++)
        {
          NODE* actnode = &(discret[fluid_idx].node[i]);
          V[i][0] = actnode->sol.a.da[icol][0];
          V[i][1] = actnode->sol.a.da[icol][1];
        }
        break;
        /*-------------------------------------------------------------------*/
      default: /* stationary streamlines */
        printf("stationary streamlines not checked yet!!!!");
        for (i=0;i<fluid_field->numnp;i++)
        {
          NODE* actnode = &(discret[fluid_idx].node[i]);
#if 0
          actgnode=actnode->gnode;
          if (actgnode->dirich!=NULL)
          {
            if (actgnode->dirich->dirich_val.a.dv[0]==0.0)
            {
              V[i][0] = actnode->sol.a.da[icol][0];
            }
            else
            {
              V[i][0] = actnode->sol.a.da[icol][0]-sstrval;
            }
            V[i][1] = actnode->sol.a.da[icol][1];
          }
          else
#endif
          {
            V[i][0] = actnode->sol.a.da[icol][0]-sstrval;
            V[i][1] = actnode->sol.a.da[icol][1];
          }
        }
        break;
      }
    }
    break;
  case structure:
    dserror("fieldtyp not implemented yet!\n");
  case ale:
    dserror("fieldtyp not implemented yet!\n");
  default:
    dserror("fieldtyp not valid!\n");
  } /* end switch(actfieldtyp) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of V2VECT */


/*----------------------------------------------------------------------*/
/*!
  \brief update time and icol

  This routine is called by VISUAL2 during the visualisation.
  Beside the actual time the actual position in the nodal sol-history
  icol is also updated

  \param  *TIME      float    (o)  actual time
  \return void

  \author genk
  \date 07/02
*/
/*----------------------------------------------------------------------*/
void v2update(float *TIME)
{

#ifdef DEBUG
  dstrc_enter("V2UPDATE");
#endif


/*-------------------------------------------------- check for stopping */
  if (icol==STSTEP)
  {
    printf("   Next Step to Stop (-1 to terminate)\n");
    scanf("%d",&STSTEP);
  }

/*------------------------------------------- check for movie creation */
  if (icol==0 && IMOVIE==2)
  {
    IMOVIE=0;
    printf("\n");
    printf("Movie creation finished\n");
  }
  if (icol==0 && IMOVIE==1) IMOVIE=2;
  if (IMOVIE==2) v2movie();


  switch (actfieldtyp)
  {
  case fluid:
    if (icol==-1 || icol+INCRE > nsteps-1)
    {
      icol=0;
      ACTSTEP=step_a.a.iv[icol];
    }
    else
    {
      icol+=INCRE;
      ACTSTEP=step_a.a.iv[icol];
    }
    *TIME = time_a.a.dv[icol];

    printf("Time: %8.4f    Step: %5d/%-5d\r",time_a.a.dv[icol],ACTSTEP,LASTSTEP);
    fflush(stdout);
    break;
  case structure:
    dserror("fieldtyp not implemented yet!\n");
  case ale:
    dserror("fieldtyp not implemented yet!\n");
  default:
    dserror("fieldtyp not valid!\n");
  } /* end switch(actfieldtyp) */


/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of V2UPDATE */


/*----------------------------------------------------------------------*/
/*!
  \brief  additional label for  VISUAL2 plots

  This routine is called by VISUAL2 during the visualisation.

  \param  STRING     char    (o)  plot on VISUAL2 window
  \return void

  \author genk
  \date 07/02
*/
/*----------------------------------------------------------------------*/
void v2string(char STRING[81])
{
  float t;
  INT i;
  INT len;

#ifdef DEBUG
  dstrc_enter("V2STRING");
#endif

  switch (actfieldtyp)
  {
  case fluid:
    t=time_a.a.dv[icol];
    sprintf(STRING,"Time: %8.4f    Step: %5d/%-5d",t,ACTSTEP,LASTSTEP);

    /* Fill the remaining chars with spaces and append a zero. */
    len = strlen(STRING);
    for (i=len; i<80; ++i) {
      STRING[i] = ' ';
    }
    /*STRING[i] = '\0';*/
    break;
  case structure:
    dserror("fieldtyp not implemented yet!\n");
  case ale:
    dserror("fieldtyp not implemented yet!\n");
  default:
    dserror("fieldtyp not valid!\n");
  } /* end switch(actfieldtyp) */


/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of V2STRING */


/*----------------------------------------------------------------------*/
/*!
  \brief compute the array WCELL for use in QAT2V2

  \param  *field     (i)  field data
  \param  *discret   (i)  element and node data
  \param  distyp     (i)  the distyp of all elements

  \warning QAT2V2 requires FORTRAN numbering, so increase node-numbers
  by one!

  \author genk
  \date 07/02
*/
/*----------------------------------------------------------------------*/
void post_v2cell(FIELD_DATA *field, POST_DISCRETIZATION* discret, DIS_TYP distyp)
{
  INT i,j;
  INT inel;

#ifdef DEBUG
  dstrc_enter("v2cell");
#endif

  inel=-1;

  switch (distyp)
  {
  case tri3: /* 3 node triangle */
    for (i=0;i<field->numele;i++) {
      ELEMENT* actele = &(discret->element[i]);
      inel++;
      for(j=0;j<3;j++) {
        PCELL[inel][j] = actele->node[j]->Id_loc+1;
      }
      PCELL[inel][3] = PCELL[inel][0];
    }
    break;
  case tri6:
    for (i=0;i<field->numele;i++) {
      ELEMENT* actele = &(discret->element[i]);
      inel++;
      /*----------------------------------------- sub-element 1 */
      PCELL[inel][0] = actele->node[0]->Id_loc+1;
      PCELL[inel][1] = actele->node[3]->Id_loc+1;
      PCELL[inel][2] = actele->node[5]->Id_loc+1;
      PCELL[inel][3] = actele->node[0]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 2 */
      PCELL[inel][0] = actele->node[3]->Id_loc+1;
      PCELL[inel][1] = actele->node[1]->Id_loc+1;
      PCELL[inel][2] = actele->node[4]->Id_loc+1;
      PCELL[inel][3] = actele->node[3]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 3 */
      PCELL[inel][0] = actele->node[3]->Id_loc+1;
      PCELL[inel][1] = actele->node[4]->Id_loc+1;
      PCELL[inel][2] = actele->node[5]->Id_loc+1;
      PCELL[inel][3] = actele->node[3]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 4 */
      PCELL[inel][0] = actele->node[4]->Id_loc+1;
      PCELL[inel][1] = actele->node[3]->Id_loc+1;
      PCELL[inel][2] = actele->node[5]->Id_loc+1;
      PCELL[inel][3] = actele->node[4]->Id_loc+1;
    }
    break;
  case quad4: /* 4 node rectangle */
    for (i=0;i<field->numele;i++) {
      ELEMENT* actele = &(discret->element[i]);
      inel++;
      for(j=0;j<4;j++) {
        PCELL[inel][j] = actele->node[j]->Id_loc+1;
      }
    }
    break;
  case quad8:
    for (i=0;i<field->numele;i++) {
      ELEMENT* actele = &(discret->element[i]);
      inel++;
      /*----------------------------------------- sub-element 1 */
      PCELL[inel][0] = actele->node[0]->Id_loc+1;
      PCELL[inel][1] = actele->node[4]->Id_loc+1;
      PCELL[inel][2] = actele->node[7]->Id_loc+1;
      PCELL[inel][3] = actele->node[0]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 2 */
      PCELL[inel][0] = actele->node[4]->Id_loc+1;
      PCELL[inel][1] = actele->node[1]->Id_loc+1;
      PCELL[inel][2] = actele->node[5]->Id_loc+1;
      PCELL[inel][3] = actele->node[4]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 3 */
      PCELL[inel][0] = actele->node[5]->Id_loc+1;
      PCELL[inel][1] = actele->node[2]->Id_loc+1;
      PCELL[inel][2] = actele->node[6]->Id_loc+1;
      PCELL[inel][3] = actele->node[5]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 4 */
      PCELL[inel][0] = actele->node[6]->Id_loc+1;
      PCELL[inel][1] = actele->node[3]->Id_loc+1;
      PCELL[inel][2] = actele->node[7]->Id_loc+1;
      PCELL[inel][3] = actele->node[6]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 5 */
      PCELL[inel][0] = actele->node[6]->Id_loc+1;
      PCELL[inel][1] = actele->node[7]->Id_loc+1;
      PCELL[inel][2] = actele->node[4]->Id_loc+1;
      PCELL[inel][3] = actele->node[5]->Id_loc+1;
    }
    break;
  case quad9:
    for (i=0;i<field->numele;i++) {
      ELEMENT* actele = &(discret->element[i]);
      inel++;
      /*----------------------------------------- sub-element 1 */
      PCELL[inel][0] = actele->node[0]->Id_loc+1;
      PCELL[inel][1] = actele->node[4]->Id_loc+1;
      PCELL[inel][2] = actele->node[8]->Id_loc+1;
      PCELL[inel][3] = actele->node[7]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 2 */
      PCELL[inel][0] = actele->node[4]->Id_loc+1;
      PCELL[inel][1] = actele->node[1]->Id_loc+1;
      PCELL[inel][2] = actele->node[5]->Id_loc+1;
      PCELL[inel][3] = actele->node[8]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 3 */
      PCELL[inel][0] = actele->node[8]->Id_loc+1;
      PCELL[inel][1] = actele->node[5]->Id_loc+1;
      PCELL[inel][2] = actele->node[2]->Id_loc+1;
      PCELL[inel][3] = actele->node[6]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 4 */
      PCELL[inel][0] = actele->node[7]->Id_loc+1;
      PCELL[inel][1] = actele->node[8]->Id_loc+1;
      PCELL[inel][2] = actele->node[6]->Id_loc+1;
      PCELL[inel][3] = actele->node[3]->Id_loc+1;
    }
    break;
  default:
    dserror("distyp %d not implemented yet!", distyp);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief  movie creation

  this function
  creates the file vis<ACTSTEP>.xwd via system command (screen-shot
  xwd -out ...)
  converts the xwd-file into a .gif file via system command (convert)
  deletes the xwd-file via system command (rm)

  see also /bau/statik/info/README.visual2_ccarat_movie

  \author genk
  \date 07/02
*/
/*----------------------------------------------------------------------*/
void v2movie()
{
  char string[100];
  char *charpointer;

  char convert[100];
  char remove[100];

#ifdef DEBUG
  dstrc_enter("v2movie");
#endif

/*--------------------------------------- add file counter to filenames */
#ifdef HPUX11
  printf("movie creation only possible under HPUX 10.20\n");
  goto end;
#endif

  sprintf(string ,"xwd -name Visual2  -out vis%04d.xwd", ACTSTEP);
  sprintf(convert ,"convert vis%04d.xwd vis%04d.png", ACTSTEP, ACTSTEP);
  sprintf(remove ,"rm vis%04d.xwd", ACTSTEP);

  printf("write vis%04d.png\n", ACTSTEP);

/*--------------------------------------------------- call UNIX-system */
  system(&string[0]);
  system(&convert[0]);
  system(&remove[0]);

/*----------------------------------------------------------------------*/
#ifdef HPUX11
end:
#endif
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Read the node displacements and put them to the node arrays.

  \param discret (o) discretization to be read
  \param result  (i) current result description
  \param place   (i) row in the node array where values go
  \param nsteps  (i) maximum number of steps

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static void post_read_displacement(POST_DISCRETIZATION* discret,
                                   RESULT_DATA* result,
                                   INT place,
                                   INT nsteps)
{
  CHUNK_DATA chunk;
  DOUBLE **array;
  INT i;
  INT j;

#ifdef DEBUG
  dstrc_enter("post_read_displacement");
#endif

  init_chunk_data(result, &chunk, "displacement");

  dsassert(chunk.value_entry_length==2, "2d problem expected");

  dsassert(place < nsteps, "nsteps exceeded");

  for (i=0; i<result->field->numnp; ++i)
  {
    if (discret->node[i].sol.Typ == cca_XX)
    {
      array = amdef("sol", &(discret->node[i].sol), nsteps, 2, "DA");
    }
    else if (place >= discret->node[i].sol.fdim)
    {
      array = amredef(&(discret->node[i].sol), nsteps, discret->node[i].sol.sdim, "DA");
    }
    else
    {
      array = discret->node[i].sol.a.da;
    }

    chunk_read_value_entry(&chunk, i);
    for (j=0; j<chunk.value_entry_length; ++j)
    {
      array[place][j] = chunk.value_buf[j];
    }
  }

  destroy_chunk_data(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Read the node velocities and pressure and put them to the
  node arrays.

  We have to read two chunks here. We do this one after the other, but
  we put all values of one step into one node array row.

  \param discret (o) discretization to be read
  \param result  (i) current result description
  \param place   (i) row in the node array where values go
  \param nsteps  (i) maximum number of steps

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static void post_read_vel_pres(POST_DISCRETIZATION* discret,
                               RESULT_DATA* result,
                               INT place,
                               INT nsteps)
{
  CHUNK_DATA chunk;
  DOUBLE **array;
  INT i;
  INT j;

#ifdef DEBUG
  dstrc_enter("post_read_vel_pres");
#endif

  /*--------------------------------------------------------------------*/
  /* read the velocity */

  init_chunk_data(result, &chunk, "velocity");

  dsassert(chunk.value_entry_length==2, "2d problem expected");

  dsassert(place < nsteps, "nsteps exceeded");

  for (i=0; i<result->field->numnp; ++i)
  {
    if (discret->node[i].sol.Typ == cca_XX)
    {
      array = amdef("sol", &(discret->node[i].sol), nsteps, 3, "DA");
    }
    else if (place >= discret->node[i].sol.fdim)
    {
      array = amredef(&(discret->node[i].sol), nsteps, discret->node[i].sol.sdim, "DA");
    }
    else
    {
      array = discret->node[i].sol.a.da;
    }

    chunk_read_value_entry(&chunk, i);
    for (j=0; j<chunk.value_entry_length; ++j)
    {
      array[place][j] = chunk.value_buf[j];
    }
  }

  destroy_chunk_data(&chunk);

  /*--------------------------------------------------------------------*/
  /* read the pressure */

  init_chunk_data(result, &chunk, "pressure");

  dsassert(chunk.value_entry_length==1, "there must be just one pressure value");

  for (i=0; i<result->field->numnp; ++i)
  {
    chunk_read_value_entry(&chunk, i);
    array = discret->node[i].sol.a.da;
    array[place][2] = chunk.value_buf[0];
  }

  destroy_chunk_data(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Helper for user interaction.

  Supplied with a question string and a list of options this function
  asks the user and returns the choice he made.

  The idea is to have the user interaction a little robust.

  \param question (i) to be asked; might contain newlines
  \param options  (i) NULL terminated list of options
  \param keys     (i) NULL terminated; one char per option
  \param default_key (i) the option choosen upon no input

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static INT user_question(CHAR* question,
                         CHAR** options,
                         CHAR* keys,
                         CHAR default_option)
{
  INT ch;
  INT c;
  INT i;
  INT len;

#ifdef DEBUG
  dstrc_enter("user_question");
#endif

  len = strlen(keys);
  for (;;)
  {

    /* prepare screen */
    printf(question);
    printf("\n");
    for (i=0; i<len; ++i)
    {
      printf("    ");
      if (keys[i] == default_option)
      {
        printf("[%c]", keys[i]);
      }
      else
      {
        printf(" %c ", keys[i]);
      }
      printf(": %s\n", options[i]);
    }
    printf("\n");

    /* handle input */
    while (isblank(ch = getchar()))
    {
      /* ignore spaces and tabs */
    }

    if ((ch == '\n') || (ch == EOF))
    {
      return default_option;
    }

    while (((c = getchar()) != '\n') && (c != EOF))
    {
      /* consume the rest of this line */
    }

    for (i=0; i<len; ++i)
    {
      if (keys[i] == ch)
      {
        /* found it */
        return ch;
      }
    }
    printf("\nDon't know how to handle '%c'.\nTry again\n\n", ch);
  }
  /* never reached */

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Prepare and run visual2

  This is called when the data has been loaded.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static void post_visual2()
{
  INT i;
  INT distype;
  INT ch;

  CHAR* color_options[] = { "colours    - black background",
                            "colours    - white background",
                            "grey scale - white background",
                            NULL };

#ifdef DEBUG
  dstrc_enter("post_visual2");
#endif

  /*---------------------------------------------------------------------*
   * input from the screen                                               *
   *---------------------------------------------------------------------*/

  /*----------------------------------------------- data structure sizes */
  MPTRI  = 1000000;
  MPPTRI = 35000;
  MFACE  = 1750000;
  MPFACE = 100000;
  MEDGE  = 100000;
  MPEDGE = 10000;
  MNODE  = fluid_field->numnp;

  for (;;) {
    printf("\n");
    printf("     Actual size of data structures:\n");
    printf("      MNODE  = %d\n",MNODE);
    printf("      MPTRI  = %d\n",MPTRI);
    printf("      MPPTRI = %d\n",MPPTRI);
    printf("      MFACE  = %d\n",MFACE);
    printf("      MPFACE = %d\n",MPFACE);
    printf("      MEDGE  = %d\n",MEDGE);
    printf("      MPEDGE = %d\n",MPEDGE);
    printf("\n");

    if (user_question("New size of data structures?",
                      yesno_options, "01", '0') == '1') {
      printf("MPTRI  =\n");
      scanf("%d",&MPTRI);
      printf("MPPTRI =\n");
      scanf("%d",&MPPTRI);
      printf("MFACE  =\n");
      scanf("%d",&MFACE);
      printf("MPFACE =\n");
      scanf("%d",&MPFACE);
      printf("MEDGE  =\n");
      scanf("%d",&MEDGE);
      printf("MPEDGE =\n");
      scanf("%d",&MPEDGE);
    }
    else {
      break;
    }
  }

  /*----------------------------------------------- calculate vorticity */
  ch = user_question("Do you want to calculate the vorticity?",
                     yesno_options, "01", '0');
  /* let's trust the user_question function and use the shortcut */
  IVORT = ch - '0';

  if (IVORT==1) {
    printf("\n");
    printf("   Calculating the vorticity...\n");

    dserror("vorticity not ported yet");

#if 0
    fdyn = alldyn[numff].fdyn;
    fdyn->nsteps=nsteps;
    *action = calc_fluid_initvort;
    calinit(actfield,actpart,action,&container);
    *action = calc_fluid_vort;
    container.actndis=0;
#ifdef D_FLUID
    container.nif=0;
    container.nii=0;
    container.nim=0;
    container.is_relax = 0;
#endif
    container.fieldtyp=fluid;
    container.dvec=NULL;
    container.dvec         = NULL;
    container.actndis  = 0;
    calelm(actfield,NULL,actpart,NULL,0,0,
           &container,action);
#endif
  } /* endif (IVORT==1) */

  for (;;)
  {
    static CHAR* run_mode_options[] =
      {
      "steady data structure, grid and variables",
      "steady data structure and grid, unsteady variables",
      "steady data structure, unsteady grid and variables",
      "unsteady data structure, grid and variables",
      NULL
    };

    ch = user_question("Please give the mode in which you want to run VISUAL2:",
                       run_mode_options, "0123", '2');
    /* let's trust the user_question function and use the shortcut */
    IOPT = ch - '0';

    if (IOPT==3)
    {
      printf("   sorry, mode not implemented yet - new input!\n");
      continue;
    }
    break;
  }

  if (IOPT==0)
    icol++;

  if ((IOPT==2) && (ale_field == NULL))
  {
    printf("Without ale there's nothing unsteady in the grid.\n"
           "Switching back to option 1.\n\n");
    IOPT = 1;
  }

  /*--------------------------------------------------------- y-scaling */
  ch = user_question("Scaling Geometry?",
                     yesno_options, "01", '0');
  /* let's trust the user_question function and use the shortcut */
  SCAL = ch - '0';

  if (SCAL==1)
  {
    printf("\n");
    printf("     Scaling parameter for y-coordinates?\n");
    scanf("%lf",&yscale);

    /* eat the final newline */
    getchar();
  }
  else
    yscale=ONE;

/*--------------------------------------------------- modiy last step */
  LASTSTEP=step_a.a.iv[nsteps-1];

  ch = user_question("Colours?", color_options, "012", '0');
  /* let's trust the user_question function and use the shortcut */
  bgcolour = ch - '0';

  /* These are hacks to make visual2 work in linux. There are hard
   * coded filenames in ccarat's fortran part (vis2_qat2v2.f) */
#ifdef SUSE73
bgcolour+=10;
#endif
#ifdef LINUX_MUENCH
bgcolour+=20;
#endif

  /* all element types are cut to 3- or 4-node elements */
  distype = discret[fluid_idx].element[0].distyp;
  switch (distype)
  {
  case tri3:
    KCELL = fluid_field->numele;
    break;
  case tri6:
    KCELL = 4*fluid_field->numele;
    break;
  case quad4:
    KCELL = fluid_field->numele;
    break;
  case quad8:
    KCELL = 5*fluid_field->numele;
    break;
  case quad9:
    KCELL = 4*fluid_field->numele;
    break;
  default:
    dserror("distyp %d not implemented yet!", distype);
  }

  /* define arrays */
  PCELL = amdef("PCELL",&PCELL_A,KCELL,4,"IA");
  WCELL = amdef("WCELL",&WCELL_A,KCELL,1,"IV");
  WNODE = amdef("WNODE",&WNODE_A,fluid_field->numnp,1,"IV");
  WFACE = amdef("WFACE",&WFACE_A,6*MFACE,1,"IV");
  CEDGE = amdef("CEDGE",&CEDGE_A,MEDGE,1,"IV");
  amzero(&PCELL_A);
  amzero(&WCELL_A);
  amzero(&WNODE_A);
  amzero(&WFACE_A);
  amzero(&CEDGE_A);

  /* check if all elements are the same distyp */
  for (i=1; i<fluid_field->numele; i++)
  {
    ELEMENT* actele = &(discret[fluid_idx].element[i]);
    if (actele->distyp!=distype)
      dserror("up to now, all elements have to be the same distyp!\n");
  }

  /* set up modified element topology */
  post_v2cell(fluid_field, &(discret[fluid_idx]), distype);

  /*--------------------------------------------------------------------*/
  /* get the x/y limits */

  /* find maximum x-coordinate */
  if (IOPT>=2)
  {
#ifdef D_FSI
    NODE* actnode = &(discret[fluid_idx].node[0]);
    minxy = actnode->x[0];
    maxxy = actnode->x[0];
    for (i=0;i<fluid_field->numnp;i++)
    {
      actnode = &(discret[fluid_idx].node[i]);
      if (fluid_ale_connect[i] == -1)
      {
        minxy = MIN(minxy, actnode->x[0]);
        maxxy = MAX(maxxy, actnode->x[0]);
      }
      else {
        INT j;
        DOUBLE xy;
        NODE* actanode = &(discret[ale_idx].node[fluid_ale_connect[i]]);
        xy = actanode->x[0];
        for (j=0;j<nsteps;j++) {
          DOUBLE dxy;
          dxy = xy+actanode->sol.a.da[j][0];
          minxy=DMIN(minxy,dxy);
          maxxy=DMAX(maxxy,dxy);
        }
      }
    }
#else
    dserror("FSI-functions not compiled in!");
#endif
  }
  else {
    NODE* actnode = &(discret[fluid_idx].node[0]);
    minxy = actnode->x[0];
    maxxy = actnode->x[0];
    for (i=1;i<fluid_field->numnp;i++) {
      actnode = &(discret[fluid_idx].node[i]);
      minxy = MIN(minxy, actnode->x[0]);
      maxxy = MAX(maxxy, actnode->x[0]);
    }
  }
  XYMIN[0] = minxy - (maxxy-minxy)*0.1;
  XYMAX[0] = maxxy + (maxxy-minxy)*0.1;

  /* find maximum y-coordinate */
  if (IOPT>=2) {
#ifdef D_FSI
    NODE* actnode = &(discret[fluid_idx].node[0]);
    minxy = actnode->x[1];
    maxxy = actnode->x[1];
    for (i=0;i<fluid_field->numnp;i++) {
      actnode = &(discret[fluid_idx].node[i]);
      if (fluid_ale_connect[i] == -1) {
        minxy = MIN(minxy, actnode->x[1]);
        maxxy = MAX(maxxy, actnode->x[1]);
      }
      else {
        INT j;
        DOUBLE xy;
        NODE* actanode = &(discret[ale_idx].node[fluid_ale_connect[i]]);
        xy = actanode->x[0];
        for (j=0;j<nsteps;j++) {
          DOUBLE dxy;
          dxy = xy+actanode->sol.a.da[j][0];
          minxy=DMIN(minxy,dxy);
          maxxy=DMAX(maxxy,dxy);
        }
      }
    }
#else
    dserror("FSI-functions not compiled in!");
#endif
  }
  else {
    NODE* actnode = &(discret[fluid_idx].node[0]);
    minxy = actnode->x[1];
    maxxy = actnode->x[1];
    for (i=1;i<fluid_field->numnp;i++) {
      actnode = &(discret[fluid_idx].node[i]);
      minxy = MIN(minxy, actnode->x[1]);
      maxxy = MAX(maxxy, actnode->x[1]);
    }
  }
  XYMIN[1] = minxy - (maxxy-minxy)*0.1;
  XYMAX[1] = maxxy + (maxxy-minxy)*0.1;

  dx = XYMAX[0] - XYMIN[0];
  dy = XYMAX[1] - XYMIN[1];
  dxy= sqrt(dx*dx+dy*dy);

  /* window size */
  /*
                         qntsc:   352x240 (NTSC quarter screen)
                         qpal:    352x288 (PAL quarter screen)
                         ntsc:    720x480 (standard NTSC)
                         pal:     720x576 (standard PAL)
                         sntsc:   640x480 (square pixel NTSC)
                         spal:    768x576 (square pixel PAL)
   */
  /*
  hsize = 720;
  vsize = 576;
  ratio = 720.0/576.0;
  */

  /* write the images (for the movie) in double size */
  hsize = 1024;
  vsize = 768;
  ratio = 1024.0/768.0;

/*------------------------ modify x/y - limits to window aspect ratio */
  centerx =XYMIN[0]+dx/TWO;
  centery =XYMIN[1]+dy/TWO;
  if (dx>=dy)
  {
    XYPIX[0] = hsize;
    XYPIX[1] = vsize;
    if (dx/dy<ratio)
    {
      dx = dy*ratio;
      XYMIN[0] = centerx - dx/TWO;
      XYMAX[0] = centerx + dx/TWO;
    }
    else
    {
      dy = dx/ratio;
      XYMIN[1] = centery - dy/TWO;
      XYMAX[1] = centery + dy/TWO;
    }
  }
  else
  {
    XYPIX[0]=vsize;
    XYPIX[1]=hsize;
    if (dy/dx<ratio)
    {
      dy = dx*ratio;
      XYMIN[1] = centery - dy/TWO;
      XYMAX[1] = centery + dy/TWO;
    }
    else
    {
      dx = dy/ratio;
      XYMIN[0] = centerx - dx/TWO;
      XYMAX[0] = centerx + dx/TWO;
    }
  }
  XYMIN[1] = XYMIN[1]/yscale;
  XYMAX[1] = XYMAX[1]/yscale;

  /* get the data limits */
  {
    NODE* actnode = &(discret[fluid_idx].node[0]);
    velx    = actnode->sol.a.da[0][0];
    vely    = actnode->sol.a.da[0][1];
    absv    = sqrt(velx*velx+vely*vely);
    pres    = actnode->sol.a.da[0][2];
    if (IVORT==1)
    {
      vort = actnode->sol.a.da[0][3];
    }
  }

  minvx   = velx;
  maxvx   = velx;
  minvy   = vely;
  maxvy   = vely;
  minpre  = pres;
  maxpre  = pres;
  minabsv = sqrt(velx*velx+vely*vely);
  maxabsv = sqrt(velx*velx+vely*vely);
  if (IVORT==1)
  {
    minvort = vort;
    maxvort = vort;
  }

  for (i=0;i<fluid_field->numnp;i++) /* loop nodes */
  {
    INT j;
    NODE* actnode = &(discret[fluid_idx].node[i]);
    for (j=0;j<nsteps;j++) /* loop columns in sol-history */
    {
      velx    = actnode->sol.a.da[j][0];
      vely    = actnode->sol.a.da[j][1];
      pres    = actnode->sol.a.da[j][2];
      absv    = sqrt(velx*velx+vely*vely);
      if(IVORT==1)
      {
        vort = actnode->sol.a.da[j][3];
      }
      minvx   = DMIN(minvx  ,velx);
      maxvx   = DMAX(maxvx  ,velx);
      minvy   = DMIN(minvy  ,vely);
      maxvy   = DMAX(maxvy  ,vely);
      minpre  = DMIN(minpre ,pres);
      maxpre  = DMAX(maxpre ,pres);
      minabsv = DMIN(minabsv,absv);
      maxabsv = DMAX(maxabsv,absv);
      if (IVORT==1)
      {
        minvort = DMIN(minvort,vort);
        maxvort = DMAX(maxvort,vort);
      }
    } /* end of loop over columns */
  } /* end of loop over nodes */

#if 0
  if (fdyn->turbu == 2 || fdyn->turbu == 3)
  {
    ITURBU = 1;
    for (i=0;i<numnp;i++) /* loop nodes */
    {
      actnode=&(actfield->dis[1].node[i]);
      for (j=0;j<nsteps;j++) /* loop columns in sol-history */
      {
        kappa      = actnode->sol.a.da[j][0];
        dissi      = actnode->sol.a.da[j][2];
        visco      = actnode->sol.a.da[j][1];

        minkappa   = DMIN(minkappa   ,kappa);
        maxkappa   = DMAX(maxkappa   ,kappa);
        mindissi   = DMIN(mindissi   ,dissi);
        maxdissi   = DMAX(maxdissi   ,dissi);
        minvisco   = DMIN(minvisco   ,visco);
        maxvisco   = DMAX(maxvisco   ,visco);
      } /* end of loop over columns */
    } /* end of loop over nodes */
  } /* end of container.turbu==2 || container.turbu==3 */
  else
#endif
  {
    ITURBU = 0;
    minkappa   = 0.0;
    maxkappa   = 1.0;
    mindissi   = 0.0;
    maxdissi   = 1.0;
    minvisco   = 0.0;
    maxvisco   = 1.0;
  }
/* grid velocity not available yet !!!! */
  mingv = 0.0;
  maxgv = 1.0;

/*--------------------------------------- store max/min data in FLIMS */
  if (maxpre==minpre) maxpre=minpre+1.0;
  FLIMS[0][0]=minpre;
  FLIMS[0][1]=maxpre;

  FLIMS[1][0]=0.0;
  FLIMS[1][1]=1.0;

  if (minvort==maxvort) maxvort=minvort+1.0;
  FLIMS[2][0]=minvort;
  FLIMS[2][1]=maxvort;

  if (minabsv==maxabsv) maxabsv=minabsv+1.0;
  FLIMS[3][0]=dxy/maxabsv;
  FLIMS[3][1]=0.0;
  FLIMS[6][0]=minabsv;
  FLIMS[6][1]=maxabsv;

  if (minvx==maxvx) maxvx=minvx+1.0;
  FLIMS[4][0]=minvx;
  FLIMS[4][1]=maxvx;

  if (minvy==maxvy) maxvy=minvy+1.0;
  FLIMS[5][0]=minvy;
  FLIMS[5][1]=maxvy;

  FLIMS[7][0]=FLIMS[0][0];
  FLIMS[7][1]=FLIMS[0][1];

  FLIMS[8][0]=mingv;
  FLIMS[8][1]=maxgv;

  FLIMS[9][0]=ZERO;
  FLIMS[9][1]=ZERO;

  FLIMS[10][0]=ZERO;
  FLIMS[10][1]=ZERO;

  if (minkappa==maxkappa) maxkappa=minkappa+1.0;
  FLIMS[11][0]=minkappa;
  FLIMS[11][1]=maxkappa;

  if (mindissi==maxdissi) maxdissi=mindissi+1.0;
  FLIMS[12][0]=mindissi;
  FLIMS[12][1]=maxdissi;

  if (minvisco==maxvisco) maxvisco=minvisco+1.0;
  FLIMS[13][0]=minvisco;
  FLIMS[13][1]=maxvisco;

  actfieldtyp = fluid;

/*-------------------------- call Fortran routine which calls VISUAL2 */
  v2call(&IOPT,&CMNCOL,&CMUNIT,
         XYPIX,XYMIN,XYMAX,
         &NKEYS,IKEYS,FKEYS,FLIMS,
         &MNODE,&MPTRI,&MPPTRI,
         &MFACE,&MPFACE,&MEDGE,&MPEDGE,&bgcolour);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief The filter's main function.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
int main(int argc, char** argv)
{
  PROBLEM_DATA problem;
  INT i;
  INT res_count;
  INT counter1;
  RESULT_DATA result;
  CHAR gimmick[]="|/-\\";
  INT gimmick_size=4;

  init_problem_data(&problem, argc, argv);

  if (!map_has_int(&(problem.control_table), "ndim", 2))
  {
    dserror("two dimensional problem expected");
  }

  /*--------------------------------------------------------------------*/
  /* Find the corresponding discretizations. We don't rely on any order. */

  for (i=0; i<problem.num_discr; ++i)
  {
    if (problem.discr[i].type == structure)
    {
      struct_field = &(problem.discr[i]);
      struct_idx = i;
    }
    else if (problem.discr[i].type == fluid)
    {
      fluid_field = &(problem.discr[i]);
      fluid_idx = i;
    }
    else if (problem.discr[i].type == ale)
    {
      ale_field = &(problem.discr[i]);
      ale_idx = i;
    }
    else
    {
      dserror("unknown field type %d", problem.discr[i].type);
    }
  }

  /*--------------------------------------------------------------------*/
  /* some tests */

  switch (problem.type)
  {
#ifdef D_FSI
  case prb_fsi:
    if (problem.num_discr != 3)
    {
      dserror("expect 3 discretizations for fsi not %d", problem.num_discr);
    }
    if ((fluid_field == NULL) || (ale_field == NULL) || (struct_field == NULL))
    {
      dserror("structure, fluid and ale field expected");
    }
    break;
#endif
  case prb_fluid:
    if (problem.num_discr == 1)
    {
      if (problem.discr[0].type != fluid)
      {
        dserror("fluid discretization expected");
      }
    }
    else if (problem.num_discr == 2)
    {
      if ((fluid_field == NULL) || (ale_field == NULL) || (struct_field != NULL))
      {
        dserror("fluid and ale field expected");
      }
    }
    else
    {
      dserror("invalid number of discretizations for fluid problem (%d)", problem.num_discr);
    }
    break;
  default:
    dserror("problem type %d not supported", problem.type);
  }

  /*--------------------------------------------------------------------*/
  /* set up the meshes */

  discret = (POST_DISCRETIZATION*)CCACALLOC(problem.num_discr,
                                            sizeof(POST_DISCRETIZATION));

  /* Iterate all discretizations. */
  for (i=0; i<problem.num_discr; ++i)
  {
    init_post_discretization(&(discret[i]), &problem, &(problem.discr[i]));
  }

  /*--------------------------------------------------------------------*/
  /* Find the number of steps. We assume that there's one result group
   * for every discretization and every step. */

  res_count = map_symbol_count(&(problem.control_table), "result");
  nsteps = res_count / problem.num_discr;
  INCRE=1;
  if ((res_count % problem.num_discr) != 0)
  {
    dserror("the number of result groups (%d) doesn't match the number of discretizations (%d)",
            res_count, problem.num_discr);
  }

  /*--------------------------------------------------------------------*/
  /* Some user interaction. We need to know what range we have to
   * visualize. */

  printf("There are %d sets of results\n", nsteps);
#if 0
  if (user_question("Do you want to see all results?",
                    yesno_options, "01", '1') == '0')
  {
    for (;;)
    {
      printf("     min = 0, max = %d\n", nsteps-1);
      printf("     At which step shall the visualisation start?\n");
      scanf("%d",&FIRSTSTEP);
      printf("     At which step shall the visualisation end?\n");
      scanf("%d",&LASTSTEP);
      if (LASTSTEP<=FIRSTSTEP || FIRSTSTEP<0 || LASTSTEP >= nsteps)
      {
        printf("   Input out of range --> try again!\n");
        continue;
      }
      break;
    }

    nsteps = LASTSTEP - FIRSTSTEP + 1;

    /* increment */
    for (;;)
    {
      printf("\n");
      printf("     Increment step number by ...? (<%d)\n",nsteps-1);
      scanf("%d",&DSTEP);
      if (DSTEP > nsteps-1 || DSTEP==0)
      {
        printf("   Increment out of range --> Try again!\n");
        continue;
      }
      break;
    }

    nsteps /= DSTEP;

    /* eat the final newline */
    getchar();
  }
  else
  {
    FIRSTSTEP  = 0;
    LASTSTEP = nsteps-1;
    DSTEP    = 1;
  }
#endif

  if (problem.start > nsteps)
  {
    printf("Start step %d requested but only %d steps available. Corrected.\n", problem.start, nsteps);
    problem.start = nsteps;
  }

  if (problem.end > nsteps)
  {
    printf("End step %d requested but only %d steps available. Corrected.\n", problem.end, nsteps);
    problem.end = nsteps;
  }

  FIRSTSTEP = problem.start;
  DSTEP     = problem.step;
  if (problem.end > -1)
  {
    LASTSTEP  = problem.end;
    nsteps    = (problem.end - problem.start) / problem.step + 1;
  }
  else
  {
    LASTSTEP  = nsteps;
    nsteps    = (nsteps - problem.start) / problem.step + 1;
  }

  printf("Visualize step %d to step %d with %d steps\n",
         FIRSTSTEP, LASTSTEP, nsteps);

  /* create time array */
  amdef("time",&time_a,nsteps,1,"DV");

  /* create step array */
  amdef("step",&step_a,nsteps,1,"IV");

  /*--------------------------------------------------------------------*/
  /* Collect all step and time information. We use the fluid field's
   * results. There must be exactly one result group for the fluid
   * field per step*/

  /* Iterate all results. */
  init_result_data(fluid_field, &result);
  for (counter1 = 0; next_result(&result); counter1++)
  {
    if (counter1>=nsteps)
      dserror("too many fluid result steps");
    printf("Find number of results: %c\r", gimmick[counter1 % gimmick_size]);
    /*printf("% 2d: pos=%d\n", counter1, result.pos);*/
    time_a.a.dv[counter1] = map_read_real(result.group, "time");
    step_a.a.iv[counter1] = map_read_int(result.group, "step");
  }
  destroy_result_data(&result);
  nsteps = counter1;
  printf("Find number of results: done.\n");

  /*--------------------------------------------------------------------*/
  /* Now read the selected steps' results. */

  if (struct_field != NULL)
  {

    /* We iterate the list of all results. Here we are interested in
     * the results of this discretization. */
    init_result_data(struct_field, &result);
    for (counter1 = 0; next_result(&result); counter1++)
    {
      if (counter1>=nsteps)
        dserror("too many structure result steps");
      if (counter1 % 10 == 0)
      {
        printf("Read structure results: %c\r", gimmick[(counter1/10) % gimmick_size]);
        fflush(stdout);
      }
      post_read_displacement(&(discret[struct_idx]), &result, counter1, nsteps);
    }
    destroy_result_data(&result);
    if (nsteps != counter1)
    {
      dserror("too few structure results");
    }
    printf("Read structure results: done.\n");
  }
  if (ale_field != NULL)
  {

    /* We iterate the list of all results. Here we are interested in
     * the results of this discretization. */
    init_result_data(ale_field, &result);
    for (counter1 = 0; next_result(&result); counter1++)
    {
      if (counter1>=nsteps)
        dserror("too many ale result steps");
      if (counter1 % 10 == 0)
      {
        printf("Read ale results: %c\r", gimmick[(counter1/10) % gimmick_size]);
        fflush(stdout);
      }
      post_read_displacement(&(discret[ale_idx]), &result, counter1, nsteps);
    }
    destroy_result_data(&result);
    if (nsteps != counter1)
    {
      dserror("too few ale results");
    }
    printf("Read ale results: done.\n");
  }

  init_result_data(fluid_field, &result);
  for (counter1 = 0; next_result(&result); counter1++)
  {
    if (counter1>=nsteps)
      dserror("too many fluid result steps: panic");
    if (counter1 % 10 == 0)
    {
      printf("Read fluid results: %c\r", gimmick[(counter1/10) % gimmick_size]);
      fflush(stdout);
    }
    post_read_vel_pres(&(discret[fluid_idx]), &result, counter1, nsteps);
  }
  destroy_result_data(&result);
  printf("Read fluid results: done.\n");

#ifdef D_FSI
  /* Find coupled nodes. If there's at least an ale field. */
  if (ale_field != NULL)
  {
    post_log(1, "Find fsi coupling...");
    fflush(stdout);
    post_find_fsi_coupling(&problem,
                           struct_field, fluid_field, ale_field,
                           &fluid_struct_connect, &fluid_ale_connect);
    post_log(1, "\n");
  }
#endif

  post_visual2();

  post_log(4, "Done.\n");
  return 0;
}


/*----------------------------------------------------------------------*/
/* Evil hack just to make it link. */

void v3_init()
{
}
void v3_init_()
{
}
void v3_init__()
{
}
