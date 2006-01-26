/*!
\file
\brief Postprocessing utility that shows the results using visual3.
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

A huge part of this file was taken from the file visual_vis3.c.

\author m.geppert / u.kue
\date 01/06

*/

#include <ctype.h>
#include "post_visual3.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INCRE 1                      /*for V3UPDATE : number of steps
                                      * forward */
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



static INT      IOPT;                 /* program mode                   */
static INT      SCAL;
static INT      ACTDIM;
static INT      NKEYS=7;              /* see VISUAL3 manual             */
static INT      IKEYS[]={120,121,122,112,102,97,109};
static INT      FKEYS[]={1,1,1,1,2,1,1};
static float    FLIMS[7][2];          /* data limits                    */
static DOUBLE  minpre,maxpre;   /*					*/
static DOUBLE  minvx,maxvx;     /*					*/
static DOUBLE  minvy,maxvy;     /*					*/
static DOUBLE  minvz,maxvz;     /*					*/
static DOUBLE  minabsv,maxabsv; /*					*/
/*static INT     FIRSTSTEP;*/
static INT     LASTSTEP;
static INT     ACTSTEP;
static INT     bgcolour;        /*background colour                     */
static DOUBLE  velx;	        /*					*/
static DOUBLE  vely;	        /*					*/
static DOUBLE  velz;	        /*					*/
static DOUBLE  pres;	        /*					*/
static DOUBLE  absv;	        /*					*/

/*static DOUBLE  vort;          not implemented node data
  static DOUBLE  kappa;
  static DOUBLE  dissi;
  static DOUBLE  visco;*/

static INT      nsteps;
static INT      ndsurf;
static INT      maxface;
static INT      KSURF;          /*number of domain surface faces*/
static INT      KCEL1;   /*number of tetraheda  */
static INT      KCEL2;   /*number of pyramids   */
static INT      KCEL3;   /*number of prisms     */
static INT      KCEL4;   /*number of hexahedra  */
static INT      KEQUIV;  /*number of non equivalency pairs    */
static INT      KPTET;   /*number of cells in poly-tetrahedra */
static INT      KNBLOCK; /*number of structured blocks        */
static INT      KSURFELE; /*number of nodes of one single surface of a
                           * element*/
static INT      BLOCKS;  /*structured blocks definitions      */
static INT      KNSURF;  /*number of domain surface groups    */
static INT      KNPTET;  /*number polztetrahedral strips      */
static INT      dnumnp;  /*number of nodes given to V3INIT    */
static INT      WIN3D;
static INT      MIRROR;
static INT      CMUNIT=99;            /* see VISUAL3 manual */

/*------------------------------- variables needed for own calculations */
static INT     numnp;          /* number of nodes of actual field*/
static INT     numnp_struct;
static INT     numele;	        /* number of elements of actual field	*/
static INT     ncols=1;         /* number of sol steps stored in sol	*/
/*static INT     nacols=1;*/        /* number of sol steps of ALE field     */
static INT     icol=-1;         /* act. num. of sol step to be visual.  */
static INT     ACTSTEP;
static INT     IMOVIE=0;        /* counter to slow down for movie creat.*/
static INT     bgcolour;        /* background colour                    */
static DOUBLE  FACX,FACY,FACZ;

static ELEMENT       *actele;          /* actual element			*/
static NODE          *actnode;         /* actual node				*/
static DIS_TYP        distyp;          /* element type  			*/
static FIELDTYP       actfieldtyp;     /* actual fieldtyp			*/
static ARRAY          time_a ;         /* time array				*/
static ARRAY          step_a ;         /* time array				*/

/* static FIELD         *actfield;         actual field  			 */
/* static FIELD         *alefield; */
/* static FIELD         *structfield; */
/* static NODE          *actanode;         actual ale-node                       */
/* static GNODE         *actgnode;         actual gnode  			 */
/* static FLUID_DYNAMIC *fdyn;             pointer to fluid dyn. inp.data    */


/* Information about the nodes belonging to the faces of a hex8
 * element, used later in V3SURFACE:*/
static INT element_structure[6][4] = {
  {0, 1, 2, 3},
  {1, 2, 6, 5},
  {2, 3, 7, 6},
  {0, 3, 7, 4},
  {4, 5, 6, 7},
  {0, 1, 5, 4}
};


void vis3caf(INT numff, INT numaf, INT numsf)
{

  INT  i,j ;

  char TITL[10]="Picture";
  char TKEYS[NKEYS][16];
  strcpy(TKEYS[0], "VELOCITY Ux     ");
  strcpy(TKEYS[1], "VELOCITY Uy     ");
  strcpy(TKEYS[2], "VELOCITY Uz     ");
  strcpy(TKEYS[3], "PRESSURE        ");
  strcpy(TKEYS[4], "FLOW VECTORS    ");
  strcpy(TKEYS[5], "ABS VEL         ");
  strcpy(TKEYS[6], "MOVIE CREATION  ");

#ifdef LINUX_MUENCH
  char CMFILE[31]="/lnm/lib/visual/spec_black.col";
#else
  /* customize for Stuttgart */
  char CMFILE[31]="/lnm/lib/visual/spec_black.col";
#endif

  INT             screen;
  INT             dummy;

#ifdef DEBUG
  dstrc_enter("vis3caf");
#endif

/*--------------------------------------------------- set some values */
  numele = discret[fluid_idx].field->numele;
  numnp  = discret[fluid_idx].field->numnp;
  actele = &(discret[fluid_idx].element[0]);
  distyp = actele->distyp;

/*------------------------- check if all elements are the same distyp */
  for (i=1;i<numele;i++) /* loop all elements */
  {
    actele=&(discret[fluid_idx].element[i]);
    if (actele->distyp!=distyp)
      dserror("up to now, all elements have to be the same distyp!\n");
  } /* end loop over all elements */

  switch (distyp)
  {
    case hex8:
      KCEL1=0;
      KCEL2=0;
      KCEL3=0;
      KCEL4=numele;
      maxface=6;
      KSURF=numele*maxface;
      KSURFELE=4;   /*number of nodes belonging to a surface of this
                     * distyp */
      break;

    case quad4:
      KCEL1=0;
      KCEL2=0;
      KCEL3=0;
      KCEL4=numele;
      KSURF=numele*6;
      break;

    default:
      printf("%d", distyp);
      dserror("distyp not implemented yet!\n");
  } /* end switch(distyp) */

/*-------------------- read solution data from flavia.res-file of proc 0 */
/*visual_readflaviares(discret->field,&ncols,&time_a,&step_a,
  &FIRSTSTEP,&LASTSTEP,&DSTEP);*/


input1:
  printf("\n");
  printf("     Please give the mode in which you want to run VISUAL3:\n");
  printf("      0 : steady data structure, grid and variables\n");
  printf("      1 : steady data structure and grid, unsteady variables\n");
  printf("     [2]: steady data structure, unsteady grid and variables\n");

  screen=getchar();
  switch(screen)
  {
    case 10: IOPT=2; break;
    case 48: IOPT=0; dummy=getchar(); break;
    case 49: IOPT=1; dummy=getchar(); break;
    case 50: IOPT=2; dummy=getchar(); break;
    default:
      printf("\nTry again!\n");
      goto input1;
  }

  if (IOPT==0) icol++;

  if (IOPT==2) /*-------------------------- read ALE field from pss-file */
  {
#if 0
/*#ifdef D_FSI*/
    if (genprob.probtyp==prb_fsi)
    {
      if (numaf<0) dserror("ALE-field does not exist!");
      if (numsf<0) dserror("STRUCTURE-field does not exist!");
      alefield=&(field[numaf]);
      structfield=&(field[numsf]);
      fsi_initcoupling(structfield,discret->field,alefield);
      visual_readflaviares(alefield,&nacols,NULL,NULL,&FIRSTSTEP,&LASTSTEP,&DSTEP);
      if (ncols!=nacols)
      {
        printf("\n");
        printf("WARNING:number of columns different in ALE and FLUID field \n");
        printf("\n");
        ncols=IMIN(ncols,nacols);
      }
    }
    if (genprob.probtyp==prb_fluid)
    {
      if (numaf<0) dserror("ALE-field does not exist!");
      alefield=&(field[numaf]);
      fluid_initmfcoupling(discret->field,alefield);
      visual_readflaviares(alefield,&nacols,NULL,NULL,&FIRSTSTEP,&LASTSTEP,&DSTEP);
      if (ncols!=nacols)
      {
        printf("\n");
        printf("WARNING:number of columns different in ALE and FLUID field \n");
        printf("\n");
        ncols=IMIN(ncols,nacols);
      }
    }
#else
    dserror("FSI functions not compiled in!\n");
#endif
  }

  /*----------------------------------------------------- get window size */
input2:
  printf("\n");
  printf("     Size of Windows:\n");
  printf("     [0]: 3D window bigger than 2D window\n");
  printf("      1 : 2D window bigger than 3D window\n");

  screen=getchar();
  switch(screen)
  {
    case 10: WIN3D=1; break;
    case 48: WIN3D=1; dummy=getchar(); break;
    case 49: WIN3D=0; dummy=getchar(); break;
    default:
      printf("\nTry again!\n");
      goto input2;
  }

#if 0
/*----------------------------------------------------- get mirror flag */
input3:

/* --------------------------------------------------- get scaling facs */
#endif

input4:
  printf("\n");
  printf("     Scaling Geometry?\n");
  printf("     [0]: no \n");
  printf("      1 : yes\n");

  screen=getchar();
  switch(screen)
  {
    case 10: SCAL=0; break;
    case 48: SCAL=0; dummy=getchar(); break;
    case 49: SCAL=1; dummy=getchar(); break;
    default:
      printf("\nTry again!\n");
      goto input4;
  }

  if (SCAL==1)
  {
    printf("\n");
    printf("     Scaling Factors:\n\n");
    printf("     factor in x-direction:\n");
    scanf("%lf",&FACX);
    printf("     factor in y-direction:\n");
    scanf("%lf",&FACY);
    printf("     factor in z-direction:\n");
    scanf("%lf",&FACZ);
    dummy=getchar();
  }
  else
    FACX=FACY=FACZ=ONE;


/*------------------------------------------------- background colour */

/*----------------------------------------------- get the data limits */
  actnode=&(discret[fluid_idx].node[0]);
  velx    = actnode->sol.a.da[0][0];
  vely    = actnode->sol.a.da[0][1];
  velz    = actnode->sol.a.da[0][2];
  absv    = sqrt(velx*velx+vely*vely+velz*velz);
  pres    = actnode->sol.a.da[0][3];

  minvx   = velx;
  maxvx   = velx;
  minvy   = vely;
  maxvy   = vely;
  minvz   = velz;
  maxvz   = velz;
  minpre  = pres;
  maxpre  = pres;
  minabsv = absv;
  maxabsv = absv;

  for (i=0;i<numnp;i++) /* loop nodes */
  {
    actnode=&(discret[fluid_idx].node[i]);
    for (j=0;j<nsteps;j++) /* loop columns in sol-history */
    {
      velx    = actnode->sol.a.da[j][0];
      vely    = actnode->sol.a.da[j][1];
      velz    = actnode->sol.a.da[j][2];
      pres    = actnode->sol.a.da[j][3];
      absv    = sqrt(velx*velx+vely*vely+velz*velz);
      minvx   = DMIN(minvx  ,velx);
      maxvx   = DMAX(maxvx  ,velx);
      minvy   = DMIN(minvy  ,vely);
      maxvy   = DMAX(maxvy  ,vely);
      minvz   = DMIN(minvz  ,velz);
      maxvz   = DMAX(maxvz  ,velz);
      minpre  = DMIN(minpre ,pres);
      maxpre  = DMAX(maxpre ,pres);
      minabsv = DMIN(minabsv,absv);
      maxabsv = DMAX(maxabsv,absv);
    } /* end of loop over columns */
  } /* end of loop over nodes */

/*--------------------------------------- store max/min data in FLIMS */
  FLIMS[0][0]=minvx;
  FLIMS[0][1]=maxvx;

  FLIMS[1][0]=minvy;
  FLIMS[1][1]=maxvy;

  FLIMS[2][0]=minvz;
  FLIMS[2][1]=maxvz;

  FLIMS[3][0]=minpre;
  FLIMS[3][1]=maxpre;

  FLIMS[5][0]=minabsv;
  FLIMS[5][1]=maxabsv;

  FLIMS[4][0]=0.000000001;
  FLIMS[4][1]=1.000000001;

  FLIMS[6][0]=0.000000001;
  FLIMS[6][1]=1.000000001;

/*------------------------------------------------------ check limits */
  for (i=0;i<NKEYS;i++)
  {
    for (j=0;j<2;j++)
    {
      if (FABS(FLIMS[i][j])<EPS9)
      {
        if (FLIMS[i][j]<ZERO) FLIMS[i][j] = -0.000000001;
        else                  FLIMS[i][j] =  0.000000001;
      }
    }
    if (FLIMS[i][0]==FLIMS[i][1]) FLIMS[i][1] += ONE;
  }

/*-------------------------------------- set real number of last step */
  LASTSTEP=step_a.a.iv[ncols-1];


  KEQUIV=0;
  KNPTET=0;
  KPTET=0;
  KNBLOCK=0;
  BLOCKS=0;
  WIN3D=1;
  MIRROR=0;
  bgcolour=0;

#ifdef SUSE73
  bgcolour+=10;
#endif
#ifdef LINUX_MUENCH
  bgcolour+=20;
#endif

  if (ACTDIM==2)
  {
    dnumnp=2*numnp;     /*we need the double amount of nodes to create a
                         * 3D block out of a 2D surface */
    KNSURF=2;
  }

  if (ACTDIM==3)
  {
    dnumnp=numnp;
    map_find_int(fluid_field->group, "ndsurf", &ndsurf);
    KNSURF=ndsurf;
  }


  V3_INIT(TITL, &IOPT,
          CMFILE,
          &CMUNIT, &WIN3D, &NKEYS, IKEYS,
          &TKEYS[0][0],
          FKEYS, &FLIMS[0][0],
          &MIRROR, &dnumnp, &KEQUIV,
          &KCEL1, &KCEL2, &KCEL3,
          &KCEL4, &KNPTET, &KPTET,
          &KNBLOCK, &BLOCKS, &KSURF,
          &KNSURF,
          10, 31, 16);


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of vis3caf */


/*---------------------------------------------------------------------------*/
/*
 * Here we read the topology information and put them in an array.
 *
 * Topology information : every node has an array, which contains the number
 * of the surfaces the node is belonging to. Counting of surfaces
 * starts with 0, -1 means theres no more information in the array.
 * */

static void post_read_topology(INT id, INT* array)
{
  CHUNK_DATA chunk;
  INT  j;

#ifdef DEBUG
  dstrc_enter("post_read_topology");
#endif

  if (map_find_int(fluid_field->group, "ndsurf", &ndsurf)==0)
  {
    dserror("'ndsurf' missing in control file");
  }

  init_chunk_data(&fluid_field->head, &chunk, "dsurf");

  chunk_read_size_entry(&chunk, id);

  for (j=0;j<chunk.size_entry_length;j++)
  {
    array[j]=chunk.size_buf[j];
  }

  destroy_chunk_data(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}
/* end of post_read_topology */


/*----------------------------------------------------------------------------*/
/*
 * A small routine to check if a certain element face belongs to a
 * surface  or not.
 *
 * The four arrays are containing topology information of the four
 * nodes of one certain element face. We start in the first array and check if theres an
 * entry !=-1. For every time that value appears in one of the other
 * arrays a counter is raised by 1.
 *
 * So if the counter reaches 3 this element face got to be visible
 * from the outside, the routine returns the surface number the
 * element belongs to. */

static INT surface_check4(INT* array1, INT* array2, INT* array3, INT* array4)
{
  INT counter, surface;
  INT j;
  INT k;

#ifdef DEBUG
  dstrc_enter("surface_check");
#endif

  for (j=0;j<3;j++)
  {
    counter=0;

    surface=array1[j];

    if (surface!=-1)
    {
      for (k=0;k<3;k++)
      {
        if (array2[k]==surface) counter++;
        if (array3[k]==surface) counter++;
        if (array4[k]==surface) counter++;
      }

      if (counter==3)
      {
#ifdef DEBUG
        dstrc_exit();
#endif
        return surface;
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return -1;
}
/* end of surface_check */


/*-------------------------------------------------------------------------*/
/*
 * V3SURFACE :
 *
 * If we got a 2D problem, we only want to see the frontside of our 3D
 * block, so the other sides are turned off.
 *
 * In case of 3D problem we use the topology information of the
 * elements.
 *
 * At the beginning we do not really know how many element faces
 * belong to a surface, so we create temporary arrays to store the
 * data in, and later transfer them to the scel[][4] array.
 *
 * We loop over all elements to check if any of the elements faces
 * belongs to a surface, using post_read_topology and
 * surface_check.
 **/

void v3surface(INT nsurf[][2], INT surf[], INT scel[][4], CHAR tsurf[][20])
{

  INT i, j, k, l;

#ifdef DEBUG
  dstrc_enter("v3surface");
#endif

  if (ACTDIM==2)
  {
    char *name = "Frontside           "
                 "Backside            "
                 "Default             ";

    for (i=0;i<numele;i++)  /* frontside */
    {
      actele=&(discret[fluid_idx].element[i]);

      for (j=0;j<4;j++)
      {
        scel[i][j]=actele->node[j]->Id_loc+1;
      }
    }

    for (i=0;i<numele;i++)  /* backside */
    {
      actele=&(discret[fluid_idx].element[i]);

      for (j=0;j<4;j++)
      {
        scel[i+numele][j]=actele->node[j]->Id_loc+1+numnp;
      }
    }

    nsurf[0][0]=numele;
    nsurf[0][1]=1;

    nsurf[1][0]=2*numele;
    nsurf[1][1]=0;

    strncpy(tsurf[0], name, 60);

  }/*end of ACTDIM == 2 -------------------------------------------------------------*/

  else if (ACTDIM==3)
  {
    INT* surface[KNSURF];                 /* temporary array for nodes
                                           * belonging to surfaces */

    INT offset;                           /* needed for writing information in the scel array*/
    INT counter_array[KNSURF];            /* array counting the number of element faces on a certain surface */
    INT actsurf=-1;
    INT face;
    INT array[KSURFELE][3];
    INT surface_size[KNSURF];
    k=0;

    for (j=0;j<KNSURF;j++)
    {
      surface[j] = (INT*)CCACALLOC(KSURFELE*(INT)(KNSURF+1)*KSURF/KNSURF, sizeof(INT));
      counter_array[j]=0;
      surface_size[j]=KSURFELE*(INT)(KNSURF+1)*KSURF/KNSURF;
    }

    for (i=0;i<numele;i++)  /* loop over all elements */
    {
      for (face=0;face<maxface;face++) /* loop over all faces of a certain element */
      {
        for (j=0;j<KSURFELE;j++)
        {
          post_read_topology(discret[fluid_idx].element[i].node[element_structure[face][j]]->Id, array[j]);
        }

        switch(KSURFELE)
        {
          case 4:
            actsurf = surface_check4(array[0], array[1], array[2], array[3]);
            break;

          default :
            printf("not implemented yet\n");
            break;
        }

        if (actsurf!=-1)
        {
          if (counter_array[actsurf]+1>=surface_size[actsurf])
          {
            surface_size[actsurf]+=KSURFELE*(INT)KSURF/KNSURF;
            surface[actsurf] = (INT*)CCAREALLOC(surface[actsurf], surface_size[actsurf]*sizeof(INT));
          }
          for (j=0;j<KSURFELE;j++)
          {
            surface[actsurf][counter_array[actsurf]*4+j]=discret[fluid_idx].element[i].node[element_structure[face][j]]->Id+1;
          }
          counter_array[actsurf]++;
        }
      }/*end of loop over all faces of a certain element*/

    }/*end of loop over all elements*/

    offset=0;
    for (j=0;j<KNSURF;j++)    /*loop over all surface groups*/
    {
      for (k=0;k<counter_array[j];k++)  /*loop over all element faces of this surface group */
      {
        for (l=0;l<KSURFELE;l++)        /*loop over all nodes */
        {
          scel[k+offset][l]=surface[j][k*4+l];
        }
      }
      offset+=counter_array[j];
      nsurf[j][0]=offset;
      nsurf[j][1]=1;
      CCAFREE(surface[j]);
    }

  }/*end of ACTDIM==3 ----------------------------------------------------------------*/

  else
  {
    dserror("Ups!");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of V3SURFACE */



void v3grid(float XYZ[][3])
{
  INT i;

#ifdef DEBUG
  dstrc_enter("v3grid");
#endif

  switch(actfieldtyp)
  {
    case fluid:
#if 0
      if (NUMA>=0 && IOPT==0)
      {
      input4:
        v3_cursor(&true);
        printf("Give the column number (min. 0; max. %d): ?\n",ncols-1);
        scanf("%d",&icol);
        v3_cursor(&false);
        if (icol<0 || icol>ncols-1)
        {
          printf("Column number out of range. New input!\n");
          printf("\n");
          goto input4;
        }
      }
#endif
      if (IOPT==2)
      {
#if 0
/*#ifdef D_FSI*/
        for (i=0;i<numnp;i++)
        {
          actnode=&(discret[fluid_idx]->node[i]);
          actgnode=actnode->gnode;
          actanode=actgnode->mfcpnode[NUMA];
          if (actanode==NULL)
          {
            XYZ[i][0] = actnode->x[0]*FACX;
            XYZ[i][1] = actnode->x[1]*FACY;
            XYZ[i][2] = actnode->x[2]*FACZ;
          }
          else
          {
            XYZ[i][0] = (actanode->x[0]+actanode->sol.a.da[icol][0])*FACX;
            XYZ[i][1] = (actanode->x[1]+actanode->sol.a.da[icol][1])*FACY;
            XYZ[i][2] = (actanode->x[2]+actanode->sol.a.da[icol][2])*FACZ;
          }
        }
#else
        dserror("FSI-functions not compiled in!\n");
#endif
      }
      else
      {
        /*---------------------------------------------------------------------
         * V3GRID in case of 2D problem :
         *
         * we create a 3D structure by doubling the number of nodes. The first
         * half is used to create the frontside of a 3D structure, so we set
         * the z-coordinate to 0. The second half is an identical copy of the
         * first ones data, only the z-coordinate is set to 1 */

        if (ACTDIM==2)
        {
          for (i=0;i<numnp;i++)
          {
            actnode=&(discret[fluid_idx].node[i]);
            XYZ[i][0] = actnode->x[0]*FACX;
            XYZ[i][1] = actnode->x[1]*FACY;
            XYZ[i][2] = 0;

            XYZ[i+numnp][0] = actnode->x[0]*FACX;
            XYZ[i+numnp][1] = actnode->x[1]*FACY;
            XYZ[i+numnp][2] = 1*FACZ;
          }
        }

        if (ACTDIM==3)
        {
          for (i=0;i<numnp;i++)
          {
            actnode=&(discret[fluid_idx].node[i]);
            XYZ[i][0] = actnode->x[0]*FACX;
            XYZ[i][1] = actnode->x[1]*FACY;
            XYZ[i][2] = actnode->x[2]*FACZ;
          }
        }
      }
      break;
    case structure:
      for (i=0;i<numnp_struct;i++)
      {
        actnode=&(discret[fluid_idx].node[i]);
        XYZ[i][0] = (actnode->x[0]+actnode->sol.a.da[icol][0])*FACX;
        XYZ[i][1] = (actnode->x[1]+actnode->sol.a.da[icol][1])*FACY;
        XYZ[i][2] = (actnode->x[2]+actnode->sol.a.da[icol][2])*FACZ;
        XYZ[i+numnp_struct][0] = (actnode->x[0]+actnode->sol.a.da[icol][0])*FACX;
        XYZ[i+numnp_struct][1] = (actnode->x[1]+actnode->sol.a.da[icol][1])*FACY;
        XYZ[i+numnp_struct][2] = (actnode->x[2]+actnode->sol.a.da[icol][2])*FACZ;
      }
      break;
    case ale:
      dserror("fieldtyp not implemented yet!");

    default : dserror("error");
  } /* end switch(actfieldtyp) */

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of V3GRID */



/*!---------------------------------------------------------------------
  \brief set VISUAL3's scalar function values

  <pre>                                                         genk 01/04

  This routine is called by VISUAL3 during the visualisation.
  It filters the scalar values out of the nodal sol-field at the actual
  position icol.
  icol is set every timestep in v3update


  </pre>
  \param  *JKEY      INT    (i)  Index of key
  \param	*S         float  (o)  Scalar function at nodes or cells
  \return void
  \sa v2update

  ------------------------------------------------------------------------*/

void v3scal(INT *JKEY, float *S)
{
  INT i;
  float vx,vy, vz;
  float dx,dy, dz;


#ifdef DEBUG
  dstrc_enter("v3scal");
#endif

  switch(actfieldtyp)
  {
    case fluid:
#if 0
      if (IOPT==0) /* read in the column number */
      {
      input3:
        v3_cursor(&true);
        printf("Give the column number (min. 0; max. %d): ?\n",ncols-1);
        scanf("%d",&icol);
        v3_cursor(&false);
        if (icol<0 || icol>ncols-1)
        {
          printf("Column number out of range. New input!\n");
          printf("\n");
          goto input3;
        }
      }
#endif
      /*-------------------------------------------------- get scalar data */

      if(ACTDIM==2)
      {
        switch(*JKEY)
        {
          /*-------------------------------------------------------------------*/
          case 1: /* Ux */
            for (i=0;i<numnp;i++)
            {
              actnode=&(discret[fluid_idx].node[i]);
              S[i]=actnode->sol.a.da[icol][0];
              S[i+numnp]=actnode->sol.a.da[icol][0];
            }
            break;
            /*-------------------------------------------------------------------*/
          case 2: /* Uy */
            for (i=0;i<numnp;i++)
            {
              actnode=&(discret[fluid_idx].node[i]);
              S[i]=actnode->sol.a.da[icol][1];
              S[i+numnp]=actnode->sol.a.da[icol][1];
            }
            break;
            /*-------------------------------------------------------------------*/
          case 3: /* Uz */
            for (i=0;i<numnp;i++)
            {
              S[i]=0;
              S[i+numnp]=0;
            }
            break;
            /*--------------------------------------------------------------------*/
          case 4: /* pressure */
            for (i=0;i<numnp;i++)
            {
              actnode=&(discret[fluid_idx].node[i]);
              S[i]=actnode->sol.a.da[icol][2];
              S[i+numnp]=actnode->sol.a.da[icol][2];
            }
            break;
            /*-------------------------------------------------------------------*/
          case 5:
            printf("VECTOR DATA IN SCALAR ROUTINE ???????\n");
            break;
            /*-------------------------------------------------------------------*/
          case 6: /* absolute velocity */
            for (i=0;i<numnp;i++)
            {
              actnode=&(discret[fluid_idx].node[i]);
              vx=actnode->sol.a.da[icol][0];
              vy=actnode->sol.a.da[icol][1];
              S[i] = sqrt(vx*vx+vy*vy);
              S[i+numnp] = sqrt(vx*vx+vy*vy);
            }
            break;
            /*-------------------------------------------------------------------*/
          case 7:
            if (IMOVIE==0)
            {
              printf("\n");
              printf("   Starting movie creation at time 0.0\n");
              printf("\n");
              IMOVIE=1;
            }
            break;
            /*-------------------------------------------------------------------*/
        } /* end switch(*JKEY) */
      } /* end of 2D part */

      if(ACTDIM==3)
      {
        switch(*JKEY)
        {
          /*-------------------------------------------------------------------*/
          case 1: /* Ux */
            for (i=0;i<numnp;i++)
            {
              actnode=&(discret[fluid_idx].node[i]);
              S[i]=actnode->sol.a.da[icol][0];
            }
            break;
            /*-------------------------------------------------------------------*/
          case 2: /* Uy */
            for (i=0;i<numnp;i++)
            {
              actnode=&(discret[fluid_idx].node[i]);
              S[i]=actnode->sol.a.da[icol][1];
            }
            break;
            /*-------------------------------------------------------------------*/
          case 3: /* Uz */
            for (i=0;i<numnp;i++)
            {
              actnode=&(discret[fluid_idx].node[i]);
              S[i]=actnode->sol.a.da[icol][2];
            }
            break;
            /*-------------------------------------------------------------------*/
          case 4: /* pressure */
            for (i=0;i<numnp;i++)
            {
              actnode=&(discret[fluid_idx].node[i]);
              S[i]=actnode->sol.a.da[icol][3];
            }
            break;
            /*-------------------------------------------------------------------*/
          case 5:
            printf("VECTOR DATA IN SCALAR ROUTINE ???????\n");
            break;
            /*-------------------------------------------------------------------*/
          case 6: /* absolute velocity */
            for (i=0;i<numnp;i++)
            {
              actnode=&(discret[fluid_idx].node[i]);
              vx=actnode->sol.a.da[icol][0];
              vy=actnode->sol.a.da[icol][1];
              vz=actnode->sol.a.da[icol][2];
              S[i] = sqrt(vx*vx+vy*vy+vz*vz);
            }
            break;
            /*-------------------------------------------------------------------*/
          case 7:
            if (IMOVIE==0)
            {
              printf("\n");
              printf("   Starting movie creation at time 0.0\n");
              printf("\n");
              IMOVIE=1;
            }
            break;
            /*-------------------------------------------------------------------*/
        } /* end switch(*JKEY) */
      } /* end of 3D part */

      break;
    case structure:
      /*-------------------------------------------------- get scalar data */
      switch(*JKEY)
      {
        /*-------------------------------------------------------------------*/
        case 1: /* Dx */
          for (i=0;i<numnp_struct;i++)
          {
            actnode=&(discret[fluid_idx].node[i]);
            S[i]=actnode->sol.a.da[icol][0];
            S[i+numnp_struct]=actnode->sol.a.da[icol][0];
          }
          break;
          /*-------------------------------------------------------------------*/
        case 2: /* Dy */
          for (i=0;i<numnp;i++)
          {
            actnode=&(discret[fluid_idx].node[i]);
            S[i]=actnode->sol.a.da[icol][1];
            S[i+numnp_struct]=actnode->sol.a.da[icol][1];
          }
          break;
          /*-------------------------------------------------------------------*/
        case 3: /* Dz */
          for (i=0;i<numnp;i++)
          {
            actnode=&(discret[fluid_idx].node[i]);
            S[i]=actnode->sol.a.da[icol][2];
            S[i+numnp_struct]=actnode->sol.a.da[icol][2];
          }
          break;
          /*-------------------------------------------------------------------*/
        case 4: /* absolute displacement */
          for (i=0;i<numnp;i++)
          {
            actnode=&(discret[fluid_idx].node[i]);
            dx=actnode->sol.a.da[icol][0];
            dy=actnode->sol.a.da[icol][1];
            dz=actnode->sol.a.da[icol][2];
            S[i] = sqrt(dx*dx+dy*dy+dz*dz);
            S[i+numnp_struct] = sqrt(dx*dx+dy*dy+dz*dz);
          }
          break;
          /*-------------------------------------------------------------------*/
        case 5:
          if (IMOVIE==0)
          {
            printf("\n");
            printf("   Starting movie creation at time 0.0\n");
            printf("\n");
            IMOVIE=1;
          }
          break;
          /*-------------------------------------------------------------------*/
      } /* end switch(*JKEY) */
      break;
    case ale:
      dserror("fieldtyp not implemented yet!\n");
    default:
      dserror("fieldtyp unknown!\n");
  } /* end switch(actfieldtyp) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of V3SCAL */



/*!---------------------------------------------------------------------
  \brief set VISUAL3's vector function values

  <pre>                                                         genk 01/04

  This routine is called by VISUAL3 during the visualisation.
  It filters the vector values out of the nodal sol-field at the actual
  position icol.
  icol is set every timestep in v2update

  </pre>
  \param  *JKEY      INT    (i)  Index of key
  \param	*V         float  (o)  Vector function at nodes or cells
  \return void
  \sa v2update

  ------------------------------------------------------------------------*/
void v3vect(INT *JKEY, float V[][3])
{

  INT i;

#ifdef DEBUG
  dstrc_enter("V3VECT");
#endif

  switch(actfieldtyp)
  {
    case fluid:
      /*-------------------------------------------------- get vector data */
      if (*JKEY==5)  /* streamlines */
      {
        if (ACTDIM==2)
        {
          for (i=0;i<numnp;i++)
          {
            actnode=&(discret[fluid_idx].node[i]);
            V[i][0] = actnode->sol.a.da[icol][0];
            V[i][1] = actnode->sol.a.da[icol][1];
            V[i][2] = 0;
          }
        }

        if (ACTDIM==3)
        {
          for (i=0;i<numnp;i++)
          {
            actnode=&(discret[fluid_idx].node[i]);
            V[i][0] = actnode->sol.a.da[icol][0];
            V[i][1] = actnode->sol.a.da[icol][1];
            V[i][2] = actnode->sol.a.da[icol][2];
          }
        }
      }
      else
        printf("SCALAR DATA IN VECTOR ROUTINE ???????\n");

      break;
    case structure:
      /*-------------------------------------------------- get vector data */
      printf("NO VECTOR DATA AVAILABLE !!!!!!\n");
      break;
    case ale:
      dserror("fieldtyp not implemented yet!\n");

    default : dserror("error");
  } /* end switch(actfieldtyp) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of V3VECT */



void v3update(float *TIME)
{
#ifdef DEBUG
  dstrc_enter("V3UPDATE");
#endif

#if 0
/*-------------------------------------------------- check for stopping */
  if (icol==STSTEP)
  {
    printf("   Next Step to Stop (-1 to terminate)\n");
    scanf("%d",&STSTEP);
  }
#endif

#if 0
/*------------------------------------------- check for movie creation */
  if (icol==0 && IMOVIE==2)
  {
    IMOVIE=0;
    printf("\n");
    printf("Movie creation finished\n");
  }
  if (icol==0 && IMOVIE==1) IMOVIE=2;
  if (IMOVIE==2) v3movie();

#endif

  switch(actfieldtyp)
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
      break;
    case structure:
      if (icol==-1 || icol+INCRE > ncols-1)
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
      break;
    case ale:
      dserror("fieldtyp not implemented yet!\n");
    default : dserror("error");
  } /* end switch(actfieldtyp) */

  /* ToDo: Wait until at least 1/25 Seconds have passed */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of V3UPDATE */


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

void v3string(char STRING[80])
{
  INT len;
  float t;

#ifdef DEBUG
  dstrc_enter("V3STRING");
#endif

  t=time_a.a.dv[icol];
  sprintf(STRING, "Time: %8.4f    Step: %5d/%-5d", t, ACTSTEP, LASTSTEP);
  len = strlen(STRING);
  memset(STRING+len, ' ', 80-len);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of V3STRING */



void v3cell(INT cel1[][4], INT cel2[][5], INT cel3[][6], INT cel4[][8],
            INT nptet[][8], INT *ptet)
{

  INT i,j;

#ifdef DEBUG
  dstrc_enter("v3cell");
#endif

  switch (actfieldtyp)
  {
    case fluid:
      for (i=0;i<numele;i++) /* loop all elements */
      {
        actele=&(discret[fluid_idx].element[i]);

        switch (actele->distyp)
        {
          case hex8:
            for (j=0;j<8;j++)
            {
              actnode = actele->node[j];
              cel4[i][j]=actnode->Id_loc+1;
            }
            break;

          case quad4: /* 4 node rectangle */

            for(j=0;j<4;j++)
            {
              cel4[i][j] = actele->node[j]->Id_loc+1;
              cel4[i][j+4] = actele->node[j]->Id_loc+1+numnp;
            }

            break;

          default:
            dserror("distyp not implemented yet!\n");

        }
      } /* end loop over all elements */
      fflush(0);
      break;

    case structure:
      for (i=0;i<numele;i++) /* loop all elements */
      {
        actele=&(discret[fluid_idx].element[i]);
        switch (actele->distyp)
        {
          case quad4:
            for (j=0;j<4;j++)
            {
              actnode = actele->node[j];
              cel4[i][j]=actnode->Id_loc+1;
              cel4[i][j+4]=actnode->Id_loc+1+numnp_struct;
            }
            break;
          default:
            dserror("distyp not implemented yet!\n");

        }
      } /* end loop over all elements */
      break;
    default:
      dserror("fieldtyp unknown!\n");
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of V3CELL */


#if 0
int strcpb(char *str1, char *str2, int len)
{
  register int i,j;

  i = 0;
  while (str2[i] != '\0') {
    str1[i] = str2[i];
    i++;
  }
  for (j=i; j < len; j++) str1[j] = ' ';

  return i;
}
#endif


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

  dsassert(chunk.value_entry_length==3, "3d problem expected");

  dsassert(place < nsteps, "nsteps exceeded");

  for (i=0; i<result->field->numnp; ++i)
  {
    if (discret->node[i].sol.Typ == cca_XX)
    {
      array = amdef("sol", &(discret->node[i].sol), nsteps,3, "DA");
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


  if (ACTDIM==2)
  {
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
        array = amredef(&(discret->node[i].sol), nsteps,3, "DA");
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
      for (j=0;j<chunk.value_entry_length;++j)
      {
        array[place+numnp][j] = chunk.value_buf[j];
      }
    }

    destroy_chunk_data(&chunk);


    /* read the pressure */

    if (map_has_map(result->group, "pressure"))
    {
      init_chunk_data(result, &chunk, "pressure");
    }
    else if (map_has_map(result->group, "average_pressure"))
    {
      printf("No pressure entry found. Only average pressure.");
      init_chunk_data(result, &chunk, "average_pressure");
    }
    else
    {
      dserror("No pressure entry found.");
    }

    dsassert(chunk.value_entry_length==1, "there must be just one pressure value");

    for (i=0; i<result->field->numnp; ++i)
    {
      chunk_read_value_entry(&chunk, i);
      array = discret->node[i].sol.a.da;
      array[place][2] = chunk.value_buf[0];
      array = discret->node[i+numnp].sol.a.da;
      array[place][2] = chunk.value_buf[0];
    }
    destroy_chunk_data(&chunk);

  } /* end of 2D part*/

  if (ACTDIM==3)
  {
    init_chunk_data(result, &chunk, "velocity");

    dsassert(chunk.value_entry_length==3, "3d problem expected");

    dsassert(place < nsteps, "nsteps exceeded");

    for (i=0; i<result->field->numnp; ++i)
    {
      if (discret->node[i].sol.Typ == cca_XX)
      {
        array = amdef("sol", &(discret->node[i].sol), nsteps, 4, "DA");
      }
      else if (place >= discret->node[i].sol.fdim)
      {
        array = amredef(&(discret->node[i].sol), nsteps,4, "DA");
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

    if (map_has_map(result->group, "pressure"))
    {
      init_chunk_data(result, &chunk, "pressure");
    }
    else if (map_has_map(result->group, "average_pressure"))
    {
      printf("No pressure entry found. Only average pressure.");
      init_chunk_data(result, &chunk, "average_pressure");
    }
    else
    {
      dserror("No pressure entry found.");
    }

    dsassert(chunk.value_entry_length==1, "there must be just one pressure value");

    for (i=0; i<result->field->numnp; ++i)
    {
      chunk_read_value_entry(&chunk, i);
      array = discret->node[i].sol.a.da;
      array[place][3] = chunk.value_buf[0];
    }

    destroy_chunk_data(&chunk);
  } /* end of 3D part */

#ifdef DEBUG
  dstrc_exit();
#endif
}/* end of post_read_vel_pres*/

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

  actfieldtyp=fluid;
  init_problem_data(&problem, argc, argv);

  if (!map_has_int(&(problem.control_table), "ndim", 2))
  {
    printf("3D problem \n");
    ACTDIM=3;
  }
  if (!map_has_int(&(problem.control_table), "ndim",3 ))
  {
    printf("2D problem \n");
    ACTDIM=2;
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
#if 0
/*#ifdef D_FSI*/
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

  if ((res_count % problem.num_discr) != 0)
  {
    dserror("the number of result groups (%d) doesn't match the number of discretizations (%d)",
            res_count, problem.num_discr);
  }

  printf("There are %d sets of results\n", nsteps);


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
  for (counter1 = 0; next_result(&result)  ; counter1++)
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

#if 0
/*#ifdef D_FSI*/
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

  vis3caf(0,0,0);

  post_log(4, "Done.\n");
  return 0;
}
