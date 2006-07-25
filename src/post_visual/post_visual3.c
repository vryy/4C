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
#include <stdio.h>
#include "post_visual3.h"
#include "post_visual3_functions.h"
#include "post_design.h"
#include "../post_common/post_common.h"
#include "../pss_full/pss_set.h"
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "../post_common/post_octtree.h"


int strcpb(char *str1, char *str2, int len);

#define INCRE 1                      /*for V3UPDATE : number of steps
                                      * forward */
/*----------------------------------------------------------------------*/
/* the result values */
/*----------------------------------------------------------------------*/
static POST_DISCRETIZATION* discret;
static POST_DESIGN fluid_design;
static FIELD_DATA* struct_field = NULL;
static FIELD_DATA* fluid_field = NULL;
static FIELD_DATA* ale_field = NULL;
static INT struct_idx=-1;
static INT fluid_idx=-1;
static INT ale_idx=-1;

#ifdef D_FSI
static INT* fluid_ale_connect;
static DOUBLE *node_director;
#endif
/*----------------------------------------------------------------------*
 |                                                        genk 07/02    |
 | all variables needed by VISUAL2 are defined extern                   |
 *----------------------------------------------------------------------*/


static INT      IOPT=0;                 /* program mode                   */
static INT      ACTDIM;
static INT      NKEYS=12;              /* see VISUAL3 manual             */
static INT      IKEYS[]={120,121,122,112,102,97,109,98, 99, 100,101, 116};
static INT      FKEYS[]={1,1,1,1,2,1,1, 1, 1, 1,1, 1};
static float    FLIMS[12][2];          /* data limits                    */
static DOUBLE  minpre;   /*					*/
static DOUBLE  minvx;     /*					*/
static DOUBLE  minvy;     /*					*/
static DOUBLE  minvz;     /*					*/
static DOUBLE  minabsv; /*					*/
static DOUBLE  mindx;    /*                                      */
static DOUBLE  mindy;    /*                                      */
static DOUBLE  mindz;    /*                                      */
static DOUBLE  minabsd;
static INT     LASTSTEP;
static INT     ACTSTEP;
static INT      nsteps;
static INT      ndsurf_tot=0; /*total number of surfaces*/
static INT      ndvol;
static INT      ndvol_tot;
static INT      ndline_tot;
static INT      ndnode_tot=0;
static INT      ndsurf_in_tot=0;
static INT      FLUID_NODE_offset=0; /*number of nodes in XYZ BEFORE
                                      * struct nodes*/
static INT      FLUID_GROUP_offset=0;
static INT      FLUID_SURF_offset=0; /*number of surfaces in scel BEFORE struct surfaces*/
static INT      FLUID_CEL3_offset=0; /*number of cells in CEL3 BEFORE struct cells*/
static INT      FLUID_CEL4_offset=0; /*number of cells in CEL4 BEFORE struct cells*/
static INT      KSURF=0;          /*number of domain surface faces*/
static INT      KCEL1=0;          /*number of tetraheda  */
static INT      KCEL2=0;          /*number of pyramids   */
static INT      KCEL3=0;          /*number of prisms     */
static INT      KCEL4=0;          /*number of hexahedra  */
static INT      KEQUIV;         /*number of non equivalency pairs    */
static INT      KPTET;          /*number of cells in poly-tetrahedra */
static INT      KNBLOCK;        /*number of structured blocks        */
static INT      KSURFELE;       /*number of nodes of one single surface of a
                                 * element*/
static INT      BLOCKS;  /*structured blocks definitions      */
static INT      KNSURF=0;  /*number of domain surface groups    */
static INT      KNSURFT=0;
static INT      KNPTET;  /*number polytetrahedral strips      */
static INT      dnumnp;  /*number of nodes given to V3INIT    */
static INT      WIN3D;
static INT      MIRROR;
static INT      CMUNIT=99;            /* see VISUAL3 manual */
static CHAR     CMFILE[80];
static CHAR     frame_folder[100];
static INT      movie_window=3;
static INT      toggle_movie_window=1;
/*------------------------------- variables needed for own calculations */
static INT     numnp=0;          /* number of nodes of actual field*/
static INT     numnp2D=0;
static INT     numnp3D=0;
static INT     numnp_fluid=0;
static INT     numnp_struct=0;
static INT     numele;	        /* number of elements of actual field*/
static INT     numele_fluid=0;
static INT     numele_struct=0;
static INT     num_discr;
static NODE    *actanode;
static INT     ncols=1;         /* number of sol steps stored in sol	*/
static INT     icol=-1;         /* act. num. of sol step to be visual.  */
static INT     ACTSTEP;
static INT     COLOUR=0;
static INT     IMOVIE=0;        /* counter to slow down for movie creat.*/
static DOUBLE  FACX,FACY,FACZ;

static ELEMENT       *actele;
static NODE          *actnode;
static DIS_TYP        distyp;          /* element type  			*/
static FIELDTYP       actfieldtyp;
static ARRAY          time_a ;         /* time array				*/
static ARRAY          step_a ;         /* time array
                                        * */
struct stat folder_check;


static DOUBLE *velocity;
static DOUBLE *pressure;
static DOUBLE *displacement;
static DOUBLE *ale_displacement;

static RESULT_DATA global_fluid_result;
static RESULT_DATA global_ale_result;
static RESULT_DATA global_struct_result;

static INT dis_dim;
static INT* numnp_tot;
static INT INPT=2;


INT  *surf_in;                      /*needed for surface in case of existing design information*/
INT  *counter_array;
INT  **tmp_array;



/*----------------------------------------------------------------------*
 |  compare the integers - qsort routine                  a.lipka 5/01  |
 |                                                                      |
 |  the call for the sorter of an INT vector is then                    |
 |                                                                      |
 |  qsort((INT*) vector, lenght, sizeof(INT), cmp_int);                 |
 |                                                                      |
 *----------------------------------------------------------------------*/
INT cmp_int(const void *a, const void *b )
{
    return *(INT *)a - * (INT *)b;
}
/*----------------------------------------------------------------------*
 |  compare the doubles - qsort routine                   a.lipka 5/01  |
 |                                                                      |
 |  the call for the sorter of a DOUBLE vector is then                  |
 |                                                                      |
 |  qsort((DOUBLE*) vector, lenght, sizeof(DOUBLE), cmp_double);        |
 |                                                                      |
 *----------------------------------------------------------------------*/
DOUBLE cmp_double(const void *a, const void *b )
{
    return *(DOUBLE *)a - * (DOUBLE *)b;
}


void vis3caf(INT numff, INT numaf, INT numsf)
{

  INT   i,j ;
  FILE* test;
  char TITL[10]="Picture";
  char TKEYS[NKEYS][16];
  strcpy(TKEYS[0], "VELOCITY Ux     ");
  strcpy(TKEYS[1], "VELOCITY Uy     ");
  strcpy(TKEYS[2], "VELOCITY Uz     ");
  strcpy(TKEYS[3], "PRESSURE        ");
  strcpy(TKEYS[4], "FLOW VECTORS    ");
  strcpy(TKEYS[5], "ABS VEL         ");
  strcpy(TKEYS[6], "MOVIE CREATION  ");
  strcpy(TKEYS[7], "DISPLACEMENT Dx ");
  strcpy(TKEYS[8], "DISPLACEMENT Dy ");
  strcpy(TKEYS[9], "DISPLACEMENT Dz ");
  strcpy(TKEYS[10],"ABS DIS         ");
  strcpy(TKEYS[11],"TOGGLE MOV 3D/2D");

#ifdef DEBUG
  dstrc_enter("vis3caf");
#endif

  /* Setting up the right colourtable.
   *
   * This table should be in the home/.Visual3/ folder.
   * If not, the right table is created.
   * The variable COLOUR depends on the command line parameters
   *       0  :  black background
   *       1  :  white background
   *       2  :  grey colourscale */

  char path[strlen(getenv("HOME"))+25];  /*   +25 = strlen(/.Visual3/spec_black.col)+1  */
  strcpy(path, getenv("HOME"));
  strcat(path, "/.Visual3/");

  /*check if .Visual3 folder exists, if not --> mkdir*/
  if (stat(path, &folder_check)==-1)
  {
    if (mkdir(path, S_IRWXU)==0) printf("\n  Created directory %s  \n", path);
    else dserror("ERROR WHILE TRYING TO CREATE FOLDER : %s", path);
  }

  switch(COLOUR)
  {
    case 0:
      strcat(path, "spec_black.col");
      test=fopen(path,"r");
      if (test==0)
      {
        write_colourtable_black();
      }
      strcpb(CMFILE,path,strlen(path));
      break;
    case 1:
      strcat(path, "spec_white.col");
      test=fopen(path,"r");
      if (test==0)
      {
        write_colourtable_white();
      }
      strcpb(CMFILE,path,strlen(path));
      break;
    case 2:
      strcat(path, "spec_grey_.col");
      test=fopen(path,"r");
      if (test==0)
      {
        write_colourtable_grey_();
      }
      strcpb(CMFILE,path,strlen(path));
      break;
    default:
      break;
  }

  /* set some values -------------*/
  numnp=0;
  numnp2D=0;
  numnp3D=0;
  numele=0;
  KCEL1=0;
  KCEL2=0;
  KCEL3=0;
  KCEL4=0;

  if (numff!=-1) /*-------------------------------init the fluid variables*/
  {
    numele_fluid = discret[numff].field->numele;
    numnp_fluid =  numnp_tot[numff];
    numnp+=numnp_fluid;
    numele+=numele_fluid;
  }
  else
  {
    numnp_fluid=0;
    numele_fluid=0;
  }

  if (numsf!=-1)/*-------------------------------init the structure variables*/
  {
    numele_struct = discret[numsf].field->numele;
    numnp_struct = discret[numsf].field->numnp;
    numnp+=numnp_struct;
    numele+=numele_struct;
  }
  else
  {
    numnp_struct=0;
    numele_struct=0;
  }

  for (j=0;j<num_discr;j++) /*loop all discretizations*/
  {

    actele = &(discret[j].element[0]);
    distyp = actele->distyp;

    /*------------------------- check if all elements are the same distyp */
    for (i=1;i<discret[j].field->numele;i++) /*loop all elements*/
    {
      actele=&(discret[j].element[i]);
      if (actele->distyp!=distyp)
        dserror("up to now, all elements have to be the same distyp!\n");
    } /* end loop over all elements */

    /*here we set important variables for V3CELL, V3SURF and V3SURFACE
     *      - numnp3D is the total amount of 3D nodes
     *      - numnp2D is the total amount of 2D nodes, we have to
     *        double them later     *
     *      - KSURFELE is the number of nodes belonging to one face of
     *        this element type     *
     *      - FLUID_CELn_offset is the number of CELn spots this fluid
     *        is going to use */

    set_CEL_values(&KCEL1, &KCEL2, &KCEL3, &KCEL4, &FLUID_CEL3_offset, &FLUID_CEL4_offset, &discret[j]);

    set_VISUAL_values(&numnp3D, &numnp2D, &KSURF, &KSURFELE, &discret[j], numnp_tot[j]);

    if (distyp==hex20)
    {
      if (INPT==2)
      { INT screen, dummy;
      input1:
        printf("\n");
        printf("     Please give the mode in which you want to interpolate missing solutions:\n");
        printf("      0 : linear interpolation\n");
        printf("     [1]: serendipity interpolation\n");

        screen=getchar();
        switch(screen)
        {
          case 10: INPT=1; break;
          case 48: INPT=0; dummy=getchar(); break;
          case 49: INPT=1; dummy=getchar(); break;
          default:
            printf("\nTry again!\n");
            goto input1;
        }
      }
    }

  }/*end of loop over all discretizations*/

    /*--------------------------------------- get the data limits */
    find_data_limits(discret, num_discr,FLIMS, ACTDIM);

#if 0
    /*--------------------------------- ...for hierarchical nodes */
    i=1;
    for (j=0; j<num_discr; j++)
    { if (discret[j].element[0].distyp==h_hex20)
        data_limits_h_hex20(&discret[j], FLIMS, numnp_tot[j], &i);
    }
#endif

    /*-------------------------------------- needed for V3SCAL later*/
    minvx=FLIMS[0][0];
    minvy=FLIMS[1][0];
    minvz=FLIMS[2][0];
    minpre=FLIMS[3][0];
    minabsv=FLIMS[5][0];
    mindx=FLIMS[7][0];
    mindy=FLIMS[8][0];
    mindz=FLIMS[9][0];
    minabsd=FLIMS[10][0];

    /*--------------------------------------- check limits*/
    for (i=0;i<NKEYS;i++)
    {
      if (FLIMS[i][0]==FLIMS[i][1]) FLIMS[i][1] += ONE;
    }

    /*-------------------------------------- set real number of last step */
    LASTSTEP=step_a.a.iv[ncols-1];

    if (numaf!=-1)
    {
      if (IOPT==0)
      {
        IOPT=2;
        printf("\n  %dD problem: steady data structure, unsteady grid and variables\n\n",ACTDIM );
      }
      if (IOPT==1) printf("\n  %dD problem : steady data structure, STEADY grid and variables\n\n",ACTDIM );
    }
    else
    {
      IOPT=1;
      printf("\n  %dD problem: steady data structure and grid, unsteady variables\n\n", ACTDIM);
    }

    KEQUIV=0;
    KNPTET=0;
    KPTET=0;
    KNBLOCK=0;
    BLOCKS=0;
    WIN3D=1;
    MIRROR=0;
    FACX=1;
    FACY=1;
    FACZ=1;

    if (ACTDIM==2)
    {
      dnumnp=2*numnp2D;     /*we need the double amount of nodes to create a
                             * 3D block out of a 2D surface */
      FLUID_NODE_offset=2*numnp_fluid;

      if (struct_idx!=-1)
      {
        KNSURFT=2;
      }
      else
      {
        KNSURFT=1;
      }

      if (numff!=-1)
      {
        velocity=(DOUBLE*)CCACALLOC(2*numnp_fluid, sizeof(DOUBLE));
        pressure=(DOUBLE*)CCACALLOC(numnp_fluid, sizeof(DOUBLE));
      }

      if (numsf!=-1)
      {
        displacement=(DOUBLE*)CCACALLOC(2*numnp_struct, sizeof(DOUBLE));
        ale_displacement=(DOUBLE*)CCACALLOC(2*numnp_fluid, sizeof(DOUBLE));
      }
    }

    if (ACTDIM==3)
    {
      dnumnp=numnp3D+2*numnp2D;
      FLUID_NODE_offset=numnp_fluid;
      KNSURF=0;
      KNSURFT=0;
      ndvol_tot=0;
      ndsurf_tot=0;
      ndline_tot=0;
      ndnode_tot=0;

#ifdef DESIGN
      /*at the moment we use the design information only for fluid
       * design bodies */
      if (map_find_int(discret[fluid_idx].field->group, "ndvol", &ndvol) &&
          discret[fluid_idx].element[0].distyp != hex20)
        /* general check if there is design info */
      {
        ndvol_tot+=fluid_design.ndvol;
        ndsurf_tot+=fluid_design.ndsurf;
        ndline_tot+=fluid_design.ndline;
        ndnode_tot+=fluid_design.ndnode;

        /*preparing the surface arrays for V3SURFACE later*/
        counter_array=(INT*)CCACALLOC(ndsurf_tot, sizeof(INT));
        tmp_array=(INT**)CCACALLOC(ndsurf_tot, sizeof(INT*));
        surf_in=(INT*)CCACALLOC(ndsurf_tot, sizeof(INT));

        for (i=0;i<ndsurf_tot;i++)
        {
          tmp_array[i]=(INT*)CCACALLOC(fluid_design.nnode, sizeof(INT));
          counter_array[i]=0;
          surf_in[i]=-1;
        }

        /*find the coupling between the discret & the design information*/
        if (map_find_int(discret[fluid_idx].field->group, "ndvol", &ndvol))
        {
          post_design_coupling(&discret[fluid_idx], &fluid_design, tmp_array, counter_array, surf_in);
        }

        ndsurf_in_tot=0;
        for (i=0;i<ndsurf_tot;i++)
        {
          if (surf_in[i]==1||counter_array[i]==0 ) ndsurf_in_tot++;
        }
        /*
         * KNSURFT is the real number of surfaces given to VISUAL :
         * its the absolute number of surfaces - the surfaces inside
         * the design body  */
        KNSURFT=ndsurf_tot-ndsurf_in_tot;
      }
#endif

      if (IOPT==2)
      {
#ifdef D_FSI
        /*if we got a shell problem we need the node director information*/
        actele=&discret[struct_idx].element[0];
        switch(actele->eltyp)
        {
          case el_shell8:
            if (struct_idx==-1) dserror("Structure field needed");
            /*at the moment in case of 3D FSI only the structure
             * surfaces are given to VISUAL --> KNSURFT += 1*/
            KNSURFT+=1;
            node_director=(DOUBLE*)CCACALLOC(3*(numnp_tot[struct_idx]),sizeof(DOUBLE));
            for (i=0;i<3*(numnp_tot[struct_idx]);i++)
            {
              node_director[i]=0;
            }
            post_read_shell8_info(discret, node_director, struct_idx);
            break;

          default:
            dserror("eltyp in vis3caf no supported yet");
            break;
        }
#else
        dserror("D_FSI not compiled in!!!");
#endif
      }

      if (numff!=-1)
      {
        velocity=(DOUBLE*)CCACALLOC(3*numnp_fluid, sizeof(DOUBLE));
        pressure=(DOUBLE*)CCACALLOC(numnp_fluid, sizeof(DOUBLE));
      }
      if (numsf!=-1)
      {
        displacement=(DOUBLE*)CCACALLOC(dis_dim*numnp_struct, sizeof(DOUBLE));
        ale_displacement=(DOUBLE*)CCACALLOC(3*numnp_fluid, sizeof(DOUBLE));
      }
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
            &KNSURFT,
            10, strlen(CMFILE), 16);

}/*end of vis3caf */

void v3movie()
{
  char string[100];
  char convert[100];
  char remove[100];

#ifdef DEBUG
  dstrc_enter("v3movie");
#endif

/*--------------------------------------- add file counter to filenames */
#ifdef HPUX11
  printf("movie creation only possible under HPUX 10.20\n");
  goto end;
#endif

  sprintf(string ,"xwd -name %d-D -silent -out %s/vis%04d.xwd",movie_window, frame_folder, ACTSTEP);
  sprintf(convert ,"convert %s/vis%04d.xwd %s/vis%04d.png",frame_folder, ACTSTEP,frame_folder, ACTSTEP);
  sprintf(remove ,"rm %s/vis%04d.xwd",frame_folder, ACTSTEP);

  printf("write %s/vis%04d.png\n",frame_folder, ACTSTEP);

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
 * (This actually happens in post_design_coupling)
 *
 * Here we only transfer this information to Visual3
 *
 * In case of 3D FSI we only give the surfaces of the structure field
 * to VISUAL.
 *
 * We always first handle the fluid surfaces, as we assume that there
 * is ALWAYS a fluid field. If there is a struct field we attach its
 * surfaces after the fluid surfaces. For that case we need FLUID_SURF_offset.
 **/

void v3surface(INT nsurf[][2], INT surf[], INT scel[][4], CHAR tsurf[][20])
{
  INT i,p, j,  inel;

#ifdef DEBUG
  dstrc_enter("v3surface");
#endif

  if (ACTDIM==2)
  {
    if (IOPT==2)
    {
      char *name = "frontside fluid     "
                   "frontside struct    "
                   "Default             ";

      strncpy((char*)tsurf, name, 60);
    }

    /* ----------------------------------------------------FLUID SURFACES 2D*/
    for (i=0;i<numele_fluid;i++) /* loop all elements */
    {
      actele=&(discret[fluid_idx].element[i]);
      switch (actele->distyp)
      {
        case quad4: /* 4 node rectangle */
          for(j=0;j<4;j++)
          {
            dsassert(actele->node[j]->Id_loc < numnp_fluid, "Ups");
            scel[i][j] = actele->node[j]->Id_loc+1;
          }
          nsurf[0][0]=numele_fluid;
          nsurf[0][1]=1;
          FLUID_SURF_offset=numele_fluid;
          FLUID_GROUP_offset=1;
          break;

        case quad9:
          inel = 0;
          p=8*i;
          /*----------------------------------------- sub-element 1 */
          scel[p+inel][0] = actele->node[0]->Id_loc+1;
          scel[p+inel][1] = actele->node[4]->Id_loc+1;
          scel[p+inel][2] = actele->node[8]->Id_loc+1;
          scel[p+inel][3] = actele->node[7]->Id_loc+1;
          inel++;
          scel[p+inel][0] = actele->node[0]->Id_loc+1+numnp_fluid;
          scel[p+inel][1] = actele->node[4]->Id_loc+1+numnp_fluid;
          scel[p+inel][2] = actele->node[8]->Id_loc+1+numnp_fluid;
          scel[p+inel][3] = actele->node[7]->Id_loc+1+numnp_fluid;
          inel++;
          /*----------------------------------------- sub-element 2 */
          scel[p+inel][0] = actele->node[4]->Id_loc+1;
          scel[p+inel][1] = actele->node[1]->Id_loc+1;
          scel[p+inel][2] = actele->node[5]->Id_loc+1;
          scel[p+inel][3] = actele->node[8]->Id_loc+1;
          inel++;
          scel[p+inel][0] = actele->node[4]->Id_loc+1+numnp_fluid;
          scel[p+inel][1] = actele->node[1]->Id_loc+1+numnp_fluid;
          scel[p+inel][2] = actele->node[5]->Id_loc+1+numnp_fluid;
          scel[p+inel][3] = actele->node[8]->Id_loc+1+numnp_fluid;
          inel++;
          /*----------------------------------------- sub-element 3 */
          scel[p+inel][0] = actele->node[8]->Id_loc+1;
          scel[p+inel][1] = actele->node[5]->Id_loc+1;
          scel[p+inel][2] = actele->node[2]->Id_loc+1;
          scel[p+inel][3] = actele->node[6]->Id_loc+1;
          inel++;
          scel[p+inel][0] = actele->node[8]->Id_loc+1+numnp_fluid;
          scel[p+inel][1] = actele->node[5]->Id_loc+1+numnp_fluid;
          scel[p+inel][2] = actele->node[2]->Id_loc+1+numnp_fluid;
          scel[p+inel][3] = actele->node[6]->Id_loc+1+numnp_fluid;
          inel++;
          /*----------------------------------------- sub-element 4 */
          scel[p+inel][0] = actele->node[7]->Id_loc+1;
          scel[p+inel][1] = actele->node[8]->Id_loc+1;
          scel[p+inel][2] = actele->node[6]->Id_loc+1;
          scel[p+inel][3] = actele->node[3]->Id_loc+1;
          inel++;
          scel[p+inel][0] = actele->node[7]->Id_loc+1+numnp_fluid;
          scel[p+inel][1] = actele->node[8]->Id_loc+1+numnp_fluid;
          scel[p+inel][2] = actele->node[6]->Id_loc+1+numnp_fluid;
          scel[p+inel][3] = actele->node[3]->Id_loc+1+numnp_fluid;

          nsurf[0][0]=8*numele_fluid;
          nsurf[0][1]=1;
          FLUID_SURF_offset=8*numele_fluid;
          FLUID_GROUP_offset=1;
          break;

        default:
          dserror("distyp not implemented yet");
          break;
      }
    }
  } /* ----------------------------------------------END OF FLUID SURFACES 2D*/

  else if (ACTDIM==3)
  {
    INT  k, n, m;
    INT offset;
    k=0;
    n=0;
    m=0;
    offset=0;

#ifdef DESIGN
    /* ----------------------------------------------------FLUID SURFACES 3D*/
    if (map_find_int(fluid_field->group, "ndvol", &ndvol) && discret[fluid_idx].element[0].distyp != hex20)
    {
      for (i=0;i<ndsurf_tot;i++)
      {
        if (surf_in[i]==-1 && counter_array[i]!=0 ) /*check if the design surface is
                                                     * inside the body*/
        {
          for (j=0;j<counter_array[i];j++)
          {
            for (k=0;k<4;k++)
            {
              scel[j+offset][k]=tmp_array[i][j*4+k]+1;
            }
          }

          offset+=counter_array[i];
          nsurf[m][0]=offset;
          nsurf[m][1]=1;
          m++;
          printf("\nKSURF : %d", KSURF);
          printf("\n OFFSET : %d \n", offset);
          FLUID_SURF_offset=offset;
          FLUID_GROUP_offset=m;

        }
        CCAFREE(tmp_array[i]);
      }/*end of loop over all surfaces*/
      CCAFREE(counter_array);
      CCAFREE(surf_in);
      CCAFREE(tmp_array);
    }/*end of "design information exists"*/

    else /*design information does not exist*/
    {
      printf(" !!!!Theres no design information in the input file!!!!\n");
      printf("      Rendering has to be activated manually.\n\n");
    }/*end of "design information does not exist"*/
    /* -------------------------------------------------END OF FLUID SURFACES 3D*/
#endif
  }

  /* -------------------------------------------------STRUCTURE SURFACES 2D & 3D*/
  if (struct_idx!=-1)
  {
    char *name = "Fluid               "
                 "Structure           ";
    strncpy((char*)tsurf, name, 40);

    INT p=0;

    for (i=0;i<numele_struct;i++) /* loop all elements */
    {
      actele=&(discret[struct_idx].element[i]);
      switch (actele->distyp)
      {

        case quad4: /* 4 node rectangle */
          inel=0;
          p=FLUID_SURF_offset+2*i;
          scel[p+inel][0] = actele->node[0]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][1] = actele->node[1]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][2] = actele->node[2]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][3] = actele->node[3]->Id_loc+1+FLUID_NODE_offset;
          inel++;
          scel[p+inel][0] = actele->node[0]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][1] = actele->node[1]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][2] = actele->node[2]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][3] = actele->node[3]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          nsurf[FLUID_GROUP_offset][0]=FLUID_SURF_offset+2*numele_struct;
          nsurf[FLUID_GROUP_offset][1]=1;
          break;

        case quad8:
          inel=0;
          p=FLUID_SURF_offset+10*i;
          /*----------------------------------------- sub element 1*/
          scel[p+inel][0] = actele->node[4]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][1] = actele->node[0]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][2] = actele->node[7]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][3] = 0;
          inel++;
          scel[p+inel][0] = actele->node[4]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][1] = actele->node[0]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][2] = actele->node[7]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][3] = 0;
          inel++;
          /*----------------------------------------- sub element 2*/
          scel[p+inel][0] = actele->node[5]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][1] = actele->node[1]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][2] = actele->node[4]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][3] = 0;
          inel++;
          scel[p+inel][0] = actele->node[5]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][1] = actele->node[1]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][2] = actele->node[4]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][3] = 0;
          inel++;
          /*----------------------------------------- sub element 3*/
          scel[p+inel][0] = actele->node[6]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][1] = actele->node[2]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][2] = actele->node[5]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][3] = 0;
          inel++;
          scel[p+inel][0] = actele->node[6]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][1] = actele->node[2]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][2] = actele->node[5]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][3] = 0;
          inel++;
          /*----------------------------------------- sub element 4*/
          scel[p+inel][0] = actele->node[7]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][1] = actele->node[3]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][2] = actele->node[6]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][3] = 0;
          inel++;
          scel[p+inel][0] = actele->node[7]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][1] = actele->node[3]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][2] = actele->node[6]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][3] = 0;
          inel++;
          /*----------------------------------------- sub element 5*/
          scel[p+inel][0] = actele->node[4]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][1] = actele->node[5]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][2] = actele->node[6]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][3] = actele->node[7]->Id_loc+1+FLUID_NODE_offset;
          inel++;
          scel[p+inel][0] = actele->node[4]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][1] = actele->node[5]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][2] = actele->node[6]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][3] = actele->node[7]->Id_loc+1+FLUID_NODE_offset+numnp_struct;

          nsurf[FLUID_GROUP_offset][0]=FLUID_SURF_offset+10*numele_struct;
          nsurf[FLUID_GROUP_offset][1]=1;
          break;

        case quad9:
          inel = 0;
          p=FLUID_SURF_offset + 8*i;
          /*----------------------------------------- sub-element 1 */
          scel[p+inel][0] = actele->node[0]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][1] = actele->node[4]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][2] = actele->node[8]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][3] = actele->node[7]->Id_loc+1+FLUID_NODE_offset;
          inel++;
          scel[p+inel][0] = actele->node[0]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][1] = actele->node[4]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][2] = actele->node[8]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][3] = actele->node[7]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          inel++;
          /*----------------------------------------- sub-element 2 */
          scel[p+inel][0] = actele->node[4]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][1] = actele->node[1]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][2] = actele->node[5]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][3] = actele->node[8]->Id_loc+1+FLUID_NODE_offset;
          inel++;
          scel[p+inel][0] = actele->node[4]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][1] = actele->node[1]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][2] = actele->node[5]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][3] = actele->node[8]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          inel++;
          /*----------------------------------------- sub-element 3 */
          scel[p+inel][0] = actele->node[8]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][1] = actele->node[5]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][2] = actele->node[2]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][3] = actele->node[6]->Id_loc+1+FLUID_NODE_offset;
          inel++;
          scel[p+inel][0] = actele->node[8]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][1] = actele->node[5]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][2] = actele->node[2]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][3] = actele->node[6]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          inel++;
          /*----------------------------------------- sub-element 4 */
          scel[p+inel][0] = actele->node[7]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][1] = actele->node[8]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][2] = actele->node[6]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][3] = actele->node[3]->Id_loc+1+FLUID_NODE_offset;
          inel++;
          scel[p+inel][0] = actele->node[7]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][1] = actele->node[8]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][2] = actele->node[6]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][3] = actele->node[3]->Id_loc+1+FLUID_NODE_offset+numnp_struct;

          nsurf[FLUID_GROUP_offset][0]=FLUID_SURF_offset+8*numele_struct;
          nsurf[FLUID_GROUP_offset][1]=1;
          break;

        case tri3:
          inel = 0;
          p= FLUID_SURF_offset + 2*i;
          scel[p+inel][0] = actele->node[0]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][1] = actele->node[1]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][2] = actele->node[2]->Id_loc+1+FLUID_NODE_offset;
          scel[p+inel][3] = 0;
          inel++;
          scel[p+inel][0] = actele->node[0]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][1] = actele->node[1]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][2] = actele->node[2]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          scel[p+inel][3] = 0;
          nsurf[FLUID_GROUP_offset][0]=FLUID_SURF_offset+2*numele_struct;
          nsurf[FLUID_GROUP_offset][1]=1;
          break;

        default:
          break;
      }
    }
  }/*end of for*/
  /* ----------------------------------------------END OF STRUCTURE SURFACES 3D*/


#ifdef DEBUG
  dstrc_exit();
#endif

} /* end of V3SURFACE */

/*---------------------------------------------------------------*/
/*
 *  V3GRID
 *
 * General information :
 * We ALWAYS store the fluid nodes before the struct nodes!!!
 *
 * ACTDIM==2 : We create a 3D block so we can visualize 2D
 * data in Visual3. To do that we double the number of
 * nodes and create a second slice parallel to the 2D data.
 * Together those two slices form a 3D block with the length
 * of 1 in z-direction (original data z-coordinate 0, copy
 * z-coordinate 1)
 *
 * ACTDIM==3 : If IOPT==2 we got to use shell information to
 * create the grid. The director vector is read by
 * post_read_shell8_information for example, and stored in
 * node_director[].*/

void v3grid(float XYZ[][3])
{
  INT i, h;

#ifdef DEBUG
  dstrc_enter("v3grid");
#endif

  for (h=0;h<num_discr;h++)
  {
    actfieldtyp=discret[h].field->type;

    if (ACTDIM==2)
    {
      if (IOPT==2)
      {
#ifdef D_FSI
        switch(actfieldtyp)
        {
          case fluid:
          for (i=0;i<numnp_fluid;i++)
          {
            actnode=&(discret[fluid_idx].node[i]);

            if (fluid_ale_connect[i] == -1)
            {
              XYZ[i][0] = actnode->x[0]*FACX;
              XYZ[i][1] = actnode->x[1]*FACY;
              XYZ[i][2] = 0;
              XYZ[i+numnp_fluid][0] = actnode->x[0]*FACX;
              XYZ[i+numnp_fluid][1] = actnode->x[1]*FACY;
              XYZ[i+numnp_fluid][2] = 1*FACZ;
            }
            else
            {
              actanode = &(discret[ale_idx].node[fluid_ale_connect[i]]);

              XYZ[i][0] = actanode->x[0] + ale_displacement[2*actanode->Id_loc + 0];
              XYZ[i][1] = actanode->x[1] + ale_displacement[2*actanode->Id_loc + 1];
              XYZ[i][2] = 0;
              XYZ[i+numnp_fluid][0] = actanode->x[0] + ale_displacement[2*actanode->Id_loc + 0];
              XYZ[i+numnp_fluid][1] = actanode->x[1] + ale_displacement[2*actanode->Id_loc + 1];
              XYZ[i+numnp_fluid][2] = 1*FACZ;
            }
          }
          break;
          case structure:
          for (i=0;i<numnp_struct;i++)
          {
            actanode = &(discret[struct_idx].node[i]);
            XYZ[i+2*numnp_fluid][0] = actanode->x[0] + displacement[2*actanode->Id_loc + 0];
            XYZ[i+2*numnp_fluid][1] = actanode->x[1] + displacement[2*actanode->Id_loc + 1];
            XYZ[i+2*numnp_fluid][2] = 0;
            XYZ[i+2*numnp_fluid+numnp_struct][0] = actanode->x[0] + displacement[2*actanode->Id_loc + 0];
            XYZ[i+2*numnp_fluid+numnp_struct][1] = actanode->x[1] + displacement[2*actanode->Id_loc + 1];
            XYZ[i+2*numnp_fluid+numnp_struct][2] = 1*FACZ;
          }
          break;

          default:
            break;
        }
#else
        dserror("FSI-functions not compiled in!\n");
#endif
      }/*end of IOPT == 2 (ACTDIM==2)*/

      if (IOPT==1)
      {
        switch(actfieldtyp)
        {
          case fluid:
            for (i=0;i<numnp_fluid;i++)
            {
              actnode=&(discret[fluid_idx].node[i]);
              XYZ[i][0] = actnode->x[0]*FACX;
              XYZ[i][1] = actnode->x[1]*FACY;
              XYZ[i][2] = 0;

              XYZ[i+numnp_fluid][0] = actnode->x[0]*FACX;
              XYZ[i+numnp_fluid][1] = actnode->x[1]*FACY;
              XYZ[i+numnp_fluid][2] = 1*FACZ;
            }
            break;
          case structure:
            for (i=0;i<numnp_struct;i++)
            {

              actnode=&(discret[struct_idx].node[i]);
              XYZ[i+2*numnp_fluid][0] = actnode->x[0]*FACX;
              XYZ[i+2*numnp_fluid][1] = actnode->x[1]*FACY;
              XYZ[i+2*numnp_fluid][2] = 0;

              XYZ[i+2*numnp_fluid+numnp_struct][0] = actnode->x[0]*FACX;
              XYZ[i+2*numnp_fluid+numnp_struct][1] = actnode->x[1]*FACY;
              XYZ[i+2*numnp_fluid+numnp_struct][2] = 1*FACZ;
            }
            break;
          case ale:
            break;
          default:
            break;

        }
      }/*end of IOPT == 1 (ACTDIM==2)*/

    }/*end of ACTDIM == 2*/

    if (ACTDIM==3)
    {
      if (IOPT==2)
      {
        DOUBLE dx, dy, dz;
        DOUBLE dvx, dvy, dvz;
        DOUBLE ddvx, ddvy, ddvz;
#ifdef D_FSI
        switch(actfieldtyp)
        {
          case fluid:
            for (i=0;i<numnp_fluid;i++)
            {
              actnode=&(discret[fluid_idx].node[i]);

              if (fluid_ale_connect[i] == -1)
              {
                XYZ[i][0] = actnode->x[0]*FACX;
                XYZ[i][1] = actnode->x[1]*FACY;
                XYZ[i][2] = actnode->x[2]*FACZ;
              }
              else
              {
                actanode = &(discret[ale_idx].node[fluid_ale_connect[i]]);

                XYZ[i][0] = actanode->x[0] + ale_displacement[3*actanode->Id_loc + 0];
                XYZ[i][1] = actanode->x[1] + ale_displacement[3*actanode->Id_loc + 1];
                XYZ[i][2] = actanode->x[2] + ale_displacement[3*actanode->Id_loc + 2];

              }
            }
            break;
          case structure:
            for (i=0;i<discret[struct_idx].field->numnp;i++)
            {
              actnode=&(discret[struct_idx].node[i]);
              dvx=node_director[3*actnode->Id_loc];
              dvy=node_director[3*actnode->Id_loc+1];
              dvz=node_director[3*actnode->Id_loc+2];

              dx=displacement[dis_dim*i];
              dy=displacement[dis_dim*i+1];
              dz=displacement[dis_dim*i+2];
              ddvx=displacement[dis_dim*i+3];
              ddvy=displacement[dis_dim*i+4];
              ddvz=displacement[dis_dim*i+5];

              XYZ[numnp_fluid+i][0] = actnode->x[0]+dx+dvx+ddvx;
              XYZ[numnp_fluid+i][1] = actnode->x[1]+dy+dvy+ddvy;
              XYZ[numnp_fluid+i][2] = actnode->x[2]+dz+dvz+ddvz;
              XYZ[numnp_fluid+numnp_struct+i][0] = actnode->x[0]+dx-dvx-ddvx;
              XYZ[numnp_fluid+numnp_struct+i][1] = actnode->x[1]+dy-dvy-ddvx;
              XYZ[numnp_fluid+numnp_struct+i][2] = actnode->x[2]+dz-dvz-ddvx;
            }
            break;
          default:
            break;
        }
#endif
      } /*end of IOPT==2 (ACTDIM==3)*/

      if (IOPT==1)
      {
        switch(actfieldtyp)
        {
          case fluid:
            for (i=0;i<numnp_fluid;i++)
            {
              actnode=&(discret[fluid_idx].node[i]);
              XYZ[i][0] = actnode->x[0]*FACX;
              XYZ[i][1] = actnode->x[1]*FACY;
              XYZ[i][2] = actnode->x[2]*FACZ;
            }
            break;
          case structure:
            for (i=0;i<numnp_struct;i++)
            {
              actnode=&(discret[struct_idx].node[i]);
              XYZ[i+numnp_fluid][0] = actnode->x[0]*FACX;
              XYZ[i+numnp_fluid][1] = actnode->x[1]*FACY;
              XYZ[i+numnp_fluid][2] = actnode->x[2]*FACZ;
            }
            break;
          case ale:
            break;
          default:
            break;
        }
      }/*end of IOPT==1 (ACTDIM==3)*/

    }/*end of ACTDIM==3*/

  }/*end of loop over discretizations*/

#ifdef DEBUG
  dstrc_exit();
#endif
  return;

} /* end of V3GRID */


/*---------------------------------------------------------------------*/
/*
  V3SCAL

  This routine is called by VISUAL3 during the visualisation.
  It filters the scalar values out of the nodal sol-field at the actual
  position icol.
  icol is set every timestep in v3update

  If data is visualized a node cannot provide (e.g. dx at a fluid node) the
  value of this node is set to the minimum of the currently visualized
  data.
  ------------------------------------------------------------------------*/

void v3scal(INT *JKEY, float *S)
{
  INT i, k;
  float vx,vy, vz;
  float dx,dy, dz;

#ifdef DEBUG
  dstrc_enter("v3scal");
#endif

  for (k=0;k<num_discr;k++)  /*go through all discretizations*/
  {
    actfieldtyp=discret[k].field->type;
        /*-------------------------------------------------- get scalar data */
        if(ACTDIM==2)  /*2D scalar data*/
        {
          switch(*JKEY)
          {
            /*-------------------------------------------------------------------*/
            case 1: /* Ux */
              if (actfieldtyp==fluid)
              {
                for (i=0;i<discret[k].field->numnp;i++)
                {
                  S[i]=velocity[2*i];
                  S[i+numnp_fluid]=velocity[2*i];
                }
              }
              if (actfieldtyp==structure)
              {
                for (i=0;i<discret[k].field->numnp;i++)
                {
                  S[i+2*numnp_fluid]=minvx;
                  S[i+2*numnp_fluid+numnp_struct]=minvx;
                }
              }
              break;
              /*-------------------------------------------------------------------*/
            case 2: /* Uy */
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_fluid;i++)
                {
                  S[i]=velocity[2*i+1];
                  S[i+numnp_fluid]=velocity[2*i+1];
                }
              }
              if (actfieldtyp==structure)
              {
                for (i=0;i<numnp_struct;i++)
                {
                  S[i+2*numnp_fluid]=minvy;
                  S[i+2*numnp_fluid+numnp_struct]=minvy;
                }
              }
              toggle_movie_window=1;
              break;
              /*-------------------------------------------------------------------*/
            case 3: /* Uz */
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_fluid;i++)
                {
                  S[i]=0;
                  S[i+numnp_fluid]=0;
                }
              }
              if (actfieldtyp==structure)
              {
                for (i=0;i<numnp_struct;i++)
                {
                  S[i+2*numnp_fluid]=0;
                  S[i+2*numnp_fluid+numnp_struct]=0;
                }
              }
              toggle_movie_window=1;
              break;
              /*--------------------------------------------------------------------*/
            case 4: /* pressure */
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_fluid;i++)
                {
                  S[i]=pressure[i];
                  S[i+numnp_fluid]=pressure[i];
                }
              }
              if (actfieldtyp==structure)
              {
                for(i=0;i<numnp_struct;i++)
                {
                  S[i+2*numnp_fluid]=minpre;
                  S[i+2*numnp_fluid+numnp_struct]=minpre;
                }
              }
              toggle_movie_window=1;
              break;
              /*-------------------------------------------------------------------*/
            case 5:
              printf("VECTOR DATA IN SCALAR ROUTINE ???????\n");
              break;
              /*-------------------------------------------------------------------*/
            case 6: /* absolute velocity */
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_fluid;i++)
                {
                  vx=velocity[2*i];
                  vy=velocity[2*i+1];
                  vz=0;
                  S[i] = sqrt(vx*vx+vy*vy+vz*vz);
                  S[i+numnp_fluid] = sqrt(vx*vx+vy*vy+vz*vz);
                }
              }
              if (actfieldtyp==structure)
              {
                for (i=0;i<numnp_struct;i++)
                {
                  S[i+2*numnp_fluid]=minabsv;
                  S[i+2*numnp_fluid+numnp_struct]=minabsv;
                }
              }
              toggle_movie_window=1;
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
            case 8: /*no Dx for fluid*/
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_fluid;i++)
                {
                  S[i]=mindx;
                  S[i+numnp_fluid]=mindx;
                }
              }
              if (actfieldtyp==structure)
              {
                for (i=0;i<numnp_struct;i++)
                {
                  S[i+2*numnp_fluid]=displacement[2*i];
                  S[i+2*numnp_fluid+numnp_struct]=displacement[2*i];
                }
              }
              toggle_movie_window=1;
              break;
              /*-------------------------------------------------------------------*/
            case 9: /* no Dy for fluid*/
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_fluid;i++)
                {
                  S[i]=mindy;
                  S[i+numnp_fluid]=mindy;
                }
              }
              if (actfieldtyp==structure)
              {
                for (i=0;i<numnp_struct;i++)
                {
                  S[i+2*numnp_fluid]=displacement[2*i+1];
                  S[i+2*numnp_fluid+numnp_struct]=displacement[2*i+1];
                }
              }
              toggle_movie_window=1;
              break;
              /*-------------------------------------------------------------------*/
            case 10: /* no Dz for fluid*/
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_fluid;i++)
                {
                  S[i]=mindz;
                  S[i+numnp_fluid]=mindz;
                }
              }
              if (actfieldtyp==structure)
              {
                for (i=0;i<numnp_struct;i++)
                {
                  S[i+2*numnp_fluid]=0;
                  S[i+2*numnp_fluid+numnp_struct]=0;
                }
              }
              toggle_movie_window=1;
              break;
              /*-------------------------------------------------------------------*/
            case 11: /*no absolute displacement for fluid*/
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_fluid;i++)
                {
                  S[i]=minabsd;
                  S[i+numnp_fluid]=minabsd;
                }
              }
              if (actfieldtyp==structure)
              {
                for (i=0;i<numnp_struct;i++)
                {
                  dx=displacement[2*i];
                  dy=displacement[2*i+1];
                  dz=0;
                  S[i+2*numnp_fluid] = sqrt(dx*dx+dy*dy+dz*dz);
                  S[i+2*numnp_fluid+numnp_struct] = sqrt(dx*dx+dy*dy+dz*dz);
                }
              }
              toggle_movie_window=1;
              break;
            case 12: /*toggle the windows for movie creation*/
            {
              if (toggle_movie_window==1)
              {
                if (movie_window==3) movie_window=2;
                else movie_window=3;
                printf("  Now recording from %dD-window\n", movie_window);
                toggle_movie_window=0;
              }
            }
            break;
              /*-------------------------------------------------------------------*/
          } /* end switch(*JKEY) */
        } /* end of 2D part */


        if(ACTDIM==3)  /*3D scalar data*/
        {
          switch(*JKEY)
          {
            /*-------------------------------------------------------------------*/
            case 1: /* Ux */
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_tot[k];i++)
                {
                  S[i]=velocity[3*i];
                }
              }
              if (actfieldtyp==structure)
              {
                for (i=0;i<numnp_tot[k];i++)
                {
                  S[i+numnp_fluid]=minvx;
                  S[i+numnp_fluid+numnp_struct]=minvx;
                }
              }
              break;
              /*-------------------------------------------------------------------*/
            case 2: /* Uy */
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_tot[k];i++)
                {
                  S[i]=velocity[3*i+1];
                }
              }
              if (actfieldtyp==structure)
              {
                for (i=0;i<numnp_tot[k];i++)
                {
                  S[i+numnp_fluid]=minvy;
                  S[i+numnp_fluid+numnp_struct]=minvy;
                }
              }
              toggle_movie_window=1;
              break;
              /*-------------------------------------------------------------------*/
            case 3: /* Uz */
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_fluid;i++)
                  {
                  S[i]=velocity[3*i+2];
                }
              }
              if (actfieldtyp==structure)
              {
                for (i=0;i<numnp_struct;i++)
                {
                  S[i+numnp_fluid]=minvz;
                  S[i+numnp_fluid+numnp_struct]=minvz;
                }
              }
              toggle_movie_window=1;
              break;
              /*--------------------------------------------------------------------*/
            case 4: /* pressure */
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_fluid;i++)
                {
                  S[i]=pressure[i];
                }
              }
              if (actfieldtyp==structure)
              {
                for(i=0;i<numnp_struct;i++)
                {
                  S[i+numnp_fluid]=minpre;
                  S[i+numnp_fluid+numnp_struct]=minpre;
                }
              }
              toggle_movie_window=1;
              break;
              /*-------------------------------------------------------------------*/
            case 5:
              printf("VECTOR DATA IN SCALAR ROUTINE ???????\n");
              break;
              /*-------------------------------------------------------------------*/
            case 6: /* absolute velocity */
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_fluid;i++)
                {
                  vx=velocity[3*i];
                  vy=velocity[3*i+1];
                  vz=velocity[3*i+2];
                  S[i] = sqrt(vx*vx+vy*vy+vz*vz);
                }
              }

              if (actfieldtyp==structure)
              {
                for (i=0;i<numnp_struct;i++)
                {
                  S[i+numnp_fluid]=minabsv;
                  S[i+numnp_fluid+numnp_struct]=minabsv;
                }
              }
              toggle_movie_window=1;
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
            case 8: /*no Dx for fluid*/
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_fluid;i++)
                {
                  S[i]=mindx;
                }
              }
              if (actfieldtyp==structure)
              {
                for (i=0;i<numnp_struct;i++)
                {
                  S[i+numnp_fluid]=displacement[dis_dim*i];
                  S[i+numnp_fluid+numnp_struct]=displacement[dis_dim*i];
                }
              }
              toggle_movie_window=1;
              break;
              /*-------------------------------------------------------------------*/
            case 9: /* no Dy for fluid*/
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_fluid;i++)
                {
                  S[i]=mindy;
                }
              }

              if (actfieldtyp==structure)
              {
                for (i=0;i<numnp_struct;i++)
                {
                  S[i+numnp_fluid]=displacement[dis_dim*i+1];
                  S[i+numnp_fluid+numnp_struct]=displacement[dis_dim*i+1];
                }
              }
              toggle_movie_window=1;
              break;
              /*-------------------------------------------------------------------*/
            case 10: /* no Dz for fluid*/
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_fluid;i++)
                {
                  S[i]=mindz;
                }
              }
              if (actfieldtyp==structure)
              {
                for (i=0;i<numnp_struct;i++)
                {
                  S[i+numnp_fluid]=displacement[dis_dim*i+2];
                  S[i+numnp_fluid+numnp_struct]=displacement[dis_dim*i+2];
                }
              }
              toggle_movie_window=1;
              break;
              /*-------------------------------------------------------------------*/
            case 11: /*no absolute displacement for fluid*/
              if (actfieldtyp==fluid)
              {
                for (i=0;i<numnp_fluid;i++)
                {
                  S[i]=minabsd;
                }
              }
              if (actfieldtyp==structure)
              {
                for (i=0;i<numnp_struct;i++)
                {
                  dx=displacement[dis_dim*i];
                  dy=displacement[dis_dim*i+1];
                  dz=displacement[dis_dim*i+2];
                  S[i+numnp_fluid] = sqrt(dx*dx+dy*dy+dz*dz);
                  S[i+numnp_fluid+numnp_struct] = sqrt(dx*dx+dy*dy+dz*dz);
                }
              }
              toggle_movie_window=1;
              break;
            case 12: /*toggle the windows for movie creation*/
            {
              if (toggle_movie_window==1)
              {
                if (movie_window==3) movie_window=2;
                else movie_window=3;
                printf("  Now recording from %dD-window\n", movie_window);
                toggle_movie_window=0;
              }
            }
            break;
              /*-------------------------------------------------------------------*/
          } /* end switch(*JKEY) */
        } /* end of 3D part */

  }/*end of loop over all discret[k]*/

return;
#ifdef DEBUG
dstrc_exit();
#endif
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

  INT i, h;

#ifdef DEBUG
  dstrc_enter("V3VECT");
#endif

  for (h=0;h<num_discr;h++)
  {
    actfieldtyp=discret[h].field->type;
    switch(actfieldtyp)
    {
      case fluid:
        /*-------------------------------------------------- get vector data */
        if (*JKEY==5)  /* streamlines */
        {
          if (ACTDIM==2)
          {
            for (i=0;i<numnp_fluid;i++)
            {
              V[i][0] = velocity[2*i];
              V[i][1] = velocity[2*i+1];
              V[i+numnp_fluid][0] = velocity[2*i];
              V[i+numnp_fluid][1] = velocity[2*i+1];

              V[i][2] = 0;
              V[i+numnp_fluid][2] = 0;
            }
          }

          if (ACTDIM==3)
          {
            for (i=0;i<numnp_fluid;i++)
            {
              V[i][0] = velocity[3*i];
              V[i][1] = velocity[3*i+1];
              V[i][2] = velocity[3*i+2];
            }
          }
        }
        else
          printf("SCALAR DATA IN VECTOR ROUTINE ???????\n");

        break;
      case structure:
        /*-------------------------------------------------- get vector data */
        break;
      case ale:
        break;

      default : dserror("error");
    } /* end switch(actfieldtyp) */
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of V3VECT */

void v3update(float *TIME)
{
  INT i;
#ifdef DEBUG
  dstrc_enter("V3UPDATE");
#endif

/*------------------------------------------- check for movie creation */
  if (icol==0 && IMOVIE==2)
  {
    IMOVIE=0;
    printf("\n");
    printf("Movie creation finished\n");
    printf("For assembling the frames to a movie enter : program %s\n", frame_folder);
  }

  if (icol==0 && IMOVIE==1)
  {
    if (stat("movie", &folder_check)==-1)
    {
      if (mkdir("movie",S_IRWXU)==-1) dserror("ERROR WHILE CREATING FOLDER /MOVIE");
      else printf("   Created directory /movie/ \n");
    }
    for (i=0;i<999;i++)
    {
      sprintf(frame_folder, "movie/frames%04d", i);
      if (stat(frame_folder, &folder_check)==-1) break;
      if (i==999) dserror("MORE THAN 999 FOLDERS");
    }
    if (mkdir(frame_folder, S_IRWXU)==-1) dserror("ERROR WHILE CREATING FOLDER : %s", frame_folder);
    else
    {
    printf("   Created directory /%s/\n", frame_folder);
    IMOVIE=2;
    }
  }
  if (IMOVIE==2) v3movie();

/*------------------------------------------ go to the next step*/
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


  CHUNK_DATA chunk;
  INT j, h;

  if (fluid_field!=NULL)
  {
    if (!next_result(&global_fluid_result))
    {
      destroy_result_data(&global_fluid_result);
      init_result_data(fluid_field, &global_fluid_result);
      if (!next_result(&global_fluid_result))
      {
        dserror("failed to reinitialize fluid");
      }
    }
  }
  if (ale_field!=NULL)
  {
    if (!next_result(&global_ale_result))
    {
      destroy_result_data(&global_ale_result);
      init_result_data(ale_field, &global_ale_result);
      if (!next_result(&global_ale_result))
      {
        dserror("failed to reinitialize ale");
      }
    }
  }
  if (struct_field!=NULL)
  {
    if (!next_result(&global_struct_result))
    {
      destroy_result_data(&global_struct_result);
      init_result_data(struct_field, &global_struct_result);
      if (!next_result(&global_struct_result))
      {
        dserror("failed to reinitialize struct");
      }
    }
  }

  /*-----------------------------------------------------*/
  for (h=0;h<num_discr;h++)
  {
    actfieldtyp=discret[h].field->type;
    switch(actfieldtyp)
    {
      case fluid:
        /*read the velocity*/
        init_chunk_data(&global_fluid_result, &chunk, "velocity");
        if (ACTDIM==2)
          dsassert(chunk.value_entry_length==2, "2d problem expected");
        if (ACTDIM==3)
          dsassert(chunk.value_entry_length==3, "3d problem expected");

        for (i=0; i<discret[h].field->numnp; ++i)
        {
          chunk_read_value_entry(&chunk, i);
          if (ACTDIM==2)
          {
            velocity[2*i+0] = chunk.value_buf[0];
            velocity[2*i+1] = chunk.value_buf[1];
          }
          if (ACTDIM==3)
          {
            velocity[3*i+0] = chunk.value_buf[0];
            velocity[3*i+1] = chunk.value_buf[1];
            velocity[3*i+2] = chunk.value_buf[2];
          }
        }
        destroy_chunk_data(&chunk);

        /*read the pressure*/
        if (map_has_map(global_fluid_result.group, "pressure"))
        {
          init_chunk_data(&global_fluid_result, &chunk, "pressure");
        }
        else if (map_has_map(global_fluid_result.group, "average_pressure"))
        {
          init_chunk_data(&global_fluid_result, &chunk, "average_pressure");
        }
        else
        {
          dserror("No pressure entry found.");
        }
        dsassert(chunk.value_entry_length==1, "there must be just one pressure value");

        for (i=0; i<discret[h].field->numnp; ++i)
        {
          chunk_read_value_entry(&chunk, i);
          pressure[i] = chunk.value_buf[0];
        }
        destroy_chunk_data(&chunk);
        if (discret[h].element[0].distyp==hex20)
          lin_interpol(&discret[h], numnp_tot[h], velocity, pressure,  INPT);
#if 0
        if (discret[h].element[0].distyp==h_hex20)
          hier_elements(&discret[h], numnp_tot[h], velocity, pressure);
#endif
        break;

      case structure:
        init_chunk_data(&global_struct_result, &chunk, "displacement");

        for (i=0; i<discret[h].field->numnp; ++i)
        {
          chunk_read_value_entry(&chunk, i);
          if (ACTDIM==2)
          {
            displacement[2*i+0] = chunk.value_buf[0];
            displacement[2*i+1] = chunk.value_buf[1];
          }
          if (ACTDIM==3)
          {
            for (j=0;j<dis_dim;j++)
            {
              displacement[dis_dim*i+j] = chunk.value_buf[j];
            }
          }
        }
        destroy_chunk_data(&chunk);
        break;

      case ale:
        init_chunk_data(&global_ale_result, &chunk, "displacement");
        if (ACTDIM==2)
          dsassert(chunk.value_entry_length==2, "2d problem expected");
        if (ACTDIM==3)
          dsassert(chunk.value_entry_length==3, "3d problem expected");

        for (i=0; i<discret[h].field->numnp; ++i)
        {
          chunk_read_value_entry(&chunk, i);
          if (ACTDIM==2)
          {
            ale_displacement[2*i+0] = chunk.value_buf[0];
            ale_displacement[2*i+1] = chunk.value_buf[1];
          }
          if (ACTDIM==3)
          {
            ale_displacement[3*i+0] = chunk.value_buf[0];
            ale_displacement[3*i+1] = chunk.value_buf[1];
            ale_displacement[3*i+2] = chunk.value_buf[2];
          }
        }
        if (discret[h].element[0].distyp==hex20)
        {
          lin_interpol(&discret[h], numnp_tot[h], ale_displacement, NULL,  INPT);

        }
#if 0
        if (discret[h].element[0].distyp==h_hex20)
          hier_elements(&discret[h], numnp_tot[h], ale_displacement, NULL);
#endif
        break;
        destroy_chunk_data(&chunk);
        break;

      default:
        break;
    }/*end of switch(actfieldtyp)*/
  }/*end of loop over discretizations*/

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



/*-----------------------------------------------------------------------
 * V3CELL
 *
 * Here we cell structure of the problem is set. Like in V3SURFACE and
 * in V3GRID we first handle the fluid and then the structure.
 *
 * FLUID_NODE_offset is an important value here. It defines how many
 * fluid nodes are in XYZ[] before the struct nodes.
 *
 * ---------------------------------------------------------------------*/
void v3cell(INT cel1[][4], INT cel2[][5], INT cel3[][6], INT cel4[][8],
            INT nptet[][8], INT *ptet)
{
  INT i,j, p, inel;

#ifdef DEBUG
  dstrc_enter("v3cell");
#endif

  /*-----------------------------------------------------------FLUID CELLS */
  for (i=0;i<numele_fluid;i++) /* loop all elements */
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

      case tet4:
        for (j=0;j<4;j++)
        {
          actnode = actele->node[j];
          cel1[i][j]=actnode->Id_loc+1;
        }
        break;

      case hex20:
      /* case h_hex20: */
      case hex27:
        /* every element is divided into 8 hex8 sub-elements */
        inel = 0;
        p=8*i;
        /*----------------------------------------- sub-element 1 */
        cel4[p+inel][0] = actele->node[ 0]->Id_loc+1;
        cel4[p+inel][1] = actele->node[ 8]->Id_loc+1;
        cel4[p+inel][2] = actele->node[20]->Id_loc+1;
        cel4[p+inel][3] = actele->node[11]->Id_loc+1;
        cel4[p+inel][4] = actele->node[12]->Id_loc+1;
        cel4[p+inel][5] = actele->node[21]->Id_loc+1;
        cel4[p+inel][6] = actele->node[26]->Id_loc+1;
        cel4[p+inel][7] = actele->node[24]->Id_loc+1;
        inel++;

        /*----------------------------------------- sub-element 2 */
        cel4[p+inel][0] = actele->node[ 8]->Id_loc+1;
        cel4[p+inel][1] = actele->node[ 1]->Id_loc+1;
        cel4[p+inel][2] = actele->node[ 9]->Id_loc+1;
        cel4[p+inel][3] = actele->node[20]->Id_loc+1;
        cel4[p+inel][4] = actele->node[21]->Id_loc+1;
        cel4[p+inel][5] = actele->node[13]->Id_loc+1;
        cel4[p+inel][6] = actele->node[22]->Id_loc+1;
        cel4[p+inel][7] = actele->node[26]->Id_loc+1;
        inel++;

        /*----------------------------------------- sub-element 3 */
        cel4[p+inel][0] = actele->node[20]->Id_loc+1;
        cel4[p+inel][1] = actele->node[ 9]->Id_loc+1;
        cel4[p+inel][2] = actele->node[ 2]->Id_loc+1;
        cel4[p+inel][3] = actele->node[10]->Id_loc+1;
        cel4[p+inel][4] = actele->node[26]->Id_loc+1;
        cel4[p+inel][5] = actele->node[22]->Id_loc+1;
        cel4[p+inel][6] = actele->node[14]->Id_loc+1;
        cel4[p+inel][7] = actele->node[23]->Id_loc+1;
        inel++;

        /*----------------------------------------- sub-element 4 */
        cel4[p+inel][0] = actele->node[11]->Id_loc+1;
        cel4[p+inel][1] = actele->node[20]->Id_loc+1;
        cel4[p+inel][2] = actele->node[10]->Id_loc+1;
        cel4[p+inel][3] = actele->node[ 3]->Id_loc+1;
        cel4[p+inel][4] = actele->node[24]->Id_loc+1;
        cel4[p+inel][5] = actele->node[26]->Id_loc+1;
        cel4[p+inel][6] = actele->node[23]->Id_loc+1;
        cel4[p+inel][7] = actele->node[15]->Id_loc+1;
        inel++;

        /*----------------------------------------- sub-element 5 */
        cel4[p+inel][0] = actele->node[12]->Id_loc+1;
        cel4[p+inel][1] = actele->node[21]->Id_loc+1;
        cel4[p+inel][2] = actele->node[26]->Id_loc+1;
        cel4[p+inel][3] = actele->node[24]->Id_loc+1;
        cel4[p+inel][4] = actele->node[ 4]->Id_loc+1;
        cel4[p+inel][5] = actele->node[16]->Id_loc+1;
        cel4[p+inel][6] = actele->node[25]->Id_loc+1;
        cel4[p+inel][7] = actele->node[19]->Id_loc+1;
        inel++;

        /*----------------------------------------- sub-element 6 */
        cel4[p+inel][0] = actele->node[21]->Id_loc+1;
        cel4[p+inel][1] = actele->node[13]->Id_loc+1;
        cel4[p+inel][2] = actele->node[22]->Id_loc+1;
        cel4[p+inel][3] = actele->node[26]->Id_loc+1;
        cel4[p+inel][4] = actele->node[16]->Id_loc+1;
        cel4[p+inel][5] = actele->node[ 5]->Id_loc+1;
        cel4[p+inel][6] = actele->node[17]->Id_loc+1;
        cel4[p+inel][7] = actele->node[25]->Id_loc+1;
        inel++;

        /*----------------------------------------- sub-element 7 */
        cel4[p+inel][0] = actele->node[26]->Id_loc+1;
        cel4[p+inel][1] = actele->node[22]->Id_loc+1;
        cel4[p+inel][2] = actele->node[14]->Id_loc+1;
        cel4[p+inel][3] = actele->node[23]->Id_loc+1;
        cel4[p+inel][4] = actele->node[25]->Id_loc+1;
        cel4[p+inel][5] = actele->node[17]->Id_loc+1;
        cel4[p+inel][6] = actele->node[ 6]->Id_loc+1;
        cel4[p+inel][7] = actele->node[18]->Id_loc+1;
        inel++;

        /*----------------------------------------- sub-element 8 */
        cel4[p+inel][0] = actele->node[24]->Id_loc+1;
        cel4[p+inel][1] = actele->node[26]->Id_loc+1;
        cel4[p+inel][2] = actele->node[23]->Id_loc+1;
        cel4[p+inel][3] = actele->node[15]->Id_loc+1;
        cel4[p+inel][4] = actele->node[19]->Id_loc+1;
        cel4[p+inel][5] = actele->node[25]->Id_loc+1;
        cel4[p+inel][6] = actele->node[18]->Id_loc+1;
        cel4[p+inel][7] = actele->node[ 7]->Id_loc+1;
        inel++;
        break;

      case quad4: /* 4 node rectangle */
        for(j=0;j<4;j++)
        {
          dsassert(actele->node[j]->Id_loc < numnp_fluid, "Ups");
          cel4[i][j] = actele->node[j]->Id_loc+1;
          cel4[i][j+4] = actele->node[j]->Id_loc+1+numnp_fluid;
        }
        break;

      case quad8:
        inel=0;
        p=4*i;
        /*----------------------------------------- sub element 1*/
        cel3[p=inel][0] = actele->node[4]->Id_loc+1;
        cel3[p+inel][1] = actele->node[4]->Id_loc+1+numnp_fluid;
        cel3[p+inel][2] = actele->node[0]->Id_loc+1+numnp_fluid;
        cel3[p+inel][3] = actele->node[0]->Id_loc+1;
        cel3[p+inel][4] = actele->node[7]->Id_loc+1+numnp_fluid;
        cel3[p+inel][5] = actele->node[7]->Id_loc+1;
        inel++;
        /*----------------------------------------- sub element 2*/
        cel3[p+inel][0] = actele->node[5]->Id_loc+1;
        cel3[p+inel][1] = actele->node[5]->Id_loc+1+numnp_fluid;
        cel3[p+inel][2] = actele->node[1]->Id_loc+1+numnp_fluid;
        cel3[p+inel][3] = actele->node[1]->Id_loc+1;
        cel3[p+inel][4] = actele->node[4]->Id_loc+1+numnp_fluid;
        cel3[p+inel][5] = actele->node[4]->Id_loc+1;
        inel++;
        /*----------------------------------------- sub element 3*/
        cel3[p+inel][0] = actele->node[6]->Id_loc+1;
        cel3[p+inel][1] = actele->node[6]->Id_loc+1+numnp_fluid;
        cel3[p+inel][2] = actele->node[2]->Id_loc+1+numnp_fluid;
        cel3[p+inel][3] = actele->node[2]->Id_loc+1;
        cel3[p+inel][4] = actele->node[5]->Id_loc+1+numnp_fluid;
        cel3[p+inel][5] = actele->node[5]->Id_loc+1;
        inel++;
        /*----------------------------------------- sub element 4*/
        cel3[p+inel][0] = actele->node[7]->Id_loc+1;
        cel3[p+inel][1] = actele->node[7]->Id_loc+1+numnp_fluid;
        cel3[p+inel][2] = actele->node[3]->Id_loc+1+numnp_fluid;
        cel3[p+inel][3] = actele->node[3]->Id_loc+1;
        cel3[p+inel][4] = actele->node[6]->Id_loc+1+numnp_fluid;
        cel3[p+inel][5] = actele->node[6]->Id_loc+1;
        inel++;
        /*----------------------------------------- sub element 5*/
        p=i;
        cel4[p+i][0] = actele->node[4]->Id_loc+1;
        cel4[p+i][1] = actele->node[5]->Id_loc+1;
        cel4[p+i][2] = actele->node[6]->Id_loc+1;
        cel4[p+i][3] = actele->node[7]->Id_loc+1;
        cel4[p+i][4] = actele->node[4]->Id_loc+1+numnp_fluid;
        cel4[p+i][5] = actele->node[5]->Id_loc+1+numnp_fluid;
        cel4[p+i][6] = actele->node[6]->Id_loc+1+numnp_fluid;
        cel4[p+i][7] = actele->node[7]->Id_loc+1+numnp_fluid;
        break;

      case quad9:
        inel = 0;
        p=4*i;
        /*----------------------------------------- sub-element 1 */
        cel4[p+inel][0] = actele->node[0]->Id_loc+1;
        cel4[p+inel][1] = actele->node[4]->Id_loc+1;
        cel4[p+inel][2] = actele->node[8]->Id_loc+1;
        cel4[p+inel][3] = actele->node[7]->Id_loc+1;
        cel4[p+inel][4] = actele->node[0]->Id_loc+1+numnp_fluid;
        cel4[p+inel][5] = actele->node[4]->Id_loc+1+numnp_fluid;
        cel4[p+inel][6] = actele->node[8]->Id_loc+1+numnp_fluid;
        cel4[p+inel][7] = actele->node[7]->Id_loc+1+numnp_fluid;
        inel++;
        /*----------------------------------------- sub-element 2 */
        cel4[p+inel][0] = actele->node[4]->Id_loc+1;
        cel4[p+inel][1] = actele->node[1]->Id_loc+1;
        cel4[p+inel][2] = actele->node[5]->Id_loc+1;
        cel4[p+inel][3] = actele->node[8]->Id_loc+1;
        cel4[p+inel][4] = actele->node[4]->Id_loc+1+numnp_fluid;
        cel4[p+inel][5] = actele->node[1]->Id_loc+1+numnp_fluid;
        cel4[p+inel][6] = actele->node[5]->Id_loc+1+numnp_fluid;
        cel4[p+inel][7] = actele->node[8]->Id_loc+1+numnp_fluid;
        inel++;
        /*----------------------------------------- sub-element 3 */
        cel4[p+inel][0] = actele->node[8]->Id_loc+1;
        cel4[p+inel][1] = actele->node[5]->Id_loc+1;
        cel4[p+inel][2] = actele->node[2]->Id_loc+1;
        cel4[p+inel][3] = actele->node[6]->Id_loc+1;
        cel4[p+inel][4] = actele->node[8]->Id_loc+1+numnp_fluid;
        cel4[p+inel][5] = actele->node[5]->Id_loc+1+numnp_fluid;
        cel4[p+inel][6] = actele->node[2]->Id_loc+1+numnp_fluid;
        cel4[p+inel][7] = actele->node[6]->Id_loc+1+numnp_fluid;
        inel++;
        /*----------------------------------------- sub-element 4 */
        cel4[p+inel][0] = actele->node[7]->Id_loc+1;
        cel4[p+inel][1] = actele->node[8]->Id_loc+1;
        cel4[p+inel][2] = actele->node[6]->Id_loc+1;
        cel4[p+inel][3] = actele->node[3]->Id_loc+1;
        cel4[p+inel][4] = actele->node[7]->Id_loc+1+numnp_fluid;
        cel4[p+inel][5] = actele->node[8]->Id_loc+1+numnp_fluid;
        cel4[p+inel][6] = actele->node[6]->Id_loc+1+numnp_fluid;
        cel4[p+inel][7] = actele->node[3]->Id_loc+1+numnp_fluid;
        break;


      default:
        dserror("fluid distyp not implemented yet : %d !\n", distyp);
    }
  } /* end loop over all elements */
  /* -----------------------------------------------END OF FLUID CELLS*/


  /* ----------------------------------------------- STRUCT CELLS*/
  if (struct_idx!=-1)
  {
    for (i=0;i<numele_struct;i++) /* loop all elements */
    {
      actele=&(discret[struct_idx].element[i]);
      switch (actele->distyp)
      {
        case hex8:
          p=FLUID_CEL4_offset;
          for (j=0;j<8;j++)
          {
            actnode = actele->node[j];
            cel4[i+p][j]=actnode->Id_loc+1;
          }
          break;

        case hex20:
        /* case h_hex20: */
        case hex27:
          /* every element is divided into 8 hex8 sub-elements */
          inel = 0;
          p=FLUID_CEL4_offset+8*i;
          /*----------------------------------------- sub-element 1 */
          cel4[p+inel][0] = actele->node[ 0]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][1] = actele->node[ 8]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][2] = actele->node[20]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][3] = actele->node[11]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][4] = actele->node[12]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][5] = actele->node[21]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][6] = actele->node[26]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][7] = actele->node[24]->Id_loc+1+FLUID_NODE_offset;
          inel++;

          /*----------------------------------------- sub-element 2 */
          cel4[p+inel][0] = actele->node[ 8]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][1] = actele->node[ 1]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][2] = actele->node[ 9]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][3] = actele->node[20]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][4] = actele->node[21]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][5] = actele->node[13]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][6] = actele->node[22]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][7] = actele->node[26]->Id_loc+1+FLUID_NODE_offset;
          inel++;

          /*----------------------------------------- sub-element 3 */
          cel4[p+inel][0] = actele->node[20]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][1] = actele->node[ 9]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][2] = actele->node[ 2]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][3] = actele->node[10]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][4] = actele->node[26]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][5] = actele->node[22]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][6] = actele->node[14]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][7] = actele->node[23]->Id_loc+1+FLUID_NODE_offset;
          inel++;

          /*----------------------------------------- sub-element 4 */
          cel4[p+inel][0] = actele->node[11]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][1] = actele->node[20]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][2] = actele->node[10]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][3] = actele->node[ 3]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][4] = actele->node[24]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][5] = actele->node[26]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][6] = actele->node[23]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][7] = actele->node[15]->Id_loc+1+FLUID_NODE_offset;
          inel++;

          /*----------------------------------------- sub-element 5 */
          cel4[p+inel][0] = actele->node[12]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][1] = actele->node[21]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][2] = actele->node[26]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][3] = actele->node[24]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][4] = actele->node[ 4]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][5] = actele->node[16]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][6] = actele->node[25]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][7] = actele->node[19]->Id_loc+1+FLUID_NODE_offset;
          inel++;

          /*----------------------------------------- sub-element 6 */
          cel4[p+inel][0] = actele->node[21]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][1] = actele->node[13]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][2] = actele->node[22]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][3] = actele->node[26]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][4] = actele->node[16]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][5] = actele->node[ 5]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][6] = actele->node[17]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][7] = actele->node[25]->Id_loc+1+FLUID_NODE_offset;
          inel++;

          /*----------------------------------------- sub-element 7 */
          cel4[p+inel][0] = actele->node[26]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][1] = actele->node[22]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][2] = actele->node[14]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][3] = actele->node[23]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][4] = actele->node[25]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][5] = actele->node[17]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][6] = actele->node[ 6]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][7] = actele->node[18]->Id_loc+1+FLUID_NODE_offset;
          inel++;

          /*----------------------------------------- sub-element 8 */
          cel4[p+inel][0] = actele->node[24]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][1] = actele->node[26]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][2] = actele->node[23]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][3] = actele->node[15]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][4] = actele->node[19]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][5] = actele->node[25]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][6] = actele->node[18]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][7] = actele->node[ 7]->Id_loc+1+FLUID_NODE_offset;
          inel++;
          break;

        case quad4: /* 4 node rectangle */
          p=FLUID_CEL4_offset;

          for(j=0;j<4;j++)
          {
            dsassert(actele->node[j]->Id_loc < numnp_struct, "Ups");
            cel4[i+p][j] = actele->node[j]->Id_loc+1+FLUID_NODE_offset;
            cel4[i+p][j+4] = actele->node[j]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          }
          break;

        case quad8:
          inel=0;
          p=FLUID_CEL3_offset + 4*i;
          /*----------------------------------------- sub element 1*/
          cel3[p+inel][0] = actele->node[4]->Id_loc+1+FLUID_NODE_offset;
          cel3[p+inel][1] = actele->node[4]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel3[p+inel][2] = actele->node[0]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel3[p+inel][3] = actele->node[0]->Id_loc+1+FLUID_NODE_offset;
          cel3[p+inel][4] = actele->node[7]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel3[p+inel][5] = actele->node[7]->Id_loc+1+FLUID_NODE_offset;
          inel++;
          /*----------------------------------------- sub element 2*/
          cel3[p+inel][0] = actele->node[5]->Id_loc+1+FLUID_NODE_offset;
          cel3[p+inel][1] = actele->node[5]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel3[p+inel][2] = actele->node[1]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel3[p+inel][3] = actele->node[1]->Id_loc+1+FLUID_NODE_offset;
          cel3[p+inel][4] = actele->node[4]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel3[p+inel][5] = actele->node[4]->Id_loc+1+FLUID_NODE_offset;
          inel++;
          /*----------------------------------------- sub element 3*/
          cel3[p+inel][0] = actele->node[6]->Id_loc+1+FLUID_NODE_offset;
          cel3[p+inel][1] = actele->node[6]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel3[p+inel][2] = actele->node[2]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel3[p+inel][3] = actele->node[2]->Id_loc+1+FLUID_NODE_offset;
          cel3[p+inel][4] = actele->node[5]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel3[p+inel][5] = actele->node[5]->Id_loc+1+FLUID_NODE_offset;
          inel++;
          /*----------------------------------------- sub element 4*/
          cel3[p+inel][0] = actele->node[7]->Id_loc+1+FLUID_NODE_offset;
          cel3[p+inel][1] = actele->node[7]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel3[p+inel][2] = actele->node[3]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel3[p+inel][3] = actele->node[3]->Id_loc+1+FLUID_NODE_offset;
          cel3[p+inel][4] = actele->node[6]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel3[p+inel][5] = actele->node[6]->Id_loc+1+FLUID_NODE_offset;
          inel++;
          /*----------------------------------------- sub element 5*/
          p=FLUID_CEL4_offset;
          cel4[p+i][0] = actele->node[4]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+i][1] = actele->node[5]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+i][2] = actele->node[6]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+i][3] = actele->node[7]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+i][4] = actele->node[4]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel4[p+i][5] = actele->node[5]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel4[p+i][6] = actele->node[6]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel4[p+i][7] = actele->node[7]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          break;

        case quad9:
          inel = 0;
          p=FLUID_CEL4_offset + 4*i;
          /*----------------------------------------- sub-element 1 */
          cel4[p+inel][0] = actele->node[0]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][1] = actele->node[4]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][2] = actele->node[8]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][3] = actele->node[7]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][4] = actele->node[0]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel4[p+inel][5] = actele->node[4]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel4[p+inel][6] = actele->node[8]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel4[p+inel][7] = actele->node[7]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          inel++;
          /*----------------------------------------- sub-element 2 */
          cel4[p+inel][0] = actele->node[4]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][1] = actele->node[1]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][2] = actele->node[5]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][3] = actele->node[8]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][4] = actele->node[4]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel4[p+inel][5] = actele->node[1]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel4[p+inel][6] = actele->node[5]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel4[p+inel][7] = actele->node[8]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          inel++;
          /*----------------------------------------- sub-element 3 */
          cel4[p+inel][0] = actele->node[8]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][1] = actele->node[5]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][2] = actele->node[2]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][3] = actele->node[6]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][4] = actele->node[8]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel4[p+inel][5] = actele->node[5]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel4[p+inel][6] = actele->node[2]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel4[p+inel][7] = actele->node[6]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          inel++;
          /*----------------------------------------- sub-element 4 */
          cel4[p+inel][0] = actele->node[7]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][1] = actele->node[8]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][2] = actele->node[6]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][3] = actele->node[3]->Id_loc+1+FLUID_NODE_offset;
          cel4[p+inel][4] = actele->node[7]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel4[p+inel][5] = actele->node[8]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel4[p+inel][6] = actele->node[6]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          cel4[p+inel][7] = actele->node[3]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
          break;

       case tri3:
         inel=0;
         p=FLUID_CEL3_offset+i;
         cel3[p+inel][0] = actele->node[0]->Id_loc+1+FLUID_NODE_offset;
         cel3[p+inel][1] = actele->node[0]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
         cel3[p+inel][2] = actele->node[1]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
         cel3[p+inel][3] = actele->node[1]->Id_loc+1+FLUID_NODE_offset;
         cel3[p+inel][4] = actele->node[2]->Id_loc+1+FLUID_NODE_offset+numnp_struct;
         cel3[p+inel][5] = actele->node[2]->Id_loc+1+FLUID_NODE_offset;
        break;

        default:
          dserror("structure distyp not implemented yet!\n");
      }
    } /* end loop over all elements */
    /* ---------------------------------------------------- END OF STRUCT CELLS*/
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of V3CELL */


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
  CHUNK_DATA chunk;
  INT numnp_temp;
  CHAR *dis_names[]=DISTYPENAMES;

  /*The .Visual3.setup file has to be deleted, otherwise VISUAL wont
   * calculate a new transformation matrix which is needed to start
   * with the right distance to the object. */

  FILE* setup=fopen(".Visual3.setup", "w");
  fclose(setup);

  CHAR gimmick[]="|/-\\";
  INT gimmick_size=4;

  printf("\n  Reading the mesh files... ");

  init_problem_data(&problem, argc, argv);

  printf("done.\n");

  /*check process command line arguments*/
  for (i=1; i<argc-1; ++i)
  {
    CHAR* arg;
    arg = argv[i];
    if (arg[0]== '-')
    {
      if (arg[1]=='w') COLOUR=1;
      if (arg[1]=='g') COLOUR=2;
      if (arg[1]=='u') IOPT=1;
    }
  }
  if (!map_has_int(&(problem.control_table), "ndim", 2))
  {
    ACTDIM=3;
  }
  if (!map_has_int(&(problem.control_table), "ndim",3 ))
  {
    ACTDIM=2;
  }
  num_discr=problem.num_discr;

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
    /*we dont need the number of ale elements and nodes*/
    else if (problem.discr[i].type == ale)
    {
      ale_field = &(problem.discr[i]);
      ale_idx = i;
    }
    else
    {
      dserror("Unknown fieldtyp");
    }
  }
  /*--------------------------------------------------------------------*/

  /* some tests */
  switch (problem.type)
  {
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


  numnp_tot = (INT*)CCACALLOC(problem.num_discr, sizeof(INT));

  /* Iterate all discretizations. */
  for (i=0; i<problem.num_discr; ++i)
  {
    /* if we got hex20 we have to create new nodes and  the field's
     * numnp is changed.we have to change it back to the original
     * amount so we can use chunk_read functions. The numnp_tot array
     * contains the number of nodes including the new created ones */
    numnp_temp = problem.discr[i].numnp;
    init_post_discretization(&(discret[i]), &problem, &(problem.discr[i]), 1);
    numnp_tot[i] = problem.discr[i].numnp;
    problem.discr[i].numnp = numnp_temp;
    printf("  distyp : %s\n",dis_names[discret[i].element[0].distyp]) ;
  }

  /*set up the design information*/
  if (ACTDIM==3)
  {
    if (map_find_int(discret[fluid_idx].field->group, "ndvol", &ndvol))
    {
      init_post_design(&fluid_design, fluid_field);
      printf("\n  FLUID_DESIGN  erfolgreich initialisiert\n");
    }
    else
      printf("  fluid has got no design information\n");

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

  printf("  There are %d sets of results.\n", nsteps);


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
  printf("  Find number of results: done.\n");

  init_result_data(fluid_field, &global_fluid_result);
  if (ale_field!= NULL)
    init_result_data(ale_field, &global_ale_result);
  if (struct_field!=NULL)
  {
    init_result_data(struct_field, &global_struct_result);

    init_result_data(struct_field, &result);
    next_result(&result);
    init_chunk_data(&result,&chunk,"displacement");
    dis_dim=chunk.value_entry_length;
    destroy_chunk_data(&chunk);
    destroy_result_data(&result);
  }

#ifdef D_FSI
  /* Find coupled nodes. If there's at least an ale field. */
  if (ale_field != NULL)
  {
    if (discret[fluid_idx].element[0].distyp == hex20)
    {
      numnp_temp=discret[fluid_idx].field->numnp;
      discret[fluid_idx].field->numnp=numnp_tot[fluid_idx];
      discret[ale_idx].field->numnp=numnp_tot[fluid_idx];
    }
    fluid_ale_connect=(INT*)CCACALLOC(numnp_tot[fluid_idx], sizeof(INT));
    printf("  Find FSI coupling : ");

    post_fsi_initcoupling(&discret[struct_idx], &discret[fluid_idx], &discret[ale_idx], fluid_ale_connect);

    if (discret[fluid_idx].element[0].distyp == hex20)
    {
      discret[fluid_idx].field->numnp=numnp_temp;
      discret[ale_idx].field->numnp=numnp_temp;
    }


    printf(" done.\n");
  }
#endif


  vis3caf(fluid_idx,ale_idx,struct_idx);

  return 0;
}
void v2_init(
    char     *titl,
    INT      *iopt,
    INT      *cmncol,
    char     *cmfile,
    INT      *cmunit,
    INT      *xypix,
    float    *xymin,
    float    *xymax,
    INT      *nkeys,
    INT      *ikeys,
    INT      *tkeys,
    INT      *fkeys,
    float   **flims,
    INT      *mnode,
    INT      *mptri,
    INT      *mpptri,
    INT      *mface,
    INT      *mpface,
    INT      *medge,
    INT      *mpedge)
{
}
void v2_init_()
{
}
void v2_init__()
{
}
