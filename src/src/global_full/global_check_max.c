/*!----------------------------------------------------------------------
\file
\brief contains functions to check the values for the max sizes

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#include "../headers/standardtypes.h"

#ifdef D_ALE
#include "../ale3/ale3.h"
#include "../ale2/ale2.h"
#endif

#ifdef D_SHELL8
#include "../shell8/shell8.h"
#endif

#ifdef D_SHELL9
#include "../shell9/shell9.h"
#endif

#ifdef D_BRICK1
#include "../brick1/brick1.h"
#endif

#ifdef D_WALL1
#include "../wall1/wall1.h"
#endif

#ifdef D_FLUID2
#include "../fluid2/fluid2.h"
#endif

#ifdef D_FLUID3
#include "../fluid3/fluid3.h"
#endif

#ifdef D_FLUID3_FAST
#include "../fluid3/fluid3.h"
#endif

#ifdef D_AXISHELL
#endif

#ifdef D_BEAM3
#include "../beam3/beam3.h"
#endif

#ifdef D_FLUID2_PRO
#include "../fluid2_pro/fluid2pro.h"
#endif

#ifdef D_FLUID3_PRO
#include "../fluid3_pro/fluid3pro.h"
#endif

#ifdef D_FLUID2_IS
#include "../fluid2_is/fluid2_is.h"
#endif

#ifdef D_FLUID3_IS
#include "../fluid3_is/fluid3_is.h"
#endif

#ifdef D_FLUID2
#include "../fluid2/fluid2_tu.h"
#endif

#ifdef D_THERM2
#include "../therm2/therm2.h"
#endif

#ifdef D_THERM3
#include "../therm3/therm3.h"
#endif

#ifdef D_SOLID3
#include "../solid3/solid3.h"
#endif




/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

/*!
\addtogroup GLOBAL
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief check the values of the max sizes

<pre>                                                              mn 04/04
This routine determines the optimal values for maxele, maxnod, maxdofpernode
and maxgauss and compares those to the given values of the respective defines.
</pre>

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ---

*----------------------------------------------------------------------*/
void check_max_sizes(
    )
{
  INT         i,j,k,l;
  INT         maxnod;
  INT         maxele;
  INT         maxdofpernode;
  INT         maxgauss;
  INT         gauss;

  FIELD      *actfield;
  ELEMENT    *actele;

  INT maxnod_axishell = -1;
  INT maxnod_beam3 = -1;
  INT maxnod_brick1 = -1;
  INT maxnod_f2 = -1;
  INT maxnod_f3 = -1;
  INT maxnod_ale3 = -1;
  INT maxnod_shell8 = -1;
  INT maxnod_shell9 = -1;
  /* INT maxnodestress_shell9 = -1; */
  INT maxnod_wall1 = -1;

#ifdef DEBUG
  dstrc_enter("check_max_sizes");
#endif

/* initialize */
maxnod        = 0;
maxele        = 0;
maxdofpernode = 0;
maxgauss      = 0;
gauss         = 0;

/* determine correct max sizes */
for(i=0; i<genprob.numfld; i++)  /* loop all fields */
{
  actfield = &(field[i]);
  for(j=0; j<actfield->ndis; j++)  /* loop all discretisations */
  {
    for (k=0; k<actfield->dis[j].numele; k++)  /* loop all elements */
    {
      actele = &(actfield->dis[j].element[k]);

      /* check maxnod */
      if (actele->numnp > maxnod)
        maxnod = actele->numnp;

      /* check maxele  and maxdofpernode */
      for (l=0; l<actele->numnp; l++)
      {
        if (actele->node[l]->numele > maxele)
          maxele = actele->node[l]->numele;

        if (actele->node[l]->numdf > maxdofpernode)
          maxdofpernode = actele->node[l]->numdf;
      }

      /* check max gauss */
      switch(actele->eltyp)
      {

#ifdef D_SHELL8
        case el_shell8:
          gauss = actele->e.s8->nGP[0]*actele->e.s8->nGP[1];
	  if (actele->numnp > maxnod_shell8)
	    maxnod_shell8 = actele->numnp;
          break;
#endif

#ifdef D_SHELL9
        case el_shell9:
          gauss=0;
          for (l=0; l<actele->e.s9->num_klay; l++)
          {
            gauss += actele->e.s9->kinlay[l].num_mlay;
          }
          gauss =gauss * actele->e.s9->nGP[0]*actele->e.s9->nGP[1]*actele->e.s9->nGP[2];

          /* check maxdofpernode for shell9 */
          if (actele->e.s9->numdf > maxdofpernode)
            maxdofpernode = actele->e.s9->numdf;
	  if (actele->numnp > maxnod_shell9)
	    maxnod_shell9 = actele->numnp;
          break;
#endif

#ifdef D_BRICK1
        case el_brick1:
          gauss = actele->e.c1->nGP[0]*actele->e.c1->nGP[1]*actele->e.c1->nGP[2];
	  if (actele->numnp > maxnod_brick1)
	    maxnod_brick1 = actele->numnp;
          break;
#endif

#ifdef D_WALL1
        case el_wall1:
          gauss = actele->e.w1->nGP[0]*actele->e.w1->nGP[1];
	  if (actele->numnp > maxnod_wall1)
	    maxnod_wall1 = actele->numnp;
          break;
#endif

#ifdef D_FLUID2
        case el_fluid2:
          gauss = actele->e.f2->nGP[0]*actele->e.f2->nGP[1];
	  if (actele->numnp > maxnod_f2)
	    maxnod_f2 = actele->numnp;
          break;
#endif

#ifdef D_FLUID3
        case el_fluid3:
          if (actele->distyp == tet4 || actele->distyp == tet10)
            gauss = actele->e.f3->nGP[0];
          if (actele->distyp == hex8 || actele->distyp == hex20 || actele->distyp == hex27)
            gauss = actele->e.f3->nGP[0]*actele->e.f3->nGP[1]*actele->e.f3->nGP[2];
	  if (actele->numnp > maxnod_f3)
	    maxnod_f3 = actele->numnp;
          break;
#endif

#ifdef D_FLUID3_F
        case el_fluid3_fast:
          if (actele->distyp == tet4 || actele->distyp == tet10)
            gauss = actele->e.f3->nGP[0];
          if (actele->distyp == hex8 || actele->distyp == hex20)
            gauss = actele->e.f3->nGP[0]*actele->e.f3->nGP[1]*actele->e.f3->nGP[2];
	  if (actele->numnp > maxnod_f3)
	    maxnod_f3 = actele->numnp;
          break;
#endif

#ifdef D_ALE
        case el_ale3:
          gauss = actele->e.ale3->nGP[0]*actele->e.ale3->nGP[1]*actele->e.ale3->nGP[2];
	  if (actele->numnp > maxnod_ale3)
	    maxnod_ale3 = actele->numnp;
          break;
#endif

#ifdef D_ALE
        case el_ale2:
          gauss = actele->e.ale2->nGP[0]*actele->e.ale2->nGP[1];
	  if (actele->numnp > maxnod_ale3)
	    maxnod_ale3 = actele->numnp;
          break;
#endif

#ifdef D_AXISHELL
        case el_axishell:
          gauss = 5;
          break;
#endif

#ifdef D_BEAM3
        case el_beam3:
          gauss = actele->e.b3->nGP[0];
	  if (actele->numnp > maxnod_beam3)
	    maxnod_beam3 = actele->numnp;
          break;
#endif

#ifdef D_FLUID2_PRO
        case el_fluid2_pro:
          gauss = actele->e.f2pro->nGP[0]*actele->e.f2pro->nGP[1];

          /* check maxdofpernode for fluid2_pro */
          if (2 > maxdofpernode)
            maxdofpernode = 2;
	  if (actele->numnp > maxnod_f2)
	    maxnod_f2 = actele->numnp;
          break;
#endif

#ifdef D_FLUID3_PRO
        case el_fluid3_pro:
          gauss = actele->e.f3pro->nGP[0]*actele->e.f3pro->nGP[1]*actele->e.f3pro->nGP[2];

          /* check maxdofpernode for fluid3_pro */
          if (3 > maxdofpernode)
            maxdofpernode = 3;
	  if (actele->numnp > maxnod_f3)
	    maxnod_f3 = actele->numnp;
          break;
#endif

#ifdef D_FLUID2_IS
        case el_fluid2_is:
          gauss = actele->e.f2is->nGP[0]*actele->e.f2is->nGP[1];

          /* check maxdofpernode for fluid2_is */
          if (3 > maxdofpernode)
            maxdofpernode = 3;
	  if (actele->numnp > maxnod_f2)
	    maxnod_f2 = actele->numnp;
          break;
#endif

#ifdef D_FLUID3_IS
        case el_fluid3_is:
          gauss = actele->e.f3is->nGP[0]*actele->e.f3is->nGP[1]*actele->e.f3is->nGP[2];

          /* check maxdofpernode for fluid3_is */
          if (4 > maxdofpernode)
            maxdofpernode = 4;
	  if (actele->numnp > maxnod_f3)
	    maxnod_f3 = actele->numnp;
          break;
#endif

#ifdef D_FLUID2
        case el_fluid2_tu:
          gauss = actele->e.f2_tu->nGP[0]*actele->e.f2_tu->nGP[1];
	  if (actele->numnp > maxnod_f2)
	    maxnod_f2 = actele->numnp;
          break;
#endif

#ifdef D_THERM2
        case el_therm2:
          gauss = actele->e.th2->nGP[0]*actele->e.th2->nGP[1];
          break;
#endif

#ifdef D_THERM3
        case el_therm3:
          gauss = actele->e.th3->gptot;
          break;
#endif

#ifdef D_SOLID3
        case el_solid3:
          gauss = actele->e.so3->gptot;
          break;
#endif

        case el_none:
          dserror("No element type given");
          break;

        default:
          dserror("Typ of element unknown");
      }
      if (gauss > maxgauss)
        maxgauss = gauss;

    } /* end of loop all elements */
  } /* end of loop all diss */
} /* end of loop all fields */


/* print result of check to screen */
printf("\n" MAGENTA_LIGHT "Checking max values ..." END_COLOR "\n");
printf("value         given  optimal\n");
printf("----------------------------\n");

#define MACRO_SIZE_OUTPUT(variable,macro) \
  if (variable > -1)						\
    printf("%-13s   %3i      %3i\n",#variable,macro,variable);

  MACRO_SIZE_OUTPUT(maxnod,MAXNOD);
  MACRO_SIZE_OUTPUT(maxele,MAXELE);
  MACRO_SIZE_OUTPUT(maxdofpernode,MAXDOFPERNODE);
  MACRO_SIZE_OUTPUT(maxgauss,MAXGAUSS);

  MACRO_SIZE_OUTPUT(maxnod_axishell,MAXNOD_AXISHELL);
  MACRO_SIZE_OUTPUT(maxnod_beam3,MAXNOD_BEAM3);
  MACRO_SIZE_OUTPUT(maxnod_brick1,MAXNOD_BRICK1);
  MACRO_SIZE_OUTPUT(maxnod_f2,MAXNOD_F2);
  MACRO_SIZE_OUTPUT(maxnod_f3,MAXNOD_F3);
  MACRO_SIZE_OUTPUT(maxnod_ale3,MAXNOD_ALE3);
  MACRO_SIZE_OUTPUT(maxnod_shell8,MAXNOD_SHELL8);
  MACRO_SIZE_OUTPUT(maxnod_shell9,MAXNOD_SHELL9);
  MACRO_SIZE_OUTPUT(maxnod_wall1,MAXNOD_WALL1);

  printf("\n");

/* checking for too small values */

#define MACRO_SIZE_CHECK(variable,macro) \
  if (variable > macro)			 \
    dserror(#macro " is too small");

  MACRO_SIZE_CHECK(maxnod,MAXNOD);
  MACRO_SIZE_CHECK(maxele,MAXELE);
  MACRO_SIZE_CHECK(maxdofpernode,MAXDOFPERNODE);
  MACRO_SIZE_CHECK(maxgauss,MAXGAUSS);

  MACRO_SIZE_CHECK(maxnod_axishell,MAXNOD_AXISHELL);
  MACRO_SIZE_CHECK(maxnod_beam3,MAXNOD_BEAM3);
  MACRO_SIZE_CHECK(maxnod_brick1,MAXNOD_BRICK1);
  MACRO_SIZE_CHECK(maxnod_f2,MAXNOD_F2);
  MACRO_SIZE_CHECK(maxnod_f3,MAXNOD_F3);
  MACRO_SIZE_CHECK(maxnod_ale3,MAXNOD_ALE3);
  MACRO_SIZE_CHECK(maxnod_shell8,MAXNOD_SHELL8);
  MACRO_SIZE_CHECK(maxnod_shell9,MAXNOD_SHELL9);
  MACRO_SIZE_CHECK(maxnod_wall1,MAXNOD_WALL1);

#if 0
/* checking for too large values */
if (maxnod < MAXNOD)
{
  printf("[37;1mMAXNOD is too large!![m\n");
  dswarning(1,8);
}
if (maxele < MAXELE)
{
  printf("[37;1mMAXELE is too large!![m\n");
  dswarning(1,9);
}
if (maxdofpernode < MAXDOFPERNODE)
{
  printf("[37;1mMAXDOFPERNODE is too large!![m\n");
  dswarning(1,10);
}
if (maxgauss < MAXGAUSS)
{
  printf("[37;1mMAXGAUSS is too large!![m\n");
  dswarning(1,11);
}
#endif

printf("\n");


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of check_max_sizes */

/*! @} (documentation module close)*/

#endif
