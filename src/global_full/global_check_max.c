#ifdef CHECK_MAX

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

#ifdef D_FLUID2
#include "../fluid2/fluid2_tu.h"
#endif

#ifdef D_THERM2
#include "../therm2/therm2.h"
#endif

#ifdef D_THERM3
#include "../therm3/therm3.h"
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
          break;
#endif

#ifdef D_BRICK1
        case el_brick1:
          gauss = actele->e.c1->nGP[0]*actele->e.c1->nGP[1]*actele->e.c1->nGP[2];
          break;
#endif

#ifdef D_WALL1
        case el_wall1:
          gauss = actele->e.w1->nGP[0]*actele->e.w1->nGP[1];
          break;
#endif

#ifdef D_FLUID2
        case el_fluid2:
          gauss = actele->e.f2->nGP[0]*actele->e.f2->nGP[1];
          break;
#endif

#ifdef D_FLUID3
        case el_fluid3:
          if (actele->distyp == tet4 || actele->distyp == tet10)
            gauss = actele->e.f3->nGP[0];
          if (actele->distyp == hex8 || actele->distyp == hex20)
            gauss = actele->e.f3->nGP[0]*actele->e.f3->nGP[1]*actele->e.f3->nGP[2];
          break;
#endif

#ifdef D_FLUID3_F
        case el_fluid3_fast:
          if (actele->distyp == tet4 || actele->distyp == tet10)
            gauss = actele->e.f3->nGP[0];
          if (actele->distyp == hex8 || actele->distyp == hex20)
            gauss = actele->e.f3->nGP[0]*actele->e.f3->nGP[1]*actele->e.f3->nGP[2];
          break;
#endif

#ifdef D_ALE
        case el_ale3:
          gauss = actele->e.ale3->nGP[0]*actele->e.ale3->nGP[1]*actele->e.ale3->nGP[2];
          break;
#endif

#ifdef D_ALE
        case el_ale2:
          gauss = actele->e.ale2->nGP[0]*actele->e.ale2->nGP[1];
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
          break;
#endif

#ifdef D_FLUID2_PRO
        case el_fluid2_pro:
          gauss = actele->e.f2pro->nGP[0]*actele->e.f2pro->nGP[1];

          /* check maxdofpernode for fluid2_pro */
          if (2 > maxdofpernode)
            maxdofpernode = 2;
          break;
#endif

#ifdef D_FLUID3_PRO
        case el_fluid3_pro:
          gauss = actele->e.f3pro->nGP[0]*actele->e.f3pro->nGP[1]*actele->e.f3pro->nGP[2];

          /* check maxdofpernode for fluid3_pro */
          if (3 > maxdofpernode)
            maxdofpernode = 3;
          break;
#endif

#ifdef D_FLUID2
        case el_fluid2_tu:
          gauss = actele->e.f2_tu->nGP[0]*actele->e.f2_tu->nGP[1];
          break;
#endif

#ifdef D_THERM2
        case el_therm2:
          gauss = actele->e.th2->nGP[0]*actele->e.th2->nGP[1];
          break;
#endif

#ifdef D_THERM3
        case el_therm3:
          gauss = actele->e.th3->gpnum[0]
	        * actele->e.th3->gpnum[1]
	        * actele->e.th3->gpnum[2];
          break;
#endif

        case el_none:
          dserror("No element typ given");
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
printf("\n[35;1mChecking max values ...[m\n");
printf("value         given  optimal\n");
printf("----------------------------\n");
printf("maxnod:         %3i      %3i\n",MAXNOD,maxnod);
printf("maxele:         %3i      %3i\n",MAXELE,maxele);
printf("maxdofpernode:  %3i      %3i\n",MAXDOFPERNODE,maxdofpernode);
printf("maxgauss:       %3i      %3i\n",MAXGAUSS,maxgauss);
printf("\n");

/* checking for too small values */
if (maxnod > MAXNOD)
  dserror("MAXNOD is too small!!");
if (maxele > MAXELE)
  dserror("MAXELE is too small!!");
if (maxdofpernode > MAXDOFPERNODE)
  dserror("MAXDOFPERNODE is too small!!");
if (maxgauss > MAXGAUSS)
  dserror("MAXGAUSS is too small!!");


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

printf("\n");


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of check_max_sizes */

/*! @} (documentation module close)*/

#endif
