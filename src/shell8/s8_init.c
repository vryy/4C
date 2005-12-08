/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/


#ifdef D_SHELL8


#include "../headers/standardtypes.h"
#include "shell8.h"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


/*----------------------------------------------------------------------*
 | initialize the element                                 m.gee 6/01    |
 *----------------------------------------------------------------------*/
void s8init(
    FIELD              *actfield,
    PARTITION          *actpart,
    INT                 disnum
    )

{

  INT          i,j,k;
  INT          nmaxw;
  INT          ngauss;
  INT         *ngp;
  ELEMENT     *actele;
  NODE        *actnode;
  S8_DATA      data;
  MATERIAL    *actmat;
  /*DOUBLE     **a3ref;*/
  INT          numa3;
  DOUBLE       a3[3];
  ARRAY        collaverdir_a;
  DOUBLE     **collaverdir;


#ifdef DEBUG
  dstrc_enter("s8init");
#endif


  for (i=0; i<actfield->dis[disnum].numele; i++)
  {
    actele = &(actfield->dis[disnum].element[i]);
    if (actele->eltyp != el_shell8) continue;
    /*------------------------------------------ init integration points */
    s8intg(actele,&data,0);/* ueberfluessig! ?*/
    /*---------------------------------------- init directors of element */
    s8a3(actele,&data,0);
    /*---------------------------------- allocate the space for stresses */
    am4def("forces",&(actele->e.s8->forces),1,18,MAXGAUSS,0,"D3");
    am4zero(&(actele->e.s8->forces));
  }
  /*--- loop elements in partition to allocate space for material history */
  for (i=0; i<actpart->pdis[disnum].numele; i++)
  {
    actele = actpart->pdis[disnum].element[i];
    actmat = &(mat[actele->mat-1]);
    if (actmat->mattyp==m_viscohyper)/* material is viscohyperelastic */
    {
      nmaxw  = actmat->m.viscohyper->nmaxw;
      ngp    = actele->e.s8->nGP;
      ngauss = ngp[0]*ngp[1]*ngp[2];
      actele->e.s8->his1 = CCACALLOC(1,sizeof(ARRAY4D));
      actele->e.s8->his2 = CCACALLOC(1,sizeof(ARRAY4D));
      am4def("mathis1",actele->e.s8->his1,ngauss,nmaxw+1,3,3,"D4");
      am4def("mathis2",actele->e.s8->his2,ngauss,nmaxw+1,3,3,"D4");
      am4zero(actele->e.s8->his1);
      am4zero(actele->e.s8->his2);
    }
  }
  /*--------------------- now do modification of directors bischoff style */
  /*------------------------------ allocate space for directors at a node */
  collaverdir = amdef("averdir",&collaverdir_a,3,MAXELE,"DA");
  /*----------------------------------------------------------------------*/
  for (i=0; i<actfield->dis[disnum].numnp; i++)
  {
    actnode = &(actfield->dis[disnum].node[i]);
    numa3=0;
    for (j=0; j<actnode->numele; j++)
    {
      actele = actnode->element[j];
      if (actele->eltyp != el_shell8) continue;
      for (k=0; k<actele->numnp; k++)
      {
        if (actele->node[k] == actnode)
        {
          collaverdir[0][numa3] = actele->e.s8->a3ref.a.da[0][k];
          collaverdir[1][numa3] = actele->e.s8->a3ref.a.da[1][k];
          collaverdir[2][numa3] = actele->e.s8->a3ref.a.da[2][k];
          numa3++;
          if (numa3 > MAXELE) dserror("Too many elements to a node, MAXELE too small");
          break;
        }
      }
    }
    /*----------- do nothing if there is only one shell8 element to actnode */
    if (numa3 <= 1) continue;
    /*------------------------------------------------ make shared director */
    s8averdir(collaverdir,numa3,a3);
    /*---------------------------------put shared director back to elements */
    for (j=0; j<actnode->numele; j++)
    {
      actele = actnode->element[j];
      if (actele->eltyp != el_shell8) continue;
      for (k=0; k<actele->numnp; k++)
      {
        if (actele->node[k] == actnode)
        {
          actele->e.s8->a3ref.a.da[0][k] = a3[0];
          actele->e.s8->a3ref.a.da[1][k] = a3[1];
          actele->e.s8->a3ref.a.da[2][k] = a3[2];
          break;
        }
      }
    }
  }/* end of loop over all nodes */
  /*----------------------------------------------------------------------*/
  amdel(&collaverdir_a);
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of s8init */
#endif

