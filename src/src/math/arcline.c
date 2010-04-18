#ifndef CCADISCRET
#ifdef D_AXISHELL

/*!-----------------------------------------------------------------------
\file
\brief collection of function for stlines and arclines

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*-----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../axishell/axishell.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;


/*!----------------------------------------------------------------------
\brief prototypes callable only in this file

*-----------------------------------------------------------------------*/
void arcline_xi(DLINE    *dline, NODE     *node, DOUBLE   *xi);
void stline_xi(DLINE    *dline, NODE     *node, DOUBLE   *xi);
void calc_arc_props(DLINE *dline);



/*!----------------------------------------------------------------------
\brief calculates xi for a node on a arcline

<pre>                                                              mn 05/03
This function calculates the relative distance xi of a node on a arcline
from the beginning of the dline
</pre>
\param *dline     DLINE   (i)   the arcline the node must be on
\param *node      NODE    (i)   the node that should be on the dline
\param *xi        DOUBLE  (o)   the relative distance

\warning There is nothing special to this routine
\return void
\sa

*----------------------------------------------------------------------*/
void arcline_xi(
    DLINE    *dline,
    NODE     *node,
    DOUBLE   *xi
    )
{

  DOUBLE   delta[2];
  DOUBLE   t,phi;
  DOUBLE   radius;
  DOUBLE   initang;
  DOUBLE   length,total_length;

#ifdef DEBUG
  dstrc_enter("arcline_xi");
#endif

  /* check whether this dline is an arcline */
  if (dline->typ != arcline)
    dserror("dline must be an arcline!!");


  delta[0]      = node->x[0] - dline->props.arcline->center[0];
  delta[1]      = node->x[1] - dline->props.arcline->center[0];
  radius        = dline->props.arcline->radius;
  initang       = dline->props.arcline->initang;
  total_length  = dline->props.arcline->total_length;

  /* check whether this node is really on this dline */
  /*if (EPS6 < sqrt((delta[0])*(delta[0]) + (delta[1])*(delta[1])) )
    dserror("Node is not on dline!!");*/


  /* calculation of the 'direction-angle' with the Theissen-formula */
  t = atan(delta[0]/delta[1]) + ( 1.0 - 0.5*SIGN(delta[0]) -
      0.5*SIGN(delta[0])*SIGN(delta[1])) * PI;

  /* calculate the angle phi,
   * defined clockwise starting at the negative y axe */
  phi = t + PI;
  if (phi>2*PI)
    phi -= 2*PI;

  /* calcute the arclength between the first dnode and the given node */
  length = radius*FABS(initang-phi);
  *xi     = length/total_length;

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of arcline_xi */



/*!----------------------------------------------------------------------
\brief calculates xi for a node on a stline

<pre>                                                              mn 05/03
This function calculates the relative distance xi of a node on a stline
from the beginning of the dline
</pre>
\param *dline     DLINE   (i)   the stline the node must be on
\param *node      NODE    (i)   the node that should be on the dline
\param *xi        DOUBLE  (o)   the relative distance

\warning There is nothing special to this routine
\return void
\sa

*----------------------------------------------------------------------*/
void stline_xi(
    DLINE    *dline,
    NODE     *node,
    DOUBLE   *xi
    )
{

  DOUBLE   length,total_length;
  DOUBLE   delta[2];

#ifdef DEBUG
  dstrc_enter("stline_xi");
#endif

  /* check whether this dline is an stline */
  if (dline->typ != stline)
    dserror("dline must be an stline!!");


  delta[0]      = node->x[0] - dline->dnode[0]->x[0];
  delta[1]      = node->x[1] - dline->dnode[0]->x[1];
  total_length  = dline->props.stline->total_length;

  /* calcute the length between the first dnode and the given node */
  length = sqrt((delta[0])*(delta[0]) + (delta[1])*(delta[1]));
  *xi     = length/total_length;

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of stline_xi */



/*!----------------------------------------------------------------------
\brief calculates geometric properties of an arcline

<pre>                                                              mn 05/03
This function calculates several geometric properties for an arcline.
These are:
 - the coordinates of the center point
</pre>
\param *dline     DLINE   (i)   the stline the node must be on

\warning There is nothing special to this routine
\return void
\sa

*----------------------------------------------------------------------*/
void calc_arc_props(
    DLINE      *dline
    )
{

  DOUBLE    radius;
  DOUBLE    initang;
  DOUBLE    t;
  DOUBLE    point[2];

#ifdef DEBUG
  dstrc_enter("calc_arc_props");
#endif

  radius = dline->props.arcline->radius;
  initang = dline->props.arcline->initang;

  /* only for 2D */
  /* dnode[0] should be the first node */
  point[0] = dline->dnode[0]->x[0];
  point[1] = dline->dnode[0]->x[1];

  /* calculate angle between center point and point[] */
  t = 1.5*PI - initang;
  if (t<0) t += 2*PI;

  dline->props.arcline->center[0] = point[0] - radius*cos(t);
  dline->props.arcline->center[1] = point[1] - radius*sin(t);

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of calc_arc_props */




/*!----------------------------------------------------------------------
\brief interpolate axishell conditions to the node

<pre>                                                              mn 05/03
This function interpolates the axishell load and thickness conditions of
a dline onto the nodes on this dline.
</pre>
\param *actdis     DISCRET   (i)   the discretization to be used

\warning There is nothing special to this routine
\return void
\sa

*----------------------------------------------------------------------*/
void interpolate_axishell_conds(
    DISCRET      *actdis
    )
{

  INT       i,j;
  DOUBLE    length;
  DOUBLE    xi;
  DOUBLE    xi2;
  DOUBLE    thick = 0.0;
  DOUBLE    pv=0.0,ph=0.0,px=0.0,pw=0.0;
  DLINE    *actdline;
  NODE     *actnode;
  GNODE    *actgnode;
  ELEMENT  *actele;

#ifdef DEBUG
  dstrc_enter("interpolate_axishell_conds");
#endif

  /* calculate the properties of the dlines */
  /* loop all dlines */
  for (i=0; i<design->ndline; i++)
  {
    actdline = &(design->dline[i]);

    if (actdline->thickness == NULL && actdline->axishellload == NULL )
      continue;

    switch (actdline->typ)
    {
      case stline:
        length = sqrt((actdline->dnode[0]->x[0]-actdline->dnode[1]->x[0])*
            (actdline->dnode[0]->x[0]-actdline->dnode[1]->x[0])
            + (actdline->dnode[0]->x[1]-actdline->dnode[1]->x[1])*
            (actdline->dnode[0]->x[1]-actdline->dnode[1]->x[1]));
        actdline->props.stline->total_length = length;
        break;
      case arcline:
        calc_arc_props(actdline);
        break;
      default:
        dserror("Unknown typ of dline for interpolation");
        break;
    }
  }

  /* now loop all elements and calculate the thickness and loads */
  for (i=0; i<actdis->numele; i++)
  {
    actele = &(actdis->element[i]);
    switch (actele->eltyp)
    {
      case el_axishell:
        for (j=0; j<2; j++)
        {
          actnode  = (actele->node[j]);
          actgnode = actnode->gnode;

          switch (actgnode->ondesigntyp)
          {
            case ondline:
              actdline = actgnode->d.dline;
              /* calculate xi for interpolation along line */
              switch (actdline->typ)
              {
                case stline:
                  stline_xi(actdline,actnode,&xi);
                  break;
                case arcline:
                  arcline_xi(actdline,actnode,&xi);
                  break;
                default:
                  dserror("Unknown typ of dline for interpolation");
                  break;
              }

              /* calculate xi2 for interpolation along vertical axis */
              xi2 = (actnode->x[1] - actdline->dnode[0]->x[1]) /
                (actdline->dnode[1]->x[1]-actdline->dnode[0]->x[1]);

              /* interpolate thickness along line */
              thick = actdline->thickness->value[0] +
                xi * (actdline->thickness->value[1] - actdline->thickness->value[0]);

              /* calculation of the loads depending on typ of interpolation */
              /* pv */
              if (actdline->axishellload->interpol_pv == 0)
              {
                /* interpolation along the curve */
                pv = actdline->axishellload->pv[0] +
                  xi * (actdline->axishellload->pv[1] - actdline->axishellload->pv[0]);
              }
              else
              {
                /* interpolation along the vertical axis */
                pv = actdline->axishellload->pv[0] +
                  xi2 * (actdline->axishellload->pv[1] - actdline->axishellload->pv[0]);
              }

              /* ph */
              if (actdline->axishellload->interpol_ph == 0)
              {
                /* interpolation along the curve */
                ph = actdline->axishellload->ph[0] +
                  xi * (actdline->axishellload->ph[1] - actdline->axishellload->ph[0]);
              }
              else
              {
                /* interpolation along the vertical axis */
                ph = actdline->axishellload->ph[0] +
                  xi2 * (actdline->axishellload->ph[1] - actdline->axishellload->ph[0]);
              }

              /* px */
              if (actdline->axishellload->interpol_px == 0)
              {
                /* interpolation along the curve */
                px = actdline->axishellload->px[0] +
                  xi * (actdline->axishellload->px[1] - actdline->axishellload->px[0]);
              }
              else
              {
                /* interpolation along the vertical axis */
                px = actdline->axishellload->px[0] +
                  xi2 * (actdline->axishellload->px[1] - actdline->axishellload->px[0]);
              }

              /* pw */
              if (actdline->axishellload->interpol_pw == 0)
              {
                /* interpolation along the curve */
                pw = actdline->axishellload->pw[0] +
                  xi * (actdline->axishellload->pw[1] - actdline->axishellload->pw[0]);
              }
              else
              {
                /* interpolation along the vertical axis */
                pw = actdline->axishellload->pw[0] +
                  xi2 * (actdline->axishellload->pw[1] - actdline->axishellload->pw[0]);
              }
              break;

            case ondnode:
              thick = actgnode->d.dnode->thickness->value[0];
              pv    = actgnode->d.dnode->axishellload->pv[0];
              ph    = actgnode->d.dnode->axishellload->ph[0];
              px    = actgnode->d.dnode->axishellload->px[0];
              pw    = actgnode->d.dnode->axishellload->pw[0];
              break;

            default:
              break;
          } /*END of switch ondesigntyp */

          actele->e.saxi->thick[j]= thick;
          actele->e.saxi->pv[j]   = pv;
          actele->e.saxi->ph[j]   = ph;
          actele->e.saxi->px[j]   = px;
          actele->e.saxi->pw[j]   = pw;
        } /* END of case axishell */

        break;

      default:
        break;
    } /* END of switch (actele->distyp) */
  } /* END of for all elements */


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of interpolate_axishell_conds */

#endif


#endif
