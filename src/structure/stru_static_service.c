/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |  routine to find the node and dof to control in field  m.gee 11/01   |
 |                                                                      |
 | actfield                             the physical field to search in |
 | control_node_global                               global node number |
 | control_dof                                number of dof to look for |
 | **node           address of pointer to hold controlled node (output) |
 | *cdof                      adress of INT to hold dof number (output) |
 |                                                                      |
 *----------------------------------------------------------------------*/
void calstatserv_findcontroldof(FIELD     *actfield,
                                INT        control_node_global,
                                INT        control_dof,
                                NODE     **node,
                                INT       *cdof)
{
INT        i;
#ifdef DEBUG
dstrc_enter("calstatserv_findcontroldof");
#endif
/*----------------------------------------------------------------------*/
*node = NULL;
for (i=0; i<actfield->dis[0].numnp; i++)
{
   if (actfield->dis[0].node[i].Id == control_node_global)
   {
      *node = &(actfield->dis[0].node[i]);
      *cdof = actfield->dis[0].node[i].dof[control_dof-1];
      break;
   }
}
if (!(*node)) dserror("Cannot find control node");
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of calstatserv_findcontroldof */



/*----------------------------------------------------------------------*
 |  routine to find the dof's for output                                |
 |                                                                      |
 | actfield                             the physical field to search in |
 | control_node_global                               global node number |
 | control_dof                                number of dof to look for |
 | **node           address of pointer to hold controlled node (output) |
 | *cdof                      adress of INT to hold dof number (output) |
 |                                                                      |
 *----------------------------------------------------------------------*/
void calstatserv_findreldofs(FIELD     *actfield,
                             INT       *reldisnode_ID,
                             INT       *reldis_dof,
                             INT        num_reldis,
                             INT       *reldof)
{
INT        i,j;
#ifdef DEBUG
dstrc_enter("calstatserv_findreldofs");
#endif
/*----------------------------------------------------------------------*/
for (j=0; j<num_reldis; j++)
{
  for (i=0; i<actfield->dis[0].numnp; i++)
  {
     if (actfield->dis[0].node[i].Id == reldisnode_ID[j])
     {
        reldof[j] = actfield->dis[0].node[i].dof[reldis_dof[j]-1];
        break;
     }
  }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of calstatserv_findreldofs */



/*----------------------------------------------------------------------*
 |  routine to get the actual stepsize in static nonlinear analysis     |
 |                                                                      |
 |                                                                      |
 *----------------------------------------------------------------------*/
void get_stepsize(INT         kstep,
                  STATIC_VAR *statvar)
{
INT        i,j;
INT        fromstep[20];
#ifdef DEBUG
dstrc_enter("get_stepsize");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<20; i++)   fromstep[i] = 0;
/*----------------------------------------------------------------------*/
/*--it is really from 1 and not from 0! -*/
for (i=1; i<statvar->numcurve; i++)
{
  for (j=0; j<i; j++)
  {
    fromstep[i] += statvar->actstep[j];
  }
}
/*----------------------------------------------------------------------*/
statvar->stepsize = 0.0;
for (i=0; i<statvar->numcurve; i++)
{
  if (kstep >= fromstep[i] && kstep <fromstep[i+1])
  {
    statvar->stepsize = statvar->actstepsize[i];
  }
}
if (kstep >= fromstep[statvar->numcurve-1])
{
  statvar->stepsize = statvar->actstepsize[statvar->numcurve-1];
}
if(statvar->stepsize ==0.0)  dserror("there must be something wrong with variable step size!");



/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* get_stepsize */
