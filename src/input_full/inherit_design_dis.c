/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/


#include "../headers/standardtypes.h"


/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01  |
  | general problem data                                               |
  | global variable GENPROB genprob is defined in global_control.c     |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;


/*----------------------------------------------------------------------*
 | functions only accessible in this file                               |
 *----------------------------------------------------------------------*/
void cal_dirich_fac(
    GNODE *gnode,
    INT index
    );



/*----------------------------------------------------------------------*
 | inherit boundary conditions from design to GNODEs         m.gee 3/02 |
 *----------------------------------------------------------------------*/
void inherit_design_dis_dirichlet(
    FIELD         *actfield,
    DISCRET       *actdis)
{

  INT     i,j;
  GNODE  *actgnode;


#ifdef DEBUG
  dstrc_enter("inherit_design_dis_dirichlet");
#endif


#ifdef D_FSI
  if (actfield->fieldtyp == ale && genprob.create_ale == 1)
  {
    /* read the ale DIRICH conditions from the fluid design */

    for (i=0; i<actdis->ngnode; i++)
    {
      actgnode = &(actdis->gnode[i]);
      switch(actgnode->ondesigntyp)
      {
        case ondnode:
          if (actgnode->d.dnode->ale_dirich != NULL)
            actgnode->dirich = actgnode->d.dnode->ale_dirich;
          break;

        case ondline:
          if (actgnode->d.dline->ale_dirich != NULL)
            actgnode->dirich = actgnode->d.dline->ale_dirich;
          break;

        case ondsurf:
          if (actgnode->d.dsurf->ale_dirich != NULL)
            actgnode->dirich = actgnode->d.dsurf->ale_dirich;
          break;

        case ondvol:
          break;

        case ondnothing:
          dserror("ALE-GNODE not owned by any design object");
          break;

        default:
          dserror("Cannot inherit dirichlet condition");
          break;
      }


      /* calculate factors for spatial functions */
      if (actgnode->dirich != NULL)
      {
        for (j=0; j<MAXDOFPERNODE; j++)
        {
          if (actgnode->dirich->funct.a.iv[j] == 0)
            actgnode->d_funct[j] = 1.0;
          else
            cal_dirich_fac(actgnode,j);

        }  /* for (j=0; j<MAXDOFPERNODE; j++) */

      }  /* if (actgnode->dirich != NULL) */

    }  /* for (i=0; i<actdis->ngnode; i++) */

  }  /* if (actfield->fieldtyp == ale && genprob.create_ale == 1) */

  else
#endif
  {

    for (i=0; i<actdis->ngnode; i++)
    {
      actgnode = &(actdis->gnode[i]);
      switch(actgnode->ondesigntyp)
      {
        case ondnode:    actgnode->dirich = actgnode->d.dnode->dirich;   break;
        case ondline:    actgnode->dirich = actgnode->d.dline->dirich;   break;
        case ondsurf:    actgnode->dirich = actgnode->d.dsurf->dirich;   break;
        case ondvol:     actgnode->dirich = actgnode->d.dvol->dirich;    break;
        case ondnothing: dserror("GNODE not owned by any design object");break;
        default: dserror("Cannot inherit dirichlet condition");          break;
      }


      /* calculate factors for spatial functions */
      if (actgnode->dirich != NULL)
      {
        for (j=0; j<MAXDOFPERNODE; j++)
        {
          if (actgnode->dirich->funct.a.iv[j] == 0)
            actgnode->d_funct[j] = 1.0;
          else
            cal_dirich_fac(actgnode,j);

        }  /* for (j=0; j<MAXDOFPERNODE; j++) */

      }  /* if (actgnode->dirich != NULL) */

    }  /* for (i=0; i<actdis->ngnode; i++) */

  }  /* else (actfield->fieldtyp == ale && genprob.create_ale == 1) */


#ifdef DEBUG
  dstrc_exit();
#endif

  return;

} /* end of inherit_design_dis_dirichlet */





/*----------------------------------------------------------------------*
 | inherit couple conditions from design to GNODEs           m.gee 3/02 |
 *----------------------------------------------------------------------*/
void inherit_design_dis_couple(DISCRET *actdis)
{
INT     i;
GNODE  *actgnode;
#ifdef DEBUG
dstrc_enter("inherit_design_dis_couple");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<actdis->ngnode; i++)
{
   actgnode = &(actdis->gnode[i]);
   switch(actgnode->ondesigntyp)
   {
   case ondnode:    actgnode->couple = actgnode->d.dnode->couple;   break;
   case ondline:    actgnode->couple = actgnode->d.dline->couple;   break;
   case ondsurf:    actgnode->couple = actgnode->d.dsurf->couple;   break;
   case ondvol:     actgnode->couple = actgnode->d.dvol->couple;    break;
   case ondnothing: dserror("GNODE not owned by any design object");break;
   default: dserror("Cannot inherit couple condition");             break;
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inherit_design_dis_couple */

#ifdef D_FSI
/*----------------------------------------------------------------------*
 | inherit boundary conditions from design to GNODEs         genk 10/02 |
 *----------------------------------------------------------------------*/
void inherit_design_dis_fsicouple(DISCRET *actdis)
{
INT     i;
GNODE  *actgnode;
#ifdef DEBUG
dstrc_enter("inherit_design_dis_fsicouple");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<actdis->ngnode; i++)
{
   actgnode = &(actdis->gnode[i]);
   switch(actgnode->ondesigntyp)
   {
   case ondnode:    actgnode->fsicouple = actgnode->d.dnode->fsicouple;   break;
   case ondline:    actgnode->fsicouple = actgnode->d.dline->fsicouple;   break;
   case ondsurf:    actgnode->fsicouple = actgnode->d.dsurf->fsicouple;   break;
   case ondvol:     actgnode->fsicouple = NULL;                           break;
   case ondnothing: dserror("GNODE not owned by any design object");      break;
   default: dserror("Cannot inherit dirichlet condition");                break;
   }
}
/* addition due to fsi mortar requirements chfoe 08/04 */
/*-------------------------------------- GLINE inherits from its DLINE */
for (i=0; i<actdis->ngline; i++)
{
   /*------------------------------- check whether gline is on a dline */
   if (actdis->gline[i].dline)
   actdis->gline[i].fsicouple = actdis->gline[i].dline->fsicouple;
   else
   actdis->gline[i].fsicouple = NULL;
}
/*-------------------------------------- GSURF inherits from its DSURF */
for (i=0; i<actdis->ngsurf; i++)
{
   /*------------------------------- check whether gline is on a dline */
   if (actdis->gsurf[i].dsurf)
   actdis->gsurf[i].fsicouple = actdis->gsurf[i].dsurf->fsicouple;
   else
   actdis->gsurf[i].fsicouple = NULL;
}
/* end addition chfoe 08/04 */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inherit_design_dis_fsicouple */

#endif


#ifdef D_SSI
/*----------------------------------------------------------------------*
 | inherit boundary conditions from design to GNODEs         genk 10/02 |
 *----------------------------------------------------------------------*/
void inherit_design_dis_ssicouple(DISCRET *actdis)
{
INT     i;
GNODE  *actgnode;
#ifdef DEBUG
dstrc_enter("inherit_design_dis_ssicouple");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<actdis->ngnode; i++)
{
   actgnode = &(actdis->gnode[i]);
   switch(actgnode->ondesigntyp)
   {
   case ondnode:    actgnode->ssicouple = actgnode->d.dnode->ssicouple;   break;
   case ondline:    actgnode->ssicouple = actgnode->d.dline->ssicouple;   break;
/* case ondsurf:    actgnode->ssicouple = actgnode->d.dsurf->ssicouple;   break; */
   case ondvol:     actgnode->ssicouple = NULL;                           break;
   case ondnothing: dserror("GNODE not owned by any design object");      break;
   default: dserror("Cannot inherit ssi-couple condition");               break;
   }
}
/*-------------------------------------- GLINE inherits from its DLINE */
for (i=0; i<actdis->ngline; i++)
{
   /*------------------------------- check whether gline is on a dline */
   if (actdis->gline[i].dline)
   actdis->gline[i].ssicouple = actdis->gline[i].dline->ssicouple;
   else
   actdis->gline[i].ssicouple = NULL;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inherit_design_dis_ssicouple */

#endif

#ifdef D_FLUID
/*----------------------------------------------------------------------*
 | inherit BCs from design to GNODEs & GLINES                genk 01/03 |
 *----------------------------------------------------------------------*/
void inherit_design_dis_freesurf(DISCRET *actdis)
{
#ifdef D_FSI
INT     i;
GNODE  *actgnode;

#ifdef DEBUG
dstrc_enter("inherit_design_dis_freesurf");
#endif

/*-------------------------------------------------------------- GNODES */
for (i=0; i<actdis->ngnode; i++)
{
   actgnode = &(actdis->gnode[i]);
   switch(actgnode->ondesigntyp)
   {
   case ondnode:    actgnode->freesurf = actgnode->d.dnode->freesurf;   break;
   case ondline:    actgnode->freesurf = actgnode->d.dline->freesurf;   break;
   case ondsurf:    actgnode->freesurf = actgnode->d.dsurf->freesurf;   break;
   case ondvol:     actgnode->freesurf = NULL;                          break;
   case ondnothing: dserror("GNODE not owned by any design object");    break;
   default: dserror("Cannot inherit dirichlet condition");              break;
   }
}
/*-------------------------------------- GLINE inherits from its DLINE */
for (i=0; i<actdis->ngline; i++)
{
   /*------------------------------- check whether gline is on a dline */
   if (actdis->gline[i].dline)
   actdis->gline[i].freesurf = actdis->gline[i].dline->freesurf;
   else
   actdis->gline[i].freesurf = NULL;
}
/*-------------------------------------- GSURF inherits from its DSURF */
for (i=0; i<actdis->ngsurf; i++)
{
   /*------------------------------- check whether gline is on a dline */
   if (actdis->gsurf[i].dsurf)
   actdis->gsurf[i].freesurf = actdis->gsurf[i].dsurf->freesurf;
   else
   actdis->gsurf[i].freesurf = NULL;
}


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
#endif
return;
} /* end of inherit_design_dis_freesurf */

/*----------------------------------------------------------------------*
 | inherit BCs from design to GNODEs & GLINES                genk 01/03 |
 *----------------------------------------------------------------------*/
void inherit_design_dis_slipdirich(DISCRET *actdis)
{
#ifdef D_FSI
INT     i;
GNODE  *actgnode;

#ifdef DEBUG
dstrc_enter("inherit_design_dis_slipdirich");
#endif

/*-------------------------------------------------------------- GNODES */
for (i=0; i<actdis->ngnode; i++)
{
   actgnode = &(actdis->gnode[i]);
   switch(actgnode->ondesigntyp)
   {
   case ondnode: break;
   case ondline:
      actgnode->slipdirich = actgnode->d.dline->slipdirich;             break;
   case ondsurf: break;
   case ondvol:  break;
   case ondnothing: dserror("GNODE not owned by any design object");    break;
   default: dserror("Cannot inherit dirichlet condition");              break;
   }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
#endif
return;
} /* end of inherit_design_dis_slipdirich */


#endif
/*----------------------------------------------------------------------*
 | inherit neumann conditions from design to discretization  m.gee 3/02 |
 *----------------------------------------------------------------------*/
void inherit_design_dis_neum(DISCRET *actdis)
{
INT     i;
GNODE  *actgnode;
#ifdef DEBUG
dstrc_enter("inherit_design_dis_neum");
#endif
/*----------------------------------------------------------------------*/
/* GNODES inherit neumann condition from DNODES if they are on a DNODE */
/*                                   otherwise, GNODE inherits nothing */
for (i=0; i<actdis->ngnode; i++)
{
  actgnode = &(actdis->gnode[i]);
  switch(actgnode->ondesigntyp)
  {
    case ondnode:
      actgnode->neum   = actgnode->d.dnode->neum;
      break;
    case ondnothing:
      dserror("GNODE not owned by any design object");
      break;
    default:
      break;
  }
}
/*-------------------------------------- GLINE inherits from its DLINE */
for (i=0; i<actdis->ngline; i++)
{
   /*------------------------------- check whether gline is on a dline */
   if (actdis->gline[i].dline)
   actdis->gline[i].neum = actdis->gline[i].dline->neum;
}
/*----------------------- GSURF inherit from DSURF, if it is on a DSURF */
for (i=0; i<actdis->ngsurf; i++)
{
   if (actdis->gsurf[i].dsurf)
   actdis->gsurf[i].neum = actdis->gsurf[i].dsurf->neum;
}
/* GVOL inherits from DVOL if it is inside a DVOL (which always should be true)*/
for (i=0; i<actdis->ngvol; i++)
{
   dsassert(actdis->gvol[i].dvol!=NULL,"volume fe-element outside any design volume");
   actdis->gvol[i].neum = actdis->gvol[i].dvol->neum;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inherit_design_dis_neum */
