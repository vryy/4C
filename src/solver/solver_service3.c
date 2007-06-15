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
#include "../solver/solver.h"
/*!----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern INT            numcurve;
extern struct _CURVE *curve;


/*----------------------------------------------------------------------*
 |                                                            m.gee 4/03|
 | init sol[place] to zero                                              |
 |                                                                      |
 | FIELD *actfield (i) active field                                     |
 | INT    disnum   (i) indize of the discretization in actfield to be used|
 | INT    arraynum (i) indize of the array, 0 = sol                     |
 |                                          1 = sol_increment           |
 |                                          2 = sol_residual            |
 |                                          3 = sol_mf                  |
 | INT    place    (i) row in ARRAY sol to be set to zero               |
 *----------------------------------------------------------------------*/
void solserv_sol_zero(FIELD *actfield, INT disnum, NODE_ARRAY arraynum, INT place)
{
INT               i,j;
INT               diff,max;
ARRAY            *array;
NODE             *actnode;
DISCRET          *actdis;
#ifdef DEBUG
dstrc_enter("solserv_sol_zero");
#endif
/*----------------------------------------------------------------------*/
actdis = &(actfield->dis[disnum]);
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   /*--------------------------------------------- select correct array */
   switch(arraynum)
   {
   case 0:
      array = &(actnode->sol);
   break;
   case 1:
      array = &(actnode->sol_increment);
   break;
   case 2:
      array = &(actnode->sol_residual);
   break;
   case 3:
      array = &(actnode->sol_mf);
   break;
   default:
      array = NULL;
      dserror("Only 0,1,2,3 allowed for arraynum to select sol, sol_increment, sol_residual, sol_mf");
   }
   /* check the size of array */
   if (place >= array->fdim)
   {
      diff = place - array->fdim;
      /*max  = IMAX(diff,5);*/
      max = diff;
      amredef(array,array->fdim+max+1,array->sdim,"DA");
   }
   for (j=0; j<array->sdim; j++)
      array->a.da[place][j] = 0.0;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of solserv_sol_zero */




/*----------------------------------------------------------------------*
 |                                                            m.gee 4/03|
 | copy values from array arrayfrom in place from      to               |
 |                  array arrayto   in place to                         |
 |                                                                      |
 | FIELD *actfield  (i) active field                                    |
 | INT    disnum    (i) indize of the discretization in actfield to be used|
 | INT    arrayfrom (i) indize of the array, 0 = sol                    |
 |                                           1 = sol_increment          |
 |                                           2 = sol_residual           |
 |                                           3 = sol_mf                 |
 | INT    arrayto   (i) indize of the array, 0 = sol                    |
 |                                           1 = sol_increment          |
 |                                           2 = sol_residual           |
 |                                           3 = sol_mf                 |
 | INT    from    (i) row in ARRAY sol to be set to zero                |
 *----------------------------------------------------------------------*/
void solserv_sol_copy(FIELD *actfield, INT disnum, NODE_ARRAY arrayfrom, NODE_ARRAY arrayto,
                      INT from, INT to)
{
INT               i,j;
INT               diff,max, min;
ARRAY            *arrayf,*arrayt;
NODE             *actnode;
DISCRET          *actdis;
#ifdef DEBUG
dstrc_enter("solserv_sol_copy");
#endif
/*----------------------------------------------------------------------*/
actdis = &(actfield->dis[disnum]);
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   /*----------------------------------------- select correct arrayfrom */
   switch(arrayfrom)
   {
   case 0:
      arrayf = &(actnode->sol);
   break;
   case 1:
      arrayf = &(actnode->sol_increment);
   break;
   case 2:
      arrayf = &(actnode->sol_residual);
   break;
   case 3:
      arrayf = &(actnode->sol_mf);
   break;
   default:
      arrayf = NULL;
      dserror("Only 0,1,2,3 allowed for arrayfrom to select sol, sol_increment, sol_residual, sol_mf");
   }
   /*----------------------------------------- select correct arrayfrom */
   switch(arrayto)
   {
   case 0:
      arrayt = &(actnode->sol);
   break;
   case 1:
      arrayt = &(actnode->sol_increment);
   break;
   case 2:
      arrayt = &(actnode->sol_residual);
   break;
   case 3:
      arrayt = &(actnode->sol_mf);
   break;
   default:
      arrayt = NULL;
      dserror("Only 0,1,2,3 allowed for arrayfrom to select sol, sol_increment, sol_residual, sol_mf");
   }
   /* check the size of arrayf */
   if (from >= arrayf->fdim)
      dserror("Cannot copy from array, because place doesn't exist");
   /* check the size of arrayt */
   if (to >= arrayt->fdim)
   {
      diff = to - arrayt->fdim;
      /*max  = IMAX(diff,5);*/
      max = diff;
      amredef(arrayt,arrayt->fdim+max+1,arrayt->sdim,"DA");
   }

   /* There are strange cases where the source array has less columns
    * than the destination. (sol_mf to sol in fluid/ale) */
   min = MIN(arrayt->sdim, arrayf->sdim);
   for (j=0; j<min; j++)
      arrayt->a.da[to][j] = arrayf->a.da[from][j];

   /* In this case zero the remaining fields. */
   for (j=0; j<arrayt->sdim-min; j++)
      arrayt->a.da[to][min+j] = 0;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of solserv_sol_copy */




#ifdef SUBDIV
void solserv_sol_trans(
    FIELD              *actfield,
    INT                 disnum,
    NODE_ARRAY          array,
    INT                 place
    )
{

  INT               i,j;
  INT               min;
  ARRAY            *arrayf,*arrayt;
  NODE             *calc_node;
  NODE             *io_node;
  DISCRET          *actdis;


#ifdef DEBUG
  dstrc_enter("solserv_sol_trans");
#endif


  actdis = &(actfield->dis[disnum]);

  /* loop calc dis */
  for (i=0; i<actdis->numnp; i++)
  {
    calc_node = &(actdis->node[i]);

    /* find corresponding io node */
    if (calc_node->master_node != NULL)
      io_node = calc_node->master_node;
    else
      continue;

    /* select correct arrayf */
    switch(array)
    {
      case 0:
        arrayf = &(calc_node->sol);
        break;

      case 1:
        arrayf = &(calc_node->sol_increment);
        break;

      case 2:
        arrayf = &(calc_node->sol_residual);
        break;

      case 3:
        arrayf = &(calc_node->sol_mf);
        break;

      default:
        arrayf = NULL;
        dserror
          ("Only 0,1,2,3 allowed for arrayfrom to select sol, sol_increment, sol_residual, sol_mf");
        break;
    }


    /* select correct arrayt */
    switch(array)
    {
      case 0:
        arrayt = &(io_node->sol);
        break;

      case 1:
        arrayt = &(io_node->sol_increment);
        break;

      case 2:
        arrayt = &(io_node->sol_residual);
        break;

      case 3:
        arrayt = &(io_node->sol_mf);
        break;

      default:
        arrayt = NULL;
        dserror
          ("Only 0,1,2,3 allowed for arrayfrom to select sol, sol_increment, sol_residual, sol_mf");
        break;
    }


    /* check the size of arrayf */
    if (place >= arrayf->fdim)
      dserror("Cannot copy from array, because place doesn't exist");

    /* check the size of arrayt */
    if (place >= arrayt->fdim)
    {
      amredef(arrayt,place+1,arrayt->sdim,"DA");
    }

    /* There are strange cases where the source array has less columns
     * than the destination. (sol_mf to sol in fluid/ale) */
    min = MIN(arrayt->sdim, arrayf->sdim);
    for (j=0; j<min; j++)
      arrayt->a.da[place][j] = arrayf->a.da[place][j];

    /* In this case zero the remaining fields. */
    for (j=0; j<arrayt->sdim-min; j++)
      arrayt->a.da[place][min+j] = 0;

  }  /* for (i=0; i<actdis->numnp; i++) */


#ifdef DEBUG
  dstrc_exit();
#endif

  return;

} /* end of solserv_sol_trans */

#endif /* ifdef SUBDIV */







/*----------------------------------------------------------------------*
 |                                                            m.gee 4/03|
 | add values from  array arrayfrom in place from      to               |
 |                  array arrayto   in place to                         |
 | scale fromvalues by fac
 |                                                                      |
 | FIELD *actfield  (i) active field                                    |
 | INT    disnum    (i) indize of the discretization in actfield to be used|
 | INT    arrayfrom (i) indize of the array, 0 = sol                    |
 |                                           1 = sol_increment          |
 |                                           2 = sol_residual           |
 |                                           3 = sol_mf                 |
 | INT    arrayto   (i) indize of the array, 0 = sol                    |
 |                                           1 = sol_increment          |
 |                                           2 = sol_residual           |
 |                                           3 = sol_mf                 |
 | INT    from    (i) row in ARRAY sol to be set to zero                |
 *----------------------------------------------------------------------*/
void solserv_sol_add(FIELD *actfield, INT disnum, NODE_ARRAY arrayfrom, NODE_ARRAY arrayto,
                      INT from, INT to, DOUBLE fac)
{
INT               i,j;
INT               diff,max;
ARRAY            *arrayf,*arrayt;
NODE             *actnode;
DISCRET          *actdis;
#ifdef DEBUG
dstrc_enter("solserv_sol_add");
#endif
/*----------------------------------------------------------------------*/
actdis = &(actfield->dis[disnum]);
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   /*----------------------------------------- select correct arrayfrom */
   switch(arrayfrom)
   {
   case 0:
      arrayf = &(actnode->sol);
   break;
   case 1:
      arrayf = &(actnode->sol_increment);
   break;
   case 2:
      arrayf = &(actnode->sol_residual);
   break;
   case 3:
      arrayf = &(actnode->sol_mf);
   break;
   default:
      arrayf = NULL;
      dserror("Only 0,1,2,3 allowed for arrayfrom to select sol, sol_increment, sol_residual, sol_mf");
   }
   /*----------------------------------------- select correct arrayfrom */
   switch(arrayto)
   {
   case 0:
      arrayt = &(actnode->sol);
   break;
   case 1:
      arrayt = &(actnode->sol_increment);
   break;
   case 2:
      arrayt = &(actnode->sol_residual);
   break;
   case 3:
      arrayt = &(actnode->sol_mf);
   break;
   default:
      arrayt = NULL;
      dserror("Only 0,1,2,3 allowed for arrayfrom to select sol, sol_increment, sol_residual, sol_mf");
   }
   /* check the size of arrayf */
   if (from >= arrayf->fdim)
      dserror("Cannot copy from array, because place doesn't exist");
   /* check the size of arrayt */
   if (to >= arrayt->fdim)
   {
      diff = to - arrayt->fdim;
      /*max  = IMAX(diff,5);*/
      max = diff;
      amredef(arrayt,arrayt->fdim+max+1,arrayt->sdim,"DA");
   }
   for (j=0; j<arrayt->sdim; j++)
      arrayt->a.da[to][j] += arrayf->a.da[from][j]*fac;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of solserv_sol_add */



/*----------------------------------------------------------------------*
 |                                                            m.gee 4/03|
 |  in the nodes in sol, sol_increment or sol_residual                  |
 |  a local vector vec is assembled                                     |
 | (it should be initialized to zero before using solserv_sol_zero)     |
 |                                                                      |
 |                                                                      |
 | ELEMENT *actele   (i)   the element to assemble the vector for       |
 | DOUBLE  *localvec (i)   a vector of element size                     |
 | INT      arraynum (i)   indize of the array, 0 = sol                 |
 |                                              1 = sol_increment       |
 |                                              2 = sol_residual        |
 | INT      place      (i) row in array that is used                    |
 *----------------------------------------------------------------------*/
void solserv_sol_localassemble(INTRA *actintra, ELEMENT *actele, DOUBLE *localvec, INT arraynum,
                              INT place)
{
INT               i,j;
INT               diff,max,numdf;
ARRAY            *array;
NODE             *actnode;
#ifdef DEBUG
dstrc_enter("solserv_sol_localassemble");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<actele->numnp; i++)
{
   actnode = actele->node[i];
   if (actnode->proc != actintra->intra_rank) continue;
   numdf   = actnode->numdf;
   /*--------------------------------------------- select correct array */
   switch(arraynum)
   {
   case 0:
      array = &(actnode->sol);
   break;
   case 1:
      array = &(actnode->sol_increment);
   break;
   case 2:
      array = &(actnode->sol_residual);
   break;
   case 3:
      array = &(actnode->sol_mf);
   break;
   default:
      array = NULL;
      dserror("Only 0,1,2,3 allowed for arraynum to select sol, sol_increment, sol_residual, sol_mf");
   }
   /*-- check for correct size of array in the adjacent nodes to actele */
   if (place >= array->fdim)
   {
      diff = place - array->fdim;
      /*max  = IMAX(diff,5);*/
      max = diff;
      amredef(array,array->fdim+max+1,array->sdim,"DA");
   }
   for (j=0; j<numdf; j++)
      array->a.da[place][j] += localvec[i*numdf+j];
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of solserv_sol_localassemble */




/*----------------------------------------------------------------------*
 |                                                            m.gee 3/02|
 | dirichlet conditions are scaled by scale and written to sol in the   |
 | given place place                                                    |
 | This routine takes values from the structure actnode->gnode->dirich  |
 | and scales them by a given factor scale. Then it writes these values |
 | the ARRAY node->sol in the given place                               |
 | actnode->sol.a.da[place][j] = actnode->gnode->dirich_val.a.dv[j]*scale|
 | Nothing is done for dofs or nodes which do not have a dirichlet      |
 | condition                                                            |
 | FIELD *actfield (i) active field                                     |
 | INT    disnum   (i) indize of the discretization in actfield to be used|
 | INT    arraynum (i) indize of the array, 0 = sol                     |
 |                                          1 = sol_increment           |
 |                                          2 = sol_residual            |
 |                                          3 = sol_mf                  |
 | DOUBLE scale    (i) scaling factor for dirichlet condition           |
 | INT place       (i) row to put values in the ARRAY sol               |
 *----------------------------------------------------------------------*/
void solserv_putdirich_to_dof(FIELD *actfield,
                              INT    disnum,
                              INT    arraynum,
                              INT    place,
                              DOUBLE T          /* actual time		        */)
{
INT               i,j;
INT               diff,max;
ARRAY            *array;
NODE             *actnode;
DISCRET          *actdis;
DIRICH_CONDITION *dirich;

INT        actcurve;	             /* actual timecurve		*/
DOUBLE     timefac[MAXTIMECURVE];    /* factors from time-curve         */
DOUBLE     acttimefac;               /* actual factor from timecurve    */
DOUBLE     initval;	             /* intial dirichlet value	        */


#ifdef DEBUG
dstrc_enter("solserv_putdirich_to_dof");
#endif
/*----------------------------------------------------------------------*/
actdis = &(actfield->dis[disnum]);

/*------------------------------------------ get values from time curve */
for (actcurve=0;actcurve<numcurve;actcurve++)
  dyn_facfromcurve(actcurve,T,&timefac[actcurve]) ;

/*------------ loop nodes and put the result back to the node structure */
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   /*--------- do nothing if there is no dirichlet condition on actnode */
   if (actnode->gnode->dirich==NULL) continue;
   /*--------------------------------------------- select correct array */
   switch(arraynum)
   {
   case 0:
      array = &(actnode->sol);
   break;
   case 1:
      array = &(actnode->sol_increment);
   break;
   case 2:
      array = &(actnode->sol_residual);
   break;
   case 3:
      array = &(actnode->sol_mf);
   break;
   default:
      array = NULL;
      dserror("Only 0,1,2,3 allowed for arraynum to select sol, sol_increment, sol_residual, sol_mf");
   }
   /*------------------------------------------ get dirichlet condition */
   dirich = actnode->gnode->dirich;
   /* if the given place is outside the dimensions of ARRAY array enlarge it */
   if (place >= array->fdim)
   {
      diff = place - array->fdim;
      /*max  = IMAX(diff,5);*/
      max = diff;
      amredef(array,array->fdim+max+1,array->sdim,"DA");
   }
   /* put values dirich->dirich_val.a.dv[j] * scale to actnode->sol.a.dv[place] */
   for (j=0; j<actnode->numdf; j++)   /* for each degree of freedom */
   {
      if (dirich->dirich_onoff.a.iv[j]==0)  /*  if no Dirichlet Cond., continue */
         continue;
      actcurve = dirich->curve.a.iv[j]-1;

      if (actcurve<0)
         acttimefac = 1.;
      else if (actcurve >= 0 && actcurve < numcurve)
         acttimefac = timefac[actcurve];
      else
      {
        acttimefac=0.0;
        dserror("Solid Dirichlet BC: actual curve > number defined curves\n");
      }

      initval = dirich->dirich_val.a.dv[j];
      /*printf("actcurve=%i\n", actcurve);
        printf("initval=%f, acttimefac=%f, node=%i, dof=%i\n", initval, acttimefac, i, j);*/
      array->a.da[place][j] = initval * acttimefac;

#ifdef D_SSI
      /* special approach for ssi - problems mfirl 03/04 */
    if (actnode->gnode->ssicouple != NULL)
      if(actnode->gnode->ssicouple->ssi_couptyp == ssi_slave)
        array->a.da[place][j] = dirich->dirich_val.a.dv[j];
#endif
/* this is special for the ortiz example m.gee*/
/*      array->a.da[place][j] *= actnode->x[j];*/
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of solserv_putdirich_to_dof */


/*----------------------------------------------------------------------*
 |                                                            m.gee 3/02|
 | entries in sol in the place  placefrom1 and placefrom2 are added     |
 | node->sol[to][..] = node->sol[from1][..] * facfrom1 +                |
 |                     node->sol[from2][..] * facfrom2                  |
 |                                                                      |
 | This is ONLY performed for dofs which have a dirichlet condition     |
 |                                                                      |
 | FIELD *actfield (i) active field                                     |
 | INT    disnum   (i) indize of the discretization in actfield to be used|
 | INT    arraynum (i) indize of the array, 0 = sol                     |
 |                                          1 = sol_increment           |
 |                                          2 = sol_residual            |
 | INT    from1    (i) place in ARRAY sol to take the values from       |
 | INT    from2    (i) place in ARRAY sol to take the values from       |
 | INT    to       (i) place in ARRAY sol to write values to            |
 | DOUBLE facfrom1 (i) scaling factor vor values from from1             |
 | DOUBLE facfrom2 (i) scaling factor vor values from from2             |
 *----------------------------------------------------------------------*/
void solserv_adddirich(FIELD *actfield, INT disnum, INT arraynum,
                              INT from1,INT from2,INT to,
                              DOUBLE facfrom1, DOUBLE facfrom2)
{
INT               i,j;
INT               diff,max;
ARRAY            *array;
NODE             *actnode;
DISCRET          *actdis;
DIRICH_CONDITION *dirich;
#ifdef DEBUG
dstrc_enter("solserv_adddirich");
#endif
/*----------------------------------------------------------------------*/
actdis = &(actfield->dis[disnum]);
/*------------ loop nodes and put the result back to the node structure */
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   if (actnode->gnode->dirich==NULL) continue;
   /*--------------------------------------------- select correct array */
   switch(arraynum)
   {
   case 0:
      array = &(actnode->sol);
   break;
   case 1:
      array = &(actnode->sol_increment);
   break;
   case 2:
      array = &(actnode->sol_residual);
   break;
   case 3:
      array = &(actnode->sol_mf);
   break;
   default:
      array = NULL;
      dserror("Only 0,1,2,3 allowed for arraynum to select sol, sol_increment, sol_residual, sol_mf");
   }
   /* get dirichlet condition */
   dirich = actnode->gnode->dirich;
   /* check for sufficient size */
   max = IMAX(from1,from2);
   max = IMAX(max,to);
   if (max >= array->fdim)
   {
      diff = max - array->fdim;
      /*max  = IMAX(diff,5);*/
      max = diff;
      amredef(array,array->fdim+max+1,array->sdim,"DA");
   }
   for (j=0; j<actnode->numdf; j++)
   {
#ifdef LOCALSYSTEMS_ST
      if (actnode->locsysId == 0)
         if (dirich->dirich_onoff.a.iv[j]==0) continue;
#else
      if (dirich->dirich_onoff.a.iv[j]==0) continue;
#endif
      array->a.da[to][j] =
      array->a.da[from1][j]*facfrom1 +
      array->a.da[from2][j]*facfrom2;
   }

}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of solserv_adddirich */





/*----------------------------------------------------------------------*
 |                                                            m.gee 3/02|
 | entries in sol in the place  placefrom1 and placefrom2 are added     |
 | node->sol[to][..] += node->sol[from1][..] * facfrom1 +               |
 |                      node->sol[from2][..] * facfrom2                 |
 |                                                                      |
 | same functionality as solserv_adddirich but adds to array in the place to|
 |                                                                      |
 *----------------------------------------------------------------------*/
void solserv_assdirich_fac(FIELD *actfield, INT disnum, INT arraynum,
                           INT from1,INT from2,INT to,
                           DOUBLE facfrom1, DOUBLE facfrom2)
{
INT               i,j;
INT               diff,max;
ARRAY            *array;
NODE             *actnode;
DISCRET          *actdis;
DIRICH_CONDITION *dirich;
#ifdef DEBUG
dstrc_enter("solserv_assdirich_fac");
#endif
/*----------------------------------------------------------------------*/
actdis = &(actfield->dis[disnum]);
/*------------ loop nodes and put the result back to the node structure */
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   if (actnode->gnode->dirich==NULL) continue;
   /*--------------------------------------------- select correct array */
   switch(arraynum)
   {
   case 0:
      array = &(actnode->sol);
   break;
   case 1:
      array = &(actnode->sol_increment);
   break;
   case 2:
      array = &(actnode->sol_residual);
   break;
   case 3:
      array = &(actnode->sol_mf);
   break;
   default:
      array = NULL;
      dserror("Only 0,1,2,3 allowed for arraynum to select sol, sol_increment, sol_residual, sol_mf");
   }
   /* get dirichlet condition */
   dirich = actnode->gnode->dirich;
   /* check the array size */
   max = IMAX(from1,from2);
   max = IMAX(max,to);
   if (max >= array->fdim)
   {
      diff = max - array->fdim;
      /*max  = IMAX(diff,5);*/
      max = diff;
      amredef(array,array->fdim+max+1,array->sdim,"DA");
   }
   for (j=0; j<actnode->numdf; j++)
   {
#ifdef LOCALSYSTEMS_ST
      if (actnode->locsysId == 0)
         if (dirich->dirich_onoff.a.iv[j]==0) continue;
#else
      if (dirich->dirich_onoff.a.iv[j]==0) continue;
#endif
      array->a.da[to][j] +=
      array->a.da[from1][j]*facfrom1 +
      array->a.da[from2][j]*facfrom2;
   }

}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of solserv_assdirich_fac */



/*----------------------------------------------------------------------*
 |                                                            m.gee 3/02|
 | entries in sol in the place  placefrom1  are copied to  place to     |
 | node->sol[to][..] = node->sol[from][..]                              |
 |                                                                      |
 | This is only performed for dofs which have a dirichlet condition on them |
 |                                                                      |
 | FIELD *actfield (i) active field                                     |
 | INT    disnum   (i) indize of the discretization in actfield to be used|
 | INT    arraynum (i) indize of the array, 0 = sol                     |
 |                                          1 = sol_increment           |
 |                                          2 = sol_residual            |
 | INT    to       (i) place in ARRAY sol to write values to            |
 | INT    from     (i) place in ARRAY sol to take the values from       |
 *----------------------------------------------------------------------*/
void solserv_cpdirich(FIELD *actfield, INT disnum, INT arraynum,
                      INT from,INT to)
{
INT               i,j;
INT               diff,max;
ARRAY            *array;
NODE             *actnode;
DISCRET          *actdis;
DIRICH_CONDITION *dirich;
#ifdef DEBUG
dstrc_enter("solserv_cpdirich");
#endif
/*----------------------------------------------------------------------*/
actdis = &(actfield->dis[disnum]);
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   if (actnode->gnode->dirich==NULL) continue;
   /*--------------------------------------------- select correct array */
   switch(arraynum)
   {
   case 0:
      array = &(actnode->sol);
   break;
   case 1:
      array = &(actnode->sol_increment);
   break;
   case 2:
      array = &(actnode->sol_residual);
   break;
   case 3:
      array = &(actnode->sol_mf);
   break;
   default:
      array = NULL;
      dserror("Only 0,1,2,3 allowed for arraynum to select sol, sol_increment, sol_residual, sol_mf");
   }
   /* get dirichlet condition */
   dirich = actnode->gnode->dirich;
   /* check dimensions */
   max = IMAX(from,to);
   if (max >= array->fdim)
   {
      diff = max - array->fdim;
      /*max  = IMAX(diff,5);*/
      max = diff;
      amredef(array,array->fdim+max+1,array->sdim,"DA");
   }
   for (j=0; j<actnode->numdf; j++)
   {
#ifdef LOCALSYSTEMS_ST
      if (actnode->locsysId == 0)
         if (dirich->dirich_onoff.a.iv[j]==0) continue;
#else
      if (dirich->dirich_onoff.a.iv[j]==0) continue;
#endif
      array->a.da[to][j] = array->a.da[from][j];
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of solserv_cpdirich */



/*----------------------------------------------------------------------*
 |                                                            m.gee 3/02|
 | init sol[place] to zero                                              |
 |                                                                      |
 | This is only performed for dofs which have a dirichlet condition on them |
 |                                                                      |
 | FIELD *actfield (i) active field                                     |
 | INT    disnum   (i) indize of the discretization in actfield to be used|
 | INT    arraynum (i) indize of the array, 0 = sol                     |
 |                                          1 = sol_increment           |
 |                                          2 = sol_residual            |
 | INT    place    (i) row in ARRAY sol to be set to zero               |
 *----------------------------------------------------------------------*/
void solserv_zerodirich(FIELD *actfield, INT disnum, INT arraynum, INT place)
{
INT               i,j;
INT               diff,max;
ARRAY            *array;
NODE             *actnode;
DISCRET          *actdis;
DIRICH_CONDITION *dirich;
#ifdef DEBUG
dstrc_enter("solserv_zerodirich");
#endif
/*----------------------------------------------------------------------*/
actdis = &(actfield->dis[disnum]);
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   if (actnode->gnode->dirich==NULL) continue;
   /*--------------------------------------------- select correct array */
   switch(arraynum)
   {
   case 0:
      array = &(actnode->sol);
   break;
   case 1:
      array = &(actnode->sol_increment);
   break;
   case 2:
      array = &(actnode->sol_residual);
   break;
   case 3:
      array = &(actnode->sol_mf);
   break;
   default:
      array = NULL;
      dserror("Only 0,1,2,3 allowed for arraynum to select sol, sol_increment, sol_residual, sol_mf");
   }
   /* get dirichlet condition */
   dirich = actnode->gnode->dirich;
   /* check the size of array */
   if (place >= array->fdim)
   {
      diff = place - array->fdim;
      /*max  = IMAX(diff,5);*/
      max = diff;
      amredef(array,array->fdim+max+1,array->sdim,"DA");
   }
   for (j=0; j<actnode->numdf; j++)
   {
#ifdef LOCALSYSTEMS_ST
      if (actnode->locsysId == 0)
         if (dirich->dirich_onoff.a.iv[j]==0) continue;
#else
      if (dirich->dirich_onoff.a.iv[j]==0) continue;
#endif
      array->a.da[place][j] = 0.0;
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of solserv_zerodirich */


