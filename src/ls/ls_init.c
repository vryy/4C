/*!----------------------------------------------------------------------
\file
\brief ls_init.c

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_LS
#include "../headers/standardtypes.h"
#include "ls_prototypes.h"
/*!
\addtogroup LEVELSET
*//*! @{ (documentation module open)*/



/*!----------------------------------------------------------------------
\brief initialize level set problem

<pre>                                                            irhan 05/04
initialize level set problem
</pre>

*----------------------------------------------------------------------*/
void ls_init(
  FIELD          *actfield,  
  LS_DYNAMIC     *lsdyn,
  INT             numr
  )      
{
  INT i;
  INT flag;
  NODE  *actnode;
  
#ifdef DEBUG 
  dstrc_enter("ls_init");
#endif
/*----------------------------------------------------------------------*/
  
  /* set control variables for element evaluation */
  lsdyn->nif = 1;			 /* evaluation of      time rhs */
  lsdyn->nii = 1;			 /* evaluation of iteration rhs */
  lsdyn->step = 0;
  /* allocate space for solution history */
  for (i=0; i<actfield->dis[0].numnp; i++)
  {
    actnode = &(actfield->dis[0].node[i]);
    dsassert(actnode->numdf==1,"inconsistency from input regarding numdf!\n");
    amredef(&(actnode->sol_increment),numr,actnode->numdf,"DA");
    amzero(&(actnode->sol_increment));
  }
  /* construct the initial profile */
  if (lsdyn->init>=1)
  {
    if (lsdyn->init==1)              /* initial data from ls_start.data */
      dserror("reading from a file not implemented!\n");
    else if (lsdyn->init==2)                /* generate initial profile */
    {
      /* initialize circles */
      if (lsdyn->lsdata->numcirc!=0)
      {
        if (lsdyn->lsdata->numcirc==1)
        {
          if (lsdyn->lsdata->is_sharp==0)
            ls_init_single_circle(actfield,lsdyn);
          else if (lsdyn->lsdata->is_sharp==1)
            ls_init_single_circle_sharp(actfield,lsdyn);
        }
        else if (lsdyn->lsdata->numcirc==2)
        {
          if (lsdyn->lsdata->is_sharp == 0)
            ls_init_double_circle(actfield,lsdyn);
          else if (lsdyn->lsdata->is_sharp == 1)
            ls_init_double_circle_sharp(actfield,lsdyn);
        }
        else if (lsdyn->lsdata->numcirc==3)
        {
          if (lsdyn->lsdata->is_sharp == 0)
            ls_init_triple_circle(actfield,lsdyn);
          else if (lsdyn->lsdata->is_sharp == 1)
            ls_init_triple_circle_sharp(actfield,lsdyn);
        }
      }
      
      /* initialize lines */
      if (lsdyn->lsdata->numline!=0)
      {
        if (lsdyn->lsdata->numline==1)
        {
          if (lsdyn->lsdata->is_sharp==0)
            ls_init_single_line(actfield,lsdyn);
          else if (lsdyn->lsdata->is_sharp==1)
            ls_init_single_line_sharp(actfield,lsdyn);
        }
        else if (lsdyn->lsdata->numline==2)
        {
          if (lsdyn->lsdata->is_sharp == 0)
            ls_init_double_line(actfield,lsdyn);
          else if (lsdyn->lsdata->is_sharp == 1)
            ls_init_double_line_sharp(actfield,lsdyn);
        }
        else if (lsdyn->lsdata->numline==3)
        {
          if (lsdyn->lsdata->is_sharp == 0)
            ls_init_triple_line(actfield,lsdyn);
          else if (lsdyn->lsdata->is_sharp == 1)
            ls_init_triple_line_sharp(actfield,lsdyn);
        }
      }
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_init */



/*!----------------------------------------------------------------------
\brief construct (smooth) level set profile for one circular interface 

<pre>                                                            irhan 05/04
construct (smooth) level set profile for one circular interface 
</pre>

*----------------------------------------------------------------------*/
void ls_init_single_circle(
  FIELD* field,
  LS_DYNAMIC* lsdyn
  )
{
  INT i;
  DOUBLE x1,x2,rad1,xc1,yc1,d,d2;
  NODE* actnode;
  
#ifdef DEBUG 
  dstrc_enter("ls_init_single_circle");
#endif
/*----------------------------------------------------------------------*/
  
  /* set the radious of the circle */
  rad1 = lsdyn->lsdata->rad1;
  /* set the center of the circle */  
  xc1 = lsdyn->lsdata->xc1;
  yc1 = lsdyn->lsdata->yc1;
  /* loop over the nodes */
  for (i=0; i<field->dis[0].numnp; i++)
  {
	actnode = &(field->dis[0].node[i]);
	/* get the coordinates of the node */
	x1 = actnode->x[0];
	x2 = actnode->x[1];
	/* calculate shortest distance from node to circle */
	d2 = (x1-xc1)*(x1-xc1) + (x2-yc1)*(x2-yc1);
	d  = sqrt(d2) - rad1;
	/* set the initial set for the node */
	actnode->sol_increment.a.da[0][0] = d;
	actnode->sol_increment.a.da[1][0] = d;			
	actnode->sol.a.da[0][0] = d;
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_init_single_circle */



/*!----------------------------------------------------------------------
\brief construct (smooth) level set profile for two circular interfaces

<pre>                                                            irhan 05/04
construct (smooth) level set profile for two circular interfaces
</pre>

*----------------------------------------------------------------------*/
void ls_init_double_circle(
  FIELD* field,
  LS_DYNAMIC* lsdyn
  )
{
  INT        i;
  DOUBLE     x1,x2,rad1,xc1,yc1,rad2,xc2,yc2,d1,d12,d2,d22,dmin;
  NODE      *actnode;
  
#ifdef DEBUG 
  dstrc_enter("ls_init_double_circle");
#endif
/*----------------------------------------------------------------------*/

  /* set the radious of the first circle  */
  rad1 = lsdyn->lsdata->rad1;
  /* set the center  of the first circle  */  
  xc1  = lsdyn->lsdata->xc1;
  yc1  = lsdyn->lsdata->yc1;

  /* set the radious of the second circle */
  rad2 = lsdyn->lsdata->rad2;
  /* set the center  of the second circle */  
  xc2  = lsdyn->lsdata->xc2;
  yc2  = lsdyn->lsdata->yc2;
  /* loop over the nodes */
  for (i=0; i<field->dis[0].numnp; i++)
  {
    actnode = &(field->dis[0].node[i]);
    /* get the coordinates of the node */
    x1 = actnode->x[0];
    x2 = actnode->x[1];
    /* calculate shortest distance from node to circles */
    d12 = (x1-xc1)*(x1-xc1) + (x2-yc1)*(x2-yc1);
    d1  = sqrt(d12) - rad1;
    d22 = (x1-xc2)*(x1-xc2) + (x2-yc2)*(x2-yc2);
    d2  = sqrt(d22) - rad2;
    dmin = d1;
    if (fabs(d2)<fabs(dmin)) dmin = d2;
    
    /* set the initial set for the node */
    actnode->sol_increment.a.da[0][0] = dmin;
    actnode->sol_increment.a.da[1][0] = dmin;			
    actnode->sol.a.da[0][0] = dmin;
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif

  return;
} /* end of ls_init_double_circle */



/*!----------------------------------------------------------------------
\brief construct (smooth) level set profile for three circular interfaces

<pre>                                                            irhan 05/04
construct (smooth) level set profile for three circular interfaces
</pre>

*----------------------------------------------------------------------*/
void ls_init_triple_circle(
  FIELD* field,
  LS_DYNAMIC* lsdyn
  )
{
  INT        i;
  DOUBLE     x1,x2,rad1,xc1,yc1,rad2,xc2,yc2,rad3,xc3,yc3;
  DOUBLE     d1,d12,d2,d22,d3,d32,dmin;
  NODE      *actnode;
  
#ifdef DEBUG 
  dstrc_enter("ls_init_triple_circle");
#endif
/*----------------------------------------------------------------------*/

  /* set the radious of the  first circle  */
  rad1 = lsdyn->lsdata->rad1;
  /* set the center  of the  first circle  */  
  xc1  = lsdyn->lsdata->xc1;
  yc1  = lsdyn->lsdata->yc1;

  /* set the radious of the second circle */
  rad2 = lsdyn->lsdata->rad2;
  /* set the center  of the second circle */  
  xc2  = lsdyn->lsdata->xc2;
  yc2  = lsdyn->lsdata->yc2;
  /* set the radious of the third  circle */
  rad3 = lsdyn->lsdata->rad3;
  /* set the center  of the third  circle */  
  xc3  = lsdyn->lsdata->xc3;
  yc3  = lsdyn->lsdata->yc3;
  /* loop over the nodes */
  for (i=0; i<field->dis[0].numnp; i++)
  {
    actnode = &(field->dis[0].node[i]);
    /* get the coordinates of the node */
    x1 = actnode->x[0];
    x2 = actnode->x[1];
    /* calculate shortest distance from node to circles */
    d12 = (x1-xc1)*(x1-xc1) + (x2-yc1)*(x2-yc1);
    d1  = sqrt(d12) - rad1;
    d22 = (x1-xc2)*(x1-xc2) + (x2-yc2)*(x2-yc2);
    d2  = sqrt(d22) - rad2;
    d32 = (x1-xc3)*(x1-xc3) + (x2-yc3)*(x2-yc3);
    d3  = sqrt(d32) - rad3;
    dmin = d1;
    if (fabs(d2)<fabs(dmin)) dmin = d2;
    if (fabs(d3)<fabs(dmin)) dmin = d3;	
    
    /* set the initial set for the node */
    actnode->sol_increment.a.da[0][0] = dmin;
    actnode->sol_increment.a.da[1][0] = dmin;			
    actnode->sol.a.da[0][0] = dmin;
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_init_triple_circle */



/*!----------------------------------------------------------------------
\brief construct (sharp) level set profile for one circular interface

<pre>                                                            irhan 05/04
construct (sharp) level set profile for one circular interface
</pre>

*----------------------------------------------------------------------*/
void ls_init_single_circle_sharp(
  FIELD *field,
  LS_DYNAMIC *lsdyn
  )	
{
  INT        i;
  DOUBLE     x1,x2,rad1,xc1,yc1,d,d2;
  DOUBLE     val;
  NODE      *actnode;
  
#ifdef DEBUG 
  dstrc_enter("ls_init_single_circle_sharp");
#endif
/*----------------------------------------------------------------------*/

  /* set the radious of the circle */
  rad1 = lsdyn->lsdata->rad1;
  /* set the center of the circle */  
  xc1 = lsdyn->lsdata->xc1;
  yc1 = lsdyn->lsdata->yc1;
  /* loop over the nodes */
  for (i=0; i<field->dis[0].numnp; i++)
  {
    actnode = &(field->dis[0].node[i]);
    /* get the coordinates of the node */
    x1 = actnode->x[0];
    x2 = actnode->x[1];
    /* calculate shortest distance from node to circle */
    d2 = (x1-xc1)*(x1-xc1) + (x2-yc1)*(x2-yc1);
    d  = sqrt(d2) - rad1;
    val = -1.0;
    if (d>0.0) val = 1.0;
    
    /* set the initial set for the node */
    actnode->sol_increment.a.da[0][0] = val;
    actnode->sol_increment.a.da[1][0] = val;			
    actnode->sol.a.da[0][0] = val;
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_init_single_circle_sharp */



/*!----------------------------------------------------------------------
\brief construct (sharp) level set profile for two circular interfaces

<pre>                                                            irhan 05/04
construct (sharp) level set profile for two circular interfaces
</pre>

*----------------------------------------------------------------------*/
void ls_init_double_circle_sharp(
  FIELD *field,
  LS_DYNAMIC *lsdyn
  )
{
  INT        i;
  DOUBLE     x1,x2,rad1,xc1,yc1,rad2,xc2,yc2,d1,d12,d2,d22,dmin;
  DOUBLE     val;
  NODE      *actnode;
  
#ifdef DEBUG 
  dstrc_enter("ls_init_double_circle_sharp");
#endif
/*----------------------------------------------------------------------*/

  /* set the radious of the first circle  */
  rad1 = lsdyn->lsdata->rad1;
  /* set the center  of the first circle  */  
  xc1  = lsdyn->lsdata->xc1;
  yc1  = lsdyn->lsdata->yc1;
  /* set the radious of the second circle */
  rad2 = lsdyn->lsdata->rad2;
  /* set the center  of the second circle */  
  xc2  = lsdyn->lsdata->xc2;
  yc2  = lsdyn->lsdata->yc2;
  /* loop over the nodes */
  for (i=0; i<field->dis[0].numnp; i++)
  {
    actnode = &(field->dis[0].node[i]);
    /* get the coordinates of the node */
    x1 = actnode->x[0];
    x2 = actnode->x[1];
    /* calculate shortest distance from node to circles */
    d12 = (x1-xc1)*(x1-xc1) + (x2-yc1)*(x2-yc1);
    d1  = sqrt(d12) - rad1;
    d22 = (x1-xc2)*(x1-xc2) + (x2-yc2)*(x2-yc2);
    d2  = sqrt(d22) - rad2;
    dmin = d1;
    if (fabs(d2)<fabs(dmin)) dmin = d2;
    val = -1.0;
    if (dmin>0.0) val = 1.0;
    
    /* set the initial set for the node */
    actnode->sol_increment.a.da[0][0] = val;
    actnode->sol_increment.a.da[1][0] = val;
    actnode->sol.a.da[0][0] = val;
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_init_double_circle_sharp */



/*!----------------------------------------------------------------------
\brief construct (sharp) level set profile for three circular interfaces

<pre>                                                            irhan 05/04
construct (sharp) level set profile for three circular interfaces
</pre>

*----------------------------------------------------------------------*/
void ls_init_triple_circle_sharp(
  FIELD *field,
  LS_DYNAMIC *lsdyn
  )
{
  INT        i;
  DOUBLE     x1,x2,rad1,xc1,yc1,rad2,xc2,yc2,rad3,xc3,yc3;
  DOUBLE     d1,d12,d2,d22,d3,d32,dmin;
  DOUBLE     val;
  NODE      *actnode;
  
#ifdef DEBUG 
  dstrc_enter("ls_init_triple_circle_sharp");
#endif
/*----------------------------------------------------------------------*/

  /* set the radious of the  first circle */
  rad1 = lsdyn->lsdata->rad1;
  /* set the center  of the  first circle */  
  xc1  = lsdyn->lsdata->xc1;
  yc1  = lsdyn->lsdata->yc1;

  /* set the radious of the second circle */
  rad2 = lsdyn->lsdata->rad2;
  /* set the center  of the second circle */  
  xc2  = lsdyn->lsdata->xc2;
  yc2  = lsdyn->lsdata->yc2;

  /* set the radious of the third  circle */
  rad3 = lsdyn->lsdata->rad3;
  /* set the center  of the third  circle */  
  xc3  = lsdyn->lsdata->xc3;
  yc3  = lsdyn->lsdata->yc3;
  /* loop over the nodes */
  for (i=0; i<field->dis[0].numnp; i++)
  {
    actnode = &(field->dis[0].node[i]);
    /* get the coordinates of the node */
    x1 = actnode->x[0];
    x2 = actnode->x[1];
    /* calculate shortest distance from node to circles */
    d12 = (x1-xc1)*(x1-xc1) + (x2-yc1)*(x2-yc1);
    d1  = sqrt(d12) - rad1;
    d22 = (x1-xc2)*(x1-xc2) + (x2-yc2)*(x2-yc2);
    d2  = sqrt(d22) - rad2;
    d32 = (x1-xc3)*(x1-xc3) + (x2-yc3)*(x2-yc3);
    d3  = sqrt(d32) - rad3;
    dmin = d1;
    if (fabs(d2)<fabs(dmin)) dmin = d2;
    if (fabs(d3)<fabs(dmin)) dmin = d3;	
    val = -1.0;
    if (dmin>0.0) val = 1.0;
    
    /* set the initial set for the node */
    actnode->sol_increment.a.da[0][0] = val;
    actnode->sol_increment.a.da[1][0] = val;
    actnode->sol.a.da[0][0] = val;
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
return;
} /* end of ls_init_triple_circle_sharp */



/*!----------------------------------------------------------------------
\brief construct (smooth) level set profile for one line interface

<pre>                                                            irhan 05/04
construct (smooth) level set profile for one line interface
</pre>

*----------------------------------------------------------------------*/
void ls_init_single_line(
  FIELD* field,
  LS_DYNAMIC* lsdyn
  )
{
  INT        i;
  DOUBLE     xs1,ys1,xe1,ye1;
  DOUBLE     xes,yes,xxs,yys;
  DOUBLE     x1,x2;
  DOUBLE     nrs,ns;
  DOUBLE     d,d_old;
  NODE      *actnode;
  
#ifdef DEBUG 
  dstrc_enter("ls_init_single_line");
#endif
/*----------------------------------------------------------------------*/

  /* set start and end points of line */  
  xs1 = lsdyn->lsdata->xs1;
  ys1 = lsdyn->lsdata->ys1;
  xe1 = lsdyn->lsdata->xe1;
  ye1 = lsdyn->lsdata->ye1;

  /* loop over the nodes */
  for (i=0; i<field->dis[0].numnp; i++)
  {
    actnode = &(field->dis[0].node[i]);
    /* get the coordinates of the node */
    x1 = actnode->x[0];
    x2 = actnode->x[1];
    /* calculate shortest distance from node to line */
    xes = xe1 - xs1;
    yes = ye1 - ys1;
    xxs = x1 - xs1;
    yys = x2 - ys1;
    /* compute norms */
    nrs = yys*xes - xxs*yes;
    ns = sqrt(xes*xes + yes*yes);
    d  = nrs/ns;

    /* check whether there exists any lines */
    if (lsdyn->lsdata->numcirc==0)
    {
      actnode->sol_increment.a.da[0][0] = d;
      actnode->sol_increment.a.da[1][0] = d;			
      actnode->sol.a.da[0][0] = d;      
    }
    else
    {
      if (d<0.0)
      {
        actnode->sol_increment.a.da[0][0] = d;
        actnode->sol_increment.a.da[1][0] = d;			
        actnode->sol.a.da[0][0] = d;
      }
      else if (d>=0.0)
      {
        d_old = actnode->sol_increment.a.da[0][0];
        if (d<d_old)
        {
          actnode->sol_increment.a.da[0][0] = d;
          actnode->sol_increment.a.da[1][0] = d;			
          actnode->sol.a.da[0][0] = d;
        }
      }
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_init_single_line */



/*!----------------------------------------------------------------------
\brief construct (smooth) level set profile for two line interfaces

<pre>                                                            irhan 05/04
construct (smooth) level set profile for two line interfaces
</pre>

*----------------------------------------------------------------------*/
void ls_init_double_line(
  FIELD* field,
  LS_DYNAMIC* lsdyn
  )
{
  INT        i;
  DOUBLE     xo,yo,x,y,xxo,yyo;
  DOUBLE     max;
  NODE      *actnode;
  
#ifdef DEBUG 
  dstrc_enter("ls_init_double_line");
#endif
/*----------------------------------------------------------------------*/

  /* set start and end points of the first line */  
  xo = lsdyn->lsdata->xs1;
  yo = lsdyn->lsdata->ys2;

  
  /* loop over the nodes */
  for (i=0; i<field->dis[0].numnp; i++)
  {
    actnode = &(field->dis[0].node[i]);
    /* get the coordinates of the node */
    x = actnode->x[0];
    y = actnode->x[1];

    /*
     * compute compoonents of the position vector with
     * respect to shifted coordinate system
     */
    xxo = x - xo;
    yyo = y - yo;

    if (xxo>=0)
    {
      if (yyo<=0)
      {
        actnode->sol_increment.a.da[0][0] = xxo;
        actnode->sol_increment.a.da[1][0] = xxo;			
        actnode->sol.a.da[0][0] = xxo;      
      }
      else
      {
        actnode->sol_increment.a.da[0][0] = sqrt(xxo*xxo+yyo*yyo);
        actnode->sol_increment.a.da[1][0] = sqrt(xxo*xxo+yyo*yyo);			
        actnode->sol.a.da[0][0] = sqrt(xxo*xxo+yyo*yyo);      
      }
    }
    else
    {
      if (yyo<=0)
      {
        if (yyo>xxo)
        {
          max = yyo;
        }
        else
        {
          max = xxo;          
        }
        actnode->sol_increment.a.da[0][0] = max;
        actnode->sol_increment.a.da[1][0] = max;			
        actnode->sol.a.da[0][0] = max;      
      }
      else
      {
        actnode->sol_increment.a.da[0][0] = yyo;
        actnode->sol_increment.a.da[1][0] = yyo;			
        actnode->sol.a.da[0][0] = yyo;      
      }
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_init_double_line */




/*!----------------------------------------------------------------------
\brief construct (smooth) level set profile for three line interfaces

<pre>                                                            irhan 05/04
construct (smooth) level set profile for three line interfaces
</pre>

*----------------------------------------------------------------------*/
void ls_init_triple_line(
  FIELD* field,
  LS_DYNAMIC* lsdyn
  )
{
  
#ifdef DEBUG 
  dstrc_enter("ls_init_triple_line");
#endif
/*----------------------------------------------------------------------*/

  /* do nothing for a while */
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_init_triple_line */



/*!----------------------------------------------------------------------
\brief construct (sharp) level set profile for one line interface

<pre>                                                            irhan 05/04
construct (sharp) level set profile for one line interface
</pre>

*----------------------------------------------------------------------*/
void ls_init_single_line_sharp(
  FIELD *field,
  LS_DYNAMIC *lsdyn
  )	
{
  
#ifdef DEBUG 
  dstrc_enter("ls_init_single_line_sharp");
#endif
/*----------------------------------------------------------------------*/

  /* do nothing for a while */
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_init_single_line_sharp */



/*!----------------------------------------------------------------------
\brief construct (sharp) level set profile for two line interfaces

<pre>                                                            irhan 05/04
construct (sharp) level set profile for two line interfaces
</pre>

*----------------------------------------------------------------------*/
void ls_init_double_line_sharp(
  FIELD *field,
  LS_DYNAMIC *lsdyn
  )
{
  
#ifdef DEBUG 
  dstrc_enter("ls_init_double_line_sharp");
#endif
/*----------------------------------------------------------------------*/

  /* do nothing for a while */
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_init_double_line_sharp */



/*!----------------------------------------------------------------------
\brief construct (sharp) level set profile for three line interfaces

<pre>                                                            irhan 05/04
construct (sharp) level set profile for three line interfaces
</pre>

*----------------------------------------------------------------------*/
void ls_init_triple_line_sharp(
  FIELD *field,
  LS_DYNAMIC *lsdyn
  )
{
  
#ifdef DEBUG 
  dstrc_enter("ls_init_triple_line_sharp");
#endif
/*----------------------------------------------------------------------*/

  /* do nothing for a while */
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
return;
} /* end of ls_init_triple_line_sharp */
/*! @} (documentation module close)*/
#endif
