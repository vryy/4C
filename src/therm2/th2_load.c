/*======================================================================*/
/*!
\file
\brief Spatial integration of loads (ie heat sources/fluxes) applied
       to element domain and its boundary.

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
*/
#ifdef D_THERM2

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm2.h"

/*!
\addtogroup THERM2
*//*! @{ (Documentation module open)*/


/*----------------------------------------------------------------------*/
/*!
\brief General problem data

global variable GENPROB genprob is defined in global_control.c

\author bborn
\date 03/06
*/
extern struct _GENPROB genprob;


/*----------------------------------------------------------------------*/
/*!
\brief Static variables: Element load vector, form functions,
derivatives of form functions, Jacobian matrix

Due to the storage class "static", these variables are accessible to all 
functions defined in this file.

\author bborn
\date 03/06
*/
static INT allocated = 0;  /* allocation falg */
static ARRAY eload_a; static DOUBLE **eload;  /* static element load vector */
static ARRAY funct_a; static DOUBLE  *funct;  /* ansatz-functions           */
static ARRAY deriv_a; static DOUBLE **deriv;  /* derivatives of ansatz-funct*/
static ARRAY xjm_a;   static DOUBLE **xjm;    /* jacobian matrix            */

/*======================================================================*/
/*!
\brief Initialise load routines

\author bborn
\date 03/06
*/
void th2_load_init()
{
#ifdef DEBUG
  dstrc_enter("th2_load_init");
#endif

  if (allocated == 0)
  {
    eload = amdef("eload", &eload_a, NUMDOF_THERM2, MAXNOD_THERM2, "DA");
    funct = amdef("funct", &funct_a, MAXNOD_THERM2, 1, "DV");
    deriv = amdef("deriv", &deriv_a, NDIM_THERM2, MAXNOD_THERM2, "DA");
    xjm = amdef("xjm_a", &xjm_a, NDIM_THERM2, NDIM_THERM2, "DA");
    /* set allocated */
    allocated = 1;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th2_load_init() */


/*======================================================================*/
/*!
\brief Finalise load routines

Deallocate previously allocated ARRAYs

\author bborn
\date 03/06
*/
void th2_load_final()
{
#ifdef DEBUG
  dstrc_enter("th2_load_final");
#endif

  if (allocated == 1)
  {
    if (eload != NULL)
    {
      amdel(&eload_a);
    }
    if (funct != NULL)
    {
      amdel(&funct_a);
    }
    if (deriv != NULL)
    {
      amdel(&deriv_a);
    }
    if (xjm != NULL)
    {
      amdel(&xjm_a);
    }
    /* reste allocation flag */
    allocated = 0;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th2_load_quit() */


/*======================================================================*/
/*!
\brief Spatial integration of heat source and boundary heat flux

The heat source is integration over the element domain (surface).
The heat flux is integratted along the elements boundaries (edges/lines).
The integration results in the external element heat load vector.

The parameter space axes are called r and s, respectively.

Quadrilateral elements can have either 4, 8 or 9 nodes, which are 
ordered as follows
          s                     s                          s
          ^                     ^                          ^                 
          |                     |                          |                
 (-1,+1)--+--(+1,+1)    1-------4--------0        1------line 0-----0       
    |     |     |       |       |        |        |        |        |       
    | parameter |       |       |        |        |        |        |       
    | space     |       |quad8/9|        |        l quad4  |        l       
    |     |     |       |       |        |        i        |        i       
    +---(0,0)---+-->  --5------(8)-------7-->   --n--------+--------n-->    
    |     |     |  r    |       |        |  r     e        |        e  r    
    |     |     |       |       |        |        1        |        3       
    |     |     |       |       |        |        |        |        |       
    |     |     |       |       |        |        |        |        |       
 (-1,-1)--+--(+1,-1)    2-------6--------3        2------line 2-----3       
                                |                          |                
                                                                            

Triangular elements can have either 3 or 6 elements, which are ordered
as follows 
    s                    s                  s
    ^                    ^                  ^
    |                    |                  |
  (0,1)                  2                  2
    |\                   |\                 |\
    | \                  | \                | \
    |  \                 l  l               |  \
    |   \                i   i              |   \
    |    \               n    n             5    4
    |     \              e     e            |     \
    | para \             2      1           |      \
    | space \            | tri3  \          | tri6  \
    |        \           |        \         |        \
 -(0,0)-----(1,0)--> r  -0--line0--1--> r  -0----3----1--> r
    |                    |                  |


Extra nautical service : compass
Useful for compass-based identification of element edges

            north 
    northwest | northeast
             \|/
        west--o--east
             /|\
    southwest | southeast
            south 


\author bborn
\date 03/06
*/
void th2_load_heat(ELEMENT *ele,  /* actual element */
                   TH2_DATA *data,
                   DOUBLE *loadvec, /* global element load vector fext */
                   INT imyrank)
{

  /* general variables */
  const INT heatminus = -1.0;  /* minus sign occuring in heat conduction */
  const INT numdf = NUMDOF_THERM2;  /* 1 DOF per node : temperature */
  INT i, j;  /* some counters */
  INT inode, idof;  /* some counters */
  INT nelenod;  /* number of element nodes */
  INT neledof;  /* element DOF */
  /* DOUBLE thickness; */
  DOUBLE det;  /* det of jacobian matrix */


  /* variables for domain integration */
  INT foundsurface;  /* flag for surfaceload present    */
  INT nir = 0;  /* number of GPs in r-direction */
  INT nis = 0;  /* number of GPs in s-direction */
  INT intc = 0;  /* "integration case" for tri-element */
  INT lr, ls;  /* integration directions */
  DOUBLE e1 = 0.0;  /* current GP r-coordinate */
  DOUBLE e2 = 0.0;  /* current GP s-coordinate */
  DOUBLE fac;  /* integration factor GP-info */
  DOUBLE facr = 0.0;  /* integration factor GP-info */
  DOUBLE facs = 0.0;  /* integration factor GP-info */


  /* variables for boundary integration */
  INT foundline;  /* flag for lineload present or not       */
  INT ngline;  /* number of geometrylines to the element */
  GLINE *gline[MAXEDG_THERM2];  /* geometry lines of the element */
  NEUM_CONDITION *lineneum[MAXEDG_THERM2]; /* short line-neum. conditions */
  INT line;  /* looper over lines                      */
  INT ngnode;  /* number of geometry-nodes on g-line     */
  INT nil;  /* number of GP's in for triangle  */
  DOUBLE glr[GLINTC_THERM2];  /* GP r-coords on edge */
  DOUBLE gls[GLINTC_THERM2];  /* GP s-coords on edge */
  DOUBLE glw[GLINTC_THERM2];  /* GP weight on edge */
  DOUBLE ds;  /* dx/dr line increment for line integration */
  DOUBLE facline;  /* integration factor for line integration */
  INT loadadd;  /* flag ==1 load is going to be added */
  DOUBLE forceline[NUMDOF_THERM2];  /* lineload value in th-direction (inp) */


  /*--------------------------------------------------------------------*/
  /*--------------------------------------------------------------------*/
  /* start */
  /* th2_load_init() has to be called before!!! */
#ifdef DEBUG
  dstrc_enter("th2_load_heat");
#endif
  /* initialize eload */
  amzero(&eload_a);
  /*------------------------------------------------------------------*/
  /* element node and DOF info */
  nelenod = ele->numnp;  /* element nodes */
  neledof = nelenod*numdf;  /* total number of element DOFs */
  /* thickness  = ele->e.th2->thick; */

  /*--------------------------------------------------------------------*/
  /*--------------------------------------------------------------------*/
  /* domain load ==> surface heat load, heat source */

  /*--------------------------------------------------------------------*/
  /* check for presence of surface loads */
  if (!(ele->g.gsurf->neum))
  {
     foundsurface = 0;
  }
  else
  {
    switch (ele->g.gsurf->neum->neum_type)
    {
      case pres_domain_load:
        foundsurface = 1;
        break;
      default:
        dserror("load case not implemented");
        foundsurface = 0;
        break;
    }
  }

  /*------------------------------------------------------------------- */
  /* integrate surface load */
  if ((foundsurface > 0) && (imyrank == ele->proc))
  {
    /*------------------------------------------------------------------*/
    /* get integraton data */ 
    switch (ele->distyp)
    {
      /* quadrilaterals */
      case quad4: case quad8: case quad9:
        nir = ele->e.th2->nGP[0];
        nis = ele->e.th2->nGP[1];
        break;
      /* triangles */
      case tri3: case tri6:
        nir = 1;
        nis = ele->e.th2->nGP[0];
        intc = ele->e.th2->gpintc;
        break;
      default:
        dserror("ele->distyp unknown! in 'w1_loadheat.c' ");
        break;
    }  /* end switch(ele->distyp) */
    /*------------------------------------------------------------------*/
    /* integration loop over GP's of actual surface */
    /* loop r-direction */
    for (lr=0; lr<nir; lr++)
    {
      /* loop s direction */
      for (ls=0; ls<nis; ls++)
      {
        /* get values of  shape functions and their derivatives */
        switch(ele->distyp)
        {
          /* --> quad - element */
          case quad4: case quad8: case quad9:
            e1   = data->gqlc[lr][nir-1];
            facr = data->gqlw[lr][nir-1];
            e2   = data->gqlc[ls][nis-1];
            facs = data->gqlw[ls][nis-1];
            break;
          /* --> tri - element */
          case tri3: case tri6:
            e1   = data->gtdcr[ls][intc];
            facr = 1.0;
            e2   = data->gtdcs[ls][intc];
            facs = data->gtdw[ls][intc];
            break;
          default:
            dserror("ele->distyp unknown!");
            break;
        }  /* end of switch(ele->distyp) */
        /*--------------------------------------------------------------*/
        /* shape functions (and (not needed)their derivatives) */
        th2_shape_deriv(funct, deriv, e1, e2, ele->distyp, 1);
        /*--------------------------------------------------------------*/
        /* Jacobian matrix */
        th2_jaco(deriv, xjm, &det, ele, nelenod);
        /*--------------------------------------------------------------*/
        /* integration factor  */
        fac = heatminus * facr * facs * det;
        /*------------------------------------------------------------- */
        /* surface-load  ==> eload modified */
        th2_load_heatsurf(ele, eload, funct, fac, nelenod);
      }  /* end of for (ls=0; ls<nis; ls++) */
    }  /* for (lr=0; lr<nir; lr++) */

    /*------------------------------------------------------------------*/
    /* WHAT ??????????????? */
    /* the surface load of this element has been done,
     * -> switch off the surface load condition */
    ele->g.gsurf->neum=NULL;  /* WHY ?????? */

  } /* end of if ((foundsurface > 0) && (imyrank==ele->proc)) */

  /*--------------------------------------------------------------------*/
  /*--------------------------------------------------------------------*/
  /* boundary load ==> heat flux on edge/line */
  /* integration of line loads on lines adjacent to this element */
  /*--------------------------------------------------------------------*/
  /* check generally for presence of line loads */
  foundline = 0;
  /* total number of lines/edges of this element */
  ngline = ele->g.gsurf->ngline;
  /* loop over lines, check for Neumann conditions on lines/edges */
  for (i=0; i<ngline; i++)
  {
    gline[i] = ele->g.gsurf->gline[i];
    lineneum[i] = gline[i]->neum;
    if (lineneum[i] != NULL)
    {
      switch (lineneum[i]->neum_type)
      {
        case neum_orthopressure:
          dserror("Load type is not available!");
          break;
        case neum_live:
          foundline = 1;
          break;
        default:
          lineneum[i] = NULL;
          foundline = 0;
          dserror("Load type is not available!");
          break;
      }
    }
  }
  /*--------------------------------------------------------------------*/
  /* integrate boundary load */
  /* at least one edge of element is subdued to a load */
  if (foundline > 0)
  {
    /* loop over lines (with Neumann conditions) */
    for (line=0; line<ngline; line++)
    {
      /* check wether current 'line' is subdued to a load */
      if (lineneum[line] != NULL)
      {
        /* check number of nodes on line */
        ngnode = gline[line]->ngnode;
        /*--------------------------------------------------------------*/
        /* get Gauss coords and weights for lines */
        /* coordinates and weights of integration points */
        switch (ele->distyp)
        {
          /*------------------------------------------------------------*/
          /* --> quad - element */
          case quad4: case quad8: case quad9:
            /* coordinates and weights of integration points */
            /* original GP-coordinates and weights for area-integration */
            /* degeneration to line-integration-info for r-s directions */
            switch (line)
            {
              /* north edge : s==+1 -- first,last,(quadratic:middle) node */
              case 0:
                nil = nir;  /* integration in r-direction */
                for (i=0; i<nil; i++)
                {
                  glr[i] = data->gqlc[i][nil-1];
                  gls[i] = +1.0;
                  glw[i] = data->gqlw[i][nil-1];
                }
                break;
              /* west edge : r==-1 */
              case 1:
                nil = nis;  /* integration along s-axis */
                for (i=0; i<nil; i++)
                {
                  glr[i] = -1.0;
                  gls[i] = data->gqlc[i][nil-1];
                  glw[i] = data->gqlw[i][nil-1];
                }
                break;
              /* south edge : s==-1 */
              case 2:
                nil = nir;  /* integration in r-direction */
                for (i=0; i<nil; i++)
                {
                  glr[i] = data->gqlc[i][nil-1];
                  gls[i] = -1.0;
                  glw[i] = data->gqlw[i][nil-1];
                }
                break;
              /* east edge : r==+1 */
              case 3:
                nil = nis;  /* integration along s-axis */
                for (i=0; i<nil; i++)
                {
                  glr[i] = +1.0;
                  gls[i] = data->gqlc[i][nil-1];
                  glw[i] = data->gqlw[i][nil-1];
                }
                break;
            }  /* end of switch (line) */
            break;
          /*------------------------------------------------------------*/
          /* --> tri-element */
          case tri3: case tri6:
            /* integration points on edge -- the same for all 3 edges*/
            nil = ele->e.th2->gpned;
            switch (line)
            {
              /* south edge : s==0 */
              case 0:
                for (i=0; i<nil; i++)
                {
                  glr[i] = data->gtlc[i][nil-1];
                  gls[i] = 0.0;
                  glw[i] = data->gtlw[i][nil-1];
                }
                break;
              /* northeast edge */
              /* quadrature along diagonal 0<=xi<=sqrt(2)
               * replaced by integration along s-axis 0<=s<=1
               * r-coords depend on r=1-s */
              case 1:
                for (i=0; i<nil; i++)
                {
                  gls[i] = data->gtlc[i][nil-1];
                  glr[i] = 1.0 - gls[i];
                  glw[i] = data->gtlw[i][nil-1];
                }
                break;
              /* west edge : r==0 */
              /* quadrature along s-axis, 0<=s<=1 */
              case 2:
                for (i=0; i<nil; i++)
                {
                  glr[i] = 0.0;
                  gls[i] = data->gtlc[i][nil-1];
                  glw[i] = data->gtlw[i][nil-1];
                }
                break;
            }  /* end of switch (line) */
            break;
          /*------------------------------------------------------------*/
          /* catch erroneous disc. types */
          default:
            dserror("Discretisation type is impossible!");
        }  /* end of switch (ele->distyp) */
        /*----------------------------------------------------------*/
        /* integration loop on actual line */
        for (i=0; i<nil; i++)
        {
          /* Gaussian point and weight at it */
          e1 = glr[i];
          e2 = gls[i];
          facr = glw[i];
          /* shape functions and their derivatives */
          th2_shape_deriv(funct, deriv, e1, e2, ele->distyp, 1);
          /* Jacobian matrix */
          th2_jaco(deriv, xjm, &det, ele, nelenod);
          /* line increment */
          ds = 0.0;
          switch (ele->distyp)
          {
            /* quadrilaterals */
            case quad4: case quad8: case quad9:
              switch (line)
              { 
                /* south and north edges */
                case 0: case 2:
                  ds = DSQR(xjm[0][0])+DSQR(xjm[0][1]);
                  ds = sqrt(ds);
                  break;
                /* west and east edges */
                case 1: case 3:
                  ds = DSQR(xjm[1][0])+DSQR(xjm[1][1]);
                  ds = sqrt(ds);
                  break;
              }  /* end of switch (line) */
              break;
            /* triangles */
            case tri3: case tri6:
              switch (line)
              {
                /* south edge */
                case 0:
                  ds = xjm[0][0]*xjm[0][0] + xjm[0][1]*xjm[0][1];
                  ds = sqrt(ds);
                  break;
                /* northeast edge */
                case 1:
                  ds = (-xjm[0][0]+xjm[1][0])*(-xjm[0][0]+xjm[1][0]) 
                    + (-xjm[0][1]+xjm[1][1])*(-xjm[0][1]+xjm[1][1]);
                  ds = sqrt(ds);
                  break;
                /* west edge */
                case 2:
                  ds = xjm[1][0]* xjm[1][0] + xjm[1][1]*xjm[1][1];
                  ds = sqrt(ds);
                  break;
              } /* end of switch(line) */
              break;
            /* catch erroneous discretisation types */
            default:
              dserror("Discretisation type is impossible!");
              break;
          }  /* end of switch (ele->distyp) */
          /* integration factor  */
          facline = heatminus * ds * facr;
          /*------------------------------------------------------*/
          /* differ load type */
          loadadd = 0;
          switch (gline[line]->neum->neum_type)
          {
            /* orthonormal pressure */
            case neum_orthopressure:
              dserror("Orthonormal pressure has not yet been implemented.");
              loadadd = 0;
              break;
            /* live -> uniform prescribed line load */
            case neum_live:
              loadadd = 1;
              for (i=0; i<numdf; i++)
              {
                forceline[i] = gline[line]->neum->neum_val.a.dv[i];
              }
              break;
            default:
              loadadd = 0;
              break;
          }  /* end of switch (gline[line]->neum->neum_type) */
          /* add load vector component to element load vector */
          if (loadadd == 1)
          {
            /* loop over nodes of element */
            for (j=0; j<nelenod; j++)
            {
              /* loop over DOFs at node */
              for (i=0; i<numdf; i++)
              {
                eload[i][j] += funct[j] * forceline[i] * facline;
              }
            }
          }  /* end if (loadadd == 1) */
        }  /* end of for (i=0; i<nil; i++) */
      }  /* end of if (lineneum[line] != NULL) */
      /*----------------------------------------------------------------*/
      /* WHY ? */
      /* What is the reason behind this move? */
      ele->g.gsurf->gline[line]->neum=NULL;  /* WHAT? */
    }  /* end of for (line=0; line<ngline; line++) */
  }  /* end of if (foundline > 0) */

  /*--------------------------------------------------------------------*/
  /*--------------------------------------------------------------------*/  
  /* add static array eload to global external element load vector */
  if (foundsurface+foundline != 0)
  {
    for (inode=0; inode<nelenod; inode++)
    {
      for (idof=0; idof<numdf; idof++)
      {
#ifdef PBB
        printf("ele %d load %d %d : %f\n", ele->Id, idof, inode, eload[idof][inode]);
#endif
        loadvec[inode*numdf+idof] += eload[idof][inode];
      }
    }
  }

  /*--------------------------------------------------------------------*/
  /*--------------------------------------------------------------------*/
  /* finish */
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th2_eleload */



/*======================================================================*/
/*!
\brief Determine load due to heat source in element domain (a surface)

\param	 *ele		 ELEMENT      (i)    actual element
\param  **eload          DOUBLE       (o)    element load vector
\param   *funct          DOUBLE       (i)    shape function
\param    fac            DOUBLE       (i)    integration factor
\param    iel            INT          (i)    number of element nodes
\return void

\author bborn
\date 03/06
*/
void th2_load_heatsurf(ELEMENT *ele,
                       DOUBLE **eload,
                       DOUBLE *funct,
                       DOUBLE fac,
                       INT iel)
{
  DOUBLE force[NUMDOF_THERM2];
  INT i, j;  /* loopers (i = loaddirection x or y)(j=node) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_load_heatsurf");
#endif
  /*--------------------------------------------------------------------*/
  /* distinguish load type */
  switch(ele->g.gsurf->neum->neum_type)
  {
    /*------------------------------------------------------------------*/
    /* uniform prescribed surface load */
    case pres_domain_load:
      for (i=0; i<NUMDOF_THERM2; i++)
      {
        force[i] = ele->g.gsurf->neum->neum_val.a.dv[i];
      }
      /* add load vector component to element load vector */
      for (j=0; j<iel; j++)
      {
        for (i=0; i<NUMDOF_THERM2; i++)
        {
          eload[i][j] += funct[j] * force[i] * fac;
        }
      }
      break;
    /*------------------------------------------------------------------*/
    case neum_live:
      dserror("load case unknown");
      break;
    /*------------------------------------------------------------------*/
    case neum_consthydro_z:
      dserror("load case unknown");
      break;
    /*------------------------------------------------------------------*/
    case neum_increhydro_z:
      dserror("load case unknown");
      break;
    /*------------------------------------------------------------------*/
    default:
      break;
  }  /* end of switch(ele->g.gsurf->neum->neum_type) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th2_load_heatsurf(...) */


/*======================================================================*/
#endif  /*end of #ifdef D_THERM2 */
/*! @} (documentation module close)*/
