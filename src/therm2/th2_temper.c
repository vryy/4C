/*======================================================================*/
/*!
\file
\brief contains the routine 'w1_funct_deriv' which calculates the shape
       functions and derivatives for a wall element

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 08/06
*/

/*----------------------------------------------------------------------*/
#ifdef D_THERM2

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm2.h"

/*----------------------------------------------------------------------*/
/*!
\brief General problem data

\author bborn
\date 03/06
*/
extern struct _GENPROB genprob;

/*----------------------------------------------------------------------*/
/*!
\brief Locally global allocatables

\author bborn
\date 03/06
*/
static INT      allocated = 0;  /* allocation flag */
static ARRAY    nodtem_a;
static DOUBLE  *nodtem;  /* current temperature at element nodes */
static ARRAY    shape_a;  /* shape functions */
static DOUBLE  *shape;
static ARRAY    deriv_a;  /* derivatives of shape functions */
static DOUBLE **deriv;


/*======================================================================*/
/*!
\brief Initialise locally globals

Allocate local globals

\author bborn
\date 03/06
*/
void th2_temper_init()
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_temper_init");
#endif

  /*--------------------------------------------------------------------*/
  /* allocation */
  if (allocated == 0)
  {
    nodtem = amdef("nodtem", &nodtem_a, MAXNOD_THERM2, 1, "DV");
    shape = amdef("shape", &shape_a, MAXNOD_THERM2, 1, "DV");
    deriv = amdef("deriv", &deriv_a, NDIM_THERM2, MAXNOD_THERM2, "DA");
    /* flag alloaction */
    allocated = 1;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th2_temper_init */


/*=====================================================================*/
/*!
\brief Finalise locally globals

Deallocate local globals

\author bborn
\date 03/06
*/
void th2_temper_final()
{

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_temper_final");
#endif

  /*--------------------------------------------------------------------*/
  /* deallocation */
  if (allocated == 1)
  {
    amdel(&nodtem_a);
    amdel(&shape_a);
    amdel(&deriv_a);
    /* reset allocation flag */
    allocated = 0;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th2_temper_final */


/*======================================================================*/
/*!
\brief Calculate temperature at point on element parameter domain

The element stiffness matrix, ie the tangent operator, is determined
for the linear planar heat conduction problem.

\param   *ele           ELEMENT     (i)   pointer to current element
\param    r             DOUBLE      (i)   r-coord of parameter (r,s)
\param    s             DOUBLE      (i)   s-coord of parameter (r,s)
\param   *tem           DOUBLE      (o)   temperature at (r,s)

\return void

\author bborn
\date 08/06
*/
void th2_temper_cal(ELEMENT *ele,
                    DOUBLE r,
                    DOUBLE s,
                    DOUBLE *tem)
{
  INT nelenod;
  INT neledof;
  INT k;  /* loop index */

  /*====================================================================*/
#ifdef DEBUG
  dstrc_enter("th2_temper_cal");
#endif

  /*--------------------------------------------------------------------*/
  /* dimensions */
  nelenod = ele->numnp;
  neledof = NUMDOF_THERM2 * nelenod;

  /*--------------------------------------------------------------------*/
  /* reset to zero nodal temper., shape funct. and derivatives */
  amzero(&nodtem_a);
  amzero(&shape_a);
  amzero(&deriv_a);

  /*--------------------------------------------------------------------*/
  /* get shape funtions at (r,s) */
  th2_shape_deriv(shape, deriv, r, s, ele->distyp, 0);

  /*--------------------------------------------------------------------*/
  /* current nodal temperature */
  for (k=0; k<nelenod; k++)
  {
    nodtem[k] = ele->node[k]->sol.a.da[0][0];
  }

  /*--------------------------------------------------------------------*/
  /* interpolate temperature */
  *tem = 0.0;
  for (k=0; k<nelenod; k++)
  {
    *tem = *tem + shape[k]*nodtem[k];
  }

  /*====================================================================*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th2_temper_cal */

#endif /* end of #ifdef D_THERM2 */
/*! @} (documentation module close)*/
