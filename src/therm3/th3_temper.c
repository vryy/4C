/*======================================================================*/
/*!
\file
\brief Calculate the temperature at given point in parameter space.
       This point is mostly a Gauss point of the BRICK1 element
       stiffness matrix / internal force vector.

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 10/06
*/

/*----------------------------------------------------------------------*/
#ifdef D_THERM3

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm3.h"

/*----------------------------------------------------------------------*/
/*!
\brief General problem data

\author bborn
\date 03/06
*/
extern GENPROB genprob;

/*----------------------------------------------------------------------*/
/*!
\brief Fields

vector of numfld FIELDs, defined in global_control.c

\author bborn
\date 10/06
*/
extern FIELD *field;

/*----------------------------------------------------------------------*/
/*!
\brief Locally global allocatables

\author bborn
\date 03/06
*/
static INT      allocated = 0;  /* allocation flag */
static ARRAY    nodtem_a;
static DOUBLE  *nodtem;  /* current temperature at element nodes */


/*======================================================================*/
/*!
\brief Initialise local globals

Allocate local globals

\author bborn
\date 10/06
*/
void th3_temper_init()
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_temper_init");
#endif

  /*--------------------------------------------------------------------*/
  /* allocation */
  if (allocated == 0)
  {
    nodtem = amdef("nodtem", &nodtem_a, MAXNOD_THERM3, 1, "DV");
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
\brief Finalise local globals

Deallocate local globals

\author bborn
\date 10/06
*/
void th3_temper_final()
{

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_temper_final");
#endif

  /*--------------------------------------------------------------------*/
  /* deallocation */
  if (allocated == 1)
  {
    amdel(&nodtem_a);
    /* reset allocation flag */
    allocated = 0;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th3_temper_final */


/*======================================================================*/
/*!
\brief Calculate temperature at point on element parameter domain

The element stiffness matrix, ie the tangent operator, is determined
for the linear planar heat conduction problem.

\param   *container     CONTAINER   (i)   container data
\param   *ele           ELEMENT     (i)   pointer to current element
\param    r             DOUBLE      (i)   r-coord of parameter (r,s,t)
\param    s             DOUBLE      (i)   s-coord of parameter (r,s,t)
\param    t             DOUBLE      (i)   t-coord of parameter (r,s,t)
\param   *tem           DOUBLE      (o)   temperature at (r,s)

\return void

\author bborn
\date 10/06
*/
void th3_temper_cal(const CONTAINER *container,
                    const ELEMENT *ele,
                    const DOUBLE r,
                    const DOUBLE s,
                    const DOUBLE t,
                    DOUBLE *tem)
{
  ARRAY_POSITION_SOL *isol 
    = &(field[genprob.numtf].dis[container->disnum_t].ipos.isol); 
  INT nelenod;
  INT neledof;
  DOUBLE rr=0.0, ss=0.0, tt=0.0;  /* Gauss coordinate in THERM3 parameter space */
  DOUBLE shape[MAXNOD_THERM3];  /* shape functions */
  DOUBLE deriv[MAXNOD_THERM3][NDIM_THERM3];  /* derivatives of shape fct */
  INT k;  /* loop index */

  /*====================================================================*/
#ifdef DEBUG
  dstrc_enter("th3_temper_cal");
#endif

  /*--------------------------------------------------------------------*/
  /* dimensions */
  nelenod = ele->numnp;
  neledof = NUMDOF_THERM3 * nelenod;

  /*--------------------------------------------------------------------*/
  /* reset to zero nodal temper., shape funct. and derivatives */
  amzero(&nodtem_a);
  memset(shape, 0, sizeof(shape));
  memset(deriv, 0, sizeof(deriv));

  /*--------------------------------------------------------------------*/
  /* transform parametric coordinates if necessary */
  switch (ele->e.th3->struct_ele->eltyp)
  {
#ifdef D_BRICK1
    case el_brick1:
      rr = s;
      ss = -r;
      tt = t;
      break;
#endif
#ifdef D_SOLID3
    case el_solid3:
      rr = r;
      ss = s;
      tt = t;
      break;
#endif
    default:
      dserror("Element type is not permissible!");
      break;
  }

  /*--------------------------------------------------------------------*/
  /* get shape funtions at (r,s,t) */
  th3_shape_deriv(ele->distyp, rr, ss, tt, 0, shape, deriv);

  /*--------------------------------------------------------------------*/
  /* current nodal temperature */
  for (k=0; k<nelenod; k++)
  {
    nodtem[k] = ele->node[k]->sol.a.da[isol->tem][0];
  }

  /*--------------------------------------------------------------------*/
  /* interpolate temperature */
  *tem = 0.0;
  for (k=0; k<nelenod; k++)
  {
    *tem += shape[k]*nodtem[k];
  }

  /*====================================================================*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th3_temper_cal */

#endif /* end of #ifdef D_THERM3 */
/*! @} (documentation module close)*/
