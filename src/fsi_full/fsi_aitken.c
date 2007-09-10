/*!----------------------------------------------------------------------
\file
\brief AITKEN iteration

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

*----------------------------------------------------------------------*/

#ifndef CCADISCRET
/*!
\addtogroup FSI
*//*! @{ (documentation module open)*/


#include "../headers/standardtypes.h"
#include "fsi_prototypes.h"


#ifdef D_FSI


static DOUBLE done = 1.0E20;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
extern struct _GENPROB       genprob;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | ranks and communicators                                              |
 | This structure struct _PAR par; is defined in main_ccarat.c
 *----------------------------------------------------------------------*/
 extern struct _PAR   par;


/*!---------------------------------------------------------------------*
\brief compute relaxation parameter via AITKEN iteration

<pre>                                                         genk 01/03

RELAX(i) = 1 - NU(i) (NU(i) at time (n+1))
NU(i) = N(i-1) + [NU(i-1) - 1)] * TOP/DEN
TOP = [DELTAd(i) - DETLAd(i+1)] * DELTAd(i+1)
DEN = [DELTAd(i) - DETLAd(i+1)]^2
DELTAd(i) = d(i-1) - d~(i)
d(i-1) ... relaxad interface displacements of last iteration step
d~(i)  ... unrelaxed interface displacements of this iteration step

see also Dissertation of D.P.MOK, chapter 6.4

</pre>

\param *structfield   FIELD	     (i)   structural field
\param  itnum         INT            (i)   actual iteration step
\return void

------------------------------------------------------------------------*/
void fsi_aitken(
  FIELD             *structfield,
  INT                disnum,
  INT                itnum,
  INT                init
  )
{
INT            i,j;           /* some counters                          */
INT            numdf;         /* number of nodal dofs                   */
INT            dof;           /* actual dof                             */
static INT     numnp_total;   /* total number of nodes                  */
static INT     numdf_total;   /* total number of coupled dofs           */
static INT    *sid;           /* flags for structural interface dofs    */
static DOUBLE  nu;            /* actual AITKEN factor                   */
DOUBLE         top=ZERO;
DOUBLE         den=ZERO;
DOUBLE         del2;
DOUBLE       **sol_mf;        /* multifield nodal solution array        */
static ARRAY   del_a;
static DOUBLE *del;           /* DELTAd(i)                              */
NODE          *actsnode;      /* actual structural node                 */
static FSI_DYNAMIC *fsidyn;
ARRAY_POSITION *ipos;

#ifdef DEBUG
dstrc_enter("fsi_aitken");
#endif


switch (init)
{
case 0: /* initialisation */
   fsidyn = alldyn[3].fsidyn;

   if (genprob.restart==0)
      nu       = ZERO;
   else
      nu       = ONE - fsidyn->relax;


   numnp_total = structfield->dis[disnum].numnp;
   sid         = fsidyn->sid.a.iv;
   numdf_total = fsidyn->sid.fdim;
   del         = amdef("del",&del_a,numdf_total,1,"DV");
break;


case 1: /* calculation */
   if (itnum==0)
      aminit(&del_a,&done);

   ipos = &(structfield->dis[disnum].ipos);

   /* the solution history sol_mf of the structfield:
      - sol_mf[0][j] holds the latest struct-displacements
      - sol_mf[1][j] holds the (relaxed) displacements of the last iteration step
    */

   /* loop structure nodes */
   for (i=0;i<numnp_total;i++)
   {
      actsnode  = &(structfield->dis[disnum].node[i]);
      numdf = actsnode->numdf;
      sol_mf = actsnode->sol_mf.a.da;

      /* loop dofs and check for coupling */
      for (j=0;j<numdf;j++)
      {
         dof = actsnode->dof[j];
         dsassert(dof<numdf_total,"dofnumber not valid!\n");
         if (sid[dof]==0) continue;
         del2     = del[dof];
         del[dof] = sol_mf[ipos->mf_reldisp][j] - sol_mf[ipos->mf_dispnp][j];
         del2     = del2 - del[dof];
         top     += del2*del[dof];
         den     += del2*del2;

      } /* end of loop over dofs */
   } /* end of loop over nodes */

   nu = nu + (nu - ONE)*top/den;

   /* finalising */
   fsidyn->relax = ONE - nu;

   /* output to the screen */
   if (par.myrank==0)
   printf("\nAITKEN ITERATION: RELAX = %.5lf\n\n",fsidyn->relax);

break;

case 2:
/* experimental post-processing due to
 * M. Krizek, L. Liu, P. Neittaanmäki: "Post-processing of Gauss-Seidel iterations",
 * Numer. Linear Algebra Appl., 6, 147-156 (1999)
 * */
{
  DOUBLE q = 0.;
  INT counter = 0;

  if (itnum==0)
    aminit(&del_a,&done);

  ipos = &(structfield->dis[disnum].ipos);

  /* the solution history sol_mf of the structfield:
     - sol_mf[0][j] holds the latest struct-displacements
     - sol_mf[1][j] holds the (relaxed) displacements of the last iteration step
   */

  /* loop structure nodes */
  for (i=0;i<numnp_total;i++)
  {
    actsnode  = &(structfield->dis[disnum].node[i]);
    numdf = actsnode->numdf;
    sol_mf = actsnode->sol_mf.a.da;

    /* loop dofs and check for coupling */
    for (j=0;j<numdf;j++)
    {
      dof = actsnode->dof[j];
      dsassert(dof<numdf_total,"dofnumber not valid!\n");
      if (sid[dof]==0) continue;
      del2     = del[dof];
      del[dof] = sol_mf[ipos->mf_reldisp][j] - sol_mf[ipos->mf_dispnp][j];

      q += del[dof] / del2;
      counter += 1;
    }
  }

  q /= counter;
  fsidyn->relax = 1. / (1 - q);

  /* output to the screen */
  if (par.myrank==0)
    printf("\nRELAX = %.5lf\n\n",fsidyn->relax);

  break;
}

default:
  dserror("parameter %d not supported",init);
}  /* switch (init) */


#ifdef DEBUG
dstrc_exit();
#endif

return;

} /* end of fsi_aitken */



void fsi_aitken_force(
  FIELD             *structfield,
  INT                sdisnum,
  FIELD             *fluidfield,
  INT                fdisnum,
  INT                itnum,
  INT                numff,
  INT                init
  )
{
  INT            i,j,k,m;           /* some counters                          */
  INT            numdf;         /* number of nodal dofs                   */
  INT            dof;           /* actual dof                             */
  static INT     numnp_total;   /* total number of nodes                  */
  static INT     numdf_total;   /* total number of coupled dofs           */
  static INT    *sid;           /* flags for structural interface dofs    */
  static DOUBLE  nu;            /* actual AITKEN factor                   */
  DOUBLE         top=ZERO;
  DOUBLE         den=ZERO;
  DOUBLE         del2;
  DOUBLE       **sol_mf;        /* multifield nodal solution array        */
  static ARRAY   del_a;
  static DOUBLE *del;           /* DELTAd(i)                              */
  NODE          *actsnode;      /* actual structural node                 */
  static FSI_DYNAMIC *fsidyn;

  ARRAY_POSITION *fluid_ipos;

#ifdef DEBUG
  dstrc_enter("fsi_aitken_force");
#endif

  switch (init)
  {
  case 0: /* initialisation */
    fsidyn = alldyn[3].fsidyn;

    if (genprob.restart==0)
      nu       = ZERO;
    else
      nu       = ONE - fsidyn->relax;


    numnp_total = structfield->dis[sdisnum].numnp;
    sid         = fsidyn->sid.a.iv;
    numdf_total = fsidyn->sid.fdim;
    del         = amdef("del",&del_a,numdf_total,1,"DV");
    break;


  case 1: /* calculation */
    if (itnum==0)
      aminit(&del_a,&done);

    fluid_ipos = &(fluidfield->dis[fdisnum].ipos);

    /* the solution history sol_mf of the structfield:
       - sol_mf[1][j] holds the latest fluid force
       - sol_mf[2][j] holds the (relaxed) forces of the last iteration step
    */


    /* loop structure nodes */
    for (i=0;i<numnp_total;i++)
    {
      actsnode  = &(structfield->dis[sdisnum].node[i]);
      numdf = actsnode->numdf;

      /* loop dofs and check for coupling */
      for (j=0;j<numdf;j++)
      {
        NODE   *actfnode;
        dof = actsnode->dof[j];
        dsassert(dof<numdf_total,"dofnumber not valid!\n");
        if (sid[dof]==0) continue;

#ifdef FSI_NONMATCH
        /* The forces are calculated by the fluid. */
        /*actfnode = actsnode->gnode->mfcpnode[numff];*/
        for (k=0;k<actsnode->numele;k++)
        {
          ELEMENT *actele;
          if (actsnode->element[k]->coupleptr==NULL) continue;

          actele=actsnode->element[k];
          /*In diesem Element sind jetzt sicher Fluid Knoten!*/
          for (m=0;m<actele->coupleptr->numnp;m++)
          {
            actfnode=actele->coupleptr->couplenode[m];
            sol_mf = actfnode->sol_mf.a.da;

            del2     = del[dof];
            del[dof] = (sol_mf[fluid_ipos->mf_forcen][j] -
                        sol_mf[fluid_ipos->mf_forcenp][j]);
            del2     = del2 - del[dof];
            top     += del2*del[dof];
            den     += del2*del2;
          }
        }
#else
        actfnode = actsnode->gnode->mfcpnode[numff];
        sol_mf = actfnode->sol_mf.a.da;

        del2     = del[dof];
        del[dof] = (sol_mf[fluid_ipos->mf_forcen][j] -
                    sol_mf[fluid_ipos->mf_forcenp][j]);
        del2     = del2 - del[dof];
        top     += del2*del[dof];
        den     += del2*del2;
#endif
      }
    }

    nu = nu + (nu - ONE)*top/den;

    /* finalising */
    fsidyn->relax = ONE - nu;

    /* output to the screen */
    if (par.myrank==0)
      printf("\nAITKEN ITERATION: RELAX = %.5lf\n\n",fsidyn->relax);

    break;

  case 2:				/* reset */
    nu = 0;
    break;

  }  /* switch (init) */


#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}



#endif

/*! @} (documentation module close)*/

#endif
