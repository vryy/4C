/*!----------------------------------------------------------------------
\file
\brief service functions for projection algorithm
------------------------------------------------------------------------*/
/*!
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fluid_prototypes.h"
#include "fluid_pm_prototypes.h"
#include "../fluid2_pro/fluid2pro.h"
/*!----------------------------------------------------------------------
\brief positions of physical values in node arrays

<pre>                                                        chfoe 11/04

This structure contains the positions of the various fluid solutions 
within the nodal array of sol_increment.a.da[ipos][dim].

extern variable defined in fluid_service.c
</pre>

------------------------------------------------------------------------*/
extern struct _FLUID_POSITION ipos;
/*----------------------------------------------------------------------*
\brief  matrix product mat1 * lumped mat2

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

<pre>                                                      basol 11/02
                                                           genk  10/03
In this routine the matrix product:
     mat1 * lumped mat2
is evaluated.
     mat1 can be any SPARSE_TYP (if someone implements it! *gg*)
     lumped mat2 is a lumped matrix, stored in a redundant vector

</pre>
\param  *m       SPARSE_ARRAY  (i/o)  matrix stored in sparse array format
\param  *m_typ	 SPARSE_TYP    (i)    sparse_typ
\param  *lmat	 DOUBLE        (i)    lumped matrix
\return void
\warning up to now only working for MSR-matrix

*-----------------------------------------------------------------------*/
void fluid_pm_matlmatmul(
                         SPARSE_ARRAY   *m,
		         SPARSE_TYP     *m_typ,
                         DOUBLE         *lmat
		        )
{
INT         i,index;
INT         numeq,nnz;
INT        *update;
INT        *bindx;
DOUBLE     *val;

#ifdef DEBUG
dstrc_enter("fluid_pm_matlmatmul");
#endif

/*----------------------------------------------------------------------*
 | REMARK:: the coming lmass is inverse lumped mass matrix              |
 *----------------------------------------------------------------------*/

switch (*m_typ)
{
case msr:
   numeq  = m->msr->numeq;
   val    = m->msr->val.a.dv;
   update = m->msr->update.a.iv;
   bindx  = m->msr->bindx.a.iv;
   nnz    = m->msr->nnz;
   /*--------------------------------------------------- loop diagonals */
   for (i=0; i<numeq; i++)
   {
      index = update[i];
      val[i] *= lmat[index];
   }
   /*----------------------------------- loop all other nonzero entries */
   for (i=numeq+1;i<=nnz;i++)
   {
      index = bindx[i];
      val[i] *= lmat[index];
   }
break;
default:
   dserror("matrix product not implemented for actual sparse typ");
break;
}

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_pm_matlmatmul*/

/*!---------------------------------------------------------------------
\brief this function does the (OLL) times a vector(redundant) multiplication

<pre>                                                        basol 11/02
                                                             genk  10/03
  option = 0: rc = P * r
  option = 1: rc = P_trans * r

</pre>
\param *rc_a       ARRAY      (o)   result of the multiplication
\param *r_a        ARRAY      (i)   vector to be multiplied
\param *P_oll      OLL        (i)   matrix in OLL format
\param *actintra   INTRA      (i)   intra communicator
\param  numeq      INT        (i)   number of equations on this proc
\param  option     INT        (i)   evaluation flag
\return void

------------------------------------------------------------------------*/
void fluid_pm_matvecmul(ARRAY *rc_a,     ARRAY *r_a, OLL *P_oll,
                        INTRA *actintra, INT numeq,  INT option)
{
INT        i;              /* simply a counter                          */
INT        cindex,rindex;  /* row/column indes                          */
DOUBLE    *r,*rc;          /* vector entries                            */
DOUBLE     sum;
MATENTRY  *actmatentry;    /* actual OLL matrix entry                   */
#ifdef PARALLEL
INT        rc_numeq;       /* number of equation of result vector       */
ARRAY      recvbuf_a;
DOUBLE    *recvbuf;        /* receive buffer                            */
#endif

#ifdef DEBUG
dstrc_enter("fluid_pm_matvecmul");
#endif

/*---------------------------------------------------------- initialise */
rc          = rc_a->a.dv;
r           = r_a->a.dv;
amzero(rc_a);

/*---------------------------------------- do the matrix vector product */
switch (option) /* rc = P * r */
{
case 0: /* rc = P * r -> OLL format */
   for (i=0; i<numeq; i++)
   {
      actmatentry = P_oll->row[i];
      rindex = actmatentry->r;
      sum = ZERO;
      while (actmatentry != NULL)
      {
	 cindex = actmatentry->c;
	 sum += actmatentry->val*r[cindex];
	 actmatentry = actmatentry->rnext;
      }
      rc[rindex] = sum;
   }
break;
case 1: /* rc = P_trans * r -> OLL format */
   for (i=0; i<numeq; i++)
   {
      actmatentry = P_oll->row[i];
      rindex = actmatentry->r;
      while (actmatentry != NULL)
      {
         cindex = actmatentry->c;
	 rc[cindex] += actmatentry->val*r[rindex];
	 actmatentry = actmatentry->rnext;
      }
   }
break;
default:
   dserror("option out of range: don't know what to do!\n");
}

#ifdef PARALLEL
rc_numeq    = rc_a->fdim;
recvbuf=amdef("recvbuf",&recvbuf_a,rc_numeq,1,"DV");
amzero(&recvbuf_a);
MPI_Allreduce(rc,recvbuf,rc_numeq,
              MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
for (i=0;i<rc_numeq;i++) rc[i] = recvbuf[i];
amdel(&recvbuf_a);
#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_pm_matvecmul */


/*!---------------------------------------------------------------------
\brief make vector by lumped mass matrix  multiplication

<pre>                                                      basol 11/02

(Ml-1)*C*phi multiplication is done here

</pre>

\param *vec          DOUBLE     (o)   vector to be multiplied (C*phi)
\param *lmat         DOUBLE     (i)   lumped mass matrix
\param  numeq_total  INT        (i)   number of equations
\return void

*----------------------------------------------------------------------*/
void fluid_pm_lmatmulvec(
                         DOUBLE   *vec,
		         DOUBLE   *lmat,
		         INT       numeq_total
			)
{
INT         i;

#ifdef DEBUG
dstrc_enter("fluid_pm_lmatmulvec");
#endif

for (i=0; i<numeq_total; i++)
vec[i] *= lmat[i];

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_pm_lmatmulvec */

/*----------------------------------------------------------------------*/
/*brief form the full velocity vector

<pre>                                                   basol 01/03

  copy fluid solution from the nodes including the dirichlet values
  to a redundant vector

<\pre>
\param *actfield FIELD  (i) actual field
\param  disnum   INT    (i) indice of the discretization to be used
\param *fullvel  DOUBLE (o) full velocity vector
\param  place	 INT    (i) row to get values in the ARRAY sol_increment
\return void

 *----------------------------------------------------------------------*/
void fluid_pm_fullvel(FIELD *actfield, INT disnum,DOUBLE *fullvel, INT place)
{
INT               i,j;     /* simply some counters                      */
INT               actdof;  /* actual dof number                         */
NODE             *actnode; /* actual node                               */
DISCRET          *actdis;  /* actual discretisation                     */

#ifdef DEBUG
dstrc_enter("fluid_pm_fullvel");
#endif

actdis = &(actfield->dis[disnum]);

/*------------ loop nodes and put the result to the full velocity vector*/
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   for (j=0;j<actnode->numdf;j++)
   {
      actdof = actnode->dof[j];
      fullvel[actdof] = actnode->sol_increment.a.da[place][j];
   }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_pm_fullvel */

/*----------------------------------------------------------------------*/
/*brief form the pressure solution on the velocity discretisaton

<pre>                                                   basol 01/03

Pressure and velocity solution exist on different discretisations.
For visualisation it's neccessary to have the solution on one discr.
So the pressure values are copied and interpolated to the velocity
discretisation.

<\pre>
\param *actfield FIELD  (i) actual field
\return void

 *----------------------------------------------------------------------*/
void fluid_pm_pretovel(FIELD *actfield,INT actpos)
{
int i,j;
double epre[MAXDOFPERELE];
ELEMENT *actvele, *actpele;
NODE *actnode;

#ifdef DEBUG
dstrc_enter("fluid_pm_pretovel");
#endif

/*------------------------------------------------------ loop elements */
for (i=0;i<actfield->dis[0].numele;i++)
{
   actpele=&(actfield->dis[1].element[i]);
   switch (actpele->e.f2pro->dm)
   {
   case dm_q2q1:
      for (j=0;j<actpele->numnp;j++)
      {
         actnode=actpele->node[j];
         epre[j] = actnode->sol.a.da[actpos][0];
      }
      /*------------------------------------- interpolate for the rest */
      epre[4]=(epre[0]+epre[1])/TWO;
      epre[5]=(epre[1]+epre[2])/TWO;
      epre[6]=(epre[2]+epre[3])/TWO;
      epre[7]=(epre[3]+epre[0])/TWO;
      epre[8]=(epre[0]+epre[1]+epre[2]+epre[3])/FOUR;
      /*------------------------------------- write values to fluiddis */
      actvele=&(actfield->dis[0].element[i]);
      for (j=0;j<actvele->numnp;j++)
      {
         actnode=actvele->node[j];
	 dsassert(actnode->sol.sdim>=3,"wrong dimension for sol-array!\n");
         actnode->sol.a.da[actpos][2]=epre[j];
      }
   break;
   default:
      dserror("Dismode unknown!\n");
   }
}


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
}
/*!---------------------------------------------------------------------
\brief init positions in sol_increment 

<pre>                                                        chfoe 01/05

This routine inits the positions in sol_increment in the case of 
projection method.
Unused flags are set to '-1' in order to detect misuse easily.

</pre>
\return void

------------------------------------------------------------------------*/
void fluid_init_pos_pm(void)
{
#ifdef DEBUG
  dstrc_enter("fluid_init_pos_pm");
#endif

/*---------------------------------------- adaptive time stepping ---*/
ipos.veln  = 0; 
ipos.velnp = 1;
ipos.velnm =-1; 
ipos.hist  =-1; 
ipos.accnm =-1;
ipos.accn  =-1;
ipos.pred  =-1;
ipos.terr  =-1;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif
/*! @} (documentation module close)*/
