/*-----------------------------------------------------------------------*/
/*!
\file
\brief contains the routine 'assemble_fast' which controls the assembling of a
       set of fast elements.


<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

 */
/*-----------------------------------------------------------------------*/

/*!
\addtogroup Fluid3_fast
*//*! @{ (documentation module open)*/

#ifdef D_FLUID3_F


#include "../headers/standardtypes.h"
#include "../solver/solver.h"


/*----------------------------------------------------------------------*
  | global dense matrices for element routines             m.gee 7/01  |
 *----------------------------------------------------------------------*/
extern struct _ARRAY etforce_global;  /* element Time RHS */
extern struct _ARRAY eiforce_global;  /* element Iteration RHS */
extern struct _ARRAY edforce_global;  /* element dirichlet RHS */

extern struct _ARRAY estif_fast;
extern struct _ARRAY emass_fast;
extern struct _ARRAY etforce_fast;
extern struct _ARRAY eiforce_fast;
extern struct _ARRAY edforce_fast;



/*-----------------------------------------------------------------------*/
/*!
  \brief routine to assemble the element matrices of a set of elements into the
         global sparse matrix

  Very detailed description what the function is doing why.

  \param sysarray         INT         (i) number of sparse system matrix
  \param actpart         *PARTITION   (i) my partition of this field
  \param actsolv         *SOLVAR      (i) active SOLVAR
  \param actintra        *INTRA       (i) my intra-communicator
  \param elevec[]        *ELEMENT     (i) vector with the elements
  \param assemble_action *ASSEMBLE_ACTION (i) assemble option
  \param iel              INT         (i) number of nodes of one element
  \param hasext[]         INT         (i) flags for external forces
  \param hasdirich[]      INT         (i) flags for dirichlet conds
  \param aloopl           INT         (i) number of elements in elevec

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void assemble_fast(
    INT                    sysarray,        /* number of first sparse system matrix */
    struct _PARTITION     *actpart,         /* my partition of the active field */
    struct _SOLVAR        *actsolv,         /* the active SOLVAR */
    struct _INTRA         *actintra,        /* the active intracommunicator */
    struct _ELEMENT       *elevec[LOOPL],   /* the element to assemble */
    enum _ASSEMBLE_ACTION  assemble_action, /* the assembly option */
    CONTAINER             *container,
    INT                    iel,
    INT                    hasext[LOOPL],
    INT                    hasdirich[LOOPL],
    INT                    aloopl)
{

static  INT        *invupd   = NULL;
static  INT        *invbindx = NULL;


  enum  _SPARSE_TYP    sysa_typ;
  union _SPARSE_ARRAY *sysa;

  DOUBLE *estif_f,*emass_f;
  DOUBLE  *etforce,*etforce_f,*eiforce,*eiforce_f,*edforce,*edforce_f;
  INT   i,j,k,l;
  INT        counter;
  ELEMENT          *actele;

  INT         nd;
  INT         numnp = MAXNOD*MAXDOFPERNODE;
  INT         loopl = LOOPL;
  INT         myrank,nprocs;
  INT         numeq_total, numeq;

  DOUBLE     *val;
  INT        *update, *bindx;
  INT       **cdofs;
  INT         ncdofs;

#if 0  /* coupling not yet implemented */
#ifdef PARALLEL
  INT         nsend;
  INT       **isend;
  DOUBLE    **dsend;
#endif
#endif

  INT         lm[MAXNOD*MAXDOFPERNODE];
  INT         owner[MAXNOD*MAXDOFPERNODE];
#ifdef SOLVE_DIRICH
  INT         dirich[MAXNOD*MAXDOFPERNODE];
#endif


  eiforce = eiforce_global.a.dv;
  etforce = etforce_global.a.dv;
  edforce = edforce_global.a.dv;

  estif_f = estif_fast.a.dv;
  emass_f = emass_fast.a.dv;
  etforce_f = etforce_fast.a.dv;
  eiforce_f = eiforce_fast.a.dv;
  edforce_f = edforce_fast.a.dv;


  /* do nothing */
  if (assemble_action==assemble_do_nothing)
    goto end;


  /* check for presence and typ of system matrices */
  if (sysarray>=0)
  {
    sysa       = &(actsolv->sysarray[sysarray]);
    sysa_typ   =   actsolv->sysarray_typ[sysarray];
  }
  else
  {
    sysa     = NULL;
    sysa_typ = sparse_none;
  }


  /* add to 1 system matrix */
  if (assemble_action==assemble_one_matrix)
  {
    switch(sysa_typ)
    {

      case msr:

        myrank     = actintra->intra_rank + 1;
        nprocs     = actintra->intra_nprocs;

        numeq_total= sysa->msr->numeq_total;
        numeq      = sysa->msr->numeq;
        update     = sysa->msr->update.a.iv;
        bindx      = sysa->msr->bindx.a.iv;
        cdofs      = actpart->pdis[0].coupledofs.a.ia;
        ncdofs     = actpart->pdis[0].coupledofs.fdim;

        val        = sysa->msr->val.a.dv;

        if (!invupd)
        {

#ifdef PERF
    perf_begin(64);
#endif
          /* allocate invupd */
          invupd = (INT*)CCACALLOC( numeq_total,sizeof(INT));
          if (!invupd) dserror("Allocation of invupd failed");

          /* initialize with minus one */
          for (k=0; k<numeq_total; k++)
          {
            invupd[k] = -1;
          }

          /* fill invupd */
          for (k=0; k<numeq; k++)
          {
            invupd[update[k]] = k+1;
          }
#ifdef PERF
    perf_end(64);
#endif
        }


#if 0  /* coupling not yet implemented */
#ifdef PARALLEL
        nsend = sysa->msr->couple_i_send->fdim;
        isend = sysa->msr->couple_i_send->a.ia;
        dsend = sysa->msr->couple_d_send->a.da;
#endif
#endif


        if (!invbindx)
        {
#ifdef PERF
    perf_begin(65);
#endif
          /* allocate invbindx */
          invbindx = (INT*)CCACALLOC( numeq_total,sizeof(INT));
          if (!invbindx) dserror("Allocation of invbindx failed");
#ifdef PERF
    perf_end(65);
#endif
        }


        /* do one element after another */
        for(l=0;l<aloopl;l++)
        {
          actele=elevec[l];



#ifdef PERF
    perf_begin(66);
#endif
          /* make location vector lm*/
          counter=0;
          for (i=0; i<actele->numnp; i++)
          {
            for (j=0; j<actele->node[i]->numdf; j++)
            {
              lm[counter]    = actele->node[i]->dof[j]+1;
#ifdef PARALLEL
              owner[counter] = actele->node[i]->proc+1;
#endif

#ifdef SOLVE_DIRICH
              if (actele->node[i]->gnode->dirich!=NULL &&
                  actele->node[i]->gnode->dirich->dirich_onoff.a.iv[j]!=0)
                dirich[counter] = 1;
              else
                dirich[counter] = 0;
#endif

              counter++;
            }
          }
          /* end of loop over element nodes */
          nd = counter;
#ifdef PERF
    perf_end(66);
#endif


          /* now start looping the dofs */
          /* loop over i (the element row) */


#ifdef PERF
    perf_begin(67);
#endif


#ifndef SOLVE_DIRICH

#ifndef PARALLEL
          faddmsr(
              estif_f,
              &lm[0],
              &owner[0],
              invupd,
              bindx,
              invbindx,
              val,
              &myrank,
              &nprocs,
              &numeq_total,
              &numeq,
              &numnp,
              &nd,
              &aloopl,
              &loopl,
              &l);
#else
          faddmsrp(
              estif_f,
              &lm[0],
              &owner[0],
              invupd,
              bindx,
              invbindx,
              val,
              &myrank,
              &nprocs,
              &numeq_total,
              &numeq,
              &numnp,
              &nd,
              &aloopl,
              &loopl,
              &l);
#endif

#else

#ifndef PARALLEL
          faddmsrd(
              estif_f,
              &lm[0],
              &owner[0],
              &dirich[0],
              invupd,
              bindx,
              invbindx,
              val,
              &myrank,
              &nprocs,
              &numeq_total,
              &numeq,
              &numnp,
              &nd,
              &aloopl,
              &loopl,
              &l);
#else
          faddmsrdp(
              estif_f,
              &lm[0],
              &owner[0],
              &dirich[0],
              invupd,
              bindx,
              invbindx,
              val,
              &myrank,
              &nprocs,
              &numeq_total,
              &numeq,
              &numnp,
              &nd,
              &aloopl,
              &loopl,
              &l);
#endif

#endif


#ifdef PERF
    perf_end(67);
#endif



#ifdef PERF
    perf_begin(68);
#endif
          for(j=0;j<iel*4;j++)
          {
            etforce[j]=etforce_f[j*LOOPL+l];
            eiforce[j]=eiforce_f[j*LOOPL+l];
            edforce[j]=edforce_f[j*LOOPL+l];
          }
#ifdef PERF
    perf_end(68);
#endif

#ifdef PERF
    perf_begin(69);
#endif
          /* assemble time rhs */
          if (container->nif!=0)
          {
            container->dvec = container->ftimerhs;
            assemble_intforce(actele,&etforce_global,container,actintra);
          }

          /* assemble iteration rhs */
          if (container->nii+hasext[l]!=0 || container->nim!=0)
          {
            container->dvec = container->fiterhs;
            assemble_intforce(actele,&eiforce_global,container,actintra);
          }

          /* assemble edforce */
          if (hasdirich[l]!=0)
          {
            container->dvec = container->fiterhs;
            assemble_intforce(actele,&edforce_global,container,actintra);
          }

          container->dvec=NULL;
#ifdef PERF
    perf_end(69);
#endif

        } /* for(i=0;i<aloopl;i++) */

        break;

      case sparse_none:
        dserror("Unspecified type of system matrix");
        break;

      default:
        dserror("Unspecified type of system matrix");
        break;
    }

    goto end;
  } /* if (assemble_action==assemble_one_matrix) */

  dserror("Unknown assemble assemble_action!!");


end:
  return;

}


#endif /* ifdef D_FLUID3_F */

/*! @} (documentation module close)*/

