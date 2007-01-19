/*!---------------------------------------------------------------------
\file
\brief controlling calculation of consistent boundary forces
       for fsi-coupling

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6572
</pre>

---------------------------------------------------------------------*/
/*!
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../fluid2/fluid2.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../fluid3/fluid3.h"
#include "../fluid3/fluid3_prototypes.h"
#include "../fluid3_fast/f3f_prototypes.h"
#include "../fluid2_is/fluid2_is.h"
#include "../fluid3_is/fluid3_is.h"
#include "../fluid2_pro/fluid2pro.h"
#include "../fluid3_pro/fluid3pro.h"
#include "../fluid2_pro/fluid2pro_prototypes.h"
#include "../fluid3_pro/fluid3pro_prototypes.h"
#include "../shell8/shell8.h"
#ifdef D_ALE
#include "../ale2/ale2.h"
#endif

#define EPS   1e-7

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 | global dense matrices for element routines             m.gee 7/01    |
 *----------------------------------------------------------------------*/
extern struct _ARRAY estif_global;
extern struct _ARRAY eforce_global;   /* element RHS                    */
extern struct _ARRAY eforce_fast;     /* element RHS(fortran)           */


/*!---------------------------------------------------------------------
\brief controls the evaluation of consistent boundary forces for for
       fsi-coupling

<pre>                                                        chfoe 05/05
This routine controls the evaluation of consistent boundary forces along
dirichlet boundary lines for fsi-coupling.
Within an init phase (init=1) all elements which participate at a fsi
coupling line are searched and labled by setting the flag
actele->e.f2->force_on on 1.
In the main part these elements are visited again, its interesting nodes
(the ones which are situated on the fsi coupling line) are found and the
elemental residual vector is calculated.
The nodal forces are written to different destinations depending on the
code version:

I.  SEQUENTIELL: -> node forces go directly into sol_mf[1] field on fluid
                -> *fcouple pointer is NULL

II. PARALLEL:    -> node forces are written into fcouple which is of
                    dimension (numdof - numeq) i.e. number of Dirichlet
                    dofs of the field.
                 -> fcouple is indicated by (dof(node) - numeq)
                 -> fcouple is allreduced after the call of fsi_cbf

III. PARALLEL with SOLVE_DIRICH defined
                 -> numdof = numeq since Dirichlet dofs are not sortet out
                 -> fcouple is a full vector indicated by dof
                 -> fcouple is allreduced after the call of fsi_cbf

</pre>
\param   *actpdis        PARTDISCRET    (i)
\param    init           INT            (i)   init flag

\warning At the moment only elements with the same number of dofs at all
         nodes can be treated.

\return void

------------------------------------------------------------------------*/
void fsi_cbf(
    PARTDISCRET    *actpdis,
    DOUBLE         *fcouple,
    ARRAY_POSITION *ipos,
    INT             numeq_total,
    INT             is_relax,
    INT             init)
{

  INT i, j;
  INT numdf;
  INT hasdirich, hasext;
  INT force_on_node[MAXNOD];
#ifdef D_FLUID3_F
  INT nfnode;  /* number of nodes of actele where forces are searched for */
  INT line;
  INT l;
#endif

#ifdef PARALLEL
  INT dofx, dofy, dofz;
#endif


  DOUBLE    rho;

  static DOUBLE  *eforce;

  NODE    *actnode;
  GNODE   *actgnode;
  ELEMENT *actele;

#ifdef D_FLUID3_F
  INT               hasext_f[LOOPL];
  INT               aloopl, coupled;
  FAST_ELES        *act_fast_eles;
#endif



#ifdef DEBUG
  dstrc_enter("fsi_cbf");
#endif



  /* initialisation */
  if (init==1)
  {

    /* preselect elements that have fsi coupling */
    switch (genprob.ndim)
    {
      case 2:
#ifdef D_FLUID2
        for(i=0; i<actpdis->numele; i++)
        {
          actele = actpdis->element[i];
          if (par.myrank!=actele->proc)
            continue; /* only my elments */

          for(j=0; j<actele->numnp; j++)
          {
            actgnode = actele->node[j]->gnode;
            /* check if there is a coupled struct node */
            if (actgnode->mfcpnode[genprob.numsf]==NULL)
              continue;

            if (actele->eltyp==el_fluid2)
              actele->e.f2->force_on = 1;
#ifdef D_FLUID2_IS
            else if (actele->eltyp==el_fluid2_is)
              actele->e.f2is->force_on = 1;
#endif
#ifdef D_FLUID2_PRO
            else if (actele->eltyp==el_fluid2_pro)
              actele->e.f2pro->force_on = 1;
#endif
            else
              dserror("element type %d unsupported",actele->eltyp);
            break;
          }
        }
#endif
        break;

      case 3:
#ifdef D_FLUID3
        for(i=0; i<actpdis->numele; i++)
        {
          actele = actpdis->element[i];
          if (par.myrank!=actele->proc)
            continue; /* only my elements */

          for(j=0; j<actele->numnp; j++)
          {
            actgnode = actele->node[j]->gnode;
#ifndef FSI_NONMATCH
            /* check if there is a coupled struct node */
            if (actgnode->mfcpnode[genprob.numsf]==NULL)
              continue;
#else
            if (actgnode->fsicouple==NULL) continue;
#endif

            if (actele->eltyp==el_fluid3 || actele->eltyp==el_fluid3_fast)
              actele->e.f3->force_on = 1;
#ifdef D_FLUID3_IS
            else if (actele->eltyp==el_fluid3_is)
              actele->e.f3is->force_on = 1;
#endif
#ifdef D_FLUID3_PRO
            else if (actele->eltyp==el_fluid3_pro)
              actele->e.f3pro->force_on = 1;
#endif
            else
              dserror("element type %d unsupported",actele->eltyp);
            break;
          }
        }
#endif
        break;

      default:
        dserror("genprob.ndim != 2 or 3!!");
    }

    goto end;
  }

  /* the elemental force vector */
  eforce  = eforce_global.a.dv;

  /* evaluate fsi coupling forces depending on dimension */
  switch (genprob.ndim)
  {
    case 2: /* problem is two-dimensional */
#ifdef D_FLUID2
      rho = mat[actpdis->element[0]->mat-1].m.fluid->density;
      numdf = 3;
      for(i=0; i<actpdis->numele; i++)
      {
        INT dof;
        actele = actpdis->element[i];

        if (actele->eltyp==el_fluid2)
        {
          if (actele->e.f2->force_on==0) continue; /* element not of interest*/
        }
#ifdef D_FLUID2_IS
        else if (actele->eltyp==el_fluid2_is)
        {
          if (actele->e.f2is->force_on==0) continue; /* element not of interest*/
        }
#endif
#ifdef D_FLUID2_PRO
        else if (actele->eltyp==el_fluid2_pro)
        {
          if (actele->e.f2pro->force_on==0) continue; /* element not of interest*/
        }
#endif
        else
          dserror("element type %d unsupported",actele->eltyp);

        for (j=0; j<actele->numnp; j++)
        {
          actgnode = actele->node[j]->gnode;
          /* check if there is a coupled struct node */
          if (actgnode->mfcpnode[genprob.numsf]!=NULL)
            force_on_node[j] = actele->node[j]->Id;
          else
            force_on_node[j] = -1;
        }

        /*--- get force vector ---*/
	if (!is_relax)
	{
	  if (actele->eltyp==el_fluid2)
	    f2_caleleres(actele,&eforce_global,ipos,&hasdirich,&hasext);
#ifdef D_FLUID2_IS
	  else if (actele->eltyp==el_fluid2_is)
	    f2is_caleleres(actele,&eforce_global,ipos,&hasdirich,&hasext);
#endif
#ifdef D_FLUID2_PRO
	  else if (actele->eltyp==el_fluid2_pro)
	    f2pro_caleleres(actele,&eforce_global,ipos,&hasdirich,&hasext);
#endif
	  else
	    dserror("element type %d unsupported",actele->eltyp);
	}
	else
	{
	  if (actele->eltyp==el_fluid2)
	    f2_caleleres_relax(actele,&estif_global,&eforce_global,ipos,&hasdirich,&hasext);
#ifdef D_FLUID2_IS
	  else if (actele->eltyp==el_fluid2_is)
	    f2is_caleleres_relax(actele,&estif_global,&eforce_global,ipos,&hasdirich,&hasext);
#endif
#ifdef D_FLUID2_PRO
	  else if (actele->eltyp==el_fluid2_pro)
	    dserror("no steepest decent with projection fsi right now");
#endif
	  else
	    dserror("element type %d unsupported",actele->eltyp);
	}

        dof = 0;
        for(j=0; j<actele->numnp; j++)
        {
          actnode = actele->node[j];

          if (force_on_node[j]!=-1)
          {
#ifdef PARALLEL
            dofx = actnode->dof[0];
            dofy = actnode->dof[1];
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
            fcouple[dofx] += eforce[dof  ]*rho;
            fcouple[dofy] += eforce[dof+1]*rho;
#else
            if(dofx<numeq_total)
              dserror("can not fsi couple free dof!\n");
            if(dofy<numeq_total)
              dserror("can not fsi couple free dof!\n");

            fcouple[dofx-numeq_total] += eforce[dof  ]*rho;
            fcouple[dofy-numeq_total] += eforce[dof+1]*rho;
#endif /* SOLVE_DIRICH */
#else
            actnode->sol_mf.a.da[ipos->mf_forcenp][0] += eforce[dof  ]*rho;
            actnode->sol_mf.a.da[ipos->mf_forcenp][1] += eforce[dof+1]*rho;
#endif  /* PARALLEL */
          }
          dof += actnode->numdf;
        }
      }
#endif
      break;



    case 3: /* problem is three-dimensional */
#ifdef D_FLUID3

      rho = mat[actpdis->element[0]->mat-1].m.fluid->density;
      numdf = 4;

      for(i=0; i<actpdis->numele; i++)
      {
        INT dof;
        actele = actpdis->element[i];

        /* fast elements will be treated later ---*/
        if (actele->eltyp == el_fluid3_fast) continue;

        if (actele->eltyp==el_fluid3)
        {
          if (actele->e.f3->force_on==0) continue; /* element not of interest*/
        }
#ifdef D_FLUID3_IS
        else if (actele->eltyp==el_fluid3_is)
        {
          if (actele->e.f3is->force_on==0) continue; /* element not of interest*/
        }
#endif
#ifdef D_FLUID3_PRO
        else if (actele->eltyp==el_fluid3_pro)
        {
          if (actele->e.f3pro->force_on==0) continue; /* element not of interest*/
        }
#endif
        else
          dserror("element type %d unsupported",actele->eltyp);

        for (j=0; j<actele->numnp; j++)
        {
          actgnode = actele->node[j]->gnode;
          /* check if there is a coupled struct node */
          if (actgnode->mfcpnode[genprob.numsf]!=NULL)
            force_on_node[j] = actele->node[j]->Id;
          else
            force_on_node[j] = -1;
        }

        /* get force vector */
	if (!is_relax)
	{
	  if (actele->eltyp==el_fluid3)
	    f3_caleleres(actele,&eforce_global,&hasdirich,&hasext,ipos);
#ifdef D_FLUID3_IS
	  else if (actele->eltyp==el_fluid3_is)
	    f3is_caleleres(actele,&eforce_global,ipos,&hasdirich,&hasext);
#endif
#ifdef D_FLUID3_PRO
	  else if (actele->eltyp==el_fluid3_pro)
	    f3pro_caleleres(actele,&eforce_global,ipos,&hasdirich,&hasext);
#endif
	  else
	    dserror("element type %d unsupported",actele->eltyp);
	}
	else
	{
	  if (actele->eltyp==el_fluid3)
	    f3_caleleres_relax(actele,&estif_global,&eforce_global,ipos,&hasdirich,&hasext);
#ifdef D_FLUID3_IS
	  else if (actele->eltyp==el_fluid3_is)
	    f3is_caleleres_relax(actele,&estif_global,&eforce_global,ipos,&hasdirich,&hasext);
#endif
#ifdef D_FLUID3_PRO
	  else if (actele->eltyp==el_fluid3_pro)
	    dserror("no steepest decent with projection fsi right now");
#endif
	  else
	    dserror("element type %d unsupported",actele->eltyp);
	}

        dof = 0;
        for(j=0; j<actele->numnp; j++)
        {
          actnode = actele->node[j];

          if (force_on_node[j]!=-1)
          {
#ifdef PARALLEL
            dofx = actnode->dof[0];
            dofy = actnode->dof[1];
            dofz = actnode->dof[2];
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
            fcouple[dofx] += eforce[dof  ]*rho;
            fcouple[dofy] += eforce[dof+1]*rho;
            fcouple[dofz] += eforce[dof+2]*rho;
#else
            if(dofx<numeq_total)
              dserror("can not fsi couple free dof!\n");
            if(dofy<numeq_total)
              dserror("can not fsi couple free dof!\n");
            if(dofz<numeq_total)
              dserror("can not fsi couple free dof!\n");

            fcouple[dofx-numeq_total] += eforce[dof  ]*rho;
            fcouple[dofy-numeq_total] += eforce[dof+1]*rho;
            fcouple[dofz-numeq_total] += eforce[dof+2]*rho;
#endif /* SOLVE_DIRICH */
#else /* the sequential case: */
            actnode->sol_mf.a.da[ipos->mf_forcenp][0] += eforce[dof  ]*rho;
            actnode->sol_mf.a.da[ipos->mf_forcenp][1] += eforce[dof+1]*rho;
            actnode->sol_mf.a.da[ipos->mf_forcenp][2] += eforce[dof+2]*rho;
#endif  /* PARALLEL */
          }
          dof += actnode->numdf;
        }
      }
#endif



#ifdef D_FLUID3_F

      rho = mat[actpdis->element[0]->mat-1].m.fluid->density;

      for (i=0; i<actpdis->num_fele; i++)
      {
        act_fast_eles = &(actpdis->fast_eles[i]);

        switch(act_fast_eles->fast_ele_typ)
        {
          case fele_f3f_hex8_e:
          case fele_f3f_hex8_a:
          case fele_f3f_hex20_e:
          case fele_f3f_hex20_a:
          case fele_f3f_tet4_e:
          case fele_f3f_tet4_a:

            aloopl = act_fast_eles->aloopl;
            coupled = 0;
            l = 0;

            /* look if there is a fsi coupled element within the set */
            while ( coupled == 0 && l < aloopl )
            {
              actele = act_fast_eles->ele_vec[l];
              if (actele->e.f3->force_on!=0)
                coupled++;
              l++;
            }

            if (!coupled)
              continue; /* goto next set if no coupling element */


            /* get force vector */
            f3fcaleleres(act_fast_eles->ele_vec,&eforce_fast,hasext_f,ipos,aloopl);


            /* element set force vector */
            eforce = eforce_fast.a.dv;


            /* loop all elements of this set */
            for(l=0; l<aloopl; l++)
            {

              actele = act_fast_eles->ele_vec[l];

              if (actele->e.f3->force_on==0)
                continue; /* element not of interest*/


              nfnode = 0;

              /* do some initialisation */
              for(j=0; j<MAXNOD_F3; j++)
                force_on_node[j] = -1;


              for(j=0; j<actele->numnp; j++)
              {
                actgnode = actele->node[j]->gnode;
                /* check if there is a coupled struct node */
                if (actgnode->mfcpnode[genprob.numsf]==NULL)
                  continue;
                force_on_node[nfnode] = j;
                nfnode++;
              }


              for(j=0; j<nfnode; j++) /* loop coupled nodes of this ele */
              {
                actnode = actele->node[force_on_node[j]];
                line = force_on_node[j] * 4;
#ifdef PARALLEL
                dofx = actnode->dof[0];
                dofy = actnode->dof[1];
                dofz = actnode->dof[2];
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
                fcouple[dofx] += eforce[line*LOOPL+l]*rho;
                fcouple[dofy] += eforce[line*LOOPL+LOOPL+l]*rho;
                fcouple[dofz] += eforce[line*LOOPL+LOOPL*2+l]*rho;
#else
                if(dofx<numeq_total)
                  dserror("can not fsi couple free dof!\n");
                if(dofy<numeq_total)
                  dserror("can not fsi couple free dof!\n");
                if(dofz<numeq_total)
                  dserror("can not fsi couple free dof!\n");

                fcouple[dofx-numeq_total] += eforce[line*LOOPL+l]*rho;
                fcouple[dofy-numeq_total] += eforce[line*LOOPL+LOOPL+l]*rho;
                fcouple[dofz-numeq_total] += eforce[line*LOOPL+LOOPL*2+l]*rho;
#endif /* SOLVE_DIRICH */
#else /* the sequential case: */
                actnode->sol_mf.a.da[ipos->mf_forcenp][0] += eforce[line*LOOPL+l]*rho;
                actnode->sol_mf.a.da[ipos->mf_forcenp][1] += eforce[line*LOOPL+LOOPL+l]*rho;
                actnode->sol_mf.a.da[ipos->mf_forcenp][2] += eforce[line*LOOPL+LOOPL*2+l]*rho;
#endif  /* PARALLEL */

              }  /* for(j=0; j<nfnode; j++) */

            }  /* for(l=0; l<aloopl; l++) */

            break;


          default:
            dserror("unknown typ of fast element");
            break;

        } /* switch(act_fast_eles->fast_ele_typ) */

      }    /* end loop over element sets */

#endif
      break;
    default:
      dserror("genprob->ndim not valid");
  }

end:


#ifdef DEBUG
  dstrc_exit();
#endif


  return;
}







/*!---------------------------------------------------------------------
\brief allreduce consistent nodal fluid forces and write it to nodes

<pre>                                                        chfoe 05/05

This routine exits in parallel only. It performs an allreduce on the
fcouple vector containing the consistent nodal fluid forces. The vector
is of different size depending on wether or not SOLVE_DIRICH is included.

I. with SOLVE_DIRICH:
                 -> numdof = numeq since Dirichlet dofs are not sortet out
                 -> fcouple is a vector of full length indicated by dof
                 -> fcouple is allreduced here

II. without SOLVE_DIRICH:
                 -> fcouple is of dimension (numdof - numeq) i.e. number
                    of Dirichlet dofs of the field.
                 -> fcouple is indicated by (dof(node) - numeq)
                 -> fcouple is allreduced here
</pre>

\param   *actpart        PARTITION    (i)   the structural partition
\param   *fsiforce       DOUBLE       (o)   full vector to be filled
\param    global_numeq   INT          (i)   the structural number of eqs


\return void
\warning This routine exists in parallel only!

------------------------------------------------------------------------*/
#ifdef PARALLEL
void fsi_allreduce_coupforce(
    DOUBLE             *fcouple,
    DOUBLE             *recvfcouple,
    INT                 numeq_total,
    INT                 numddof,
    INTRA              *actintra,
    FIELD              *actfield,
    INT                 disnum
    )

{

  INT  i, j;
  INT  dof;

  NODE  *actnode;
  GNODE *actgnode;

  ARRAY_POSITION *ipos;

#ifdef DEBUG
  dstrc_enter("fsi_allreduce_coupforce");
#endif

  ipos = &(actfield->dis[disnum].ipos);

  /* I. global solution includes Dirchlet dofs */
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
  /* fcouple and recvfcouple are of length 'numeq_total' */
  MPI_Allreduce(fcouple,recvfcouple,numeq_total,MPI_DOUBLE,MPI_SUM,
      actintra->MPI_INTRA_COMM);
  for (i=0; i<actfield->dis[disnum].numnp; i++)
  {
    actnode  = &(actfield->dis[disnum].node[i]);
    actgnode = actnode->gnode;
    if(actgnode->fsicouple == NULL) continue; /* do nothing for uncoupled nodes */
    for (j=0; j<genprob.ndim; j++)
    {
      dof = actnode->dof[j];
      /* recvfcouple indicated by dof */
      actnode->sol_mf.a.da[ipos->mf_forcenp][j] = recvfcouple[dof];
    }
  }
  /* II. global solution without Dirchlet dofs */
#else
  /* fcouple and recvfcouple are of length 'numddof' = numdf - numeq_total */
  MPI_Allreduce(fcouple,recvfcouple,numddof,MPI_DOUBLE,MPI_SUM,
      actintra->MPI_INTRA_COMM);
  for (i=0; i<actfield->dis[disnum].numnp; i++)
  {
    actnode  = &(actfield->dis[disnum].node[i]);
    actgnode = actnode->gnode;
    if(actgnode->fsicouple == NULL) continue; /* do nothing for uncoupled nodes */
    for (j=0; j<genprob.ndim; j++)
    {
      dof = actnode->dof[j];
      if (dof < numeq_total)
        dserror("something is wrong!\n");
      /* recvfcouple indicated by dof-numeq_total */
      actnode->sol_mf.a.da[ipos->mf_forcenp][j] = recvfcouple[dof-numeq_total];
    }
  }
#endif /* SOLVE_DIRICH */


#ifdef DEBUG
  dstrc_exit();
#endif

  return;

}
#endif /* PARALLEL */



/*!---------------------------------------------------------------------
\brief put fsi coupling forces on fullvector on struct. field

<pre>                                                        chfoe 05/05

This routine serves to sort the consistent nodal fluid forces which have
been evaluated on the fluid field into a full vector on the structural
field.
The forces are sorted into fsiforce by their global dof number.

</pre>

\param   *actpart        PARTITION    (i)   the structural partition
\param   *fsiforce       DOUBLE       (o)   full vector to be filled
\param    global_numeq   INT          (i)   the structural number of eqs


\return void

------------------------------------------------------------------------*/
void fsi_load(
    PARTITION          *actpart,
    INT                 disnum,
    FIELD              *fluidfield,
    INT                 fdisnum,
    DOUBLE             *fsiforce,
    INT                 global_numeq
    )
{
  INT  i,j;
  INT  numnp;  /* number of nodes of this pdis */
  INT  dof;

  NODE  *actfnode;
  GNODE *actsgnode;

  ARRAY_POSITION *fluid_ipos;

#ifdef FSI_NONMATCH
  ARRAY        funct_a;
  DOUBLE       *funct;
  DIS_TYP      typ;
  INT          k,m,n;       /*counters*/
  GNODE        *actfgnode;
  DOUBLE       r,s,t;
  NODE         *actsnode;
  INT          g,h;
  DOUBLE       structcoord[3];
  ELEMENT      *hostele;
#endif

#ifdef DEBUG
  dstrc_enter("fsi_load");
#endif

#ifdef FSI_NONMATCH

  /* initialisation */
  /*number of nodes on this processor*/
  numnp = actpart->pdis[disnum].numnp;

  fluid_ipos = &(fluidfield->dis[fdisnum].ipos);

  /*allocate space for the evaluation of the shape functions*/
  funct = amdef("funct"  ,&funct_a  ,MAXNOD_F3,1,"DV");

  /* To project the fluid nodal forces on the structure, a loop over
   * all fluid nodes is performed. Since the position of all FSI fluid
   * nodes in relation to their structure host elements is known
   * (GNODE->coupleptr), the shape functions are used to interpolate
   * the fluid loads to the structure nodes.  The sum of all
   * contributions from fluid nodes delivers the structure nodal
   * loads.
   *
   * _               _          _
   * fsiforce = Sum( N(r,s,t) * ffluid )
   */

  /* loop all fluid nodes */
  /* We rely on the allreduce of the fluid interface forces in a
   * parallel setting. */
  for (i=0;i<fluidfield->dis[fdisnum].numnp;i++)
  {
    actfnode = &(fluidfield->dis[fdisnum].node[i]);
    actfgnode = actfnode->gnode;

    /* node is on FSI interface */
    if (actfgnode->fsicouple!=NULL)
    {
      /*access the coupled structure host element*/
      hostele=actfgnode->coupleptr->hostele;
      typ=hostele->distyp;

      switch (hostele->eltyp)
      {
#ifdef D_BRICK1
      case el_brick1:

	/*access the local coordinates of the fluid node*/
	r=actfgnode->coupleptr->x[0];
	s=actfgnode->coupleptr->x[1];
	t=actfgnode->coupleptr->x[2];

	/*evaluate shape functions at the position of the fluid node*/
	f3_hex(funct,NULL,NULL,r,s,t,typ,0);

	/*loop all struct nodes belonging to this element*/
	for (k=0;k<hostele->numnp;k++)
	{
	  actsnode=hostele->node[k];
	  actsgnode=actsnode->gnode;

	  /* loop nodal dofs */
	  for(m=0;m<actsnode->numdf;m++)
	  {
	    /*get the global dof number*/
	    dof = actsnode->dof[m];

	    if (dof>= global_numeq) continue; /* falls es sich um eine
					       * Dirichlet RB handelt
					       * muß keine Kraft
					       * berechnet werden*/
	    if (par.myrank == actsnode->proc) /* falls der Knoten auf
					       * meinem Prozessor
					       * liegt*/
	    {
	      /* calculate the struct nodal force by interpolation via
	       * the shape functions*/
	      fsiforce[dof] += (actfnode->sol_mf.a.da[fluid_ipos->mf_forcenp][m])*(funct[k]);
	    }
	  }
	}
	break;
#endif
#ifdef D_SHELL8
      case el_shell8:
      {
	/*access the local coordinates of the fluid node*/
	r=actfgnode->coupleptr->x[0];
	s=actfgnode->coupleptr->x[1];

	/*evaluate shape functions at the position of the fluid node*/
	s8_funct_deriv(funct,NULL,r,s,typ,0);

	/*loop all struct nodes belonging to this element*/
	for (k=0;k<hostele->numnp;k++)
	{
	  actsnode=hostele->node[k];
	  actsgnode=actsnode->gnode;

	  /* loop nodal dofs */
	  /* for(m=0;m<actsnode->numdf;m++) */
	  for(m=0;m<3;m++)
	  {
	    /*get the global dof number*/
	    dof = actsnode->dof[m];

	    if (dof>= global_numeq) continue; /* falls es sich um eine
					       * Dirichlet RB handelt
					       * muß keine Kraft
					       * berechnet werden*/
	    if (par.myrank == actsnode->proc) /* falls der Knoten auf
					       * meinem Prozessor
					       * liegt*/
	    {
	      /* calculate the struct nodal force by interpolation via
	       * the shape functions*/
	      fsiforce[dof] += (actfnode->sol_mf.a.da[fluid_ipos->mf_forcenp][m])*(funct[k]);
	    }
	  }
	}
        break;
      }
#endif
      default:
        dserror("unknown element type %d", hostele->eltyp);
      }
    }
  }/*end of loop over fluid nodes*/

  amdel(&funct_a);

#else

  /* initialisation */
  numnp = actpart->pdis[disnum].numnp;
  fluid_ipos = &(fluidfield->dis[fdisnum].ipos);

  /* loop over all nodes on this proc */
  for (i=0; i<numnp; i++)
  {
    NODE* actnode;
    actnode = actpart->pdis[disnum].node[i];
    actsgnode = actnode->gnode;
    if (actsgnode->fsicouple) /* node is on FSI interface */
    {
      actfnode  = actsgnode->mfcpnode[genprob.numff];
      /*for(j=0; j<actnode->numdf && j<actfnode->numdf-1; j++)*/ /* loop nodal dofs */
      for(j=0; j<actnode->numdf; j++) /* loop nodal dofs */
      {
        dof = actnode->dof[j];
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
        if (actnode->gnode->dirich && actnode->gnode->dirich->dirich_onoff.a.iv[j])
          continue;
#else
        if ( dof>= global_numeq) continue;
#endif
        if (par.myrank == actnode->proc)
          fsiforce[dof] = actfnode->sol_mf.a.da[fluid_ipos->mf_forcenp][j];
      }
    }
  } /* end loop over nodes on proc*/

#endif

#ifdef DEBUG
  dstrc_exit();
#endif

  return;

}
#endif
