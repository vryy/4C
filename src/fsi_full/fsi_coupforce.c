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
/*#include "fluid_prototypes.h"*/
#include "../fluid2/fluid2.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../fluid3/fluid3.h"
#include "../fluid3/fluid3_prototypes.h"
#include "../fluid3_fast/f3f_prototypes.h"
#ifdef D_ALE
#include "../ale2/ale2.h"
#endif

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
void fsi_cbf(PARTDISCRET    *actpdis,
	     DOUBLE         *fcouple,
             ARRAY_POSITION *ipos,
             INT             numeq_total,
             INT             init)
{
INT i, j, l;
INT line;
INT numdf;
INT hasdirich, hasext;
INT force_on_node[MAXNOD];
INT nfnode;  /* number of nodes of actele where forces are searched for */
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

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("fsi_cbf");
#endif

/*--------------------------------------------------- initialisation ---*/
if (init==1)
{
   for(j=0; j<MAXNOD; j++) /* ensure that MAXNOD >= MAX(MAXNOD_F2, MAXNOD_F3)! */
      force_on_node[j] = -1;
   /*--- preselect elements that have fsi coupling ---*/
   switch (genprob.ndim)
   {
   case 2:
#ifdef D_FLUID2
      for(i=0; i<actpdis->numele; i++)
      {
         actele = actpdis->element[i];
         if (par.myrank!=actele->proc) continue; /* only my elments */
         for(j=0; j<actele->numnp; j++)
         {
            actgnode = actele->node[j]->gnode;
            /* check if there is a coupled struct node */
            if (actgnode->mfcpnode[genprob.numsf]==NULL)
               continue;
            actele->e.f2->force_on = 1;
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
         if (par.myrank!=actele->proc) continue; /* only my elments */
         for(j=0; j<actele->numnp; j++)
         {
            actgnode = actele->node[j]->gnode;
            /* check if there is a coupled struct node */
            if (actgnode->mfcpnode[genprob.numsf]==NULL)
               continue;
            actele->e.f3->force_on = 1;
            break;
         }
      }      
#endif
   break;
   }
   goto end;
}

/*--------------------------------------- the elemental force vector ---*/
eforce  = eforce_global.a.dv;

/*-------------- evaluate fsi coupling forces depending on dimension ---*/
switch (genprob.ndim)
{
case 2: /* problem is two-dimensional */
#ifdef D_FLUID2
   rho = mat[actpdis->element[0]->mat-1].m.fluid->density;
   numdf = 3;
   for(i=0; i<actpdis->numele; i++)
   { 
      actele = actpdis->element[i];

     /*------------------------------------------------------------------*/
      if (actele->e.f2->force_on==0) continue; /* element not of interest*/
      nfnode = 0;
      for(j=0; j<MAXNOD_F2; j++)
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
      /*--- get force vector ---*/
      f2_caleleres(actele,&eforce_global,ipos,&hasdirich,&hasext);

      for(j=0; j<nfnode; j++)
      {
         actnode = actele->node[force_on_node[j]];
         line = force_on_node[j] * 3;
#ifdef PARALLEL
         dofx = actnode->dof[0];
         dofy = actnode->dof[1];
 #if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
         fcouple[dofx] += eforce[line]*rho;
         fcouple[dofy] += eforce[line+1]*rho;
 #else
         if(dofx<numeq_total)
            dserror("can not fsi couple free dof!\n");
         if(dofy<numeq_total)
            dserror("can not fsi couple free dof!\n");

         fcouple[dofx-numeq_total] += eforce[line]*rho;
         fcouple[dofy-numeq_total] += eforce[line+1]*rho;
 #endif /* SOLVE_DIRICH */
#else
         actnode->sol_mf.a.da[1][0] += eforce[line]*rho;
         actnode->sol_mf.a.da[1][1] += eforce[line+1]*rho;
#endif  /* PARALLEL */
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
      actele = actpdis->element[i];
      /*------------------------ fast elements will be treated later ---*/
      if (actele->eltyp == el_fluid3_fast) continue;

      /*----------------------------------------------------------------*/
      if (actele->e.f3->force_on==0) continue; /* element not of interest*/
      nfnode = 0;
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

      /*--- get force vector ---*/
      f3_caleleres(actele,&eforce_global,&hasdirich,&hasext,ipos);

      for(j=0; j<nfnode; j++)
      {
         actnode = actele->node[force_on_node[j]];
         line = force_on_node[j] * 4;
#ifdef PARALLEL
         dofx = actnode->dof[0];
         dofy = actnode->dof[1];
         dofz = actnode->dof[2];
 #if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
         fcouple[dofx] += eforce[line]*rho;
         fcouple[dofy] += eforce[line+1]*rho;
         fcouple[dofz] += eforce[line+2]*rho;
 #else
         if(dofx<numeq_total)
            dserror("can not fsi couple free dof!\n");
         if(dofy<numeq_total)
            dserror("can not fsi couple free dof!\n");
         if(dofz<numeq_total)
            dserror("can not fsi couple free dof!\n");

         fcouple[dofx-numeq_total] += eforce[line]*rho;
         fcouple[dofy-numeq_total] += eforce[line+1]*rho;
         fcouple[dofz-numeq_total] += eforce[line+2]*rho;
 #endif /* SOLVE_DIRICH */
#else /* the sequential case: */
         actnode->sol_mf.a.da[1][0] += eforce[line]*rho;
         actnode->sol_mf.a.da[1][1] += eforce[line+1]*rho;
         actnode->sol_mf.a.da[1][2] += eforce[line+2]*rho;
#endif  /* PARALLEL */
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
         if (!coupled) continue; /* goto next set if no coupling element */


         /*--- get force vector ---*/
         f3fcaleleres(act_fast_eles->ele_vec,&eforce_fast,hasext_f,ipos,aloopl);
         /* element set force vector */
         eforce = eforce_fast.a.dv;

         /* loop all elements of this set */
         for(l=0; l<aloopl; l++)
         {
            actele = act_fast_eles->ele_vec[l];
            if (actele->e.f3->force_on==0) continue; /* element not of interest*/
            nfnode = 0;

            /* do some initialisation */
            for(j=0; j<MAXNOD_F3; j++)
               force_on_node[j] = -1;

            for(i=0; i<actele->numnp; i++)
            {
               actgnode = actele->node[i]->gnode;
               /* check if there is a coupled struct node */
               if (actgnode->mfcpnode[genprob.numsf]==NULL)
                  continue;
               force_on_node[nfnode] = i;
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
               actnode->sol_mf.a.da[1][0] += eforce[line*LOOPL+l]*rho;
               actnode->sol_mf.a.da[1][1] += eforce[line*LOOPL+LOOPL+l]*rho;
               actnode->sol_mf.a.da[1][2] += eforce[line*LOOPL+LOOPL*2+l]*rho;
#endif  /* PARALLEL */
            }
         } /* end the loop */
      break;
      default:
         dserror("unknown typ of fast element");
      } /* switch(act_fast_eles->fast_ele_typ) */
   }    /* end loop over element sets */
#endif
break;
default:
   dserror("genprob->ndim not valid");
}

end:
/*----------------------------------------------------------------------*/
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
void fsi_allreduce_coupforce( DOUBLE *fcouple,
                              DOUBLE *recvfcouple,
                              INT     numeq_total,
                              INT     numddof,
                              INTRA  *actintra,
                              FIELD  *actfield
                            )
{
INT  i, j;
INT  dof;

NODE  *actnode;
GNODE *actgnode;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("fsi_allreduce_coupforce");
#endif

/* I. global solution includes Dirchlet dofs */
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
/* fcouple and recvfcouple are of length 'numeq_total' */
MPI_Allreduce(fcouple,recvfcouple,numeq_total,MPI_DOUBLE,MPI_SUM,
              actintra->MPI_INTRA_COMM);
for (i=0; i<actfield->dis[0].numnp; i++)
{
   actnode  = &(actfield->dis[0].node[i]);
   actgnode = actnode->gnode;
   if(actgnode->fsicouple == NULL) continue; /* do nothing for uncoupled nodes */
   for (j=0; j<genprob.ndim; j++)
   {
      dof = actnode->dof[j];
      /* recvfcouple indicated by dof */
      actnode->sol_mf.a.da[1][j] = recvfcouple[dof];
   }
}
/* II. global solution without Dirchlet dofs */
#else
/* fcouple and recvfcouple are of length 'numddof' = numdf - numeq_total */
MPI_Allreduce(fcouple,recvfcouple,numddof,MPI_DOUBLE,MPI_SUM,
              actintra->MPI_INTRA_COMM);
for (i=0; i<actfield->dis[0].numnp; i++)
{
   actnode  = &(actfield->dis[0].node[i]);
   actgnode = actnode->gnode;
   if(actgnode->fsicouple == NULL) continue; /* do nothing for uncoupled nodes */
   for (j=0; j<genprob.ndim; j++)
   {
      dof = actnode->dof[j];
      if (dof < numeq_total)
         dserror("something is wrong!\n");
      /* recvfcouple indicated by dof-numeq_total */
      actnode->sol_mf.a.da[1][j] = recvfcouple[dof-numeq_total];
   }
}
#endif /* SOLVE_DIRICH */

/*----------------------------------------------------------------------*/
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
void fsi_load(PARTITION *actpart, DOUBLE *fsiforce, INT global_numeq)
{
INT  i,j;
INT  numnp;  /* number of nodes of this pdis */
INT  dof;

NODE  *actnode, *actfnode;
GNODE *actsgnode;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("fsi_load");
#endif

/*--------------------------------------------------- initialisation ---*/
numnp = actpart->pdis->numnp;

/*--------------------------------- loop over all nodes on this proc ---*/
for (i=0; i<numnp; i++)
{
   actnode = actpart->pdis->node[i];
   actsgnode = actnode->gnode;
   if (actsgnode->fsicouple) /* node is on FSI interface */
   {
      for(j=0; j<actnode->numdf; j++) /* loop nodal dofs */
      {
         actfnode  = actsgnode->mfcpnode[genprob.numff];
         dof = actnode->dof[j];
         if ( dof>= global_numeq) continue;
         if (par.myrank == actnode->proc)
            fsiforce[dof] = actfnode->sol_mf.a.da[1][j];
      }
   }
} /* end loop over nodes on proc*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
return;
}

#endif
