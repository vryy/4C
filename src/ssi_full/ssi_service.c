/*!----------------------------------------------------------------------
\file
\brief service routines for ssi algorithms

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*! 
\addtogroup SSI
*//*! @{ (documentation module open)*/
#ifdef D_SSI
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "ssi_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par; 
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
 | a dense matrix                                        m.gee 11/01    |
 | this structure holds a dense matrix to be solved with lapack         |
 *----------------------------------------------------------------------*/
extern struct _DENSE lhs_dens;

/*!----------------------------------------------------------------------*
\brief coupling algorithm for ssi problems

<pre>
                                                   firl / genk 10/03    
                                                                       
  routine to transfer displacements at time t from the master field to 
  the slave field in ssi problems, the displacements are directly written
  on the nodes, only possible for conforming meshes 
  and two field problems;
  WARNING! With this approach Dirichlet b.c. on coupling nodes are over-
           written, this is not physically consistent !       
  
  The routine is placed in ssi_service.c
                                                                       
\param  *actfield  FIELD      (i) pointer to the field under consideration                   

</pre>                                                                       
                                                                       
------------------------------------------------------------------------*/
void ssiserv_put_disp2slave(FIELD *actfield, CONTAINER *container, 
                            DOUBLE *rldfac, SSI_DYNAMIC *ssidyn)
{
DISCRET *actdis;
NODE *actsnode;                     /* used for the actual slave node */
NODE *actmnode;                     /* used for the actual master node */
DIRICH_CONDITION *dirich;
INT i, j;                           /* counters */

#ifdef DEBUG 
dstrc_enter("ssiserv_put_disp2node");
#endif

actdis = &(actfield->dis[0]);
/* loop nodes of slave field */
for (i=0; i<actdis->numnp; i++)
{
  actsnode = &(actdis->node[i]);
  /* do nothing if the node is not a coupling node */
  if (actsnode->gnode->ssicouple == NULL) continue;
  if(actsnode->gnode->ssicouple->ssi_couptyp == ssi_slave)
  {
    if(container->coupl_typ == 0)
    {
    /* set pointer to the corresponding master node */
    actmnode = actsnode->gnode->mfcpnode[0];
    dsassert(actmnode != NULL,"Pointer to master node not properly set!");

    /* set pointer for dirichlet boundary condition */
    dirich = actsnode->gnode->dirich;
    dsassert(dirich != NULL,"Pointer to dirichlet condition for slave node \
             not properly set!");

     /* put values from actmnode->sol_mf.a.da[0][j] to 
        dirich->dirich_val.a.dv[j] */
     for (j=0; j<actsnode->numdf; j++)
       dirich->dirich_val.a.dv[j] =((-1.0 * container->dirichfacs[6] * 
                                   actmnode->sol.a.da[9][j]) + 
                                   ((1.0 + container->dirichfacs[6]) * 
                                   (ssidyn->relax * actmnode->sol_mf.a.da[0][j] + 
                                   (1.0 - ssidyn->relax) * actmnode->sol_mf.a.da[1][j])));
                                   /* end of loop over numdf */
  } /* end of if clause (if(container->coupl_typ == 0) )*/
  else /*nonconforming discretization */
  {
    /* set pointer for dirichlet boundary condition */
    dirich = actsnode->gnode->dirich;
    dsassert(dirich != NULL,"Pointer to dirichlet condition for slave node \
             not properly set!");

     /* put values from actsnode->sol_mf.a.da[0][j] to 
        dirich->dirich_val.a.dv[j] */
     for (j=0; j<actsnode->numdf; j++)
     {
     /* new approach 08/03/04 the parameter rldfac was deleted in subroutine 
         solserv_putdirich_to_dof; hence the dirichlet condition is not divided
         by rldfac anymore  */

      dirich->dirich_val.a.dv[j] =((-1.0 * container->dirichfacs[6] * 
                                 actsnode->sol.a.da[9][j]) + 
                                 ((1.0 + container->dirichfacs[6]) * 
                                 actsnode->sol_mf.a.da[0][j]));
    } /* end of loop over numdf */
  }
  } /* ens of if clause ssi_slave*/
} /* end of loop over numnp */
#ifdef DEBUG 
dstrc_exit();
#endif
} /* end of ssiserv_put_disp2slave */



/*!----------------------------------------------------------------------*
\brief put the displacements at time t to the coupling nodes on the slave 
       side

<pre>
                                                   firl / genk 11/03    
                                                                       
  routine to transfer displacements at time t from the master field to 
  the slave field in ssi problems, the displacements are directly written
  on the nodes, only possible for conforming meshes 
  and two field problems;
  WARNING! With this approach Dirichlet b.c. on coupling nodes are over-
           written, this is not physically consistent !       
  This is necessary to obtain correct stresses and displacments in the 
  output file.
  
  
  The routine is placed in ssi_service.c
                                                                       
\param  *actfield  FIELD      (i) pointer to the field under consideration                   

</pre>                                                                       
                                                                       
------------------------------------------------------------------------*/
void ssiserv_put_true_disp2slave(FIELD *actfield, CONTAINER *container)
{
  DISCRET *actdis;
  NODE *actsnode;                     /* used for the actual slave node */
  NODE *actmnode;                     /* used for the actual master node */
  INT i, j;                           /* counters */

  #ifdef DEBUG 
  dstrc_enter("ssiserv_put_true_disp2node");
  #endif
  
  actdis = &(actfield->dis[0]);
  /* loop nodes of slave field */
  for (i=0; i<actdis->numnp; i++)
  {
    actsnode = &(actdis->node[i]);
    /* do nothing if the node is not a coupling node */ 
    if(actsnode->gnode->ssicouple == NULL) continue; 
    if(actsnode->gnode->ssicouple->ssi_couptyp == ssi_slave)
    {
      if(container->coupl_typ == 0) /* conforming discretization */
      {
        /* set pointer to the corresponding master node */
        actmnode = actsnode->gnode->mfcpnode[0];
        dsassert(actmnode != NULL,"Pointer to master node not properly set!");

        /* put values from actmnode->sol_mf.a.da[0][j] to 
           actsnode->sol.a.da[0][j] */
        for (j=0; j<actsnode->numdf; j++)
        {
	  actsnode->sol.a.da[0][j] = actmnode->sol.a.da[9][j]; 
        }
      }
      else /* non-conforming discretization */
      {
        /* put values from actsnode->sol_mf.a.da[0][j] to 
           actsnode->sol.a.da[0][j] */
        for (j=0; j<actsnode->numdf; j++)
        {
	  actsnode->sol.a.da[0][j] = actsnode->sol_mf.a.da[6][j];
	  actsnode->sol_mf.a.da[0][j] = actsnode->sol_mf.a.da[6][j]; 
	  /*actsnode->sol.a.da[0][j] = actsnode->sol_mf.a.da[1][j];*/ 
        }
      }
    }
  }
#ifdef DEBUG 
dstrc_exit();
#endif
} /* end of ssiserv_put_true_disp2slave */


/*!---------------------------------------------------------------------                                         
\brief compute internal forces for slave nodes on the interface

<pre>                                                  firl / genk 2/04

Here are computed the internal forces in the slave elements at the 
interface. These values are transferred to the master nodes on the 
interface. This subroutine is used for conforming desretizations as well
as for non-conforming discretizations.

</pre>
\param *actele        ELEMENT	     (i)   actual slave element
\param *estif_global  ARRAY          (i)   stiffness matrix of act. elem.
\param *emass_global  ARRAY          (i)   mass matrix of act. ele.
\param *container     CONTAINER      (i)   the container
\return 0                                                                             

------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 |                                                      firl / genk 02/04 |
 |  routine based on assemble_dirich_dyn written by m.gee 3/02               |
 |  dirichlet conditions to an element vector elevec_a                       |
 |  and then assembles this element vector of cond. dirich.conditions to the |
 |  global vector container->coup_force                                      |
 |  this assembly is a special design for structural dynamics with           |
 |  generalized alfa method                                                  |
 |                                                                           |
 |  facs[0] = -(1.0-alpham)*(1.0/beta)/(DSQR(dt))                            |
 |  facs[1] =  (1.0-alpham)*(1.0/beta)/dt                                    |
 |  facs[2] =  (1.0-alpham)/(2*beta) - 1                                     |
 |  facs[3] = -(1.0-alphaf)*(gamma/beta)/dt                                  |
 |  facs[4] =  (1.0-alphaf)*gamma/beta - 1                                   |
 |  facs[5] =  (gamma/(2*beta)-1)*(1.0-alphaf)*dt                               |
 |  facs[6] = -(1.0-alphaf)                                              |
 |  facs[7] =  raleigh damping factor for mass                               |
 |  facs[8] =  raleigh damping factor for stiffness                          |
 |  facs[9] =  dt                                                            |
 |                                                                           |
 | in the nodes the results are stored the following way:                    |
 |                                                                           |
 | results from sol_increment.a.da[0][0..numdf-1]                            |
 | place 0 holds converged incremental displacements (without prescribed dofs) 
 |                                                                           |
 | in ARRAY sol.a.da[place][0..numdf-1]:                                     |
 | place 0 holds total displacements of free dofs at time t                  |
 | place 1 holds velocities at time t                                        |
 | place 2 holds accels at time t                                            |
 | place 3 holds prescribed displacements at time t-dt                       |
 | place 4 holds prescribed displacements at time t                          |
 | place 5 holds place 4 - place 3                                           |
 | place 6 holds the  velocities of prescribed dofs                          |
 | place 7 holds the  accels of prescribed dofs                              |
 | place 8 is working space                                                  |
 |                                                                           |
 |  container->coup_forces contains the rhs mass, damping and stiffness      |
 |  parts due to prescribed displacements from the master field              |
 |  (dirichlet conditions !=0)                                               |
 |                                                                           |
 | see PhD theses Mok page 165                                               |
 *---------------------------------------------------------------------------*/
void calc_ssi_coupforce_mod(ELEMENT *actele, ARRAY *estif_global, 
                         ARRAY *emass_global, CONTAINER *container)
{
INT                   i,j;
INT                   counter,hasdirich;
INT                   numdf;
INT                   nd=0;
INT                   idamp=0;
DOUBLE                mdamp = 0.0;
DOUBLE                kdamp = 0.0;
DOUBLE              **estif;
DOUBLE              **emass;
DOUBLE                dirich[MAXDOFPERELE];
DOUBLE                dforces[MAXDOFPERELE];
DOUBLE                disp_incre[MAXDOFPERELE];
INT                   dirich_onoff[MAXDOFPERELE];
INT                   lm[MAXDOFPERELE];
GNODE                *actgnode;
NODE                 *actmnode, *actsnode;

#ifdef DEBUG 
dstrc_enter("calc_ssi_coupforce");
#endif
    
/*----------------------------------------------------------------------*/
/*-------- check presence of any ssi coupling condition to this element */
hasdirich=0;
for (i=0; i<actele->numnp; i++)
{
   if (actele->node[i]->gnode->dirich == NULL) continue;
   if (actele->node[i]->gnode->dirich->dirich_type == dirich_SSI) 
   {
      hasdirich=1;
      break;
   }
}
/*--------------------- there are no dirichlet conditions here so leave */
if (hasdirich==0) goto end;
/*--------------------------------------- check for presence of damping */

if (ABS(container->dirichfacs[7]) > EPS13 || ABS(container->dirichfacs[8]) > EPS13) 
{
   idamp=1;
   mdamp = container->dirichfacs[7];
   kdamp = container->dirichfacs[8];
}

/*----------------------------------------------------------------------*/
estif  = estif_global->a.da;
emass  = emass_global->a.da;
/*---------------------------------- set number of dofs on this element */
for (i=0; i<actele->numnp; i++) nd += actele->node[i]->numdf;
/*---------------------------- init the vectors dirich and dirich_onoff */
for (i=0; i<nd; i++)
{
   dirich[i] = 0.0;
   dforces[i] = 0.0;
   disp_incre[i] = 0.0;
   dirich_onoff[i] = 0;
}
/*-------------------------------- fill vectors dirich and dirich_onoff */
for (i=0; i<actele->numnp; i++)
{
   numdf    = actele->node[i]->numdf;
   actgnode = actele->node[i]->gnode;
   for (j=0; j<numdf; j++)
   {
      lm[i*numdf+j] = actele->node[i]->dof[j];
      if (actgnode->dirich==NULL) continue;
      dirich_onoff[i*numdf+j] = actgnode->dirich->dirich_onoff.a.iv[j];
   }
}


/*----------------------------------------------------------------------*/
/*------------------- make the entry (1 - alpha_f) * K * (u(t)-u(t-dt)) */
/*------------------------------------ = -facs[6] * K * (sol[0]-sol[9]) */
/*----------------------------------------------------------------------*/
/*-----  actual displacements u(t) are in node->sol.a.da[0][0..numdf-1] */
/*------------------------- old displacements are stored in sol.a.da[9] */
counter=0;
/*---------------------------------- loop over all nodes of the element */
for (i=0; i<actele->numnp; i++)
{
  actsnode = actele->node[i];
  /*---------- check if the node under consideration is a coupling node */
  if (actsnode->gnode->ssicouple == NULL)  /* free node */ 
  {
    if (actsnode->gnode->dirich == NULL) 
    {
      for (j=0; j<actsnode->numdf; j++)
      {
        dirich[counter] = actsnode->sol.a.da[0][j] - 
                          actsnode->sol.a.da[9][j];
        counter++;
      }
    }
    /* does the following ever happen??? */
    else if (actsnode->gnode->dirich->dirich_onoff.a.iv[counter] == 0)
    {
      for (j=0; j<actsnode->numdf; j++)
      {
        dirich[counter] = actsnode->sol.a.da[0][j] - 
                          actsnode->sol.a.da[9][j];
        counter++;
      }
    }
  }
  else if (actsnode->gnode->ssicouple->ssi_couptyp == ssi_slave) /* slave node */
  {
    if(container->coupl_typ == 0) /* conforming discretization */
    {  
      /*------------------------------ set the corresponing master node */
      actmnode = actsnode->gnode->mfcpnode[0];
      for (j=0; j<actsnode->numdf; j++)
      {
        /* the displacements at the master nodes are not scaled with the */
        /* relaxation factor omega up to now; Hence this is done here to */
        /* obtain consistent forces; it is assumed the there are the actual */
        /* displacements in sol_mf.a.da[0][] and the old ones in sol_mf.a.da[1][]*/ 
        dirich[counter] = ((container->relax_param * 
                           actmnode->sol_mf.a.da[0][j]) +
                          (1-container->relax_param) * 
                           actmnode->sol_mf.a.da[1][j]) - 
                           actmnode->sol.a.da[9][j];
        counter++;
      }
    }
    else  /* non-conforming discretization */
    {
      for (j=0; j<actsnode->numdf; j++)
      {
        /*dirich[counter] = actsnode->sol_mf.a.da[0][j];*/
        dirich[counter] = actsnode->sol_mf.a.da[6][j]-
                          actsnode->sol.a.da[9][j];
        counter++;
      }
    }
  }
  /*-------------- check if the node under consideration is a free node */
  else if (actsnode->gnode->ssicouple->ssi_couptyp == ssi_none && 
           actsnode->gnode->dirich->dirich_onoff.a.iv[counter] == 0)
  {
    for (j=0; j<actsnode->numdf; j++)
    {
      dirich[counter] = actsnode->sol.a.da[0][j] - 
                        actsnode->sol.a.da[9][j];
      counter++;
    }
  }
}
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      dforces[i] += estif[i][j] * dirich[j] * 
                    (-1.0 * container->dirichfacs[6]);
   }/* loop j over columns */
}/* loop i over rows */

/*----------------------------------------------------------------------*/
/*--------------------------- make the entries M  * h(du,udot,udotdot) */
/* 
this consists of three parts:
   M * -facs[0] * sol_increment.a.da[0][0...numdf-1]  
   M * -facs[1] * sol[1][0..numdf-1]       
   M * -facs[2] * sol[2][0..numdf-1]       

   sol_increment.a.da[0][0...numdf-1]    contains    (u(t)-u(t-dt))
   sol[1][0..numdf-1]                    contains    (udot(t-dt))
   sol[2][0..numdf-1]                    contains    (uddot(t-dt))
*/   
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   M * -facs[0] * sol_increment.a.da[0][0..numdf-1] 
   sol_increment.a.da[0] holds (u(t) - u(t-dt))
*/
for (i=0; i<nd; i++) dirich[i]=0.0;
counter=0;
/*----------------------------------- loop over all node of the element */
for (i=0; i<actele->numnp; i++)
{
  actsnode = actele->node[i];
  /*---------- check if the node under consideration is a coupling node */
  if (actsnode->gnode->ssicouple == NULL) continue; 
  if (actsnode->gnode->ssicouple->ssi_couptyp == ssi_slave)
  {
    if(container->coupl_typ == 0) /* conforming discretization */
    {  
      /*----------------------------- set the corresponding master node */
      actmnode = actsnode->gnode->mfcpnode[0];
      for (j=0; j<actsnode->numdf; j++)
      {
        dirich[counter] = actmnode->sol_increment.a.da[0][j];
        counter++;
      }
    }
    else  /* non-conforming discretization */
    {
      for (j=0; j<actsnode->numdf; j++)
      {
        dirich[counter] = actsnode->sol.a.da[5][j];
        counter++;
      }
    }
  }
  /*-------------- check if the node under consideration is a free node */
  if (actsnode->gnode->ssicouple->ssi_couptyp == ssi_none && 
      actsnode->gnode->dirich->dirich_onoff.a.iv[counter] == 0)
  {
    for (j=0; j<actsnode->numdf; j++)
    {
      dirich[counter] = actsnode->sol_increment.a.da[0][j];
      counter++;
    }
  }
}
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      dforces[i] += emass[i][j] * dirich[j] * 
                    container->dirichfacs[0] * (-1.0);
   }/* loop j over columns */
}/* loop i over rows */
/*----------------------------------------------------------------------*/
for (i=0; i<nd; i++) dirich[i]=0.0;
counter=0;
/*----------------------------------- loop over all node of the element */
for (i=0; i<actele->numnp; i++)
{
  actsnode = actele->node[i];
  /*---------- check if the node under consideration is a coupling node */
  if (actsnode->gnode->ssicouple == NULL) continue; 
  if (actsnode->gnode->ssicouple->ssi_couptyp == ssi_slave)
  {
    if(container->coupl_typ == 0) /* conforming discretization */
    {  
      /*----------------------------- set the corresponding master node */
      actmnode = actsnode->gnode->mfcpnode[0];
      for (j=0; j<actsnode->numdf; j++)
      {
        dirich[counter] = actmnode->sol.a.da[1][j];
        counter++;
      }
    }
    else  /* non-conforming discretization */
    {
      for (j=0; j<actsnode->numdf; j++)
      {
        dirich[counter] = actsnode->sol.a.da[6][j];
        counter++;
      }
    }
  }
  /*-------------- check if the node under consideration is a free node */
  if (actsnode->gnode->ssicouple->ssi_couptyp == ssi_none && 
      actsnode->gnode->dirich->dirich_onoff.a.iv[counter] == 0)
  {
    for (j=0; j<actsnode->numdf; j++)
    {
      dirich[counter] = actsnode->sol.a.da[1][j];
      counter++;
    }
  }
}
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      dforces[i] += emass[i][j] * dirich[j] * 
                    container->dirichfacs[1] * (-1.0);
   }/* loop j over columns */
}/* loop i over rows */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   M * -facs[2] * sol[2][0..numdf-1]       sol[2] holds uddot(t-dt)
   for the nodes with a dirichlet condition, the values of sol[6] are used
   M * -facs[2] * sol[7][0..numdf-1]       sol[7] holds uddot(t-dt)

*/
for (i=0; i<nd; i++) dirich[i]=0.0;
counter=0;
/*----------------------------------- loop over all node of the element */
for (i=0; i<actele->numnp; i++)
{
  actsnode = actele->node[i];
  /*---------- check if the node under consideration is a coupling node */
  if (actsnode->gnode->ssicouple == NULL) continue; 
  if (actsnode->gnode->ssicouple->ssi_couptyp == ssi_slave)
  {
    if(container->coupl_typ == 0) /* conforming discretization */
    {  
      /*----------------------------- set the corresponding master node */
      actmnode = actsnode->gnode->mfcpnode[0];
      for (j=0; j<actsnode->numdf; j++)
      {
        dirich[counter] = actmnode->sol.a.da[2][j];
        counter++;
      }
    }
    else  /* non-conforming discretization */
    {
      for (j=0; j<actsnode->numdf; j++)
      {
        dirich[counter] = actsnode->sol.a.da[7][j];
        counter++;
      }
    }
  }
  /*-------------- check if the node under consideration is a free node */
  if (actsnode->gnode->ssicouple->ssi_couptyp == ssi_none && 
      actsnode->gnode->dirich->dirich_onoff.a.iv[counter] == 0)
  {
    for (j=0; j<actsnode->numdf; j++)
    {
      dirich[counter] = actsnode->sol.a.da[2][j];
      counter++;
    }
  }
}
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      dforces[i] += emass[i][j] * dirich[j] * 
                    container->dirichfacs[2] * (-1.0);
   }/* loop j over columns */
}/* loop i over rows */



if (idamp)
{

/*----------------------------------------------------------------------*/
/*--------------------------- make the entries C  * e(-du,udot,udotdot) */
/* 
this consists of three parts:
   C * -facs[3] * sol_increment.a.da[0][0..numdf-1]       
   sol_increment.a.da[0] holds (u(t) - u(t-dt))
   C * -facs[4] * sol[1][0..numdf-1]       sol[1] holds udot(t)
   C * -facs[5] * sol[2][0..numdf-1]       sol[2] holds uddot(t)
*/   
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   C * facs[3] * sol_increment.a.da[0][0..numdf-1]       
   sol_increment.a.da[0] holds (u(t) - u(t-dt))
*/
for (i=0; i<nd; i++) dirich[i]=0.0;
counter=0;
/*---------------------------------- loop over all nodes of the element */
for (i=0; i<actele->numnp; i++)
{
  actsnode = actele->node[i];
  /*---------- check if the node under consideration is a coupling node */
  if (actsnode->gnode->ssicouple->ssi_couptyp == ssi_slave)
  {
    if(container->coupl_typ == 0) /* conforming discretization */
    {  
      /*------------------------------ set the corresponing master node */
      actmnode = actsnode->gnode->mfcpnode[0];
      for (j=0; j<actsnode->numdf; j++)
      {
        dirich[counter] = actmnode->sol_increment.a.da[0][j];
        counter++;
      }
    }
    else  /* non-conforming discretization */
    {
      for (j=0; j<actsnode->numdf; j++)
      {
        dirich[counter] = actsnode->sol.a.da[5][j];
        counter++;
      }
    }
  }
  /*-------------- check if the node under consideration is a free node */
  if (actsnode->gnode->ssicouple->ssi_couptyp == ssi_none && 
      actsnode->gnode->dirich->dirich_onoff.a.iv[counter] == 0)
  {
    for (j=0; j<actsnode->numdf; j++)
    {
      dirich[counter] = actsnode->sol_increment.a.da[0][j];
      counter++;
    }
  }
}
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      dforces[i] += (mdamp*emass[i][j]+kdamp*estif[i][j]) * dirich[j] * 
                     container->dirichfacs[3] * (-1);
   }/* loop j over columns */
}/* loop i over rows */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   C * facs[4] * sol[1][0..numdf-1]       sol[1] holds ddot(t)
*/
for (i=0; i<nd; i++) dirich[i]=0.0;
counter=0;
/*----------------------------------- loop over all node of the element */
for (i=0; i<actele->numnp; i++)
{
  actsnode = actele->node[i];
  /*---------- check if the node under consideration is a coupling node */
  if (actsnode->gnode->ssicouple->ssi_couptyp == ssi_slave)
  {
    if(container->coupl_typ == 0) /* conforming discretization */
    {  
      /*------------------------------ set the corresponing master node */
      actmnode = actsnode->gnode->mfcpnode[0];
      for (j=0; j<actsnode->numdf; j++)
      {
        dirich[counter] = actmnode->sol.a.da[1][j];
        counter++;
      }
    }
    else  /* non-conforming discretization */
    {
      for (j=0; j<actsnode->numdf; j++)
      {
        dirich[counter] = actsnode->sol.a.da[1][j];
        counter++;
      }
    }
  }
  /*-------------- check if the node under consideration is a free node */
  if (actsnode->gnode->ssicouple->ssi_couptyp == ssi_none && 
      actsnode->gnode->dirich->dirich_onoff.a.iv[counter] == 0)
  {
    for (j=0; j<actsnode->numdf; j++)
    {
      dirich[counter] = actsnode->sol.a.da[1][j];
      counter++;
    }
  }
}
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      dforces[i] += (mdamp*emass[i][j]+kdamp*estif[i][j]) * dirich[j] * 
                     container->dirichfacs[4] * (-1);
   }/* loop j over columns */
}/* loop i over rows */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   C * facs[5] * sol[2][0..numdf-1]       sol[2] holds ddotdot(t)
*/
for (i=0; i<nd; i++) dirich[i]=0.0;
counter=0;
/*----------------------------------- loop over all node of the element */
for (i=0; i<actele->numnp; i++)
{
  actsnode = actele->node[i];
  /*---------- check if the node under consideration is a coupling node */
  if (actsnode->gnode->ssicouple->ssi_couptyp == ssi_slave)
  {
    if(container->coupl_typ == 0) /* conforming discretization */
    {  
      /*------------------------------ set the corresponing master node */
      actmnode = actsnode->gnode->mfcpnode[0];
      for (j=0; j<actsnode->numdf; j++)
      {
        dirich[counter] = actmnode->sol.a.da[2][j];
        counter++;
      }
    }
    else  /* non-conforming discretization */
    {
      for (j=0; j<actsnode->numdf; j++)
      {
        dirich[counter] = actsnode->sol.a.da[2][j];
        counter++;
      }
    }
  }
  /*-------------- check if the node under consideration is a free node */
  if (actsnode->gnode->ssicouple->ssi_couptyp == ssi_none && 
      actsnode->gnode->dirich->dirich_onoff.a.iv[counter] == 0)
  {
    for (j=0; j<actsnode->numdf; j++)
    {
      dirich[counter] = actsnode->sol.a.da[2][j];
      counter++;
    }
  }
}
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      dforces[i] += (mdamp*emass[i][j]+kdamp*estif[i][j]) * dirich[j] * 
                     container->dirichfacs[5] * (-1);
   }/* loop j over columns */
}/* loop i over rows */
/*----------------------------------------------------------------------*/

} /* end of if (idamp) */



/*------------------- multiply dforces with -1 to obtain support forces */
/* in this routine the displacements, the velocities and the accelerations*/
/* of the time n+1 are used to compute the internal forces. Hence, the */
/* forces are divided by (1-alpha_f). The factor (1-alpha_f) is included */
/* in the dirichfacs. */
for (i=0; i<nd; i++)
{
  dforces[i] *= ((1.0) / (container->dirichfacs[6]));
}

/*-------- now assemble the vector dforces to the global vector coup_forces */
counter = 0;
for (i=0; i<actele->numnp; i++)
{
  actsnode = actele->node[i];
  /*---------- check if the node under consideration is a coupling node */
  if (actsnode->gnode->ssicouple == NULL) 
  {
    counter += 2;
    continue; 
  }
  if (actsnode->gnode->ssicouple->ssi_couptyp == ssi_slave)
  {
    if(container->coupl_typ == 0) /* conforming discretization */
    {  
      /*------------------------------ set the corresponing master node */
      actmnode = actsnode->gnode->mfcpnode[0];
      for (j=0; j<actsnode->numdf; j++)
      {
        actmnode->sol_mf.a.da[4][j] += dforces[counter];
        counter++;
      }
    }
    else  /* non-conforming discretization */
    {
      for (j=0; j<actsnode->numdf; j++)
      {
        actsnode->sol_mf.a.da[4][j] += dforces[counter];
        counter++;
      }
    }
  }
  else
  {
    counter += 2;
  }
}
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* calc_ssi_coupforce_mod */




/*!----------------------------------------------------------------------*
\brief enlarge node.sol.a.da to [11][numdf-1] and init to zero

<pre>
                                                   firl / genk 10/03    
                                                                       
  The routine is placed in ssi_service.c
                                                                       
\param  *actfield  FIELD      (i) pointer to the field under consideration                   

</pre>                                                                       
                                                                       
------------------------------------------------------------------------*/        
/*----------------------------------------------------------------------*
 |  Resize node.sol.a.da to [11][numdf-1] and set      firl / genk 10/03|
 |  values to zero                                                      |
 |  necessary due to ssi-coupling                                       |
 |  The basis for this routine is solserv_result_total written by m.gee |
 |  FIELD *actfield (i) the active field                                |
 |                                                                      |
 *----------------------------------------------------------------------*/
void init_sol_a_tozero(FIELD *actfield)
{
INT      i;

INT      numeq_total;
NODE    *actnode;

#ifdef DEBUG 
dstrc_enter("init_sol_a_tozero");
#endif
/*----------------------------------------------------------------------*/
numeq_total = actfield->dis->numeq;
/*------------------------ loop nodes and set sol.a.da[0..1][i] to zero */
for (i=0; i<actfield->dis[0].numnp; i++)
{
   actnode = &(actfield->dis[0].node[i]);
   /*--------------------------------------------------- enlarge sol ---*/
   
   amredef(&(actnode->sol),11,actnode->sol.sdim,"DA");
   amzero(&(actnode->sol));
#if 0
   for (j=0; j<actnode->numdf; j++)
   {
      dof = actnode->dof[j];
      if (dof>=numeq_total) continue;
      actnode->sol.a.da[0][j] = 0.0;
      actnode->sol.a.da[1][j] = 0.0;
      actnode->sol.a.da[2][j] = 0.0;
   }   
#endif
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of init_sol_a_tozero */

/*!----------------------------------------------------------------------*
\brief enlarge node.sol_mf.a.da to [7][numdf-1] and init to zero

<pre>
                                                   firl / genk 10/03    
                                                                       
  The routine is placed in ssi_service.c
                                                                       
\param  *actfield  FIELD      (i) pointer to the field under consideration                   

</pre>                                                                       
                                                                       
------------------------------------------------------------------------*/    
/*----------------------------------------------------------------------*
 |  Resize node.sol_mf.a.da to [7][numdf-1] and set    firl / genk 10/03|
 |  values to zero                                                      |
 |  necessary due to ssi-coupling                                       |
 |  The basis for this routine is solserv_result_total written by m.gee |
 |  FIELD *actfield (i) the active field                                |
 |                                                                      |
 *----------------------------------------------------------------------*/
void init_sol_mf_a_tozero(FIELD *actfield)
{
INT      i;
INT      numeq_total;
NODE    *actnode;

#ifdef DEBUG 
dstrc_enter("init_sol_mf_a_tozero");
#endif
/*----------------------------------------------------------------------*/
numeq_total = actfield->dis->numeq;
/*---------------------- loop nodes and set sol_mf.a.da[0..1][j] to zero */
for (i=0; i<actfield->dis[0].numnp; i++)
{
   actnode = &(actfield->dis[0].node[i]);
   /*---------------------------------------- enlarge sol, if necessary */
   amredef(&(actnode->sol_mf),7,actnode->sol_mf.sdim,"DA");
   amzero(&(actnode->sol_mf));
#if 0
   for (j=0; j<actnode->numdf; j++)
   {
      dof = actnode->dof[j];
      if (dof>=numeq_total) continue;
      actnode->sol_mf.a.da[0][j] = 0.0;
      actnode->sol_mf.a.da[1][j] = 0.0;
      actnode->sol_mf.a.da[2][j] = 0.0;
      actnode->sol_mf.a.da[3][j] = 0.0;
      actnode->sol_mf.a.da[4][j] = 0.0;
      actnode->sol_mf.a.da[5][j] = 0.0;
      actnode->sol_mf.a.da[6][j] = 0.0;
   }   
#endif
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of init_sol_a_tozero */




/*!----------------------------------------------------------------------*
\brief point neumann conditions

<pre>
                                                   firl / genk 10/03    
  The routine is placed in ssi_service.c
                                                                       
\param  *rhs      DOUBLE      (o) rhs for master field
         dimrhs   INT         (i) length of rhs
        *actpart  PARTITION   (i) pointer to actual partition

</pre>                                                                       
                                                                       
------------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |  point neumann conditions                      firl / genk 10/03     |
 | this routine is based on rhs_point_neum written by m.gee in 3/02     |
 | Attention! This assembly of nodal forces works only correct with     |
 | partitioning in the "Cut_Elements" style !!!!                        |
 |                                                                      |
 | Here we set up the rhs for the master field                          |
 | The forces assembled in this routine are resulting from ssi coupling |
 *----------------------------------------------------------------------*/
void ssiserv_rhs_point_neum(DOUBLE *rhs, INT dimrhs, PARTITION *actpart)     
{
INT             i,j;
NODE           *actnode;
INT             actdof;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("ssiserv_rhs_point_neum");
#endif
/*------------------------ loop all nodes of actpart (all master nodes) */
for (i=0; i<actpart->pdis[0].numnp; i++)
{
  /*--------------------------------------------------- set active node */
  actnode = actpart->pdis[0].node[i];
  /*---------------- check if the node is a coupling node (master node) */
  if (actnode->gnode->ssicouple == NULL) continue;
  if (actnode->gnode->ssicouple->ssi_couptyp == ssi_master)
  {
    /*----------------------------------------------- loop dofs of node */
    for (j=0; j<actnode->numdf; j++)
    {
      /*------------------------------------------------ set active dof */
      actdof = actnode->dof[j];
               /*---------- check whether this dof is inside free dofs (<dimrhs) */
      if (actdof < dimrhs)
        /*------------------------------------------ assemble the value */
        rhs[actdof] += actnode->sol_mf.a.da[4][j];
    }
  }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ssiserv_rhs_point_neum */




/*!---------------------------------------------------------------------                                         
\brief convergence check for SSI-iteration

<pre>                                                  firl / genk 10/03

in this routine the convergence ratios for the SSI iteration over the
fields is calculated.
There are two possibilities to check the convergence:
(see dissertation of D.P. MOK chapter 6.2)

- scaled_2-norm_of_residual (inrmfsi=1):
   || g(i) || / sqrt(neq) <= TOL

-  2-norm_of_residual_of_1st_iter (inrmfsi=2):  
   || g(i) || / || g(0) || <= TOL

where g(i) = d~(i+1) - d(i);
      neq ... number of structural inteface dofs
		     
</pre>
\param *structfield   FIELD	     (i)   structural field
\param *ssidyn 	      SSI_DYNAMIC    (i)   
\return INT                                                                             

------------------------------------------------------------------------*/
INT ssi_convcheck(FIELD *structfield, SSI_DYNAMIC *ssidyn, INT itnum)
{
INT     i,j;           /* some counters                                 */
INT     converged=0;   /* flag for convergence                          */
INT     numnp_total;   /* total number of nodes                         */
INT     numdf_total;   /* total number of dofs                          */
INT     numdf,dof;     /* actual number of dofs, actual dof             */
INT    *sid;           /* structural interface dofs                     */
DOUBLE  fac; 
DOUBLE  gnorm=ZERO;    
DOUBLE  g;
DOUBLE  grat=10.0;
NODE   *actsnode;      /* actual struct node                            */
static DOUBLE g0norm;  /* norm of first iteration                       */

#ifdef DEBUG 
dstrc_enter("ssi_convcheck");
#endif

if (itnum==0)
{
   grat=ONE;
   goto check;
}
    
/*----------------------------------------------------- set some values */
numnp_total = structfield->dis[0].numnp;
sid         = ssidyn->sid.a.iv;
numdf_total = ssidyn->sid.fdim;

/*---------------------- loop nodes and calculate norm at the interface */
for (i=0;i<numnp_total;i++)
{
   actsnode  = &(structfield->dis[0].node[i]);
   /*----------------------------------------- check for coupling nodes */
   numdf = actsnode->numdf; 
   /*-------------------------------------------------------- loop dofs */
   for (j=0;j<numdf;j++)
   {   
      dof = actsnode->dof[j];
      dsassert(dof<numdf_total,"dofnumber not valid!\n");
      if (sid[dof]==0) continue;
      g = actsnode->sol_mf.a.da[0][j] - actsnode->sol_mf.a.da[1][j];
      gnorm += g*g;
   } /* end of loop over dofs */  
} /* end of loop over nodes */   

/*------------------------------------- determine the convergence ratio */
gnorm = sqrt(gnorm);
switch (ssidyn->inrmfsi)
{
case 1: /* scaled_2-norm_of_residual */
   fac  = sqrt(ssidyn->numsid);
   grat = gnorm/fac;
break;
case 2: /* 2-norm_of_residual_of_1st_iter */
   if (itnum==1) 
   {
      g0norm = gnorm;
      if (g0norm<EPS5) g0norm=ONE;
   }
   grat = gnorm/g0norm;   
break;
default:
   dserror("parameter out of range: inrmfsi\n");
} /* end switch (fsidyn->inrmfsi) */

/*------------------------------ check for convergence over the fields */
check:
if (grat<ssidyn->convtol)
   converged=2;
if (itnum==ssidyn->itemax)
   converged++;

/*----------------------------------------------- output to the screen */
if (par.myrank==0)
{
printf("CONVERGENCE CHECK FOR ITERATION OVER FIELDS (ITNUM = %4d/%4d):\n",
        itnum,ssidyn->itemax);
switch (ssidyn->inrmfsi)
{
case 1:
   if (converged==0)
   {
      printf("|| g(i) || / sqrt(neq) = %10.3E >= TOL = %10.3E \n",
              grat,ssidyn->convtol);
      printf("NO CONVERGENCE OF ITERATION OVER FIELDS!  STEP %3i \n\n", ssidyn->step);
   }
   if (converged==1)
   {
      printf("|| g(i) || / sqrt(neq) = %10.3E >= TOL = %10.3E \n",
              grat,ssidyn->convtol);
      printf("NO CONVERGENCE OF ITERATION OVER FIELDS AFTER ITEMAX STEPS!\n");
      printf("                ***** CONTINUING ****\n\n");
   }    	
   if (converged>=2)
   {
      printf("|| g(i) || / sqrt(neq) = %10.3E < TOL = %10.3E \n",
              grat,ssidyn->convtol);
      printf("CONVERGENCE OF ITERATION OVER FIELDS!  STEP %3i \n\n", ssidyn->step);
   }
break;
case 2:   
   if (converged==0)
   {
      printf("|| g(i) || / || g(0) || = %10.3E >= TOL = %10.3E \n",
              grat,ssidyn->convtol);
      printf("NO CONVERGENCE OF ITERATION OVER FIELDS!\n\n");
   }
   if (converged==1)
   {
      printf("|| g(i) || / || g(0) || = %10.3E >= TOL = %10.3E \n",
              grat,ssidyn->convtol);
      printf("NO CONVERGENCE OF ITERATION OVER FIELDS AFTER ITEMAX STEPS!\n");
      printf("                ***** CONTINUING ****\n\n");
   }    	
   if (converged>=2)
   {
      printf("|| g(i) || / || g(0) || = %10.3E < TOL = %10.3E \n",
              grat,ssidyn->convtol);
      printf("CONVERGENCE OF ITERATION OVER FIELDS!\n\n");
   }
break;
default:
   dserror("parameter out of range: inrmfsi\n");   
}/* end switch (ssidyn->inrmfsi) */
}/* end if (par.myrank==0)*/

#ifdef DEBUG 
dstrc_exit();
#endif
return (converged);
} /* end of  ssi_convcheck*/




/*!---------------------------------------------------------------------                                         
\brief service routine to update the coupling forces in sol_mf[4]

<pre>                                                  firl / genk 10/03
The old coupling forces are written on sol_mf[5] these where scaled and 
added to the actual coupling forces written in sol_mf[4]
</pre>

\param *actfield              FIELD          (i)   the actual field
\param *container             CONTAINER	     (i)   a pointer on the container
------------------------------------------------------------------------
 facs[0] = -(1.0-alpham)*(1.0/beta)/(DSQR(dt))  		
 facs[1] =  (1.0-alpham)*(1.0/beta)/dt  			
 facs[2] =  (1.0-alpham)/(2*beta) - 1				
 facs[3] = -(1.0-alphaf)*(gamma/beta)/dt			
 facs[4] =  (1.0-alphaf)*gamma/beta - 1 			
 facs[5] =  (gamma/(2*beta)-1)*(1.0-alphaf)*dt  		
 facs[6] = -(1.0-alphaf)					
 facs[7] =  raleigh damping factor for mass			
 facs[8] =  raleigh damping factor for stiffness		
 facs[9] =  dt  	
------------------------------------------------------------------------*/
 					
void ssiserv_update_coupforc(FIELD *actfield, CONTAINER *container)
{
INT i,j;
NODE *actsnode, *actmnode;
for (i=0; i<actfield->dis->numnp; i++)
{
  actmnode = &(actfield->dis[0].node[i]);
  if (actmnode->gnode->ssicouple->ssi_couptyp == ssi_master)
  {
    actsnode = actmnode->gnode->mfcpnode[1];
    for (j=0; j<actmnode->numdf; j++)
    {
      actmnode->sol_mf.a.da[5][j] *= (((container->dirichfacs[6])+1.0) / 
                                     (container->dirichfacs[6] * (-1.0)));
      actmnode->sol_mf.a.da[4][j] *= (1.0 / 
                                     (container->dirichfacs[6] * (-1.0)));
      actmnode->sol_mf.a.da[4][j] += actmnode->sol_mf.a.da[5][j];
    }
    
  }
}
}

void ssiserv_erase_coupforc(FIELD *actfield, CONTAINER *container)
{
INT i,j;
NODE *actsnode, *actmnode;
for (i=0; i<actfield->dis->numnp; i++)
{
  actsnode = &(actfield->dis[0].node[i]);
  if (actsnode->gnode->ssicouple->ssi_couptyp == ssi_slave)
  {
    actmnode = actsnode->gnode->mfcpnode[0];
    for (j=0; j<actsnode->numdf; j++)
    {
      actmnode->sol_mf.a.da[4][j] = 0.0; 
    }
    
  }
}
}


/*!---------------------------------------------------------------------                                         
\brief detect number of nodes on slave and master interface

<pre>                                                  firl / genk 3/04

In this subroutine the interface data are stored in the structure 
interfaces.  

</pre>
\param *masterfield     FIELD	      (i)   master field
\param *slavefield      FIELD         (i)   slave field
\param *int_faces       INTERFACES    (i)   the interfaces of the model

\return 0                                                                             

------------------------------------------------------------------------*/
void ssi_detect_nmb_nods_int(FIELD *masterfield, FIELD *slavefield,
                             INTERFACES *int_faces)
{
INT i,j,a,l,b;                      /* counters */
NODE *actnode;                      /* actual node under consideration */
NODE *actnode1;                     /* actual node under consideration */
NODE *node_work;
ELEMENT *actele;                    /* actual element under consideration */
INT k;                              /* number of dline */
GLINE *actgline;                    /* actual gline under consideration */
GNODE *actgnode;                    /* actual gnode under consideration */
/*ARRAY slv_nods_onint;*/               /* ARRAY for the slave nodes on the */ 
                                    /* interface */
/*ARRAY mst_nods_onint;*/               /* ARRAY for the master nodes on the*/ 
                                    /* interface */
/*INT work;*/                           /* working integer */
INTERFACE *actinterface, *actinterface2;/* pointer to the actual interface */
DOUBLE work;
DOUBLE dist1, dist1x, dist1y, dist2, dist2x, dist2y; /* nodal distances */
#ifdef DEBUG 
dstrc_enter("ssi_detect_nmb_nods_int");
#endif


/* assign the Id to the interfaces */
for(i=0; i<int_faces->numint;i++)
{
  int_faces->interface[i].Id = int_faces->int_ids[i];
}
/* loop the interfaces */
for(i=0; i<int_faces->numint;i++)
{
  actinterface = &int_faces->interface[i];
  /*--------------------------------------------------------------------*/
  /*                            MASTER DOMAIN                           */
  /*--------------------------------------------------------------------*/
  /* -I------ detect the number of glines at the corresponding interface*/
  a = 0;
  for(j=0; j<masterfield->dis->ngline; j++)
  {
     actgline = &masterfield->dis->gline[j];
     if(actgline->ssicouple == NULL) continue;
     
     if((actgline->ssicouple->ssi_couptyp == ssi_master) &&
        (actgline->ssicouple->ssi_coupleId == actinterface->Id))
        a++;
  }
  /* -------------assign the number of master nodes at the interface to */
  /* -----------------------------------------------actinterface->numnpm*/
  actinterface->numnpm = a+1;
  /* ----------assign the number of master elements at the interface to */
  /* ----------------------------------------------actinterface->numelem*/
  actinterface->numelem = a;

  /* -II--allocate memory for the pointer vector to the master nodes and*/
  /* -------------------------------the master elements of actinterface */
  actinterface->nodem = 
     (NODE**)CCACALLOC(actinterface->numnpm,sizeof(NODE));
  actinterface->elementm = 
     (ELEMENT**)CCACALLOC(actinterface->numelem,sizeof(ELEMENT));

  /* -III define the pointers in actinterface->elementm to the elements */
  /* --------------------------------------------------at the interface */
  a = 0;
  for(j=0; j<masterfield->dis->numele; j++)
  {
     actele = &masterfield->dis->element[j];
     /* loop the glines of this element */
     for(k=0; k<actele->numnp; k++)
     {
       actgline = actele->g.gsurf->gline[k];
       if(actgline->ssicouple == NULL) continue;
       
       if((actgline->ssicouple->ssi_couptyp == ssi_master) &&
          (actgline->ssicouple->ssi_coupleId == actinterface->Id))
       {
          actinterface->elementm[a] = actele;
          a++;
          goto end_of_this_ele;
       }
     }     
    end_of_this_ele: ;
  }
  /* -IV- define the pointers in actinterface->nodem to the nodes at the*/
  /* ---------------------------------------------------------interface */
  a = 0;
  for(j=0; j<masterfield->dis->ngline; j++)
  {
    actgline = &masterfield->dis->gline[j];
    if(actgline->ssicouple == NULL) continue;
    
    if((actgline->ssicouple->ssi_couptyp == ssi_master) &&
       (actgline->ssicouple->ssi_coupleId == actinterface->Id))
    {
       /* loop gnodes of this gline */
       for(k=0; k<2; k++)
       {
         actgnode = actgline->gnode[k];
         actnode = actgnode->node;
         /* loop vector of actinterface->nodem[] to detect if actnode is */
         /* already there; this happens, because two glines share a node */
         for(l=0; l<actinterface->numnpm; l++)
         {
           if(actinterface->nodem[l] == NULL) continue;
           if(actinterface->nodem[l]->Id == actnode->Id)
           { 
             goto end_of_this_node;
           }
         }
         actinterface->nodem[a] = actnode;
         a++;
         end_of_this_node: ;
       }
    }       
  }


  /*--------------------------------------------------------------------*/
  /*                             SLAVE DOMAIN                           */
  /*--------------------------------------------------------------------*/
  /* -I------ detect the number of glines at the corresponding interface*/
  a = 0;
  for(j=0; j<slavefield->dis->ngline; j++)
  {
     actgline = &slavefield->dis->gline[j];
     if(actgline->ssicouple == NULL) continue;
     
     if((actgline->ssicouple->ssi_couptyp == ssi_slave) &&
        (actgline->ssicouple->ssi_coupleId == actinterface->Id))
        a++;
  }
  /* --------------assign the number of slave nodes at the interface to */
  /* -----------------------------------------------actinterface->numnps*/
  actinterface->numnps = a+1;
  /* -----------assign the number of slave elements at the interface to */
  /* ----------------------------------------------actinterface->numeles*/
  actinterface->numeles = a;

  /* -II---allocate memory for the pointer vector to the slave nodes and*/
  /* ------------------------------- the slave elements of actinterface */
  actinterface->nodes = 
     (NODE**)CCACALLOC(actinterface->numnps,sizeof(NODE));
  actinterface->elements = 
     (ELEMENT**)CCACALLOC(actinterface->numeles,sizeof(ELEMENT));

  /* -III define the pointers in actinterface->elements to the elements */
  /* --------------------------------------------------at the interface */
  a = 0;
  for(j=0; j<slavefield->dis->numele; j++)
  {
     actele = &slavefield->dis->element[j];
     /* loop the glines of this element */
     for(k=0; k<actele->numnp; k++)
     {
       actgline = actele->g.gsurf->gline[k];
       if(actgline->ssicouple == NULL) continue;
       if((actgline->ssicouple->ssi_couptyp == ssi_slave) &&
          (actgline->ssicouple->ssi_coupleId == actinterface->Id))
       {
          actinterface->elements[a] = actele;
          a++;
          goto end_of_this_ele_s;
       }
     }     
    end_of_this_ele_s: ;
  }
  /* -IV- define the pointers in actinterface->nodes to the nodes at the*/
  /* ---------------------------------------------------------interface */
  a = 0;
  for(j=0; j<slavefield->dis->ngline; j++)
  {
    actgline = &slavefield->dis->gline[j];
    if(actgline->ssicouple == NULL) continue;
       
    if((actgline->ssicouple->ssi_couptyp == ssi_slave) &&
       (actgline->ssicouple->ssi_coupleId == actinterface->Id))
    {
       /* loop gnodes of this gline */
       for(k=0; k<2; k++)
       {
         actgnode = actgline->gnode[k];
         actnode = actgnode->node;
         /* loop vector of actinterface->nodem[] to detect if actnode is */
         /* already there; this happens, because two glines share a node */
         for(l=0; l<actinterface->numnps; l++)
         {
           if(actinterface->nodes[l] == NULL) continue;
           if(actinterface->nodes[l]->Id == actnode->Id)
           { 
             goto end_of_this_node_s;
           }
         }
         actinterface->nodes[a] = actnode;
         a++;
         end_of_this_node_s: ;
       }
    }       
  }
  /* allocate memory for the gnodes which bound actinterface */
  actinterface->gnode_bound1s = (GNODE*)CCACALLOC(1,sizeof(GNODE));
  actinterface->gnode_bound2s = (GNODE*)CCACALLOC(1,sizeof(GNODE));
  actinterface->gnode_bound1m = (GNODE*)CCACALLOC(1,sizeof(GNODE));
  actinterface->gnode_bound2m = (GNODE*)CCACALLOC(1,sizeof(GNODE));

  /* look for the slave gnodes which bound the interface */
  /* loop the slave elements of actinterface */
  b=0;
  for(j=0; j<actinterface->numeles; j++)
  {
    actele = actinterface->elements[j];
    /* loop gnodes of actele */
    for(k=0; k<actele->numnp; k++)
    {
      actgnode = actele->node[k]->gnode;
      /* loop glines of actgnode */
      a=0;
      for(l=0; l<actgnode->ngline; l++)
      {
        actgline = actgnode->gline[l];
        if(actgline->ssicouple == NULL) continue;
        if(actgline->ssicouple->ssi_coupleId == actinterface->Id)
          a++;
      }
      if(a==1 && b==1)
        actinterface->gnode_bound2s = actgnode;
      if(a==1 && b==0)
      {
        actinterface->gnode_bound1s = actgnode;
        b++;
      }
    }
  }
  /* look for the master gnodes which bound the interface */
  /* loop the master elements of actinterface */
  b=0;
  for(j=0; j<actinterface->numelem; j++)
  {
    actele = actinterface->elementm[j];
    /* loop gnodes of actele */
    for(k=0; k<actele->numnp; k++)
    {
      actgnode = actele->node[k]->gnode;
      /* loop glines of actgnode */
      a=0;
      for(l=0; l<actgnode->ngline; l++)
      {
        actgline = actgnode->gline[l];
        if(actgline->ssicouple == NULL) continue;
        
        if(actgline->ssicouple->ssi_coupleId == actinterface->Id)
          a++;
      }
      if(a==1 && b==1)
        actinterface->gnode_bound2m = actgnode;
      if(a==1 && b==0)
      {
        actinterface->gnode_bound1m = actgnode;
        b++;
      }
    }
  }
  /* define the positive interface direction */
  actnode = actinterface->gnode_bound1s->node;
  actnode1 = actinterface->gnode_bound2s->node;
  /* allocate an vector to actinterface->int_vec (2-d case) */
  /*amdef("int_vec",actinterface->int_vec, 2, 1, "DV");*/
  actinterface->int_vec = (DOUBLE*)CCACALLOC(2,sizeof(DOUBLE));

  actinterface->int_vec[0] = actnode1->x[0] - actnode->x[0];  
  actinterface->int_vec[1] = actnode1->x[1] - actnode->x[1];
  work = sqrt(actinterface->int_vec[0] * actinterface->int_vec[0] +
              actinterface->int_vec[1] * actinterface->int_vec[1]);  
  actinterface->int_vec[0] *=(1.0/work);
  actinterface->int_vec[1] *=(1.0/work);

/* allocate memory for the continuity equation */
  actinterface->continuity_eq = (DENSE*)CCACALLOC(1,sizeof(DENSE));
  amdef("LHS",&(actinterface->continuity_eq->A),actinterface->numnps,
         actinterface->numnps,"DA");
  amdef("ipiv",&(actinterface->continuity_eq->ipiv),actinterface->numnps,1,"IV");
  amdef("work",&(actinterface->continuity_eq->work),2*actinterface->numnps,1,"DV");

/* -------allocate memory for the continuity equation which is used in  */
/* ----------------------------------------------------ssi_calc_disp4slv*/
  actinterface->conti_eq_save = (DENSE*)CCACALLOC(1,sizeof(DENSE));
  amdef("LHS",&(actinterface->conti_eq_save->A),actinterface->numnps,
         actinterface->numnps,"DA");
  amdef("ipiv",&(actinterface->conti_eq_save->ipiv),actinterface->numnps,1,"IV");

  /* sort the vectors of nodes at the interface */
  /* loop the slave nodes at the interface */
  for(j=0; j<actinterface->numnps; j++)
  {
    actnode = actinterface->nodes[j];
    dist1x = actnode->x[0] - actinterface->gnode_bound1s->node->x[0];
    dist1y = actnode->x[1] - actinterface->gnode_bound1s->node->x[1];
    dist1 = sqrt(dist1x * dist1x + dist1y * dist1y);
    for(k=0; (k+j)<actinterface->numnps; k++)
    {
      actnode1 = actinterface->nodes[k+j];
      dist2x = actnode1->x[0] - actinterface->gnode_bound1s->node->x[0];
      dist2y = actnode1->x[1] - actinterface->gnode_bound1s->node->x[1];
      dist2 = sqrt(dist2x * dist2x + dist2y * dist2y);
      if(dist2<dist1)
      {
        /* switch the pointers in the vector */
        node_work = actinterface->nodes[j];
        actinterface->nodes[j] = actinterface->nodes[k+j];
        actinterface->nodes[k+j] = node_work;
        /* switch the distances */
        work = dist1;
        dist1 = dist2;
        dist2 = work;
      }
    }
  }
  /* loop the master nodes at the interface */
  for(j=0; j<actinterface->numnpm; j++)
  {
    actnode = actinterface->nodem[j];
    dist1x = actnode->x[0] - actinterface->gnode_bound1s->node->x[0];
    dist1y = actnode->x[1] - actinterface->gnode_bound1s->node->x[1];
    dist1 = sqrt(dist1x * dist1x + dist1y * dist1y);
    for(k=0; (k+j)<actinterface->numnpm; k++)
    {
      actnode1 = actinterface->nodem[k+j];
      dist2x = actnode1->x[0] - actinterface->gnode_bound1s->node->x[0];
      dist2y = actnode1->x[1] - actinterface->gnode_bound1s->node->x[1];
      dist2 = sqrt(dist2x * dist2x + dist2y * dist2y);
      if(dist2<dist1)
      {
        /* switch the pointers in the vector */
        node_work = actinterface->nodem[j];
        actinterface->nodem[j] = actinterface->nodem[k+j];
        actinterface->nodem[k+j] = node_work;
        /* switch the distances */
        work = dist1;
        dist1 = dist2;
        dist2 = work;
      }
    }
  }

  /* store a 1 in vector actinterface->conti_eq_save->ipiv, if actnode1_1 is on the */
  /* boundary of the interface (if it is on a dnode) */
  actinterface->conti_eq_save->ipiv.a.iv[0] = 1;
  actinterface->conti_eq_save->ipiv.a.iv[actinterface->numnps-1] = 1;
  

} /* end of loop over the interfaces (i) */
/* loop over the interfaces */
for(i=0; i<int_faces->numint; i++)
{
  actinterface = &int_faces->interface[i];
  /* detect the gnodes which connect the several interfaces */

  /* slave domain */
  for(j=0; j<actinterface->numnps; j++)
  {
    actnode = actinterface->nodes[j];
    if(actnode->gnode->node->Id == actinterface->gnode_bound1s->node->Id || 
       actnode->gnode->node->Id == actinterface->gnode_bound2s->node->Id)
    {
      /* loop the other interfaces */
      for(k=0; k<int_faces->numint; k++)
      {
        actinterface2 = &int_faces->interface[k];
        if(actinterface->Id == actinterface2->Id)
          goto end_of_this_interface;
        else
        {
          if(actnode->gnode->node->Id == actinterface2->gnode_bound1s->node->Id ||
             actnode->gnode->node->Id == actinterface2->gnode_bound2s->node->Id)
          {
            if(actinterface->gnode_bound1sc == NULL)
              actinterface->gnode_bound1sc = actnode->gnode;
            else 
              actinterface->gnode_bound2sc = actnode->gnode;
            break;
          }
        }
        end_of_this_interface: ;
      }
      
    }    
  }

  /* master domain */
  for(j=0; j<actinterface->numnpm; j++)
  {
    actnode = actinterface->nodem[j];
    if(actnode->gnode->node->Id == actinterface->gnode_bound1m->node->Id || 
       actnode->gnode->node->Id == actinterface->gnode_bound2m->node->Id)
    {
      /* loop the other interfaces */
      for(k=0; k<int_faces->numint; k++)
      {
        actinterface2 = &int_faces->interface[k];
        if(actinterface->Id == actinterface2->Id)
          goto end_of_this_interface1;
        else
        {
          if(actnode->gnode->node->Id == actinterface2->gnode_bound1m->node->Id ||
             actnode->gnode->node->Id == actinterface2->gnode_bound2m->node->Id)
          {
            if(actinterface->gnode_bound1mc == NULL)
              actinterface->gnode_bound1mc = actnode->gnode;
            else 
              actinterface->gnode_bound2mc = actnode->gnode;
            break;
          }
        }
        end_of_this_interface1: ;
      }
      
    }    
  }

} /* end og loop over the interfaces (i) */

#ifdef DEBUG 
dstrc_exit();
#endif
} /* end of ssi_detect_nmb_nods_int */
/*! @} (documentation module close)*/


/*!---------------------------------------------------------------------                                         
\brief compute coefficients for mortar method

<pre>                                                  firl / genk 1/04

Here we assemble a coefficient matrix and a right hand side. One obtains 
the coefficients (called zeta_j^r) as solution of this system of 
equations. We refer to the master thesis of Matthias Firl for more 
informations.
Later the factorized coefficient matrix is needed again. Thus, the pointer
to the dense structure is a parameter of the function call.

</pre>
\param *masterfield   FIELD	     (i)   master field
\param *slavefield    FIELD          (i)   slave field
\param nmb_nods_mstint INT           (i)   nmb of nods on master interface
\param nmb_nods_slvint INT           (i)   nmb of nods on slave interface
\param *slv_nods_onint ARRAY         (i)   slave nodes in the interface
\param *mst_nods_onint ARRAY         (i)   master nodes in the interface
\param step            INT           (i)   the actual iteration step
\param *lhs_dens       DENSE         (i)   pointer to a dense structure
\return 0                                                                             

------------------------------------------------------------------------*/
void ssi_mortar_coeff(FIELD *masterfield, FIELD *slavefield, 
                      INT step, SSI_DYNAMIC *ssidyn, INTERFACES *int_faces)
{
DOUBLE gauss2[2][2];                /* gauss points and weights (2)*/
INT intersection;                   /* flag for intersection of nodal*/
                                    /* basis fct's 0=no, 1=yes*/
INT info;                           /* flag for lapack solver */
INT ione;                           /* number of right hand sides*/
INT a,b,i,j,k,l,m,n,p,q,r,act;    /* counters */
INT numnods, switch_act;            /* number of elements of the actual side*/
                                    /* of actinterface */            
char uplo[1];                       /* character for lapack solver */
NODE *actnode1_1, *actnode1_2, *actnode2_1, *actnode2_2;
                                    /* pointer to actual nodes */
NODE **actnodes;
DOUBLE d1_1[2], d1_2[2], d2_1[2], d2_2[2];
                                    /* position vectors of nodes on the  */
                                    /* interface */
DOUBLE nrmdr1_2[2];                 /* normed vector r1_2 */
DOUBLE nr1_2, nr2_1, nr2_2;         /* norm of vectors r1_2, r2_1, r2_2  */
DOUBLE r1_2[2], r2_1[2], r2_2[2];   /* position vectors relative to node */
                                    /* d1_1 */
DOUBLE lambdr1_2, lambdr2_1, lambdr2_2; /* e.g lambdr1_2 =               */
                                        /*         r1_2[0] / nrmdr1_2[0] */
DOUBLE b1, b2;                       /* integration bounds */
DOUBLE m1, n1, m2, n2;               /* parameters for linear functions */
INT flag;

ARRAY rhs1;                          /* right hand side of system of eq.  */
DOUBLE *rhs;                         /* to compute the mortar coeff.      */

                                     /* first loop (i)                    */
INTERFACE *actinterface;             /* the actual interface */


#ifdef DEBUG 
dstrc_enter("mortar_coefficients");
#endif

/* -------------------------------------------------values to variables */
a=0; i=0; j=0; k=0; l=0; m=0; n=0; p=0; switch_act = 0;

gauss2[0][0] = 0.577350269189626;
gauss2[1][0] = -0.577350269189626;
gauss2[0][1] = 1.0;
gauss2[1][1] = 1.0;

/*----------------------------------------------------------------------*/
/* loop over the interfaces */
/*----------------------------------------------------------------------*/
for(r=0; r<int_faces->numint; r++)
{
  actinterface = &int_faces->interface[r];
  /* reset the parameters in actinterface->continuity_eq */
  uplo[0] = 'L';
  /* nmb of rhs's */
  ione = actinterface->numnpm;
  actinterface->continuity_eq->numeq_total = actinterface->numnps ;
  actinterface->continuity_eq->lwork=
    2*(actinterface->continuity_eq->numeq_total);

  amzero(&(actinterface->continuity_eq->A));

  amzero(&(actinterface->continuity_eq->ipiv));

  amzero(&(actinterface->continuity_eq->work));

  rhs = amdef("RHS",&rhs1,ione,actinterface->continuity_eq->numeq_total,"DA");
  amzero(&rhs1);

/* the DENSE structure conti_eq_save is used to store the original coefficient */
/* matrix; this matrix is needed in the routine ssi_calc_disp4slv again */
/* furthermore the integer vector ipiv is used to store the number 1 if the */
/* actual node is on a dnode; if the actual node is on a dline a zero is stored*/
/* this is also needed in ssi_calc_disp4slv to indicate the boundary of the */
/* interface */

  amzero(&(actinterface->conti_eq_save->A));

  /* ---In the first loop we compute the left hand side and in the second */
  /* the right hand side; hence, act == 0 -> lhs, act == 1 -> rhs */
  for(act=0;act<2;act++)
  {
    if(act==0) /* slave nodes */
    {
      actnodes = 
      (NODE**)CCACALLOC(actinterface->numnps,sizeof(NODE));
      for(i=0; i<actinterface->numnps; i++)
      {
        actnodes[i] = actinterface->nodes[i];
      }
      /*field1 = slavefield;*/
      numnods = actinterface->numnps;
    }
    else /* master nodes */
    {
      actnodes = 
      (NODE**)CCACALLOC(actinterface->numnpm,sizeof(NODE));
      for(i=0; i<actinterface->numnpm; i++)
      {
        actnodes[i] = actinterface->nodem[i];
      }
      /*field1 = slavefield;*/
      numnods = actinterface->numnpm;
    }
    /* ------------------------------------------------- loop actnodes */
    for(i=0;i<numnods-1;i++)
    {
      for(j=0;j<2;j++)
      {
        actnode1_1 = actnodes[i];
        actnode1_2 = actnodes[i+1];
        /* decide which vector of the array sol_mf should be used */
        /*if(field1 == slavefield)*/
        if(act == 0)
        {
          /* the deformed geometry is used to set up the functions*/
          d1_1[0] = actnode1_1->x[0]+actnode1_1->sol_mf.a.da[6][0];
          d1_1[1] = actnode1_1->x[1]+actnode1_1->sol_mf.a.da[6][1];
          d1_2[0] = actnode1_2->x[0]+actnode1_2->sol_mf.a.da[6][0];
          d1_2[1] = actnode1_2->x[1]+actnode1_2->sol_mf.a.da[6][1];
        }
        else /* actual elements = master elements */
        {
          /* the deformed geometry is used to set up the functions*/
          actnode1_1->sol_mf.a.da[6][0] = (ssidyn->relax * 
                                           actnode1_1->sol_mf.a.da[0][0])+
                                          (1-ssidyn->relax) * 
                                           actnode1_1->sol_mf.a.da[1][0];
          actnode1_1->sol_mf.a.da[6][1] = (ssidyn->relax * 
                                           actnode1_1->sol_mf.a.da[0][1])+
                                          (1-ssidyn->relax) * 
                                           actnode1_1->sol_mf.a.da[1][1];
          actnode1_2->sol_mf.a.da[6][0] = (ssidyn->relax * 
                                           actnode1_2->sol_mf.a.da[0][0])+
                                          (1-ssidyn->relax) * 
                                           actnode1_2->sol_mf.a.da[1][0];
          actnode1_2->sol_mf.a.da[6][1] = (ssidyn->relax * 
                                           actnode1_2->sol_mf.a.da[0][1])+
                                          (1-ssidyn->relax) * 
                                           actnode1_2->sol_mf.a.da[1][1];

          d1_1[0] = actnode1_1->x[0] + actnode1_1->sol_mf.a.da[6][0];
          d1_1[1] = actnode1_1->x[1] + actnode1_1->sol_mf.a.da[6][1];
          d1_2[0] = actnode1_2->x[0] + actnode1_2->sol_mf.a.da[6][0];
          d1_2[1] = actnode1_2->x[1] + actnode1_2->sol_mf.a.da[6][1];
        }
        /* -----------loop elements of slave domain, biorthogonal shp. fcts */
        b=0;
        for(l=0;l<actinterface->numnps-1;l++)
        {
          intersection = 0;
          actnode2_1 = actinterface->nodes[l];
          actnode2_2 = actinterface->nodes[l+1];
          for(m=0; m<2; m++)
          {
            /* the deformed geometry is used to set up the functions*/
            d2_1[0] = actnode2_1->x[0]+actnode2_1->sol_mf.a.da[6][0];
            d2_1[1] = actnode2_1->x[1]+actnode2_1->sol_mf.a.da[6][1];
            d2_2[0] = actnode2_2->x[0]+actnode2_2->sol_mf.a.da[6][0];
            d2_2[1] = actnode2_2->x[1]+actnode2_2->sol_mf.a.da[6][1];

            /* near the boundary of the interface we reduce the order of the */
            /* nodal basis function; hence the bounds in which this function */
            /* ----is defined also change, this is adjusted in the following */
            if(actnode2_1->gnode->node->Id == actinterface->gnode_bound1s->node->Id ||
               actnode2_1->gnode->node->Id == actinterface->gnode_bound2s->node->Id)
            {
              if(m == 0)
              {
                d2_2[0]=d2_1[0]+0.5*(d2_2[0]-d2_1[0]);
                d2_2[1]=d2_1[1]+0.5*(d2_2[1]-d2_1[1]);
              }
              if(m == 1)
              {
                d2_1[0]=d2_1[0]+0.5*(d2_2[0]-d2_1[0]);
                d2_1[1]=d2_1[1]+0.5*(d2_2[1]-d2_1[1]);
              }
            }
            if(actnode2_2->gnode->node->Id == actinterface->gnode_bound1s->node->Id ||
               actnode2_2->gnode->node->Id == actinterface->gnode_bound2s->node->Id) 
            {
              if(m == 0)
              {
                d2_2[0]=d2_1[0]+0.5*(d2_2[0]-d2_1[0]);
                d2_2[1]=d2_1[1]+0.5*(d2_2[1]-d2_1[1]);
              }
              if(m == 1)
              {
                d2_1[0]=d2_1[0]+0.5*(d2_2[0]-d2_1[0]);
                d2_1[1]=d2_1[1]+0.5*(d2_2[1]-d2_1[1]);
              }
            }
            /* -------------------------compute vectors relative to node 1_1 */
            r1_2[0] = d1_2[0] - d1_1[0];
            r1_2[1] = d1_2[1] - d1_1[1];
            r2_1[0] = d2_1[0] - d1_1[0];
            r2_1[1] = d2_1[1] - d1_1[1];
            r2_2[0] = d2_2[0] - d1_1[0];
            r2_2[1] = d2_2[1] - d1_1[1];
            /* ---------------------compute norm of vectors r1_2, r2_1, r2_2 */
            nr1_2 = sqrt((r1_2[0])*(r1_2[0])+(r1_2[1])*(r1_2[1]));
            nr2_1 = sqrt((r2_1[0])*(r2_1[0])+(r2_1[1])*(r2_1[1]));
            nr2_2 = sqrt((r2_2[0])*(r2_2[0])+(r2_2[1])*(r2_2[1]));
            /* ------------------------compute normed vector r1_2 = nrmdr1_2 */
            nrmdr1_2[0] = r1_2[0] / nr1_2;
            nrmdr1_2[1] = r1_2[1] / nr1_2;
            /* --------compute lambda values of the vectors r1_2, r2_1, r2_2 */
            if((r1_2[0] > sqrt(r1_2[1]*r1_2[1])) || (r1_2[0]<-sqrt(r1_2[1]*r1_2[1])))
            {
              lambdr1_2 = r1_2[0] / nrmdr1_2[0];
              lambdr2_1 = r2_1[0] / nrmdr1_2[0];
              lambdr2_2 = r2_2[0] / nrmdr1_2[0];
            }
            else
            {
              lambdr1_2 = r1_2[1] / nrmdr1_2[1];
              lambdr2_1 = r2_1[1] / nrmdr1_2[1];
              lambdr2_2 = r2_2[1] / nrmdr1_2[1];
            }
            /* Check for a possible nonzero intersection, goto end if such an*/
            /* ---------------------------------intersection is not possible */        
            if((lambdr2_1 >= nr1_2) && (lambdr2_2 >= nr1_2)) goto end_of_node2;
            if((lambdr2_1 < 0) && (lambdr2_2 < 0)) goto end_of_node2;
            
            /* ---------------------here we determine the integration bounds */
            ssi_detect_intersection(lambdr1_2, lambdr2_1, lambdr2_2, nr1_2, 
                                    nr2_1, nr2_2, &b1, &b2, &intersection);
    
            if(intersection == 1)
            {
              /* if we have a nonzero inters. we compute the nodal basis fcts*/
              /* ---calculate slope and value n of the nodal basis function, */ 
              /* here its the nodal basis function of the trace space of the */ 
              /* ------slave elements on the interface. Hence, these are the */ 
              /* -------------------------------------standard hat functions */
              if(j == 0)
              {
                m1=-1.0/lambdr1_2;
                n1=1.0;
              }
              else
              {
                m1=1.0/lambdr1_2;
                n1=0;
              }
              /* -now we compute the nodal basis functions for the Lagrange */
              /* ---multiplier space, these functions fulfill the so called */
              /* biorthogonality relation, we refer to the Master Thesis of */
              /* ------------------------Matthias Firl for more information */
              if(m == 0)
              {
                m2 = -3.0/(lambdr2_2 - lambdr2_1);
                n2 = 2.0 - m2 * lambdr2_1;
              }
              else
              {
                m2 = 3.0/(lambdr2_2 - lambdr2_1);
                n2 = -1.0 - m2 * lambdr2_1;
              }
    
              if(actnode2_1->gnode->node->Id == actinterface->gnode_bound1s->node->Id ||
                 actnode2_1->gnode->node->Id == actinterface->gnode_bound2s->node->Id ||
                 actnode2_2->gnode->node->Id == actinterface->gnode_bound1s->node->Id ||
                 actnode2_2->gnode->node->Id == actinterface->gnode_bound2s->node->Id) 
              {
                m2=0;
                n2=1;
              }
              if(act == 0)
              /* act == 0 -> lhs */
              {
                /* detect position of integral value in the coefficient matrix */
                for(q=0;q<actinterface->numnps;q++)
                {
                  if(j == 0 && actnode1_1->Id == actinterface->nodes[q]->Id)
                  /*if(out_node1->Id == actinterface->nodes[q]->Id)*/ 
                  { 
                    a=q;
                    break;
                  }
                  if(j == 1 && actnode1_2->Id == actinterface->nodes[q]->Id)
                  /*if(out_node1->Id == actinterface->nodes[q]->Id)*/ 
                  { 
                    a=q;
                    break;
                  }
                }
                for(q=0;q<actinterface->numnps;q++)
                {
                  if(m==0 && actnode2_1->Id == actinterface->nodes[q]->Id)
                  /*if(out_node2->Id == actinterface->nodes[q]->Id) */
                  { 
                    b=q;
                    break;
                  }
                  if(m==1 && actnode2_2->Id == actinterface->nodes[q]->Id)
                  /*if(out_node2->Id == actinterface->nodes[q]->Id) */
                  { 
                    b=q;
                    break;
                  }
                }
    
                for(p=0;p<2;p++)
                {
                  actinterface->continuity_eq->A.a.da[a][b] += 
                                     (m1 * ((b1+b2)/2 + (b1-b2)/2*gauss2[p][0])+n1) *
                                     (m2 * ((b1+b2)/2 + (b1-b2)/2*gauss2[p][0])+n2) *
                                     ((b2-b1)/2) * gauss2[p][1];
                  actinterface->conti_eq_save->A.a.da[a][b] += 
                                     (m1 * ((b1+b2)/2 + (b1-b2)/2*gauss2[p][0])+n1) *
                                     (m2 * ((b1+b2)/2 + (b1-b2)/2*gauss2[p][0])+n2) *
                                     ((b2-b1)/2) * gauss2[p][1];
    
                }
              } 
              else
              /* act == 1 -> rhs */
              {
                /* detect position of integral value the rhs matrix */
                for(q=0;q<actinterface->numnpm;q++)
                {
                  if(j==0 && actnode1_1->Id == actinterface->nodem[q]->Id) 
                  { 
                    b=q;
                    break;
                  }
                  if(j==1 && actnode1_2->Id == actinterface->nodem[q]->Id) 
                  { 
                    b=q;
                    break;
                  }
                }
                for(q=0;q<actinterface->numnps;q++)
                {
                  if(m==0 && actnode2_1->Id == actinterface->nodes[q]->Id) 
                  { 
                    a=q;
                    break;
                  }
                  if(m==1 && actnode2_2->Id == actinterface->nodes[q]->Id) 
                  { 
                    a=q;
                    break;
                  }
                }
    
                for(p=0;p<2;p++)
                {
                  /* --Here we save the transposed of the rhs solution matrix */
                  /* ------because the solver is in fortran and fortran has a */
                  /* -------------------------different memory storage system */
                  rhs1.a.da[b][a] += (m1 * ((b1+b2)/2 + (b1-b2)/2*gauss2[p][0])+n1) *
                                     (m2 * ((b1+b2)/2 + (b1-b2)/2*gauss2[p][0])+n2) *
                                     ((b2-b1)/2) * gauss2[p][1];
                }
              } /* end if if clause (act)*/
            }
            end_of_node2: ;
          } /* end of loop over the nodes of this slave element (m) */
        } /* end of loop over the slave elements (l)*/
      } /* end of loop over the nodes of this slave element (j) */
    } /* end of loop over the slave elements (i)*/
  } /* end of general loop lhs--rhs (act)*/
  
  /*-----------------------------------------------------------------------*/
  /*                     SOLUTION USING LAPACK LU-DECOMPOSITION            */
  /*-----------------------------------------------------------------------*/
  
  /* decomposition of lhs */
  
  dsytrf(
          uplo,              
          &(actinterface->continuity_eq->numeq_total),  
          actinterface->continuity_eq->A.a.da[0],  
          &(actinterface->continuity_eq->numeq_total),  
          actinterface->continuity_eq->ipiv.a.iv,  
          actinterface->continuity_eq->work.a.dv,  
          &(actinterface->continuity_eq->lwork),   
          &info              
        );
  
  if (info!=0) dserror("Lapack factorization failed");
  
  /* solution for different rhs's*/
  
  info=1;
  dsytrs(
           uplo,
           &(actinterface->continuity_eq->numeq_total),
           &ione,
           actinterface->continuity_eq->A.a.da[0],
           &(actinterface->continuity_eq->numeq_total),
           actinterface->continuity_eq->ipiv.a.iv,
           rhs1.a.da[0],
           &(actinterface->continuity_eq->numeq_total),
           &info
              );
  if (info!=0) dserror("Lapack solve failed");

  /*-----------------------------------------------------------------------*/
  /*            ALLOCATE MEMORY TO STORE THE SOLUTION AT THE NODES         */
  /*-----------------------------------------------------------------------*/
  
  /* -------------------memory allocation only in the first iteration step */
  if(step == 1)
  {
    a=0;
    b=0;
    /* ---------with this loops one obtains the maximum number of non-zero */
    /* ----------------------------------------cefficients per master node */
    for(i=0;i<actinterface->numnpm;i++)
    {
      a=0;
      for(j=0;j<actinterface->numnps;j++)
      {
        if((rhs1.a.da[i][j] > EPS10) || (rhs1.a.da[i][j] < -EPS10))
          a++;
      }
      if(a > b) 
        b=a;
    }
  
    /* --------------allocate memory for the master nodes on the interface */
    for(i=0; i<actinterface->numnpm; i++)
    {
      flag = 0;
      actnode1_1 = actinterface->nodem[i];
      if(actinterface->gnode_bound1mc != NULL)
         if(actnode1_1->gnode->node->Id == actinterface->gnode_bound1mc->node->Id)
            flag++;
      if(actinterface->gnode_bound2mc != NULL)
         if(actnode1_1->gnode->node->Id == actinterface->gnode_bound2mc->node->Id)
            flag++;

      if(flag>0)
      {
        if(actnode1_1->mtr_coeff.fdim < 15) 
        /* allocate memory to the node if the pointer mtr_coeff points to zero */
        {
          /*
          amdef("mrt coeff.",&(actnode1_1->mtr_coeff),15,int_faces->numint,"DA");
          */
          amdef("mrt coeff.",&(actnode1_1->mtr_coeff),actinterface->numnps,int_faces->numint,"DA");
          amzero(&(actnode1_1->mtr_coeff));
        }
      }
      else
      {
        /*
        amdef("mrt coeff.",&(actnode1_1->mtr_coeff),b,1,"DV");
        */
        amdef("mrt coeff.",&(actnode1_1->mtr_coeff),actinterface->numnps,1,"DV");
        amzero(&(actnode1_1->mtr_coeff));
      }
    }
  } /* end of if clause (step)*/
  
  /*-----------------------------------------------------------------------*/
  /*                      STORE THE SOLUTION AT THE NODES                  */
  /*-----------------------------------------------------------------------*/
  /* loop nodes of actinterface->nodem */
  for(i=0; i<actinterface->numnpm; i++)
  { 
    actnode1_1 = actinterface->nodem[i];
    /*amzero(&(actnode1_1->mtr_coeff));*/

    /* check if the node is at the intersection point between two interfaces */
    flag = 0;
    if(actinterface->gnode_bound1mc != NULL)
       if(actnode1_1->gnode->node->Id == actinterface->gnode_bound1mc->node->Id)
          flag++;
    if(actinterface->gnode_bound2mc != NULL)
       if(actnode1_1->gnode->node->Id == actinterface->gnode_bound2mc->node->Id)
          flag++;

    if(flag>0)
    {
      /* store the corresponding value of rhs1[][] at the node */
      b=0;
      for(k=0;k<actinterface->numnps;k++)
      {
        actnode1_1->mtr_coeff.a.da[k][r] = rhs1.a.da[i][k];
      }
      
    }
    else
    {
      /* store the corresponding value of rhs1[][] at the node */
      b=0;
      for(k=0;k<actinterface->numnps;k++)
      {
          actnode1_1->mtr_coeff.a.dv[k] = rhs1.a.da[i][k];
      }
      
    }
    /* update actnode1_1->nmb_zeros */
    a=0;
  } /* end of loop over the nodes of this element */
  
} /* end of loop over the interfaces (r) */

#ifdef DEBUG 
dstrc_exit();
#endif
} /*end of ssi_mortar_coeff */
  



/*!---------------------------------------------------------------------                                         
\brief detect intersection of 1-d nodal basis functions

<pre>                                                  firl / genk 1/04

With this subroutine we check the given bounds for a nonzero inter-
section. Resulting from this data we estimate the integration bounds.

</pre>
\param x1_1           DOUBLE         (i)   x coord. of first node of fct 1
\param x2_1           DOUBLE         (i)   x coord. of first node of fct 2
\param length1        DOUBLE         (i)   length of element one
\param length2        DOUBLE         (i)   length of element two
\param *b1            DOUBLE         (i)   lower integration bound
\param *b2            DOUBLE         (i)   upper integration bound
\param *intersection  INT            (i)   flag to detect intersection
\return 0                                                                             

------------------------------------------------------------------------*/
void ssi_detect_intersection(DOUBLE lambdr1_2, DOUBLE lambdr2_1, 
     DOUBLE lambdr2_2, DOUBLE nr1_2, DOUBLE nr2_1, DOUBLE nr2_2, 
     DOUBLE *b1, DOUBLE *b2, INT *intersection)
{

#ifdef DEBUG 
dstrc_enter("ssi_detect_intersection");
#endif

if(lambdr2_1 < lambdr2_2)
{
                                      /* Intersection Scenarios 1-d     */
                                      /*      1          2              */
  if ((lambdr2_1 < EPS8)&&            /*  1   o                         */
      (lambdr2_1 > -EPS8))            /*  2   o                         */
  {                                   
    *b1 = 0.0;
    if ((lambdr1_2<lambdr2_2) && (nr1_2 < nr2_2))          
    {                                 /*  1   o--------o                */
      *b2 = nr1_2;                    /*  2   o------------o            */
      *intersection = 1;
      goto end_of_if;
    }
    if ((lambdr1_2 < lambdr2_2+EPS10 && lambdr1_2 > lambdr2_2-EPS10) && 
        (nr1_2 < nr2_2+EPS10 && nr1_2 > nr2_2-EPS10 ))       
    {                                 /*  1   o------------o            */
      *b2 = nr2_2;                    /*  2   o------------o            */
      *intersection = 1;
      goto end_of_if;
    }
    if ((lambdr1_2 > lambdr2_2) && (nr1_2 > nr2_2))           
    {                                 /*  1   o------------o            */
      *b2 = nr2_2;                    /*  2   o--------o                */
      *intersection = 1;
      goto end_of_if;
    }
  }
  if (lambdr2_1 < -EPS8)              /*  1      o                      */
  {                                   /*  2   o            o            */
    *b1 = 0.0;
    if ((lambdr1_2 < lambdr2_2) && (nr1_2 < nr2_2))           
    {                                 /*  1      o------o               */
      *b2 = nr1_2;                    /*  2   o------------o            */
      *intersection = 1;
      goto end_of_if;
    }
    if ((lambdr1_2 < lambdr2_2+EPS10 && lambdr1_2 > lambdr2_2-EPS10) && 
        (nr1_2 < nr2_2+EPS10 && nr1_2 > nr2_2-EPS10 ))       
    {                                 /*  1      o---------o            */
      *b2 = nr1_2;                    /*  2   o------------o            */
      *intersection = 1;
      goto end_of_if;
    }
    if ((lambdr1_2 > lambdr2_2) && (nr1_2 > nr2_2))           
    {                                 /*  1      o------------o         */
      *b2 = nr2_2;                    /*  2   o------------o            */
      *intersection = 1;
      goto end_of_if;
    }
  }
  if (lambdr2_1 > EPS8)               /*  1   o            o            */
  {                                   /*  2      o                      */
    *b1 = nr2_1; 
    if ((lambdr1_2 < lambdr2_2) && (nr1_2 < nr2_2))           
    {                                 /*  1   o------------o            */
      *b2 = nr1_2;                    /*  2      o------------o         */
      *intersection = 1;
      goto end_of_if;
    }
    if ((lambdr1_2 < lambdr2_2+EPS10 && lambdr1_2 > lambdr2_2-EPS10) && 
        (nr1_2 < nr2_2+EPS10 && nr1_2 > nr2_2-EPS10 ))       
    {                                 /*  1   o------------o            */
      *b2 = nr1_2;                    /*  2      o---------o            */
      *intersection = 1;
      goto end_of_if;
    }
    if ((lambdr1_2 > lambdr2_2) && (nr1_2 > nr2_2))           
    {                                 /*  1   o------------o            */
      *b2 = nr2_2;                    /*  2      o------o               */
      *intersection = 1;
      goto end_of_if;
    }
  }
} /* end of if (lambdr2_1 < lambdr2_2)*/
if(lambdr2_1 > lambdr2_2)
{
                                      /* Intersection Scenarios 1-d     */
                                      /*      1          2              */
  if ((lambdr2_2 < EPS8)&&            /*  1   o                         */
      (lambdr2_2 > -EPS8))            /*  2   o                         */
  {                                   
    *b1 = 0.0;
    if ((lambdr1_2<lambdr2_1) && (nr1_2 < nr2_1))          
    {                                 /*  1   o--------o                */
      *b2 = nr1_2;                    /*  2   o------------o            */
      *intersection = 1;
      goto end_of_if;
    }
    if ((lambdr1_2 == lambdr2_1) && (nr1_2 == nr2_1))       
    {                                 /*  1   o------------o            */
      *b2 = nr2_1;                    /*  2   o------------o            */
      *intersection = 1;
      goto end_of_if;
    }
    if ((lambdr1_2 > lambdr2_1) && (nr1_2 > nr2_1))           
    {                                 /*  1   o------------o            */
      *b2 = nr2_1;                    /*  2   o--------o                */
      *intersection = 1;
      goto end_of_if;
    }
  }
  if (lambdr2_2 < -EPS8)              /*  1      o                      */
  {                                   /*  2   o            o            */
    *b1 = 0.0;
    if ((lambdr1_2 < lambdr2_1) && (nr1_2 < nr2_1))           
    {                                 /*  1      o------o               */
      *b2 = nr1_2;                    /*  2   o------------o            */
      *intersection = 1;
      goto end_of_if;
    }
    if ((lambdr1_2 == lambdr2_1) && (nr1_2 == nr2_1))      
    {                                 /*  1      o---------o            */
      *b2 = nr1_2;                    /*  2   o------------o            */
      *intersection = 1;
      goto end_of_if;
    }
    if ((lambdr1_2 > lambdr2_1) && (nr1_2 > nr2_1))           
    {                                 /*  1      o------------o         */
      *b2 = nr2_1;                    /*  2   o------------o            */
      *intersection = 1;
      goto end_of_if;
    }
  }
  if (lambdr2_2 > EPS8)               /*  1   o            o            */
  {                                   /*  2      o                      */
    *b1 = nr2_2; 
    if ((lambdr1_2 < lambdr2_1) && (nr1_2 < nr2_1))           
    {                                 /*  1   o------------o            */
      *b2 = nr1_2;                    /*  2      o------------o         */
      *intersection = 1;
      goto end_of_if;
    }
    if ((lambdr1_2 == lambdr2_1) && (nr1_2 == nr2_1))         
    {                                 /*  1   o------------o            */
      *b2 = nr1_2;                    /*  2      o---------o            */
      *intersection = 1;
      goto end_of_if;
    }
    if ((lambdr1_2 > lambdr2_1) && (nr1_2 > nr2_1))           
    {                                 /*  1   o------------o            */
      *b2 = nr2_1;                    /*  2      o------o               */
      *intersection = 1;
      goto end_of_if;
    }
  }
} /* end of if (lambdr2_1 < lambdr2_2)*/

end_of_if: ;

#ifdef DEBUG 
dstrc_exit();
#endif
} /*end of ssi_detect_intersection*/

  
  
/*!---------------------------------------------------------------------                                         
\brief put coupling forces from slave nodes to master nodes, mortar method

<pre>                                                  firl / genk 02/04

The internal forces of the slave nodes  are stored at 
sol_mf.a.da[4][..]. This subroutine assembles a array 
coupforc.a.da [actinterface->numnps][2]. 
In a second step the corresponding forces on the master nodes are 
computed. This forces depend on the mortar coefficients stored at the 
master nodes on the inteface in node->mtr_coeff.a.dv[]. The coupling 
forces on the master nodes are also stored at sol_mf.a.da[4][..].

</pre>
\param *masterfield   FIELD	     (i)   master field
\param *slavefield    FIELD          (i)   slave field
\param *int_faces    INTERFACES      (i)   the interfaces of the model
\return 0                                                                             

------------------------------------------------------------------------*/
void ssi_put_coupforc2mst(FIELD *masterfield, FIELD *slavefield, 
                          INTERFACES *int_faces)
{
ARRAY coupforc;                         /* array of coupling forces     */
DOUBLE *coupforc_a;
NODE *actnode;                          /* pointer to the actual node   */
INT  a, i, j, k;                        /* counters*/
DOUBLE b;                               /* parameter for nodal position */
                                        /* on the interface             */
INTERFACE *actinterface;                /* the actual interface         */

#ifdef DEBUG 
dstrc_enter("ssi_put_coupforc2mst");
#endif
/* erase sol_mf[4] for all master nodes */
solserv_sol_zero(masterfield, 0, 3, 4);

/* loop over the interfaces */
for(i=0; i<int_faces->numint; i++)
{
  actinterface = &int_faces->interface[i];
  /*------------------------------------------------- values to variables */
  coupforc_a = amdef("coupforc", &coupforc, actinterface->numnps, 2, "DA");
  amzero(&coupforc);
  
  /*----------------------------------------------------------------------*/
  /*                               PART I                                 */
  /*----------------------------------------------------------------------*/
  /* construct array of coupling forces according to the order of nodes in*/
  
  /* loop slave nodes at the interface */
  for(j=0; j<actinterface->numnps; j++)
  {
    actnode = actinterface->nodes[j];
    /* the factor b is necessary to scale the nodal force if the actual */
    /* node is an interface node at the intersection between two interfaces */
    b=1.0;
    if(actinterface->gnode_bound1sc != NULL &&
       actnode->gnode->node->Id == actinterface->gnode_bound1sc->node->Id)
      b=0.5;
    if(actinterface->gnode_bound2sc != NULL &&
       actnode->gnode->node->Id == actinterface->gnode_bound2sc->node->Id)
      b=0.5;

    coupforc.a.da[j][0] = actnode->sol_mf.a.da[4][0] * b;
    coupforc.a.da[j][1] = actnode->sol_mf.a.da[4][1] * b;
  }
  
  /*----------------------------------------------------------------------*/
  /*                               PART II                                */
  /*----------------------------------------------------------------------*/
  
  /* loop master nodes at the interface */
  for(j=0; j<actinterface->numnpm; j++)
  {
    actnode = actinterface->nodem[j];
    /*
    a = actnode->nmb_zeros;
    */
    a=0;
    /* loop the mortar coefficients */
    for(k=0;k<actnode->mtr_coeff.fdim;k++)
    {
      if((a+k)<coupforc.fdim)
      {
        /* check if the actual node is on the intersection between two interfaces*/
        if((actinterface->gnode_bound1mc != NULL && 
           actnode->gnode->node->Id == actinterface->gnode_bound1mc->node->Id) || 
          (actinterface->gnode_bound2mc != NULL && 
           actnode->gnode->node->Id == actinterface->gnode_bound2mc->node->Id))
        {    
          /* the index i runs with the interfaces */
          actnode->sol_mf.a.da[4][0] += actnode->mtr_coeff.a.da[k][i] * 
                                        coupforc.a.da[a+k][0];
          actnode->sol_mf.a.da[4][1] += actnode->mtr_coeff.a.da[k][i] * 
                                          coupforc.a.da[a+k][1];
        }
        else
        {    
          actnode->sol_mf.a.da[4][0] += actnode->mtr_coeff.a.dv[k] * 
                                         coupforc.a.da[a+k][0];
          actnode->sol_mf.a.da[4][1] += actnode->mtr_coeff.a.dv[k] * 
                                        coupforc.a.da[a+k][1];
        }
      }
    } /* end of loop over the mortar coefficients (k)*/     
  } /* end of loop over the master nodes at the interface */
} /* end of loop over the interfaces (i)*/

#ifdef DEBUG 
dstrc_exit();
#endif
} /* end of ssi_put_coupforc2mst */

  
  
/*!---------------------------------------------------------------------                                         
\brief put displacements from master nodes ot slave nodes, interpol. method

<pre>                                                  firl / genk 02/04

The actual displacements are relaxed by ssidym-->relax. Then the 
displacements of the slave nodes are linearly interpolated between the 
relaxed displacements of the master nodes. 

</pre>
\param *field_to      FIELD	     (i)   field which gets values
\param *field_from    FIELD          (i)   field with given values
\param place          INT            (i)   position in sol_mf[place][0..1]
\param *ssidyn        SSI_DYNAMIC    (i)   pointer to structure ssi_dynamic
\return 0                                                                             

------------------------------------------------------------------------*/
void ssi_interpolation_disp(FIELD *field_to, FIELD *field_from, INT place,
                       SSI_DYNAMIC *ssidyn)
{
ELEMENT *actele_from;                   /* pointer to the actual element */
NODE *actnode_to;                       /* pointer to the actual node */
                                        /* which gets values */
NODE *actnode_from1;                    /* pointer to the actual node 1*/
                                        /* which supplies values */
NODE *actnode_from2;                    /* pointer to the actual node 2*/
                                        /* which supplies values */
INT i, j, k, l;                         /* counters*/
DOUBLE from1from2[2];                   /* vec betw nods from1 & from 2*/
DOUBLE to_from1[2];                     /* vec betw nods to & from 1*/
DOUBLE to_from2[2];                     /* vec betw nods to & from 2*/
DOUBLE nfrom1from2;                     /* norm of vector from1from2 */
DOUBLE nto_from1;                       /* norm of vector to_from1 */
DOUBLE nto_from2;                       /* norm of vector to_from2 */
DOUBLE fac1, fac2;                      /* working values */


#ifdef DEBUG 
dstrc_enter("ssi_interpolation_disp");
#endif
/* clean vector sol_mf[place] of the nodes of field_to */
solserv_sol_zero(field_to, 0, 3, place);

/* loop nodes of field_to */
for(i=0; i<field_to->dis->numnp; i++)
{
  actnode_to = &(field_to->dis->node[i]);
  /*check if actnode_to is on the interface and not on the boundary of the interface*/
  if((actnode_to->gnode->ssicouple->ssi_couptyp != ssi_none) && 
    (actnode_to->gnode->ondesigntyp == ondline))
  {
    /* loop elements of field_from */
    for(j=0; j<field_from->dis->numele; j++)
    {
      actele_from = &field_from->dis->element[j];
      /* loop nodes of actele_from */
      for(k=0; k<actele_from->numnp; k++)
      {
        actnode_from1 = actele_from->node[k];
        /* check if actnode_from1 is a coupling node */
        if(actnode_from1->gnode->ssicouple->ssi_couptyp != ssi_none)
        {
          /* loop all other nodes of this element */
          for(l=0; l<actele_from->numnp; l++)
          {
            actnode_from2 = actele_from->node[l];
            /* check if actnode_from2 is a coupling node */
            if(actnode_from2->gnode->ssicouple->ssi_couptyp != ssi_none)
            {
              /*check if actnode_from1 != actnode_from2 */
              if(actnode_from1->Id != actnode_from2->Id)
              {
                /* compute relaxed displacement if place == 6 */
                if(place==6)
                {
                  actnode_from1->sol_mf.a.da[place][0] = 
                    ssidyn->relax * actnode_from1->sol_mf.a.da[0][0] +
                    (1-ssidyn->relax) * actnode_from1->sol_mf.a.da[1][0];
                  actnode_from1->sol_mf.a.da[place][1] = 
                    ssidyn->relax * actnode_from1->sol_mf.a.da[0][1] +
                    (1-ssidyn->relax) * actnode_from1->sol_mf.a.da[1][1];
                  actnode_from2->sol_mf.a.da[place][0] = 
                    ssidyn->relax * actnode_from2->sol_mf.a.da[0][0] +
                    (1-ssidyn->relax) * actnode_from2->sol_mf.a.da[1][0];
                  actnode_from2->sol_mf.a.da[place][1] = 
                    ssidyn->relax * actnode_from2->sol_mf.a.da[0][1] +
                    (1-ssidyn->relax) * actnode_from2->sol_mf.a.da[1][1];
                }
                /* compute the vectors between the nodes */
                /* the parameters of the vector depend on the actual position*/
                /* between this nodes; thus, the actual relaxed displacement */
                /* is added to the initial geometry */
                from1from2[0] = (actnode_from1->x[0]+actnode_from1->sol_mf.a.da[6][0]) -
                                (actnode_from2->x[0]+actnode_from2->sol_mf.a.da[6][0]); 
                from1from2[1] = (actnode_from1->x[1]+actnode_from1->sol_mf.a.da[6][1]) -
                                (actnode_from2->x[1]+actnode_from2->sol_mf.a.da[6][1]); 
                to_from1[0] = (actnode_to->x[0]+actnode_to->sol_mf.a.da[6][0]) -
                              (actnode_from1->x[0]+actnode_from1->sol_mf.a.da[6][0]); 
                to_from1[1] = (actnode_to->x[1]+actnode_to->sol_mf.a.da[6][1]) -
                              (actnode_from1->x[1]+actnode_from1->sol_mf.a.da[6][1]); 
                to_from2[0] = (actnode_to->x[0]+actnode_to->sol_mf.a.da[6][0]) -
                              (actnode_from2->x[0]+actnode_from2->sol_mf.a.da[6][0]); 
                to_from2[1] = (actnode_to->x[1]+actnode_to->sol_mf.a.da[6][1]) -
                              (actnode_from2->x[1]+actnode_from2->sol_mf.a.da[6][1]); 
                /* compute norms of the three vectors */
                nfrom1from2 = sqrt((from1from2[0]*from1from2[0]) + 
                                   (from1from2[1]*from1from2[1]));
                nto_from1 = sqrt((to_from1[0]*to_from1[0]) +
                                 (to_from1[1]*to_from1[1]));
                nto_from2 = sqrt((to_from2[0]*to_from2[0]) +
                                 (to_from2[1]*to_from2[1]));
                /* check if actnode_to is between actnode_from1 & actnode_from2*/
                if((nto_from1 < nfrom1from2)&&(nto_from2 < nfrom1from2))
                {
                  /* compute working values slope, shift */
                  fac1 = (actnode_from2->sol_mf.a.da[place][0]-
                          actnode_from1->sol_mf.a.da[place][0]) / nfrom1from2;
                  fac2 = actnode_from1->sol_mf.a.da[place][0];
                  /* apply interpolated value to actnode_to->sol_mf[place][0]*/
                  actnode_to->sol_mf.a.da[place][0] = 
                    fac1 * nto_from1 + fac2;
                  /* compute working values slope, shift */
                  fac1 = (actnode_from2->sol_mf.a.da[place][1]-
                          actnode_from1->sol_mf.a.da[place][1]) / nfrom1from2;
                  fac2 = actnode_from1->sol_mf.a.da[place][1];
                  actnode_to->sol_mf.a.da[place][1] = 
                    fac1 * nto_from1 + fac2;
                  if(place==6)
                  {
                    actnode_to->sol_mf.a.da[0][0] = 
                    actnode_to->sol_mf.a.da[place][0];
                    actnode_to->sol_mf.a.da[0][1] = 
                    actnode_to->sol_mf.a.da[place][1];
                  }
                  goto end_of_loop;
                }
              } /* end of if clause Id1 != Id2 */
            } /* end f if clause, if actnode_from2 is on the interface */
          } /* end of loop over other nodes of element actele */
        } /* end f if clause, if actnode_from1 is on the interface */
      } /* end of loop over nodes of actele */
    } /* end of loop over elements of field_from */
  } /* end of if clause if actnode_to is on the interface */
  /*check if actnode_to is on the interface and on the boundary of the interface*/
  if((actnode_to->gnode->ssicouple->ssi_couptyp != ssi_none) && 
    (actnode_to->gnode->ondesigntyp == ondnode))
  {
    /* loop nodes of field_from */
    for(j=0; j<field_from->dis->numnp; j++)
    {
      actnode_from1 = &field_from->dis->node[j];
      /* check if actnode_from1 is also on a dnode */
      if((actnode_from1->gnode->ssicouple->ssi_couptyp != ssi_none) &&
         (actnode_from1->gnode->ondesigntyp == ondnode))
      {
        /* compute relaxed displacement if place == 6 */
        if(place==6)
        {
          actnode_from1->sol_mf.a.da[place][0] = 
            ssidyn->relax * actnode_from1->sol_mf.a.da[0][0] +
            (1-ssidyn->relax) * actnode_from1->sol_mf.a.da[1][0];
          actnode_from1->sol_mf.a.da[place][1] = 
            ssidyn->relax * actnode_from1->sol_mf.a.da[0][1] +
            (1-ssidyn->relax) * actnode_from1->sol_mf.a.da[1][1];
        }
        /* compute vector between this two nodes */
        to_from1[0] = (actnode_to->x[0]+actnode_to->sol_mf.a.da[6][0]) -
                      (actnode_from1->x[0]+actnode_from1->sol_mf.a.da[6][0]); 
        to_from1[1] = (actnode_to->x[1]+actnode_to->sol_mf.a.da[6][1]) -
                      (actnode_from1->x[1]+actnode_from1->sol_mf.a.da[6][1]); 
        /* compute norm of this vector */
        nto_from1 = sqrt((to_from1[0]*to_from1[0]) +
                         (to_from1[1]*to_from1[1]));
        /* check if the norm is smaller than 0.01, this indicates */
        /* coincident nodes */
        if(nto_from1 < 0.01)
        {
         /* apply value from actnode_from1->sol_mf[place][0..1] directly */
         /* to actnode_to->sol_mf[place][0..1] */
          actnode_to->sol_mf.a.da[place][0] = 
          actnode_from1->sol_mf.a.da[place][0];
          actnode_to->sol_mf.a.da[place][1] = 
          actnode_from1->sol_mf.a.da[place][1];
          if(place==6)
          {
            actnode_to->sol_mf.a.da[0][0] = 
            actnode_to->sol_mf.a.da[place][0];
            actnode_to->sol_mf.a.da[0][1] = 
            actnode_to->sol_mf.a.da[place][1];
          }
          goto end_of_loop;
        }
      } /* end of if clause */
    } /* end of loop over the nodes of field_from (j) */
  } /* end of if clause */
  end_of_loop: ;
} /* end of loop over master nodes */

#ifdef DEBUG 
dstrc_exit();
#endif
} /* end of ssi_interpolation */



/*!---------------------------------------------------------------------                                         
\brief transfer forces between subdomains, interpol. method

<pre>                                                  firl / genk 02/04

The internal forces of field_from are stored in sol_mf[4][]. These forces
are transfered to the corresponding nodes of field_to. 

</pre>
\param *field_to      FIELD	     (i)   field which gets values
\param *field_from    FIELD          (i)   field with given values
\param place          INT            (i)   position in sol_mf[place][0..1]
\param *ssidyn        SSI_DYNAMIC    (i)   pointer to structure ssi_dynamic
\return 0                                                                             

------------------------------------------------------------------------*/
void ssi_interpolation_frc(FIELD *field_to, FIELD *field_from, INT place,
                       SSI_DYNAMIC *ssidyn)
{
ELEMENT *actele_to;                     /* pointer to the actual element */
NODE *actnode_to1;                      /* pointer to the actual node 1*/
                                        /* which gets values */
NODE *actnode_to2;                      /* pointer to the actual node 2*/
                                        /* which gets values */
NODE *actnode_from;                     /* pointer to the actual node */
                                        /* which supplies values */
INT i, j, k, l;                         /* counters*/
DOUBLE to1to2[2];                       /* vec betw nods from1 & from 2*/
DOUBLE from_to1[2];                     /* vec betw nods to & from 1*/
DOUBLE from_to2[2];                     /* vec betw nods to & from 2*/
DOUBLE nto1to2;                         /* norm of vector from1from2 */
DOUBLE nfrom_to1;                       /* norm of vector to_from1 */
DOUBLE nfrom_to2;                       /* norm of vector to_from2 */
DOUBLE fac1, fac2;                      /* working values */



#ifdef DEBUG 
dstrc_enter("ssi_interpolation_frc");
#endif

/* clean vector sol_mf[place] of the nodes of field_to */
solserv_sol_zero(field_to, 0, 3, place);

/* loop elements of field_to */
for(i=0; i<field_to->dis->numele; i++)
{
  actele_to = &(field_to->dis->element[i]);
  /*loop nodes of actele_to */
  for(j=0; j<actele_to->numnp; j++)
  {
    actnode_to1 = actele_to->node[j];
    /* check if actnode_to1 is a coupling node */
    if(actnode_to1->gnode->ssicouple->ssi_couptyp != ssi_none)
    {
      /* loop all other nodes of this element */
      for(k=0; k<actele_to->numnp; k++)
      {
        actnode_to2 = actele_to->node[k];
        /* check if actnode_to2 is a coupling node */
        if(actnode_to2->gnode->ssicouple->ssi_couptyp != ssi_none)
        {
          /*check if actnode_to1 != actnode_to2 and that j < k*/
          if((actnode_to1->Id != actnode_to2->Id)&&(j < k))
          {
            /* loop nodes of field_from */
            for(l=0; l<field_from->dis->numnp;l++)
            {
              actnode_from = &(field_from->dis->node[l]);
              /* check if actnode_from is a coupling node and check that 
                 actnode_from is on a dline */
              if((actnode_from->gnode->ssicouple->ssi_couptyp != ssi_none) && 
                 (actnode_from->gnode->ondesigntyp == ondline))
              {
                /* compute the vectors between the nodes */
                /* the parameters of the vector depend on the actual position*/
                /* between this nodes; thus, the actual relaxed displacement */
                /* is added to the initial geometry */
                to1to2[0] = 
                  (actnode_to2->x[0] + actnode_to2->sol_mf.a.da[6][0]) - 
                  (actnode_to1->x[0] + actnode_to1->sol_mf.a.da[6][0]);                
                to1to2[1] = 
                  (actnode_to2->x[1] + actnode_to2->sol_mf.a.da[6][1]) - 
                  (actnode_to1->x[1] + actnode_to1->sol_mf.a.da[6][1]);                
                from_to1[0] = 
                  (actnode_to1->x[0] + actnode_to1->sol_mf.a.da[6][0]) - 
                  (actnode_from->x[0] + actnode_from->sol_mf.a.da[6][0]);
                from_to1[1] = 
                  (actnode_to1->x[1] + actnode_to1->sol_mf.a.da[6][1]) - 
                  (actnode_from->x[1] + actnode_from->sol_mf.a.da[6][1]);
                from_to2[0] = 
                  (actnode_to2->x[0] + actnode_to2->sol_mf.a.da[6][0]) - 
                  (actnode_from->x[0] + actnode_from->sol_mf.a.da[6][0]);
                from_to2[1] = 
                  (actnode_to2->x[1] + actnode_to2->sol_mf.a.da[6][1]) - 
                  (actnode_from->x[1] + actnode_from->sol_mf.a.da[6][1]);
                /* compute norms of the three vectors */
                nto1to2 = sqrt((to1to2[0]*to1to2[0])+
                               (to1to2[1]*to1to2[1]));
                nfrom_to1 = sqrt((from_to1[0]*from_to1[0])+
                                 (from_to1[1]*from_to1[1]));
                nfrom_to2 = sqrt((from_to2[0]*from_to2[0])+
                                 (from_to2[1]*from_to2[1]));
                /* check if actnode_from is between actnode_to1 & actnode_to2*/
                if((nfrom_to1 < nto1to2)&&(nfrom_to2 < nto1to2))
                {
                  /* compute working factors */
                  fac1 = (nto1to2 - nfrom_to1) / nto1to2;
                  fac2 = (nto1to2 - nfrom_to2) / nto1to2;
                  /* compute nodal forces for the nodes of actele_to */
                  actnode_to1->sol_mf.a.da[4][0] += fac1 * actnode_from->sol_mf.a.da[4][0];
                  actnode_to1->sol_mf.a.da[4][1] += fac1 * actnode_from->sol_mf.a.da[4][1];
                  actnode_to2->sol_mf.a.da[4][0] += fac2 * actnode_from->sol_mf.a.da[4][0];
                  actnode_to2->sol_mf.a.da[4][1] += fac2 * actnode_from->sol_mf.a.da[4][1];
                  goto end_of_node_from; 
                }
              } /* end of if-clause */
              /* check if actnode_from is a coupling node and check that 
                 actnode_from is on a dnode */
              if((actnode_from->gnode->ssicouple->ssi_couptyp != ssi_none) && 
                 (actnode_from->gnode->ondesigntyp == ondnode))
              {
                /* compute the vectors between the nodes */
                /* the parameters of the vector depend on the actual position*/
                /* between this nodes; thus, the actual relaxed displacement */
                /* is added to the initial geometry */
                from_to1[0] = 
                  (actnode_to1->x[0] + actnode_to1->sol_mf.a.da[6][0]) - 
                  (actnode_from->x[0] + actnode_from->sol_mf.a.da[6][0]);
                from_to1[1] = 
                  (actnode_to1->x[1] + actnode_to1->sol_mf.a.da[6][1]) - 
                  (actnode_from->x[1] + actnode_from->sol_mf.a.da[6][1]);
                from_to2[0] = 
                  (actnode_to2->x[0] + actnode_to2->sol_mf.a.da[6][0]) - 
                  (actnode_from->x[0] + actnode_from->sol_mf.a.da[6][0]);
                from_to2[1] = 
                  (actnode_to2->x[1] + actnode_to2->sol_mf.a.da[6][1]) - 
                  (actnode_from->x[1] + actnode_from->sol_mf.a.da[6][1]);
                /* compute norms of the three vectors */
                nfrom_to1 = sqrt((from_to1[0]*from_to1[0])+
                                 (from_to1[1]*from_to1[1]));
                nfrom_to2 = sqrt((from_to2[0]*from_to2[0])+
                                 (from_to2[1]*from_to2[1]));
                /* check if actnode_to1 and actnode_from are coincident */
                if(nfrom_to1 < EPS8)
                {
                  /* compute nodal forces for the nodes of actele_to */
                  actnode_to1->sol_mf.a.da[4][0] += actnode_from->sol_mf.a.da[4][0];
                  actnode_to1->sol_mf.a.da[4][1] += actnode_from->sol_mf.a.da[4][1];
                  goto end_of_node_from; 
                }
                if(nfrom_to2 < EPS8)
                {
                  /* compute nodal forces for the nodes of actele_to */
                  actnode_to2->sol_mf.a.da[4][0] += actnode_from->sol_mf.a.da[4][0];
                  actnode_to2->sol_mf.a.da[4][1] += actnode_from->sol_mf.a.da[4][1];
                  goto end_of_node_from; 
                }
              } /* end of if-clause */
              end_of_node_from: ;
            } /* end of loop over nodes of field_from (l)*/
          } /* end of if actnode_to1 != actnode_to2 */
        } /* end of if actnode_to2 is a coupling node */
      } /* end of loop all other nodes of this element (k) */
    } /* end of check if actnode_to1 is a coupling node */
  } /*end of loop nodes of actele_to (j) */
} /* end of loop elements of field_to (i) */


#ifdef DEBUG 
dstrc_exit();
#endif
} /* end of ssi_interpolation_frc */






/*!---------------------------------------------------------------------                                         
\brief compute displacements for slave nodes, mortar method

<pre>                                                  firl / genk 02/04

The consistency condition is applied to evaluate the dirichlet b.c. for
the slave nodes. Here is assumed that the displacements of the master 
nodes are stored in sol_mf.a.da[0][..]. After the computation of the 
corresponding displacements for the slave nodes these values are also 
stored in sol_mf[0][..]. Here of course at the slave nodes. 
The algorithm is similar to the subroutine ssi_mortar_coeff.

</pre>
\param *masterfield   FIELD	     (i)   master field
\param *slavefield    FIELD          (i)   slave field
\param nmb_nods_mstint INT           (i)   nmb of nods on master interface
\param nmb_nods_slvint INT           (i)   nmb of nods on slave interface
\return 0                                                                             

------------------------------------------------------------------------*/
void ssi_calc_disp4slv(FIELD *masterfield, FIELD *slavefield, 
                       SSI_DYNAMIC *ssidyn, INTERFACES *int_faces)
{
DOUBLE gauss2[2][2];                /* gauss points and weights (2)     */
INT intersection;                   /* flag for intersection of nodal   */
                                    /* basis fct's 0=no, 1=yes          */
INT a,b,i,j,k,l,m,n,p,q,r;          /* counters                         */
NODE *actnode1_1, *actnode1_2, *actnode2_1, *actnode2_2;
                                    /* pointer to actual nodes          */
DOUBLE d1_1[2], d1_2[2], d2_1[2], d2_2[2];
                                    /* position vectors of nodes on the */
                                    /* interface */
DOUBLE nrmdr1_2[2];                 /* normed vector r1_2 */
DOUBLE nr1_2, nr2_1, nr2_2;         /* norm of vectors r1_2, r2_1, r2_2 */
DOUBLE r1_2[2], r2_1[2], r2_2[2];   /* position vectors relative to node*/
                                    /* d1_1 */
DOUBLE lambdr1_2, lambdr2_1, lambdr2_2; /* e.g lambdr1_2 =              */
                                        /*         r1_2[0] / nrmdr1_2[0]*/
DOUBLE b1, b2;                       /* integration bounds              */
DOUBLE m1, n1, m2, n2;               /* parameters for linear functions */

ARRAY rhs;                           /* right hand side of equation     */
DOUBLE *rhs_a;                       /* to compute the displacements    */

ARRAY rhs_d;                         /* additional rhs due to D.b.c.    */
DOUBLE *rhs_d_a;                     /*                                 */

ARRAY slv_nods_onint_work;           /* additional rhs due to D.b.c.    */
DOUBLE *slv_nods_onint_work_a;       /*                                 */

char uplo[1];                        /* character for lapack solver     */
INT info;                            /* flag for lapack solver          */
INT ione;                            /* number of right hand sides      */

ARRAY solu_arr;                      /* nodal displacement array        */
DOUBLE *solu_arr_a;

ARRAY hasdline;                      /* entry == 1 if node is on a dline*/
DOUBLE *hasdline_a;

ARRAY work_vec;                      /* working vector to rearrange the lhs*/
DOUBLE *work_vec_a;
INTERFACE *actinterface;             /* the actual interface */

ARRAY rel_d_2_1;                     /* the relaxed displacements of the */
DOUBLE *rel_d_2_1_a;                 /* master node 2_1 */ 
ARRAY rel_d_2_2;                     /* the relaxed displacements of the */
DOUBLE *rel_d_2_2_a;                 /* master node 2_2 */ 
                                     
#ifdef DEBUG 
dstrc_enter("ssi_calc_disp4slv");
#endif


/* -------------------------------------------------values to variables */
a=0; i=0; j=0; k=0; l=0; m=0; n=0; p=0;

gauss2[0][0] = 0.577350269189626;
gauss2[1][0] = -0.577350269189626;
gauss2[0][1] = 1.0;
gauss2[1][1] = 1.0;

/*----------------------------------------------------------------------*/
/* loop over the interfaces */
/*----------------------------------------------------------------------*/
for(r=0; r<int_faces->numint; r++)
{
  actinterface = &int_faces->interface[r];

  /* the rhs values are stored in a transposed matrix, this is neccesary */
  /* because the LAPACK solver is a FORTRAN routine, the memory storage */
  /* in FORTRAN is different with respect to C*/
  rhs_a = amdef("RHS",&rhs,2,actinterface->numnps,"DA");
  amzero(&rhs);
  slv_nods_onint_work_a = amdef("slv_nods_int_work",&slv_nods_onint_work,actinterface->numnps,1,"IV");
  amzero(&slv_nods_onint_work);
  rhs_d_a = amdef("add RHS",&rhs_d,actinterface->numnps,2,"DA");
  amzero(&rhs_d);
  work_vec_a = amdef("work",&work_vec,actinterface->numnps,1,"DV");
  amzero(&work_vec);
  solu_arr_a = amdef("solu",&solu_arr,actinterface->numnps,2,"DA");
  amzero(&solu_arr);
  hasdline_a = amdef("hasdline",&hasdline,actinterface->numnps,1,"IV");
  amzero(&hasdline);
  rel_d_2_1_a = amdef("rel_d_2_1",&rel_d_2_1,2,1,"DV");
  amzero(&rel_d_2_1);
  rel_d_2_2_a = amdef("rel_d_2_2",&rel_d_2_2,2,1,"DV");
  amzero(&rel_d_2_2);
  /* ---the pointer lhs_dens points on factorized matrix, this matrix has */
  /* -------------------------------their coefficients below the diagonal */
  uplo[0] = 'L'; 
  /* the unmber of right hand sides corresponds to the number of dof of the*/
  /* ---------------------------------------------------------master nodes */
  ione = masterfield->dis->node[0].numdf;
  
  
  /*loop elements of slavefield, biorthogonal shape functions, test space */
  for(i=0;i<actinterface->numnps-1;i++)
  {
    actnode1_1 = actinterface->nodes[i];
    actnode1_2 = actinterface->nodes[i+1];

    /* ----------------------------------------loop nodes of this element */
    for(j=0;j<2;j++)
    {
      /* the deformed geometry is used to set up the functions*/
      /* for the slave nodes the actual displacements are used */
      d1_1[0] = actnode1_1->x[0]+actnode1_1->sol_mf.a.da[6][0];
      d1_1[1] = actnode1_1->x[1]+actnode1_1->sol_mf.a.da[6][1];
      d1_2[0] = actnode1_2->x[0]+actnode1_2->sol_mf.a.da[6][0];
      d1_2[1] = actnode1_2->x[1]+actnode1_2->sol_mf.a.da[6][1];
      /* near the boundary of the interface we reduce the order of the */
      /* nodal basis function; hence the bounds in which this function */
      /* ----is defined also change, this is adjusted in the following */
      if(actnode1_1->gnode->node->Id == actinterface->gnode_bound1s->node->Id ||
         actnode1_1->gnode->node->Id == actinterface->gnode_bound2s->node->Id)
      {
        if(j==0)
        {
          d1_2[0]=d1_1[0]+0.5*(d1_2[0]-d1_1[0]);
          d1_2[1]=d1_1[1]+0.5*(d1_2[1]-d1_1[1]);
        }
        if(j==1)
        {
          d1_1[0]=d1_1[0]+0.5*(d1_2[0]-d1_1[0]);
          d1_1[1]=d1_1[1]+0.5*(d1_2[1]-d1_1[1]);
        }
      }
      if(actnode1_2->gnode->node->Id == actinterface->gnode_bound1s->node->Id ||
         actnode1_2->gnode->node->Id == actinterface->gnode_bound2s->node->Id)
      {
        if(j==0)
        {
          d1_2[0]=d1_1[0]+0.5*(d1_2[0]-d1_1[0]);
          d1_2[1]=d1_1[1]+0.5*(d1_2[1]-d1_1[1]);
        }
        if(j==1)
        {
          d1_1[0]=d1_1[0]+0.5*(d1_2[0]-d1_1[0]);
          d1_1[1]=d1_1[1]+0.5*(d1_2[1]-d1_1[1]);
        }
      }
      /* detect position of integral value in the vectors lhs and rhs */
      for(q=0;q<actinterface->numnps;q++)
      {
        if(j==0 && actnode1_1->Id == actinterface->nodes[q]->Id) 
        { 
          a=q;
          break;
        }
        if(j==1 && actnode1_2->Id == actinterface->nodes[q]->Id) 
        { 
          a=q;
          break;
        }
      }
      /* --loop elements of master domain, standard hat fcts, trial space */
      b=0;
      for(l=0;l<actinterface->numnpm-1;l++)
      {
        actnode2_1 = actinterface->nodem[l];
        actnode2_2 = actinterface->nodem[l+1];

        rel_d_2_1.a.dv[0] = actnode2_1->sol_mf.a.da[0][0];
        rel_d_2_1.a.dv[1] = actnode2_1->sol_mf.a.da[0][1];
        rel_d_2_2.a.dv[0] = actnode2_2->sol_mf.a.da[0][0];
        rel_d_2_2.a.dv[1] = actnode2_2->sol_mf.a.da[0][1];
        /* ------------------------------loop nodes of this slave element */
        for(m=0;m<2;m++)
        {
          intersection = 0;
          /* the deformed geometry is used to set up the functions*/
          /* for the master nodes the displacements of the last iteration */
          /* are used */
          d2_1[0] = actnode2_1->x[0]+actnode2_1->sol_mf.a.da[6][0];
          d2_1[1] = actnode2_1->x[1]+actnode2_1->sol_mf.a.da[6][1];
          d2_2[0] = actnode2_2->x[0]+actnode2_2->sol_mf.a.da[6][0];
          d2_2[1] = actnode2_2->x[1]+actnode2_2->sol_mf.a.da[6][1];
          /* -------------------------compute vectors relative to node 1_1 */
          r1_2[0] = d1_2[0] - d1_1[0];
          r1_2[1] = d1_2[1] - d1_1[1];
          r2_1[0] = d2_1[0] - d1_1[0];
          r2_1[1] = d2_1[1] - d1_1[1];
          r2_2[0] = d2_2[0] - d1_1[0];
          r2_2[1] = d2_2[1] - d1_1[1];
          /* ---------------------compute norm of vectors r1_2, r2_1, r2_2 */
          nr1_2 = sqrt((r1_2[0])*(r1_2[0])+(r1_2[1])*(r1_2[1]));
          nr2_1 = sqrt((r2_1[0])*(r2_1[0])+(r2_1[1])*(r2_1[1]));
          nr2_2 = sqrt((r2_2[0])*(r2_2[0])+(r2_2[1])*(r2_2[1]));
          /* ------------------------compute normed vector r1_2 = nrmdr1_2 */
          nrmdr1_2[0] = r1_2[0] / nr1_2;
          nrmdr1_2[1] = r1_2[1] / nr1_2;
          /* --------compute lambda values of the vectors r1_2, r2_1, r2_2 */
          /* the component with the biggest norm is taken to compute the */
          /* lambdas, thus the error due to floating point division is */ 
          /* minimized*/
          if(sqrt(r1_2[0]*r1_2[0]) > sqrt(r1_2[1]*r1_2[1])) 
          {
            lambdr1_2 = r1_2[0] / nrmdr1_2[0];
            lambdr2_1 = r2_1[0] / nrmdr1_2[0];
            lambdr2_2 = r2_2[0] / nrmdr1_2[0];
          }
          else
          {
            lambdr1_2 = r1_2[1] / nrmdr1_2[1];
            lambdr2_1 = r2_1[1] / nrmdr1_2[1];
            lambdr2_2 = r2_2[1] / nrmdr1_2[1];
          }
          /* Check for a possible nonzero intersection, goto end if such an*/
          /* ---------------------------------intersection is not possible */        
          if((lambdr2_1 >= nr1_2) && (lambdr2_2 >= nr1_2)) goto end_of_node2;
          if((lambdr2_1 < 0) && (lambdr2_2 < 0)) goto end_of_node2;
          
          /* ---------------------here we determine the integration bounds */
          ssi_detect_intersection(lambdr1_2, lambdr2_1, lambdr2_2, nr1_2, 
                                  nr2_1, nr2_2, &b1, &b2, &intersection);
  
          if(intersection == 1)
          {
            /* if we have a nonzero inters. we compute the nodal basis fcts*/
            /* -now we compute the nodal basis functions for the Lagrange */
            /* ---multiplier space, these functions fulfill the so called */
            /* biorthogonality relation, we refer to the Master Thesis of */
            /* ------------------------Matthias Firl for more information */
            if(j==0)
            {
              m1 = -3.0/lambdr1_2;
              n1 = 2.0;
            }
            else
            {
              m1 = 3.0/lambdr1_2;
              n1 = -1.0;
            }
            if(actnode1_1->gnode->node->Id == actinterface->gnode_bound1s->node->Id ||
               actnode1_1->gnode->node->Id == actinterface->gnode_bound2s->node->Id ||
               actnode1_2->gnode->node->Id == actinterface->gnode_bound1s->node->Id ||
               actnode1_2->gnode->node->Id == actinterface->gnode_bound2s->node->Id)
            {
              m1=0;
              n1=1;
            }
            /* loop over numdf */
            for(q=0;q<actnode2_1->numdf;q++)
            {
              /* ---calculate slope and value n of the nodal basis function, */ 
              /* here its the nodal basis function of the trace space of the */ 
              /* -----master elements on the interface. Hence, these are the */ 
              /* ---standard hat functions multiplied with the corresponding */
              /* ----------------------------------------displacement values */
              if(m==0)   
              {
                m2 = (rel_d_2_2.a.dv[q] - rel_d_2_1.a.dv[q]) / 
                     (lambdr2_2 - lambdr2_1);
                n2 = rel_d_2_1.a.dv[q] - m2 * lambdr2_1;
              }
              else
              {
                m2 = (rel_d_2_1.a.dv[q] - rel_d_2_2.a.dv[q]) / 
                     (lambdr2_2 - lambdr2_1);
                n2 = rel_d_2_1.a.dv[q] - m2 * lambdr2_1;
              }
              for(p=0;p<2;p++)
              {
                rhs.a.da[q][a] += (m1 * ((b1+b2)/2 + (b1-b2)/2*gauss2[p][0])+n1) *
                                  (m2 * ((b1+b2)/2 + (b1-b2)/2*gauss2[p][0])+n2) *
                                  ((b2-b1)/2) * gauss2[p][1];
              }
            } /* end of loop over numdf (q)*/
          } /* end of if clause (if intersection == 1)*/
          break;
          end_of_node2: ;
        } /* end of loop over the nodes of the master element (m) */
      } /* end of loop over the master elements (l)*/

    } /* end of loop over the nodes of this slave element (j) */
  } /* end of loop over the slave elements (i)*/
  

  /***********************************************************************/
  /*                                                                     */
  /*****************************  SOLUTION  ******************************/
  /*                                                                     */
  /***********************************************************************/

  /* In general two different approaches are applicable. The first one en-
     forces strong continuity at the boundary of the interface. Hence, the
     displacements of the master nodes at the boundary of the interface are
     directly transferred to the slave nodes at the boundary of the inter-
     face. 
     In the second approach the continuity equation is applied over the 
     whole interface. Hence, the continuity at the boundary of the interface 
     is also fulfilled in a weak sense. */   

  /*---------------------------------------------------------------------*/
  /* COMPUTATION OF NODAL DISPLACEMENTS WITH STRONG CONTINUITY AT THE    */
  /*                     BOUNDARY OF THE INTERFACE                       */
  /*---------------------------------------------------------------------*/
  
  /* fill vector solu with known displacements, the displacements of the */
  /* master node on the corresponding dnode */

#if 0

  /* loop over vector actinterface->conti_eq_save->ipiv */
  for(i=0; i<actinterface->numnps; i++)
  {
    actnode1_1 = actinterface->nodes[i];
    if(actnode1_1->gnode->Id == actinterface->gnode_bound1s->Id ||
       actnode1_1->gnode->Id == actinterface->gnode_bound2s->Id)
    {
      actnode2_1 = actnode1_1->gnode->mfcpnode[0];
      rel_d_2_1.a.dv[0] = ssidyn->relax * actnode2_1->sol_mf.a.da[0][0] +
                         (1-ssidyn->relax) * actnode2_1->sol_mf.a.da[1][0];
      rel_d_2_1.a.dv[1] = ssidyn->relax * actnode2_1->sol_mf.a.da[0][1] +
                         (1-ssidyn->relax) * actnode2_1->sol_mf.a.da[1][1];
      solu_arr.a.da[i][0] = rel_d_2_1.a.dv[0];
      solu_arr.a.da[i][1] = rel_d_2_1.a.dv[1];

      actinterface->conti_eq_save->ipiv.a.iv[i] = 1;
    }
    else
    {
      hasdline.a.iv[i] = 1;
      actinterface->conti_eq_save->ipiv.a.iv[i] = 0;
    }
    /* store a copy of vector actinterface->conti_eq_save->ipiv in slv_nods_onint_work */
    slv_nods_onint_work.a.iv[i] = actinterface->nodes[i]->Id;
  }
  
  /* rearrange the system of equations */
  /* loop over the vector actinterface->conti_eq_save->ipiv */
  for(i=0; i<actinterface->numnps; i++)
  {
    if(actinterface->conti_eq_save->ipiv.a.iv[i] == 0)
      goto end_01;
    /* loop other elements of vector actinterface->conti_eq_save->ipiv*/
    for(j=0; (j+i)<actinterface->numnps; j++)
    {
      if(actinterface->conti_eq_save->ipiv.a.iv[i+j] == 0)
      {
        /* change column in lhs matrix actinterface->conti_eq_save->A.a.da */
        for(k=0; k<actinterface->numnps; k++)
        {
          work_vec.a.dv[k] = actinterface->conti_eq_save->A.a.da[k][i];
          actinterface->conti_eq_save->A.a.da[k][i] = actinterface->conti_eq_save->A.a.da[k][i+j];
          actinterface->conti_eq_save->A.a.da[k][i+j] = work_vec.a.dv[k];
        }
        /* change entries in array solu_arr */
        work_doub = solu_arr.a.da[i][0];
        solu_arr.a.da[i][0] = solu_arr.a.da[i+j][0];
        solu_arr.a.da[i+j][0] = work_doub;
        work_doub = solu_arr.a.da[i][1];
        solu_arr.a.da[i][1] = solu_arr.a.da[i+j][1];
        solu_arr.a.da[i+j][1] = work_doub;
        /* change entry in vector hasdnode */
        work_int = actinterface->conti_eq_save->ipiv.a.iv[i];
        actinterface->conti_eq_save->ipiv.a.iv[i] = actinterface->conti_eq_save->ipiv.a.iv[i+j];
        actinterface->conti_eq_save->ipiv.a.iv[i+j] = work_int;
        goto end_01;
      }
    }
    end_01: ;
  }

  /* rearrange the system of equations */
  /* loop over the vector actinterface->conti_eq_save->ipiv */
  for(i=0; i<actinterface->numnps; i++)
  {
    if(hasdline.a.iv[i] == 1)
      goto end_02;
    /* loop other elements of vector actinterface->conti_eq_save->ipiv*/
    for(j=0; (j+i)<actinterface->numnps; j++)
    {
      if(hasdline.a.iv[i+j] == 1)
      {
        /* change row in lhs matrix actinterface->conti_eq_save->A.a.da */
        for(k=0; k<actinterface->numnps; k++)
        {
          work_vec.a.dv[k] = actinterface->conti_eq_save->A.a.da[i][k];
          actinterface->conti_eq_save->A.a.da[i][k] = actinterface->conti_eq_save->A.a.da[i+j][k];
          actinterface->conti_eq_save->A.a.da[i+j][k] = work_vec.a.dv[k];
        }
        /* change entries in array rhs */
        work_doub = rhs.a.da[0][i];
        rhs.a.da[0][i] = rhs.a.da[0][i+j];
        rhs.a.da[0][i+j] = work_doub;
        work_doub = rhs.a.da[1][i];
        rhs.a.da[1][i] = rhs.a.da[1][i+j];
        rhs.a.da[1][i+j] = work_doub;
        /* change entry in vector hasdline */
        work_int = hasdline.a.iv[i];
        hasdline.a.iv[i] = hasdline.a.iv[i+j];
        hasdline.a.iv[i+j] = work_int;
        /* change entry in vector slv_nods_onint_work */
        work_int = slv_nods_onint_work.a.iv[i];
        slv_nods_onint_work.a.iv[i] = slv_nods_onint_work.a.iv[i+j];
        slv_nods_onint_work.a.iv[i+j] = work_int;
        goto end_02;
      }
    }
    end_02: ;
  }
  
  /* loop the rearranged vector hasdnode to obtain the number of free nodes */
  for(i=0; i<actinterface->numnps; i++)
  {
    if(actinterface->conti_eq_save->ipiv.a.iv[i] == 1)
      work_int = i-1;
  }
  /*
    printf("\n");
    printf("Number of free dofs %2i", work_int);
    printf("\n");
  */
  /* calculate the additional values for the rhs vectors */
  for(i=0; i<work_int; i++)
  {
    for(j=work_int; j<actinterface->numnps; j++)
    {
      rhs_d.a.da[i][0] += actinterface->conti_eq_save->A.a.da[i][j] * solu_arr.a.da[j][0];
      rhs_d.a.da[i][1] += actinterface->conti_eq_save->A.a.da[i][j] * solu_arr.a.da[j][1];
    }
  }
  
  /* sum up the two rhs's */
  for(i=0; i<work_int; i++)
  {
    rhs_d.a.da[i][0] = -rhs_d.a.da[i][0] + rhs.a.da[0][i];
    rhs_d.a.da[i][1] = -rhs_d.a.da[i][1] + rhs.a.da[1][i];
  }
  
  /* compute the nodal displacements for the slave nodes at the interface */
  for(i=0; i<work_int; i++)
  {
    solu_arr.a.da[i][0] = rhs_d.a.da[i][0] / actinterface->conti_eq_save->A.a.da[i][i];
    solu_arr.a.da[i][1] = rhs_d.a.da[i][1] / actinterface->conti_eq_save->A.a.da[i][i];
  }
  
  /*-----------------------------------------------------------------------*/
  /*                      STORE THE SOLUTION AT THE NODES                  */
  /*-----------------------------------------------------------------------*/
  for(i=0; i<actinterface->numnps; i++)
  {
    actnode1_1 = actinterface->nodes[i];
    /* look for the position of this node in the array solu_arr */
    a=0;
    for(j=0; j<actinterface->numnps; j++)
    {
      if(actnode1_1->Id == slv_nods_onint_work.a.iv[j])
      {
        a=j;
        break;
      }
    }

    actnode1_1->sol_mf.a.da[0][0] = (ssidyn->relax * solu_arr.a.da[a][0]) + 
                                    (1 - ssidyn->relax) * 
                                     actnode1_1->sol_mf.a.da[6][0];
    actnode1_1->sol_mf.a.da[0][1] = (ssidyn->relax * solu_arr.a.da[a][1]) + 
                                    (1 - ssidyn->relax) * 
                                     actnode1_1->sol_mf.a.da[6][1];

    /* copy the displacements to field sol_mf[6], in the */ 
    /* subroutine ssi_mortar_coeff this values are needed */
    actnode1_1->sol_mf.a.da[6][0] = actnode1_1->sol_mf.a.da[0][0];
    actnode1_1->sol_mf.a.da[6][1] = actnode1_1->sol_mf.a.da[0][1];
  }

#endif

  /*---------------------------------------------------------------------*/
  /* COMPUTATION OF NODAL DISPLACEMENTS WITH WEAK CONTINUITY AT THE      */
  /*                     BOUNDARY OF THE INTERFACE                       */
  /*---------------------------------------------------------------------*/
  /*---------------------------------------------------------------------*/
  /*                     SOLUTION USING LAPACK LU-DECOMPOSITION          */
  /*---------------------------------------------------------------------*/

  /* solution for different rhs's*/
  
  /* here the second part of the LAPACK solver is used. In actinterface->
     continuity_eq->A the factorized lhs matrix is stored. The 
     factorization is computed in the subroutine ssi_mortar_coeff */
  info=1;
  dsytrs(
           uplo,
           &(actinterface->continuity_eq->numeq_total),
           &ione,
           actinterface->continuity_eq->A.a.da[0],
           &(actinterface->continuity_eq->numeq_total),
           actinterface->continuity_eq->ipiv.a.iv,
           rhs.a.da[0],
           &(actinterface->continuity_eq->numeq_total),
           &info
              );
  if (info!=0) dserror("Lapack solve failed");
  
  /*-----------------------------------------------------------------------*/
  /*                      STORE THE SOLUTION AT THE NODES                  */
  /*-----------------------------------------------------------------------*/
  for(i=0; i<actinterface->numnps; i++)
  {
    actnode1_1 = actinterface->nodes[i];

    actnode1_1->sol_mf.a.da[0][0] = ssidyn->relax * rhs.a.da[0][i] + 
                                    (1-ssidyn->relax) * 
                                    actnode1_1->sol_mf.a.da[6][0];
    actnode1_1->sol_mf.a.da[0][1] = ssidyn->relax * rhs.a.da[1][i] + 
                                    (1-ssidyn->relax) * 
                                    actnode1_1->sol_mf.a.da[6][1];

    /* copy the displacements to field sol_mf[6], in the */ 
    /* subroutine ssi_mortar_coeff this values are needed */
    actnode1_1->sol_mf.a.da[6][0] = actnode1_1->sol_mf.a.da[0][0];
    actnode1_1->sol_mf.a.da[6][1] = actnode1_1->sol_mf.a.da[0][1];
  }


} /* end of loop over the interfaces (r) */

#ifdef DEBUG 
dstrc_exit();
#endif
} /*end of ssi_calc_disp4slv */

                                                                                                                                        



/*!---------------------------------------------------------------------                                         
\brief allocate dirich condition for each slave gnode at the interface

<pre>                                                  firl / genk 02/04

In general each slave node at the interface has its own dirichlet 
condition. The value of this condition results from the coupling 
algorithm. Per default all nodes on a dirich line obtain their Dirichlet
boundary condition from the dline. Hence, all nodes at this dline get the 
same Dirichlet boundary condition.

</pre>
\param *slavefield    FIELD          (i)   slave field
\return 0                                                                             

------------------------------------------------------------------------*/
void ssi_alloc_dirich(FIELD *slavefield)
{
INT i    ;                            /* counter */
GNODE *actgnode;                      /* actual gnode  */
NODE *actnode;                        /* actual node  */


#ifdef DEBUG 
dstrc_enter("ssi_alloc_dirich");
#endif
/*loop all nodes of slavefield*/
for(i=0; i<slavefield->dis->numnp; i++)
{
  actnode=&(slavefield->dis->node[i]);
  actgnode = actnode->gnode;  
  /*check if actgnode is a coupling node*/
  if(actgnode->ssicouple == NULL) continue;
  if(actgnode->ssicouple->ssi_couptyp == ssi_slave)
  {
    actgnode->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
    if (!actgnode->dirich) dserror("Allocation of memory failed");  
    amdef("onoff",&(actgnode->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
    amzero(&(actgnode->dirich->dirich_onoff));
    amdef("val",&(actgnode->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
    amdef("curve",&(actgnode->dirich->curve),MAXDOFPERNODE,1,"IV"); 
    amdef("function",&(actgnode->dirich->funct),MAXDOFPERNODE,1,"IV");
    amzero(&(actgnode->dirich->dirich_val));
    amzero(&(actgnode->dirich->curve));
    amzero(&(actgnode->dirich->funct));
    /*----------------------------------- initialise for ssi-coupling */
    actgnode->dirich->dirich_onoff.a.iv[0] = 1;   
    actgnode->dirich->dirich_onoff.a.iv[1] = 1;       
    actgnode->dirich->dirich_type=dirich_SSI;
  }
}

#ifdef DEBUG 
dstrc_exit();
#endif
}


#endif
/*! @} (documentation module close)*/
