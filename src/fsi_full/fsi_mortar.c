/*!----------------------------------------------------------------------
\file
\brief service routines for fsi algorithms

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*! 
\addtogroup FSI
*//*! @{ (documentation module open)*/
#ifdef D_FSI
#ifdef D_MORTAR
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fsi_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"
#include "../wall1/wall1.h"
#include "../wall1/wall1_prototypes.h"
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
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
 extern ALLDYNA      *alldyn;
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



/*!--------------------------------------------------------------------- 
\brief initialise interface for fsi coupling

<pre>                                                         firl 04/04
                                                             chfoe 07/04

 here the number of interfaces and the number of dlines per interface
 is estimated
 			     
</pre>   

\param *masterfield   FIELD         (i)      structure field 
\param *slavefield    FIELD         (i)      fluid field     
\param *int_faces     INTERFACES    (i)      pointer to interface strut.

\return void 

------------------------------------------------------------------------*/
void fsi_initcoupling_intfaces( 
                                FIELD       *masterfield,
                                FIELD       *slavefield, 
                                INTERFACES  *int_faces
		              )
{
INT     i,j,a;                   /* simply some counters                */
GLINE  *actgline;                /* actual glines                       */
ARRAY   int_ids;                 /* a vector with the interface Ids     */
DOUBLE  *int_ids_a;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("fsi_initcoupling_intfaces");
#endif

/*------------ detect the number of glines with a coupling condition ---*/
a = 0;
for(i=0; i<masterfield->dis->ngline; i++)/* loop glines of master field */
{
  actgline = &(masterfield->dis->gline[i]);
  if(actgline->fsicouple == NULL) continue;
  if(actgline->fsicouple->fieldtyp == structure)
    a++;
}
/*--- store the couplingId's of the glines with a coupling condition ...
  ----------------------------------- ... in the vector int_ids.a.iv ---*/
int_ids_a = amdef("int_ids", &int_ids,a,1,"IV");
amzero(&int_ids);
a=0;
for(i=0; i<masterfield->dis->ngline; i++)
{
  actgline = &(masterfield->dis->gline[i]);
  if(actgline->fsicouple == NULL) continue;  
  if(actgline->fsicouple->fieldtyp == structure)
  {
    int_ids.a.iv[a] = actgline->fsicouple->fsi_coupleId;
    a++;
  }
}
/*------------------------------------------ sort the vector int_ids ---*/
a=0;
for(i=0; i<int_ids.fdim; i++)
{
  for(j=0; j<int_ids.fdim-1; j++)
  {
    if(int_ids.a.iv[j]>int_ids.a.iv[j+1])
    {
      a = int_ids.a.iv[j];
      int_ids.a.iv[j] = int_ids.a.iv[j+1];
      int_ids.a.iv[j+1] = a;
    }  
  }
}
/*---------------------------------- detect the number of interfaces ---*/
a=1;
for(i=0; i<int_ids.fdim-1; i++)
{
  if(int_ids.a.iv[i] != int_ids.a.iv[i+1])
    a++;
}

/*-------------- store the number of interfaces in int_faces->numint ---*/
int_faces->numint = a;

/*----------------- allocate memory for the vector of interface Id's ---*/
int_faces->int_ids = (INT*)CCACALLOC(int_faces->numint,
                     sizeof(INT));

/* store the interface Id's in the vector int_faces->int_ids */
/* loop over the vector int_ids.a.iv which contains the interface Id's*/
a=1;
for(i=0; i<int_ids.fdim-1; i++)
{
  if(i==0)
    int_faces->int_ids[i] = int_ids.a.iv[i];
  if(int_ids.a.iv[i] != int_ids.a.iv[i+1] && int_ids.a.iv[i+1] != 0)
  {
    int_faces->int_ids[a] = int_ids.a.iv[i+1];
    a++;
  }
}

/*--------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fsi_initcoupling_intfaces */


/*!---------------------------------------------------------------------                                         
\brief detect number of nodes on slave and master interface

<pre>                                                  firl / genk 3/04

In this subroutine the interface data are stored in the structure 
interfaces.  

</pre>
\param *masterfield     FIELD	      (i)   structure field
\param *slavefield      FIELD         (i)   fluid field
\param *int_faces       INTERFACES    (i)   the interfaces of the model

\return 0                                                                             

------------------------------------------------------------------------*/
void fsi_init_interfaces(FIELD *masterfield, FIELD *slavefield,
                         INTERFACES *int_faces)
{
INT i,j,a,l,b;                      /* counters */
INT dof;                            /* the dof of actnode */
NODE *actnode;                      /* actual node under consideration */
NODE *actnode1;                     /* actual node under consideration */
NODE *node_work;
ELEMENT *actele;                    /* actual element under consideration */
INT k;                              /* number of dline */
GLINE *actgline;                    /* actual gline under consideration */
GNODE *actgnode;                    /* actual gnode under consideration */
INTERFACE *actinterface, *actinterface2;/* pointer to the actual interface */
DOUBLE work;
DOUBLE dist1, dist1x, dist1y, dist2, dist2x, dist2y; /* nodal distances */

FSI_DYNAMIC *fsidyn;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("fsi_init_interfaces");
#endif

/*--------------------------------------------------- initialisation ---*/
fsidyn = alldyn[genprob.numaf+1].fsidyn;

/*---------------------------------- assign the Id to the interfaces ---*/
for(i=0; i<int_faces->numint;i++)
{
  int_faces->interface[i].Id = int_faces->int_ids[i];
}

/*----------------------------------------------------------------------*/
for(i=0; i<int_faces->numint;i++) /* loop the interfaces */
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
     if(actgline->fsicouple == NULL) continue;
     
     if((actgline->fsicouple->fieldtyp == structure) &&
        (actgline->fsicouple->fsi_coupleId == actinterface->Id))
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
       if(actgline->fsicouple == NULL) continue;
       if((actgline->fsicouple->fieldtyp == structure) &&
          (actgline->fsicouple->fsi_coupleId == actinterface->Id))
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
    if(actgline->fsicouple == NULL) continue;
    if((actgline->fsicouple->fieldtyp == structure) &&
       (actgline->fsicouple->fsi_coupleId == actinterface->Id))
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
           if(actinterface->nodem[l] == actnode)
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
     if (actgline->fsicouple == NULL) continue;
     if((actgline->fsicouple->fieldtyp == fluid) &&
        (actgline->fsicouple->fsi_coupleId == actinterface->Id))
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

  /* ----III define the pointers in actinterface->elements to the fluid */
  /* ---------------------------------------- elements at the interface */
  a = 0;
  for(j=0; j<slavefield->dis->numele; j++)
  {
     actele = &slavefield->dis->element[j];
     /* loop the glines of this element */
     for(k=0; k<actele->g.gsurf->ngline; k++)
     {
       actgline = actele->g.gsurf->gline[k];
       if(actgline->fsicouple == NULL) continue;
       if((actgline->fsicouple->fieldtyp == fluid) &&
          (actgline->fsicouple->fsi_coupleId == actinterface->Id))
       {
          actinterface->elements[a] = actele;
          a++;
          break;
       }
     }     
  }
  /* -IV- define the pointers in actinterface->nodes to the fluid nodes */
  /* ------------------------------------------------- at the interface */
  a = 0;
  for(j=0; j<slavefield->dis->ngline; j++)
  {
    actgline = &slavefield->dis->gline[j];
   if(actgline->fsicouple == NULL) continue;
   if((actgline->fsicouple->fieldtyp == fluid) &&
       (actgline->fsicouple->fsi_coupleId == actinterface->Id))
    /* if both gnodes at gline have fsicouple condition */
    {
       /* loop gnodes of this gline */
       for(k=0; k<2; k++)
       {
         actgnode = actgline->gnode[k];
         actnode = actgnode->node;
         /* loop vector of actinterface->nodes[] to detect if actnode is */
         /* already there; this happens, because two glines share a node */
         for(l=0; l<actinterface->numnps; l++)
         {
           if(actinterface->nodes[l] == actnode)
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
      actgnode = actele->node[k]->gnode->mfcpnode[2]->gnode;

      /* loop glines of actgnode */
      a=0;
      for(l=0; l<actgnode->ngline; l++)
      {
        actgline = actgnode->gline[l];
        if(actgline->fsicouple == NULL) continue;
        if(actgline->fsicouple->fsi_coupleId == actinterface->Id)
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
        if(actgline->fsicouple == NULL) continue;
        if(actgline->fsicouple->fsi_coupleId == actinterface->Id)
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
    if(actnode->Id == actinterface->gnode_bound1s->node->Id || 
       actnode->Id == actinterface->gnode_bound2s->node->Id)
    {
      /* loop the other interfaces */
      for(k=0; k<int_faces->numint; k++)
      {
        actinterface2 = &int_faces->interface[k];
        if(actinterface->Id == actinterface2->Id)
          goto end_of_this_interface;
        else
        {
          if(actnode->Id == actinterface2->gnode_bound1s->node->Id ||
             actnode->Id == actinterface2->gnode_bound2s->node->Id)
          {
            if(actinterface->gnode_bound1sc == NULL)
              actinterface->gnode_bound1sc = actnode->gnode->mfcpnode[2]->gnode;
            else 
              actinterface->gnode_bound2sc = actnode->gnode->mfcpnode[2]->gnode;
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
    if(actnode->Id == actinterface->gnode_bound1m->node->Id || 
       actnode->Id == actinterface->gnode_bound2m->node->Id)
    {
      /* loop the other interfaces */
      for(k=0; k<int_faces->numint; k++)
      {
        actinterface2 = &int_faces->interface[k];
        if(actinterface->Id == actinterface2->Id)
          goto end_of_this_interface1;
        else
        {
          if(actnode->Id == actinterface2->gnode_bound1m->node->Id ||
             actnode->Id == actinterface2->gnode_bound2m->node->Id)
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

  /* modify the entries in fsidyn->sid.a.iv, they are originally stored in 
     fsi_struct_intdofs(); but this original algorithm does not work 
     correctly for non-conforming discretizations                          */
  
  /* loop over the structure nodes of the interface */
  b=0;
  for(j=0; j<actinterface->numnpm; j++)
  {
    actnode = actinterface->nodem[j];
    actnode1 = actnode->gnode->mfcpnode[2];
    /* loop over the dofs of actnode */
    for(k=0; k<actnode->numdf; k++)
    {
      if(actnode->gnode->dirich == NULL)
      {
        dof = actnode->dof[k];
        fsidyn->sid.a.iv[dof] = 1;
        b++;
      }
      else if(actnode->gnode->dirich->dirich_onoff.a.iv[k] != 1)
      {
        dof = actnode->dof[k];
        fsidyn->sid.a.iv[dof] = 1;
        b++;
      }
    }
  }
  fsidyn->numsid = b;
} /* end og loop over the interfaces (i) */
} /* end of fsi_init_interfaces */


/*!---------------------------------------------------------------------                                         
\brief compute coefficients for mortar method

<pre>                                                  firl / genk 1/04

Here is assembled a coefficient matrix and a right hand side. One obtains 
the coefficients (called zeta_j^r) as solution of this system of 
equations. It is referred to the master thesis of Matthias Firl for more 
information.
Later the factorized coefficient matrix is needed again. Thus, the pointer
to the dense structure is a parameter of the function call.

</pre>
\param *fsidyn      FSI_DYNAMIC      (i)   pointer to fsi dynamic
\param int_faces    INTERFACES       (i)   the interfaces of the model
\return 0                                                                             

------------------------------------------------------------------------*/
void fsi_mortar_coeff(FSI_DYNAMIC *fsidyn, INTERFACES *int_faces)
{
DOUBLE gauss2[2][2];                /* gauss points and weights (2)*/
INT intersection;                   /* flag for intersection of nodal*/
                                    /* basis fct's 0=no, 1=yes*/
INT info;                           /* flag for lapack solver */
INT ione;                           /* number of right hand sides*/
INT a,b,i,j,k,l,m,n,p,q,r,act;      /* counters */
INT numnods, switch_act;            /* number of elements of the actual side*/
                                    /* of actinterface */            
INT flag;

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
DOUBLE lambdr1_2_xy[2];              /* lambdr1_2[0] = r1_2[0]/nrmdr1_2[0]*/
                                     /* lambdr1_2[1] = r1_2[1]/nrmdr1_2[1]*/
DOUBLE lambdr2_1_xy[2];              /* lambdr2_1[0] = r2_1[0]/nrmdr2_1[0]*/
                                     /* lambdr2_1[1] = r2_1[1]/nrmdr2_1[1]*/
DOUBLE lambdr2_2_xy[2];              /* lambdr2_2[0] = r2_2[0]/nrmdr2_2[0]*/
                                     /* lambdr2_2[1] = r2_2[1]/nrmdr2_2[1]*/
DOUBLE b1, b2;                       /* integration bounds */
DOUBLE m1, n1, m2, n2;               /* parameters for linear functions */

ARRAY rhs1;                          /* right hand side of system of eq.  */
DOUBLE *rhs;                         /* to compute the mortar coeff.      */

INTERFACE *actinterface;             /* the actual interface */

#ifdef DEBUG 
dstrc_enter("fsi_mortar_coeff");
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
/* matrix; this matrix is needed in the routine fsi_calc_disp4slv again */
/* furthermore the integer vector ipiv is used to store the number 1 if the */
/* actual node is on a dnode; if the actual node is on a dline a zero is stored*/
/* this is also needed in fsi_calc_disp4slv to indicate the boundary of the */
/* interface */

  amzero(&(actinterface->conti_eq_save->A));

  /* ---In the first loop is computed the left hand side and in the second */
  /* the right hand side; hence, act == 0 -> lhs, act == 1 -> rhs */
  for(act=0;act<2;act++)
  {
    if(act==0) /* slave nodes */
    {
      actnodes = 
      (NODE**)CCACALLOC(actinterface->numnps,sizeof(NODE));
      for(i=0; i<actinterface->numnps; i++)
      {
        actnodes[i] = actinterface->nodes[i]->gnode->mfcpnode[2];
        dsassert(actnodes[i]!=NULL,"cannot read from NULL-pointer!\n");
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
        if(act == 0) /* slave elements */
        {
          /* the deformed geometry is used to set up the functions*/
          d1_1[0] = actnode1_1->x[0]+actnode1_1->sol_mf.a.da[1][0];
          d1_1[1] = actnode1_1->x[1]+actnode1_1->sol_mf.a.da[1][1];
          d1_2[0] = actnode1_2->x[0]+actnode1_2->sol_mf.a.da[1][0];
          d1_2[1] = actnode1_2->x[1]+actnode1_2->sol_mf.a.da[1][1];
        }
        else /* actual elements = master elements */
        {
          /* the deformed geometry is used to set up the functions*/
          /* computation of the relaxed interface displacements */
          d1_1[0] = actnode1_1->x[0] + actnode1_1->sol_mf.a.da[1][0];
          d1_1[1] = actnode1_1->x[1] + actnode1_1->sol_mf.a.da[1][1];
          d1_2[0] = actnode1_2->x[0] + actnode1_2->sol_mf.a.da[1][0];
          d1_2[1] = actnode1_2->x[1] + actnode1_2->sol_mf.a.da[1][1];
        }
        /* -----------loop elements of slave domain, biorthogonal shp. fcts */
        b=0;
        for(l=0;l<actinterface->numnps-1;l++)
        {
          intersection = 0;
          actnode2_1 = actinterface->nodes[l]->gnode->mfcpnode[2];
          dsassert(actnode2_1!=NULL,"cannot read from NULL-pointer!\n");
          actnode2_2 = actinterface->nodes[l+1]->gnode->mfcpnode[2];
          dsassert(actnode2_2!=NULL,"cannot read from NULL-pointer!\n");
          for(m=0; m<2; m++)
          {
            /* the deformed geometry is used to set up the functions*/
            /* d(n+1) are stored in sol_mf.a.da[1][.] */
            d2_1[0] = actnode2_1->x[0]+actnode2_1->sol_mf.a.da[1][0];
            d2_1[1] = actnode2_1->x[1]+actnode2_1->sol_mf.a.da[1][1];
            d2_2[0] = actnode2_2->x[0]+actnode2_2->sol_mf.a.da[1][0];
            d2_2[1] = actnode2_2->x[1]+actnode2_2->sol_mf.a.da[1][1];
            /* near the boundary of the interface we reduce the order of the */
            /* nodal basis function; hence the bounds in which this function */
            /* ----is defined also change, this is adjusted in the following */
            if(actnode2_1->Id == actinterface->gnode_bound1s->node->Id ||
               actnode2_1->Id == actinterface->gnode_bound2s->node->Id)
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
            if(actnode2_2->Id == actinterface->gnode_bound1s->node->Id ||
               actnode2_2->Id == actinterface->gnode_bound2s->node->Id) 
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
            /* x - component */
            lambdr1_2_xy[0] = r1_2[0] / nrmdr1_2[0];
            lambdr2_1_xy[0] = r2_1[0] / nrmdr1_2[0];
            lambdr2_2_xy[0] = r2_2[0] / nrmdr1_2[0];
            /* y - component */
            lambdr1_2_xy[1] = r1_2[1] / nrmdr1_2[1];
            lambdr2_1_xy[1] = r2_1[1] / nrmdr1_2[1];
            lambdr2_2_xy[1] = r2_2[1] / nrmdr1_2[1];
            /* check for a possible intersection */
            if(lambdr2_1_xy[0]/lambdr1_2_xy[0] <= EPS8 && 
               lambdr2_2_xy[0]/lambdr1_2_xy[0] <= EPS8) 
                 goto end_of_node2;
            if(nr2_1 > nr1_2+EPS6 && lambdr2_1_xy[0] > lambdr1_2_xy[0] &&
               nr2_2 > nr1_2+EPS6 && lambdr2_2_xy[0] > lambdr1_2_xy[0])
               goto end_of_node2;
            /* look for the best lambda for fsi_detect_intersection(),
               the goal of this choice is to minimize the error due to 
               floating point division */
            if((r1_2[0] > sqrt(r1_2[1]*r1_2[1])) || (r1_2[0]<-sqrt(r1_2[1]*r1_2[1])))
            {
              lambdr1_2 = lambdr1_2_xy[0];
              lambdr2_1 = lambdr2_1_xy[0];
              lambdr2_2 = lambdr2_2_xy[0];
            }
            else
            {
              lambdr1_2 = lambdr1_2_xy[1];
              lambdr2_1 = lambdr2_1_xy[1];
              lambdr2_2 = lambdr2_2_xy[1];
            }
            /* Check for a possible nonzero intersection, goto end if such an*/
            /* ---------------------------------intersection is not possible */        
            if((lambdr2_1 >= nr1_2) && (lambdr2_2 >= nr1_2)) goto end_of_node2;
            if((lambdr2_1 < 0) && (lambdr2_2 < 0)) goto end_of_node2;
            
            /* ---------------------here we determine the integration bounds */
            fsi_detect_intersection(lambdr1_2, lambdr2_1, lambdr2_2, nr1_2, 
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
    
              if(actnode2_1->Id == actinterface->gnode_bound1s->node->Id ||
                 actnode2_1->Id == actinterface->gnode_bound2s->node->Id ||
                 actnode2_2->Id == actinterface->gnode_bound1s->node->Id ||
                 actnode2_2->Id == actinterface->gnode_bound2s->node->Id) 
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
                  if(j == 0 && actnode1_1->Id == actinterface->nodes[q]->gnode->mfcpnode[2]->Id)
                  /*if(out_node1->Id == actinterface->nodes[q]->Id)*/ 
                  { 
                    a=q;
                    break;
                  }
                  if(j == 1 && actnode1_2->Id == actinterface->nodes[q]->gnode->mfcpnode[2]->Id)
                  /*if(out_node1->Id == actinterface->nodes[q]->Id)*/ 
                  { 
                    a=q;
                    break;
                  }
                }
                for(q=0;q<actinterface->numnps;q++)
                {
                  if(m==0 && actnode2_1->Id == actinterface->nodes[q]->gnode->mfcpnode[2]->Id)
                  /*if(out_node2->Id == actinterface->nodes[q]->Id) */
                  { 
                    b=q;
                    break;
                  }
                  if(m==1 && actnode2_2->Id == actinterface->nodes[q]->gnode->mfcpnode[2]->Id)
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
                  if(m==0 && actnode2_1->Id == actinterface->nodes[q]->gnode->mfcpnode[2]->Id) 
                  { 
                    a=q;
                    break;
                  }
                  if(m==1 && actnode2_2->Id == actinterface->nodes[q]->gnode->mfcpnode[2]->Id) 
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
  if(fsidyn->step == 1)
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
      actnode1_1 = actinterface->nodem[i];
      flag=0;
      if(actinterface->gnode_bound1mc != NULL)      
        if (actnode1_1->gnode->node->Id == actinterface->gnode_bound1mc->node->Id)
        flag++;

      if(actinterface->gnode_bound2mc != NULL)      
        if (actnode1_1->gnode->node->Id == actinterface->gnode_bound2mc->node->Id)
        flag++;
    
      if(flag>0)
      {
        if(actnode1_1->mtr_coeff.fdim < 15) 
        /* allocate memory to the node if the pointer mtr_coeff points to zero */
        {
          amdef("mrt coeff.",&(actnode1_1->mtr_coeff),actinterface->numnps,int_faces->numint,"DA");
          amzero(&(actnode1_1->mtr_coeff));
        }
      }
      else
      {
        amdef("mrt coeff.",&(actnode1_1->mtr_coeff),actinterface->numnps,1,"DV");
        amzero(&(actnode1_1->mtr_coeff));
      }
    }
  } /* end of if clause (fsidyn->step)*/
  
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
    a=0;
  } /* end of loop over the nodes of this element */
  
} /* end of loop over the interfaces (r) */
#ifdef DEBUG 
dstrc_exit();
#endif
} /*end of fsi_mortar_coeff */



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
void fsi_detect_intersection(DOUBLE lambdr1_2, DOUBLE lambdr2_1, 
     DOUBLE lambdr2_2, DOUBLE nr1_2, DOUBLE nr2_1, DOUBLE nr2_2, 
     DOUBLE *b1, DOUBLE *b2, INT *intersection)
{

#ifdef DEBUG 
dstrc_enter("fsi_detect_intersection");
#endif

if(lambdr2_1 < lambdr2_2)
{
                                      /* Intersection Scenarios 1-d     */
                                      /*      1          2              */
  if ((lambdr2_1 < EPS8)&&            /*  1   o                         */
      (lambdr2_1 > -EPS8))            /*  2   o                         */
  {                                   
    *b1 = 0.0;
    if (lambdr1_2<lambdr2_2)          
    {                                 /*  1   o--------o                */
      *b2 = nr1_2;                    /*  2   o------------o            */
      *intersection = 1;
      goto end_of_if;
    }
    if ((lambdr1_2 < lambdr2_2+EPS6 && lambdr1_2 > lambdr2_2-EPS6) && 
        (nr1_2 < nr2_2+EPS6 && nr1_2 > nr2_2-EPS6 ))       
    {                                 /*  1   o------------o            */
      *b2 = nr2_2;                    /*  2   o------------o            */
      *intersection = 1;
      goto end_of_if;
    }
    if (lambdr1_2 > lambdr2_2)           
    {                                 /*  1   o------------o            */
      *b2 = nr2_2;                    /*  2   o--------o                */
      *intersection = 1;
      goto end_of_if;
    }
  }
  if (lambdr2_1 < -EPS8)              /*  1      o                      */
  {                                   /*  2   o            o            */
    *b1 = 0.0;
    if (lambdr1_2 < lambdr2_2)            
    {                                 /*  1      o------o               */
      *b2 = nr1_2;                    /*  2   o------------o            */
      *intersection = 1;
      goto end_of_if;
    }
    if ((lambdr1_2 < lambdr2_2+EPS6 && lambdr1_2 > lambdr2_2-EPS6) && 
        (nr1_2 < nr2_2+EPS6 && nr1_2 > nr2_2-EPS6 ))       
    {                                 /*  1      o---------o            */
      *b2 = nr1_2;                    /*  2   o------------o            */
      *intersection = 1;
      goto end_of_if;
    }
    if (lambdr1_2 > lambdr2_2)           
    {                                 /*  1      o------------o         */
      *b2 = nr2_2;                    /*  2   o--------o                */
      *intersection = 1;
      goto end_of_if;
    }
  }
  if (lambdr2_1 > EPS8)               /*  1   o            o            */
  {                                   /*  2      o                      */
    *b1 = nr2_1; 
    if (lambdr1_2 < lambdr2_2)           
    {                                 /*  1   o------------o            */
      *b2 = nr1_2;                    /*  2      o------------o         */
      *intersection = 1;
      goto end_of_if;
    }
    if ((lambdr1_2 < lambdr2_2+EPS6 && lambdr1_2 > lambdr2_2-EPS6) && 
        (nr1_2 < nr2_2+EPS6 && nr1_2 > nr2_2-EPS6 ))       
    {                                 /*  1   o------------o            */
      *b2 = nr1_2;                    /*  2      o---------o            */
      *intersection = 1;
      goto end_of_if;
    }
    if (lambdr1_2 > lambdr2_2)         
    {                                 /*  1   o------------o            */
      *b2 = nr2_2;                    /*  2      o------o               */
      *intersection = 1;
      goto end_of_if;
    }
  }
} /* end of if (lambdr2_1 < lambdr2_2)*/
if(lambdr2_1 > lambdr2_2)
{
  dserror("Hey, lambdr2_1 > lambdr2_2, check your algorithm! ");
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

end_of_if:
;
#ifdef DEBUG 
dstrc_exit();
#endif
} /*end of fsi_detect_intersection*/



/*!---------------------------------------------------------------------                                         
\brief compute displacements for ale nodes, mortar method

<pre>                                                  firl / genk 02/04

The consistency condition is applied to evaluate the dirichlet b.c. for
the ale nodes. Here is assumed that the displacements of the structure 
nodes are stored in sol_mf.a.da[0][..]. After the computation of the 
corresponding displacements for the ale nodes these values are  
stored in sol_increment[0][..] and sol[0].  
The algorithm is nearly similar to the subroutine fsi_mortar_coeff.

</pre>
\param    *fsidyn        FSI_DYNAMIC   (i)   pointer to fsi_dynamic 
\param    *int_faces     INTERFACES    (i)   pointer to the interfaces of the model
\return 0                                                                             

------------------------------------------------------------------------*/
void fsi_calc_disp4ale(FSI_DYNAMIC *fsidyn, INTERFACES *int_faces)
{
DOUBLE gauss2[2][2];                /* gauss points and weights (2)*/
INT intersection;                   /* flag for intersection of nodal*/
                                    /* basis fct's 0=no, 1=yes*/
INT a,b,i,j,k,l,m,n,p,q,r;          /* counters */            
NODE *actnode1_1, *actnode1_2, *actnode2_1, *actnode2_2;
                                    /* pointer to actual nodes */
DOUBLE d1_1[2], d1_2[2], d2_1[2], d2_2[2];
                                    /* position vectors of nodes on the  */
                                    /* interface */
DOUBLE nrmdr1_2[2];                 /* normed vector r1_2 */
DOUBLE nr1_2, nr2_1, nr2_2;         /* norm of vectors r1_2, r2_1, r2_2  */
DOUBLE r1_2[2], r2_1[2], r2_2[2];   /* position vectors relative to node */
                                    /* d1_1 */
DOUBLE lambdr1_2, lambdr2_1, lambdr2_2; /* e.g lambdr1_2 =               */
                                        /*         r1_2[0] / nrmdr1_2[0] */
DOUBLE lambdr1_2_xy[2];              /* lambdr1_2[0] = r1_2[0]/nrmdr1_2[0]*/
                                     /* lambdr1_2[1] = r1_2[1]/nrmdr1_2[1]*/
DOUBLE lambdr2_1_xy[2];              /* lambdr2_1[0] = r2_1[0]/nrmdr2_1[0]*/
                                     /* lambdr2_1[1] = r2_1[1]/nrmdr2_1[1]*/
DOUBLE lambdr2_2_xy[2];              /* lambdr2_2[0] = r2_2[0]/nrmdr2_2[0]*/
                                     /* lambdr2_2[1] = r2_2[1]/nrmdr2_2[1]*/
DOUBLE b1, b2;                       /* integration bounds */
DOUBLE m1, n1, m2, n2;               /* parameters for linear functions */

ARRAY rhs;                           /* right hand side of equation       */
DOUBLE *rhs_a;                       /* to compute the displacements      */

ARRAY rhs_d;                         /* additional rhs due to D.b.c.      */
DOUBLE *rhs_d_a;                     /*                                   */

ARRAY slv_nods_onint_work;           /* additional rhs due to D.b.c.      */
DOUBLE *slv_nods_onint_work_a;       /*                                   */

char uplo[1];                        /* character for lapack solver */
INT info;                            /* flag for lapack solver */
INT ione;                            /* number of right hand sides*/

ARRAY solu_arr;                      /* nodal displacement array */
DOUBLE *solu_arr_a;

ARRAY hasdline;                      /* entry == 1 if node is on a dline */
DOUBLE *hasdline_a;

ARRAY work_vec;                      /* working vector to rearrange the lhs*/
DOUBLE *work_vec_a;
INTERFACE *actinterface;             /* the actual interface */

ARRAY rel_d_2_1;                     /* the relaxed displacements of the */
DOUBLE *rel_d_2_1_a;                 /* master node 2_1 */ 
ARRAY rel_d_2_2;                     /* the relaxed displacements of the */
DOUBLE *rel_d_2_2_a;                 /* master node 2_2 */ 
                                     
#ifdef DEBUG 
dstrc_enter("fsi_calc_disp4ale");
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
  ione = actinterface->nodem[0]->numdf;
  
  
  /*loop elements of slavefield, biorthogonal shape functions, test space */
  for(i=0;i<actinterface->numnps-1;i++)
  {
    actnode1_1 = actinterface->nodes[i]->gnode->mfcpnode[2];
    dsassert(actnode1_1!=NULL,"cannot read from NULL-pointer!\n");
    actnode1_2 = actinterface->nodes[i+1]->gnode->mfcpnode[2];
    dsassert(actnode1_2!=NULL,"cannot read from NULL-pointer!\n");

    /* ----------------------------------------loop nodes of this element */
    for(j=0;j<2;j++)
    {
      /* the deformed geometry is used to set up the functions*/
      /* for the slave nodes the actual displacements are used */
      /* d(n+1) are stored in sol_mf.a.da[1][.] */
      d1_1[0] = actnode1_1->x[0]+actnode1_1->sol_mf.a.da[1][0];
      d1_1[1] = actnode1_1->x[1]+actnode1_1->sol_mf.a.da[1][1];
      d1_2[0] = actnode1_2->x[0]+actnode1_2->sol_mf.a.da[1][0];
      d1_2[1] = actnode1_2->x[1]+actnode1_2->sol_mf.a.da[1][1];
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
        if(j==0 && actnode1_1->gnode->mfcpnode[1]->Id == actinterface->nodes[q]->Id) 
        { 
          a=q;
          break;
        }
        if(j==1 && actnode1_2->gnode->mfcpnode[1]->Id == actinterface->nodes[q]->Id) 
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
        /* d(n+1) are stored in sol_mf.a.da[0][.] */
        rel_d_2_1.a.dv[0] = actnode2_1->sol_mf.a.da[0][0];
        rel_d_2_1.a.dv[1] = actnode2_1->sol_mf.a.da[0][1];
        rel_d_2_2.a.dv[0] = actnode2_2->sol_mf.a.da[0][0];
        rel_d_2_2.a.dv[1] = actnode2_2->sol_mf.a.da[0][1];
        /* -----------------------------------loop nodes of this element */
        for(m=0;m<2;m++)
        {
          intersection = 0;
          /* the deformed geometry is used to set up the functions*/
          /* for the master nodes the displacements of the last iteration */
          /* are used */
          /* d(n) are stored in sol_mf.a.da[1][.] */
          d2_1[0] = actnode2_1->x[0]+actnode2_1->sol_mf.a.da[1][0];
          d2_1[1] = actnode2_1->x[1]+actnode2_1->sol_mf.a.da[1][1];
          d2_2[0] = actnode2_2->x[0]+actnode2_2->sol_mf.a.da[1][0];
          d2_2[1] = actnode2_2->x[1]+actnode2_2->sol_mf.a.da[1][1]; 
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
          /* x - component */
          lambdr1_2_xy[0] = r1_2[0] / nrmdr1_2[0];
          lambdr2_1_xy[0] = r2_1[0] / nrmdr1_2[0];
          lambdr2_2_xy[0] = r2_2[0] / nrmdr1_2[0];
          /* y - component */
          lambdr1_2_xy[1] = r1_2[1] / nrmdr1_2[1];
          lambdr2_1_xy[1] = r2_1[1] / nrmdr1_2[1];
          lambdr2_2_xy[1] = r2_2[1] / nrmdr1_2[1];
            /* check for a possible intersection */
            if(lambdr2_1_xy[0]/lambdr1_2_xy[0] <= EPS8 && 
               lambdr2_2_xy[0]/lambdr1_2_xy[0] <= EPS8) 
                 goto end_of_node2;
            if(nr2_1 > nr1_2+EPS6 && lambdr2_1_xy[0] > lambdr1_2_xy[0] &&
               nr2_2 > nr1_2+EPS6 && lambdr2_2_xy[0] > lambdr1_2_xy[0])
               goto end_of_node2;
          /* look for the best lambda for fsi_detect_intersection(),
             the goal of this choice is to minimize the error due to 
             floating point division */
          if((r1_2[0] > sqrt(r1_2[1]*r1_2[1])) || (r1_2[0]<-sqrt(r1_2[1]*r1_2[1])))
          {
            lambdr1_2 = lambdr1_2_xy[0];
            lambdr2_1 = lambdr2_1_xy[0];
            lambdr2_2 = lambdr2_2_xy[0];
          }
          else
          {
            lambdr1_2 = lambdr1_2_xy[1];
            lambdr2_1 = lambdr2_1_xy[1];
            lambdr2_2 = lambdr2_2_xy[1];
          }
          /* Check for a possible nonzero intersection, goto end if such an*/
          /* ---------------------------------intersection is not possible */        
          if((lambdr2_1 >= nr1_2) && (lambdr2_2 >= nr1_2)) goto end_of_node2;
          if((lambdr2_1 < 0) && (lambdr2_2 < 0)) goto end_of_node2;
          
          /* ---------------------here we determine the integration bounds */
          fsi_detect_intersection(lambdr1_2, lambdr2_1, lambdr2_2, nr1_2, 
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
      actnode2_1 = actnode1_1->gnode->mfcpnode[0]; /* struct node */
      rel_d_2_1.a.dv[0] = actnode2_1->sol_mf.a.da[0][0];
      rel_d_2_1.a.dv[1] = actnode2_1->sol_mf.a.da[0][1];

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
    slv_nods_onint_work.a.iv[i] = actinterface->nodes[i]->gnode->
                                  mfcpnode[2]->Id; /* ALE node */
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
    actnode1_1 = actinterface->nodes[i]->gnode->mfcpnode[2]; /* ALE node */
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
    actnode1_1->sol_mf.a.da[3][0] = solu_arr.a.da[a][0];
    actnode1_1->sol_mf.a.da[3][1] = solu_arr.a.da[a][1];

    /* copy the displacements to field sol_mf[6], in the */ 
    /* subroutine fsi_mortar_coeff this values are needed */
  }


#endif

  /*---------------------------------------------------------------------*/
  /* COMPUTATION OF NODAL DISPLACEMENTS WITH WEAK CONTINUITY AT THE      */
  /*                     BOUNDARY OF THE INTERFACE                       */
  /*---------------------------------------------------------------------*/
  /*---------------------------------------------------------------------*/
  /*                     SOLUTION USING LAPACK LU-DECOMPOSITION          */
  /*---------------------------------------------------------------------*/
  
  /* here the second part of the LAPACK solver is used. In actinterface->
     continuity_eq->A the factorized lhs matrix is stored. The 
     factorization is computed in the subroutine fsi_mortar_coeff */
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
    /* storage scheme in gnode->mfcpnode: (if the mesh is conforming)*/
    /* mfcpnode[0] = structure node */
    /* mfcpnode[1] = fluid node */
    /* mfcpnode[2] = ale node */
    /* storage scheme in gnode->mfcpnode: (if the mesh is non-conforming)*/
    /* mfcpnode[0] = 0 */
    /* mfcpnode[1] = fluid node */
    /* mfcpnode[2] = ale node */
    actnode1_1 = actinterface->nodes[i]->gnode->mfcpnode[2];

    actnode1_1->sol_mf.a.da[3][0] = rhs.a.da[0][i] ;
    actnode1_1->sol_mf.a.da[3][1] = rhs.a.da[1][i];
  }


} /* end of loop over the interfaces (r) */

} /*end of fsi_calc_disp4ale */


/*!---------------------------------------------------------------------
\brief fsi loads, mortar method

<pre>                                                         firl 04/04

in this function the function f2_fsiload is called to calculate the 
internal forces of the fluid elements at the interface

</pre>
\param	*int_faces       INTERFACES   (i)    the interfaces of the model
\return void

------------------------------------------------------------------------*/
void fsi_calc_intforces(INTERFACES *int_faces)
{
ELEMENT *actfele;                   /* the actual fluid element */
INT i,j;                            /* counter */
INTERFACE *actinterface;            /* the actual interface */

/* loop over the interfaces */
for(i=0; i<int_faces->numint; i++)
{
  actinterface = &int_faces->interface[i];
  /*loop over all fluid elements of actinterface, called slave elements */
  for(j=0; j<actinterface->numeles; j++)
  {
    actfele = actinterface->elements[j];
    f2_fsiload(actfele);
  }
}
#ifdef DEBUG 
dstrc_exit();
#endif
} /* end of fsi_calc_intforces*/


/*!---------------------------------------------------------------------
\brief fsi loads, mortar method

<pre>                                                         firl 04/04

in this function the nodal forces on the element edges due to the fluid
are calculated. the nodal forces are stored at sol_mf[3][0..1]

</pre>
\param	*ele		 ELEMENT      (i)    actual element
\return void

------------------------------------------------------------------------*/
void f2_fsiload(ELEMENT  *ele)
{
INT          lr;                 /* integration directions          */
INT          i,j,jj;             /* some loopers                    */
INT          nir,nis;            /* number of GP's in r-s direction */
INT          iel;                /* number of element nodes         */
INT          nd;                 /* element DOF                     */

DOUBLE       e1;              /* GP-koordinates in r-s-system   */
DOUBLE       facr;  /* integration factor  GP-info    */

/*--------------------- variables needed for integration of line loads */
INT             foundline;   /* flag for lineload present or not       */
INT             ngline;      /* number of geometrylines to the element */
GNODE          *actgnode;    /* the actual gnode under consideration   */
GLINE          *gline[4];    /* geometrylines of the element           */   
FSI_COUPLE_CONDITION *linefsi[4]; /* short grep on line-neum. conditions    */
NODE           *actfnode;    /* actual fluid node                      */
NODE           *actanode;    /* actual ale node                      */
INT             line;        /* looper over lines                      */
INT             ngnode;      /* number of geometry-nodes on g-line     */
INT             ngr;         /* number of GP'e for line-integration    */
INT             iegnod[4];
DOUBLE          xgp[3];      /* Coordinates of GP'es                   */
DOUBLE          wgx[3];      /* weights at GP'es                       */
DOUBLE          vnx,vny;     /* comp, of normal vector at INT point    */ 

DOUBLE          forceline[2];/* lineload value in x and y direct.(inp) */
DOUBLE          sigmaint[3]; /* fluid stresses at integration point    */
DOUBLE          nsigma[3][4];/* nodal fluid stresses       */
DOUBLE          xyzl[2][4];  /* nodal coordinates            */
DOUBLE          funct[2];    /* shape function values */
DOUBLE          deriv[2][2]; /* derivative of shape functions */
ARRAY eload_a; 
DOUBLE **eload;              /* static element load vector */
#ifdef DEBUG 
dstrc_enter("f2_fsiload");
#endif

#ifdef D_FSI
/*----------------------------------------------------------------------*/
/* calculation phase        (init=0)                                    */
/*----------------------------------------------------------------------*/
eload   = amdef("eload"  ,&eload_a,2,4,"DA");
/*---------------------------------------------------- initialize eload */
amzero(&eload_a);
/*-------------------------------------------- init the gaussian points */
/*w1intg(ele,data,1);*/
nir = 2;           /* number of gauss points in r direction */
nis = 2;           /* number of gauss points in s direction */
iel = 4;           /* number of element nodes         */
nd  = 8;           /* element DOF (wall element)    */

/*----------------------------------------------------------------------*/
/*--------- integration of line loads on lines adjacent to this element */
/*------------------------------------ check for presence of line loads */
foundline=0;
/*------------------------------------- number of lines to this element */
ngline=ele->g.gsurf->ngline;
/*-------------- loop over lines, check for neumann conditions on lines */
for (i=0; i<ngline; i++)
{
   gline[i] = ele->g.gsurf->gline[i];
   linefsi[i] = gline[i]->fsicouple;
   if(linefsi[i]==NULL) continue;
   if (linefsi[i]->fieldtyp==fluid) foundline=1;
   else linefsi[i]=NULL ;
}
if (foundline==0) goto endline;
/*------------------------------- loop over lines (with fsi conditions) */
for (line=0; line<ngline; line++)
{
   if (linefsi[line]==NULL) continue;
   /*------------------------------------ check number of nodes on line */
   ngnode = gline[line]->ngnode;
   /*--------------------------------------------------- get edge nodes */
   w1_iedg(iegnod,ele,line,0);
   /*------------------------------------------------ get the fsi-loads */
   for (j=0;j<ngnode;j++)
   {
      /*--------------------------------- find corresponding fluid node */
      /*actsgnode = gline[line]->gnode[j];*/
      /* storage scheme in gnode->mfcpnode: */
      /* mfcpnode[0] = structure node (if present) */
      /* mfcpnode[1] = fluid node (if present) */
      /* mfcpnode[2] = ale node (if present) */
      actgnode = gline[line]->gnode[j];
      /*actsnode  = actsgnode->node;*/
      actanode  = actgnode->mfcpnode[2];
      actfnode  = actgnode->mfcpnode[1];
      /*------------- get the coordinates in the deformed configuration */
      xyzl[0][j] = actanode->x[0]+actanode->sol.a.da[0][0];
      xyzl[1][j] = actanode->x[1]+actanode->sol.a.da[0][1];      
      /*-- loop the 3 stresses and get values from sol_mf of fluid node */
      for (i=0;i<3;i++)
      {
         nsigma[i][j] = actfnode->sol_mf.a.da[1][i];
      }
   }
   /*-------------------- coordinates and weights of integration points */
   /*--------- original GP-coordinates and weights for area-integration */
   ngr = nir; 
   /* _MF_
   for (i=0; i<ngr; i++) { xgp[i] = data->xgrr[i];
   			   wgx[i] = data->wgtr[i]; } */
   xgp[0] = -0.5773502691896;
   xgp[1] = 0.5773502691896;
   wgx[0] = 1.0;
   wgx[1] = 1.0;
   /*----------------------------------- integration loop on actual line */

   for (lr=0; lr<ngr; lr++)/*-------------------------- loop r-direction */
   {
      /*============================= gaussian point and weight at it ===*/
      e1   = xgp[lr];
      facr = wgx[lr]; 
      /*------------------- get shape function values on the actual edge */
      /*w1_degfuncderiv(funct,deriv,e1,ele->distyp,1);*/
      funct[0] = 1.0/2.0 * (1-e1);
      funct[1] = 1.0/2.0 * (1+e1);
      deriv[0][0]= -1.0/2.0;
      deriv[0][1]=  1.0/2.0;
      /*------------------------------------ compute load at gauss point */
      for (i=0;i<3;i++) sigmaint[i]=ZERO;
      for (i=0;i<3;i++)
      {
   	 for (j=0;j<ngnode;j++) sigmaint[i]+=funct[j]*nsigma[i][j];
      } 	 
      /*-------------------------- compute normal vector at gauss point *
       | see Mok et al in Engineering Computations (1999)               |
       *----------------------------------------------------------------*/
      vnx=ZERO;
      vny=ZERO;
      for(i=0;i<ngnode;i++) vnx+=deriv[0][i]*xyzl[1][i];
      for(i=0;i<ngnode;i++) vny-=deriv[0][i]*xyzl[0][i];
      /*-------------------------- compute stress vector at gauss point *
       |  force = sigma * n  (Cauchy's law)                             |
       *----------------------------------------------------------------*/
      forceline[0] = sigmaint[0]*vnx + sigmaint[2]*vny;
      forceline[1] = sigmaint[2]*vnx + sigmaint[1]*vny;
      /*-------------- add load vector component to element load vector */
      /* jacobian determinant cancels with length of normal vector      */
      for (j=0; j<ngnode; j++)
      {
         jj=iegnod[j];
	 eload[0][jj] += funct[j] * forceline[0] * facr;
         eload[1][jj] += funct[j] * forceline[1] * facr;
      }
   }/*============================================= end of loop over lr */
   /* line number lie has been done,switch of the neumann pointer of it */
   /*ele->g.gsurf->gline[line]->neum=NULL;*/
}/* end loop line over lines */
/*----------------------------------------------------------------------*/
endline:
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/* store the solutions of eload at the fluid nodes                      */
/*----------------------------------------------------------------------*/
/* loop over the nodes of ele */
for(i=0; i<ele->numnp; i++)
{ 
  actfnode = ele->node[i];
  /* loop over the dofs of actfnode */ 
  for(j=0; j<2; j++)
  {
    actfnode->sol_mf.a.da[3][j] += eload[j][i];
  }
}

#else
dserror("FSI-functions not compiled in\n");
#endif
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of f2_fsiload */



/*!---------------------------------------------------------------------                                         
\brief put coupling forces from fluid nodes to structure nodes, mortar method

<pre>                                                        firl 04/04

The internal forces of the fluid nodes  are stored at 
sol_mf.a.da[3][..]. This subroutine assembles a array 
coupforc.a.da [actinterface->numnps][2]. 
In a second step the corresponding forces on the structure nodes are 
computed. This forces depend on the mortar coefficients stored at the 
structure nodes on the interface in node->mtr_coeff.a.dv[]. The coupling 
forces on the structure nodes are stored at sol_mf.a.da[4][..].

</pre>
\param *masterfield   FIELD	     (i)   struct field
\param *int_faces    INTERFACES      (i)   the interfaces of the model
\return 0                                                                             

------------------------------------------------------------------------*/
void fsi_put_coupforc2struct(FIELD *masterfield, INTERFACES *int_faces)
{
ARRAY coupforc;                    /* array of coupling forces          */
DOUBLE *coupforc_a;
NODE *actnode;                     /* pointer to the actual node        */
INT a, i, j, k;                    /* counters                          */
DOUBLE b;                          /* parameter for nodal position      */
                                   /* on the interface                  */
INTERFACE *actinterface;           /* the actual interface              */

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
    
    /*------ the internal forces of sol_mf.a.da[3][.] are scaled by -1 to */
    /*-------------------------------------------- obtain external forces */
    coupforc.a.da[j][0] = -actnode->sol_mf.a.da[3][0] * b;
    coupforc.a.da[j][1] = -actnode->sol_mf.a.da[3][1] * b;
  }
  
  /*----------------------------------------------------------------------*/
  /*                               PART II                                */
  /*----------------------------------------------------------------------*/
  /* ---put forces from coupforc[..][0..1] to the respective master nodes */
  
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
  } /* end of loop over the master nodes at the interface (i) */
} /* end of loop over the interfaces (i)*/
} /* end of fsi_put_coupforc2struct */


/*!----------------------------------------------------------------------
\brief  point neumann conditions                            firl  04/04     |

 <pre> 

 this routine is based on rhs_point_neum written by m.gee in 3/02     
 Attention! This assembly of nodal forces works only correct with    
 partitioning in the "Cut_Elements" style !!!!                       
                                                                      
 Here we set up the rhs[4] for the structure field                    
 The forces assembled in this routine are resulting from fsi coupling 

</pre>
\param *rhs          DOUBLE	     (i)   pointer to rhs
\param dimrhs        INT             (i)   the dimension of the rhs
\param *actpart      PARTITION       (i)   the actual part
\return 0                                                 
                            
------------------------------------------------------------------------*/
void fsiserv_rhs_point_neum(DOUBLE *rhs, INT dimrhs, PARTITION *actpart)     
{
INT             i,j;
NODE           *actnode;
INT             actdof;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("fsiserv_rhs_point_neum");
#endif
/*--------------------- loop all nodes of actpart (all structure nodes) */
for (i=0; i<actpart->pdis[0].numnp; i++)
{
  /*--------------------------------------------------- set active node */
  actnode = actpart->pdis[0].node[i];
  /*------------- check if the node is a coupling node (structure node) */
  if (actnode->gnode->fsicouple == NULL) continue;
  if (actnode->gnode->fsicouple->fieldtyp == structure)
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
} /* end of fsiserv_rhs_point_neum */

/* end of file "fsi_mortar.c" */
#endif /* end ifdef D_FSI */
#endif /* end ifdef D_MORTAR */
/*! @} (documentation module close)*/
