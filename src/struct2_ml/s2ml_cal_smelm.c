/*!----------------------------------------------------------------------
\file
\brief contains the routine 

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>
*----------------------------------------------------------------------*/
#ifdef D_MLSTRUCT

#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../wall1/wall1.h"
#include "s2ml.h"
#include "s2ml_prototypes.h"

/*! 
\addtogroup MLSTRUCT 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief  routine for calling the submesh elements to calculate the
stiffnesses and internal forces; these are then assembled

\param  *actmaele        ELEMENT     (I)    actual macro element
\param  *actsmfield      FIELD       (I)    submesh field 
\param  *actsmsolv       SOLVAR      (I/O)  submesh SOLVAR->sm assembled stif mi_mi
\param  *actsmintra      INTRA       (I)    the sm pseudo intra-communicator
\param  *smintforcemi    DOUBLE      (O)    submesh SOLVAR->sm assembled stif mi_mi
\param  *smstiffmima_csr DBCSR       (O)    sm assebled stiffness_mi_ma
\param  *smstiffmami_csr DBCSR       (O)    sm assebled stiffness_ma_mi
\param  *smintforcema    DOUBLE      (O)    sm "assembled" int force macro
\param **smstiffmama     DOUBLE      (O)    sm "assembled" stiffness macro-macro
\param   istore          INT         (I)    is it update step?

\return void                                               

*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *sm_mat;
/*!----------------------------------------------------------------------
proc's info about his partition, global variable defined in input_submesh.c
*----------------------------------------------------------------------*/
extern struct _PARTITION  *sm_part;
/*----------------------------------------------------------------------*
 | global dense matrices for submesh-element routines                   |
 *----------------------------------------------------------------------*/
struct _ARRAY estif_ma_ma;    /* element stiffness matrix (macro-macro) */
struct _ARRAY estif_ma_mi;    /* element stiffness matrix (macro-micro) */
struct _ARRAY estif_mi_ma;    /* element stiffness matrix (micro-macro) */
struct _ARRAY estif_mi_mi;    /* element stiffness matrix (micro-micro) */
struct _ARRAY intforce_ma;    /* element internal force   (micro)       */  
struct _ARRAY intforce_mi;    /* element internal force   (macro)       */  
/*----------------------------------------------------------------------*/
void s2ml_cal_smelm(ELEMENT     *actmaele,   /*  actual macro element    */
                    FIELD       *actsmfield, /*  submesh field           */
                    SOLVAR      *actsmsolv,  /*  submesh SOLVAR->sm assembled stif mi_mi*/
                    INTRA       *actsmintra, /* the sm pseudo intra-communicator */
                    DOUBLE      *smintforcemi, /*  sm assembled int force micro */
                    DBCSR       *smstiffmima_csr, /* sm assebled stiffness_mi_ma*/
                    DBCSR       *smstiffmami_csr, /* sm assebled stiffness_ma_mi*/
                    DOUBLE      *smintforcema, /*  sm "assembled" int force macro */
                    DOUBLE     **smstiffmama, /*  sm "assembled" stiffness macro-macro */
                    INT          istore)     /*  is it update step?      */
{
INT              i,j,k,a,b; /* loopers */
ELEMENT         *actsmele;
MATERIAL        *actsmmat;
PARTITION       *actsmpart;       /* my partition of the active submeshfield */
ASSEMBLE_ACTION  assemble_action;
INT              actsysarray = 0; 

INT              nnz_guess;       /* number of nonzero guess for opening CSR-matrix */
INT              numeq_sm;        /* number of submesh DOF's (without the dirichlet ones) */
INT              firstdof=0;      /* first dofnumber on this processor    */
INT              lastdof;         /* last dofnumber on this proc          */
INT              lm[18];          /* location matrix of submesh          */
DOUBLE           val;             /* value i.e. elemententity to assemble          */
INT              rowindex;        /* row number in global matrix, to place elemententity   */  
INT              columindex;      /* column index in global matrix, to place elemententity*/
INT              counter;         /* */
NODE            *actsmnode;       /* actual node                            */
INT              numdof_maele;    /* */
INT              numdof_smele;    /* */
INT              dof;             /* global number of DOF */
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s2ml_cal_smelm");
#endif
/*----------------------------------------------------------------------*/

actsmpart  = &(sm_part[0]);
/*----------------------------------------------------------------------*/
/*      allocate assembled stifness micro-macro (not quadratic!)        */
/*            csr-format: -> doesn't need a mask beforehand             */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
numeq_sm = actsmfield->dis[0].numeq;
lastdof  = numeq_sm - 1;
/*------------------------------------------ number of element DOF's ---*/
numdof_maele = actmaele->numnp * 2;
numdof_smele = actsmfield->dis[0].element[0].numnp * 2;
/*---------------- guess number of nonzero entities of this matrices ---*/
nnz_guess = 4 * numeq_sm;
/*-------------------------- open assebled stiffness_mi_ma and ma_mi ---*/
mlpcg_csr_open(smstiffmima_csr,firstdof,lastdof,
               smstiffmima_csr->numeq_total,2*nnz_guess,actsmintra);
mlpcg_csr_zero(smstiffmima_csr,actsmintra);
/*------------------------------------------------------------------ ---*/
mlpcg_csr_open(smstiffmami_csr,firstdof,(numdof_maele-1),
               numdof_maele,2*nnz_guess,actsmintra);
mlpcg_csr_zero(smstiffmami_csr,actsmintra);
/* ========================================= loop over submesh elements */
for (i=0; i<actsmfield->dis[0].numele; i++)
{
   /*---------------------------- set pointer to active submesh-element */
   actsmele = &(actsmfield->dis[0].element[i]);
   actsmmat = &(sm_mat[actsmele->mat - 1]);
   /*----------------------------------- call submesh element routines  */
   switch(actsmele->eltyp)
   {
   case el_wall1:
      s2ml_stiff_wall(actsmmat,actmaele,actsmele,
                      &estif_ma_ma,&estif_ma_mi,&estif_mi_ma,&estif_mi_mi,
                      &intforce_ma,&intforce_mi,istore,0); 
   break;
   case el_interf:
      s2ml_stiff_interf(actsmmat,actmaele,actsmele,
                       &estif_ma_ma,&estif_ma_mi,&estif_mi_ma,&estif_mi_mi,
                       &intforce_ma,&intforce_mi,istore,0);
   break;
   default:
     dserror("element unknown for submesh");
   break;
   }/* end of switch elements */
/*----------------------------------------------------------------------*/
/*   assemble stiffness micro-micro in ccf-format ( <- for umpfpack)    */
/*       quadr. matrix -> standard assemley, with submesh-mask          */
/*----------------------------------------------------------------------*/
   assemble_action = assemble_one_matrix;
   assemble(actsysarray,     /*  INT==0 */
            &estif_mi_mi,    /*  element stiffness micro-micro form actual smele */
            -1,              /*  no second matrix to assemble with this ccf-mask */
            NULL,            /*  no second matrix to assemble with this ccf-mask */
            actsmpart,
            actsmsolv,
            actsmintra,
            actsmele,
            assemble_action, /* =assemble-one-matrix */
            NULL);           /* container is never used in routine "assemble" */
/*----------------------------------------------------------------------*/
/*   assemble stifness micro-macro and macro-micro (not quadratic!)     */
/*            csr-format: -> doesn't need a mask beforehand             */
/*----------------------------------------------------------------------*/

/*------------------------------------- location matrix for sm-DOF's ---*/
   counter = 0;
   for (j=0; j<actsmele->numnp; j++)
   {
      actsmnode = actsmele->node[j];
      for (k=0;k<actsmnode->numdf;k++)
      {
         lm[counter] = actsmnode->dof[k];
         counter ++;
      }
   }
/*----------------------------------------------------------------------*/
/*   assemble stiffness micro-macro: allembley only rows,               */
/*                   colums are just summed up                          */
/*----------------------------------------------------------------------*/
   /*-------------------------- loop over the sm-element DOF's row a ---*/
   for (a=0;a<numdof_smele;a++)
   {
      rowindex = lm[a];
      /*------------------ Dirichlet-bounded SM-DOF: do not assemble ---*/
      if (rowindex >= numeq_sm) continue;
      /*----------------- loop over the macro-element DOF's column b ---*/
      for (b=0;b<numdof_maele;b++)
      {
      /*- Dir-bound Macro-DOF:(not global assemb)but assemble here in SM*/
         columindex = b;
         val = estif_mi_ma.a.da[a][b];
         mlpcg_csr_addentry(smstiffmima_csr, val, rowindex, columindex, actsmintra);   
      }
   }
/*----------------------------------------------------------------------*/
/*   assemble stiffness macro-micro: allembley only colums,             */
/*                   rows are just summed up                            */
/*----------------------------------------------------------------------*/
   /*----------------------- loop over the macor-element DOF's row a ---*/
   for (a=0;a<numdof_maele;a++)
   {
      rowindex = a;
      /*-------------------- loop over the sm-element DOF's column b ---*/
      for (b=0;b<numdof_smele;b++)
      {
         columindex = lm[b];
         if (columindex >= numeq_sm) continue;
         val = estif_ma_mi.a.da[a][b];
         mlpcg_csr_addentry(smstiffmami_csr, val, rowindex, columindex, actsmintra);   
      }
   }
/*----------------------------------------------------------------------*/
/*                    assemble internal forces micro                    */
/*----------------------------------------------------------------------*/
   counter = -1;
   for (j=0; j<actsmele->numnp; j++)
   {
      for (k=0; k<actsmele->node[j]->numdf; k++)
      {
         counter++;
         dof = actsmele->node[j]->dof[k];
         if (dof >= numeq_sm) continue;
         smintforcemi[dof] += intforce_mi.a.dv[counter];
      }
   }
   
/*----------------------------------------------------------------------*/
/*          add stiffness macro-macro and internal force macro          */
/*----------------------------------------------------------------------*/
   for (a=0;a<numdof_maele;a++)
   {
     smintforcema[a] += intforce_ma.a.dv[a];
     for (b=0;b<numdof_maele;b++)
     {
       smstiffmama[a][b] += estif_ma_ma.a.da[a][b];
     }
   }
   
}/* end of loop over elements */
   
/*------------------------- close assebled stiffness_mi_ma and ma_mi ---*/
mlpcg_csr_close(smstiffmima_csr);
mlpcg_csr_close(smstiffmami_csr);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of w1ml_cal_smelm */
/*----------------------------------------------------------------------*/



/*!----------------------------------------------------------------------
\brief  routine for initialisation of structural multiscale if it
is the first time for an macroelement to do multiscale

<pre>                                                           ah 07/04 
This routine 

</pre>
\param  *actmaele        ELEMENT     (I)    actual macro element
\param   numsmnodes      INT         (I)    number of submeshnodes
\param   numsmele        INT         (I)    number of submeshelements

\return void                                               

*----------------------------------------------------------------------*/
void s2ml_init(ELEMENT  *actmaele,       /*  actual macro element       */
               INT       numsmnodes,     /*  number of submeshnodes     */
               INT       numsmele)       /*  number of submeshelements  */
{
INT    i,j;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s2ml_init");
#endif
/*----------------------------------------------------------------------*/
if (estif_ma_ma.Typ != cca_DA)/*-> nur einmalige allocierung von Speicher, nicht fuer jedes Makroelement */
{
  amdef("k_ma_ma",&estif_ma_ma,(MAXNOD*3),(MAXNOD*3),"DA");
  amdef("k_ma_mi",&estif_ma_mi,(MAXNOD*3),(MAXNOD*3),"DA");
  amdef("k_mi_ma",&estif_mi_ma,(MAXNOD*3),(MAXNOD*3),"DA");
  amdef("k_mi_mi",&estif_mi_mi,(MAXNOD*3),(MAXNOD*3),"DA");
  amdef("intf_ma",&intforce_ma,(MAXNOD*3),1,"DV");
  amdef("intf_mi",&intforce_mi,(MAXNOD*3),1,"DV");
}
/*----------------- allocate and initialise macro info at submesh nodes */      
actmaele->e.w1->sm_nodaldata = (SM_NODAL_DATA*)CCACALLOC(numsmnodes,sizeof(SM_NODAL_DATA));
for (i=0; i<numsmnodes; i++)
{
   actmaele->e.w1->sm_nodaldata[i].displ_mi[0] = 0.0;
   actmaele->e.w1->sm_nodaldata[i].displ_mi[1] = 0.0;
   actmaele->e.w1->sm_nodaldata[i].incre_displ_mi[0] = 0.0;
   actmaele->e.w1->sm_nodaldata[i].incre_displ_mi[1] = 0.0;
   actmaele->e.w1->sm_nodaldata[i].store_displ_mi[0] = 0.0;
   actmaele->e.w1->sm_nodaldata[i].store_displ_mi[1] = 0.0;
}
/*----------------------------------------------------------------------*/
/*---- allocate and initialise history variables at submesh elementGP's */      
actmaele->e.w1->sm_eledata = (SM_ELEMENT_DATA*)CCACALLOC(numsmele,sizeof(SM_ELEMENT_DATA));
for(i=0; i<numsmele; i++)
{
  actmaele->e.w1->sm_eledata[i].sm_GPdata = (SM_GP_DATA*)CCACALLOC(9,sizeof(SM_GP_DATA));
  for(j=0; j<9; j++)
  {
    actmaele->e.w1->sm_eledata[i].sm_GPdata[j].kappa    = 0.0;
    actmaele->e.w1->sm_eledata[i].sm_GPdata[j].dt       = 0.0;
    actmaele->e.w1->sm_eledata[i].sm_GPdata[j].dn       = 0.0;
    actmaele->e.w1->sm_eledata[i].sm_GPdata[j].utpl     = 0.0;
    actmaele->e.w1->sm_eledata[i].sm_GPdata[j].D[0][0]  = 0.0;
    actmaele->e.w1->sm_eledata[i].sm_GPdata[j].D[0][1]  = 0.0;
    actmaele->e.w1->sm_eledata[i].sm_GPdata[j].D[1][0]  = 0.0;
    actmaele->e.w1->sm_eledata[i].sm_GPdata[j].D[1][1]  = 0.0;
    actmaele->e.w1->sm_eledata[i].sm_GPdata[j].T[0]     = 0.0;
    actmaele->e.w1->sm_eledata[i].sm_GPdata[j].T[1]     = 0.0;
    actmaele->e.w1->sm_eledata[i].sm_GPdata[j].yip      = -1;
    actmaele->e.w1->sm_eledata[i].sm_GPdata[j].kappa_n  = 0.0;
    actmaele->e.w1->sm_eledata[i].sm_GPdata[j].kappa_t  = 0.0;
    
  }
}
/*---- allocate stifness matrices which have to be stored at macroelement */      
actmaele->e.w1->stiff_mi_ma_csr = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
actmaele->e.w1->stiff_mi_mi_ccf = (SPARSE_ARRAY*)CCACALLOC(1,sizeof(SPARSE_ARRAY));
actmaele->e.w1->stiff_mi_mi_ccf[0].ccf = (CCF*)CCACALLOC(1,sizeof(CCF));
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of w1ml_init */
/*----------------------------------------------------------------------*/

/*! @} (documentation module close)*/

#endif /* D_MLSTRUCT */
