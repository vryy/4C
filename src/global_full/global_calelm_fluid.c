/*!----------------------------------------------------------------------
\file
\brief initialise and control element calculations 

------------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "../shell8/shell8.h"
#include "../wall1/wall1.h"
#include "../brick1/brick1.h"
#include "../fluid3/fluid3.h"
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
enum _CALC_ACTION calc_action[MAXFIELD];
/*----------------------------------------------------------------------*
 | global dense matrices for element routines             genk 04/02    |
 *----------------------------------------------------------------------*/
struct _ARRAY estif_global;    /* element stiffness matrix              */
struct _ARRAY emass_global;    /* element mass matrix                   */  
struct _ARRAY etforce_global;  /* element Time RHS                      */
struct _ARRAY eiforce_global;  /* element Iteration RHS                 */
struct _ARRAY edforce_global;  /* element dirichlet RHS                 */
struct _ARRAY intforce_global;
/*!---------------------------------------------------------------------                                         
\brief call fluid elements

<pre>                                                         genk 04/02       
</pre>
\param *actfield      FIELD	  (i)	 active field		 
\param *actsolv       SOLVAR	  (i)    active SOLVAR  	 
\param *actpart       PARTITION   (i)	 my part. of this field  
\param *actintra      INTRA	  (i)	 my intra-communicator   
\param  sysarray1     int	  (i)	 num. of 1st sparse system matrix  
\param  sysarray2     int	  (i)	 number of secnd system matrix, if present, else -1  
\param *ftimerhs      double	  (o)	 time RHS		 
\param *fiterhs       double	  (o)	 iteration RHS  	 
\param  global_numeq  int	  (i)	 total number of equations on all procs 
\param  nii	      int	  (i)	 flag for iteration RHS  
\param  nif	      int	  (i)	 flag for time RHS	 
\param  kstep	      int	  (-)
\param *action	      CALC_ACTION (i)	 calculation option passed to element routines  
\return void                                                                       

------------------------------------------------------------------------*/
void calelm_fluid(FIELD        *actfield,      
                  SOLVAR       *actsolv,      
                  PARTITION    *actpart,      
                  INTRA        *actintra,     
                  int           sysarray1,    
                  int           sysarray2,    
                  double       *ftimerhs,     
		  double       *fiterhs,      
                  int           global_numeq, 
                  int           nii,          
		  int           nif,          
                  int           kstep,
                  CALC_ACTION  *action)       
/*----------------------------------------------------------------------*/
{
int               i;
int               hasdirich=0;      /* flag                             */
ELEMENT          *actele;           /* actual element                   */
SPARSE_TYP        sysarray1_typ;     
SPARSE_TYP        sysarray2_typ;
ASSEMBLE_ACTION   assemble_action;

#ifdef DEBUG 
dstrc_enter("calelm_fluid");
#endif

/*-------------- zero the parallel coupling exchange buffers if present */  
#ifdef PARALLEL 
/*------------------------ check the send & recv buffers from sysarray1 */
if (sysarray1 != -1)
{
   switch(actsolv->sysarray_typ[sysarray1])
   {
   case msr:
      if (actsolv->sysarray[sysarray1].msr->couple_d_send)
         amzero(actsolv->sysarray[sysarray1].msr->couple_d_send);
      if (actsolv->sysarray[sysarray1].msr->couple_d_recv)
         amzero(actsolv->sysarray[sysarray1].msr->couple_d_recv);
   break;
   case parcsr:
      if (actsolv->sysarray[sysarray1].parcsr->couple_d_send)
         amzero(actsolv->sysarray[sysarray1].parcsr->couple_d_send);
      if (actsolv->sysarray[sysarray1].parcsr->couple_d_recv)
         amzero(actsolv->sysarray[sysarray1].parcsr->couple_d_recv);
   break;
   case ucchb:
      if (actsolv->sysarray[sysarray1].ucchb->couple_d_send)
         amzero(actsolv->sysarray[sysarray1].ucchb->couple_d_send);
      if (actsolv->sysarray[sysarray1].ucchb->couple_d_recv)
         amzero(actsolv->sysarray[sysarray1].ucchb->couple_d_recv);
   break;
   case dense:
      if (actsolv->sysarray[sysarray1].dense->couple_d_send)
         amzero(actsolv->sysarray[sysarray1].dense->couple_d_send);
      if (actsolv->sysarray[sysarray1].dense->couple_d_recv)
         amzero(actsolv->sysarray[sysarray1].dense->couple_d_recv);
   break;
   case rc_ptr:
      if (actsolv->sysarray[sysarray1].rc_ptr->couple_d_send)
         amzero(actsolv->sysarray[sysarray1].rc_ptr->couple_d_send);
      if (actsolv->sysarray[sysarray1].rc_ptr->couple_d_recv)
         amzero(actsolv->sysarray[sysarray1].rc_ptr->couple_d_recv);
   break;
   case ccf:
      if (actsolv->sysarray[sysarray1].ccf->couple_d_send)
         amzero(actsolv->sysarray[sysarray1].ccf->couple_d_send);
      if (actsolv->sysarray[sysarray1].ccf->couple_d_recv)
         amzero(actsolv->sysarray[sysarray1].ccf->couple_d_recv);
   break;
   case skymatrix:
      if (actsolv->sysarray[sysarray1].sky->couple_d_send)
         amzero(actsolv->sysarray[sysarray1].sky->couple_d_send);
      if (actsolv->sysarray[sysarray1].sky->couple_d_recv)
         amzero(actsolv->sysarray[sysarray1].sky->couple_d_recv);
   break;
   default:
      dserror("Unknown typ of system matrix");
   break;
   }
}
/*------------------------ check the send & recv buffers from sysarray2 */
if (sysarray2 != -1)
{
   switch(actsolv->sysarray_typ[sysarray2])
   {
   case msr:
      if (actsolv->sysarray[sysarray2].msr->couple_d_send)
         amzero(actsolv->sysarray[sysarray2].msr->couple_d_send);
      if (actsolv->sysarray[sysarray2].msr->couple_d_recv)
         amzero(actsolv->sysarray[sysarray2].msr->couple_d_send);
   break;
   case parcsr:
      if (actsolv->sysarray[sysarray2].parcsr->couple_d_send)
         amzero(actsolv->sysarray[sysarray2].parcsr->couple_d_send);
      if (actsolv->sysarray[sysarray2].parcsr->couple_d_recv)
         amzero(actsolv->sysarray[sysarray2].parcsr->couple_d_send);
   break;
   case ucchb:
      if (actsolv->sysarray[sysarray2].ucchb->couple_d_send)
         amzero(actsolv->sysarray[sysarray2].ucchb->couple_d_send);
      if (actsolv->sysarray[sysarray2].ucchb->couple_d_recv)
         amzero(actsolv->sysarray[sysarray2].ucchb->couple_d_send);
   break;
   case dense:
      if (actsolv->sysarray[sysarray2].dense->couple_d_send)
         amzero(actsolv->sysarray[sysarray2].dense->couple_d_send);
      if (actsolv->sysarray[sysarray2].dense->couple_d_recv)
         amzero(actsolv->sysarray[sysarray2].dense->couple_d_send);
   break;
   case rc_ptr:
      if (actsolv->sysarray[sysarray2].rc_ptr->couple_d_send)
         amzero(actsolv->sysarray[sysarray2].rc_ptr->couple_d_send);
      if (actsolv->sysarray[sysarray2].rc_ptr->couple_d_recv)
         amzero(actsolv->sysarray[sysarray2].rc_ptr->couple_d_recv);
   break;
   case ccf:
      if (actsolv->sysarray[sysarray2].ccf->couple_d_send)
         amzero(actsolv->sysarray[sysarray2].ccf->couple_d_send);
      if (actsolv->sysarray[sysarray2].ccf->couple_d_recv)
         amzero(actsolv->sysarray[sysarray2].ccf->couple_d_recv);
   break;
   case skymatrix:
      if (actsolv->sysarray[sysarray2].sky->couple_d_send)
         amzero(actsolv->sysarray[sysarray2].sky->couple_d_send);
      if (actsolv->sysarray[sysarray2].sky->couple_d_recv)
         amzero(actsolv->sysarray[sysarray2].sky->couple_d_recv);
   break;
   default:
      dserror("Unknown typ of system matrix");
   break;
   }
}
#endif

/* =======================================================call elements */
/*---------------------------------------------- loop over all elements */
for (i=0; i<actpart->pdis[0].numele; i++)
{
   /*------------------------------------ set pointer to active element */
   actele = actpart->pdis[0].element[i];
   /* if present, init the element vectors intforce_global and dirich_global */
/*   if (dvec) amzero(&intforce_global);*/
   switch(actele->eltyp)/*======================= call element routines */
   {
   case el_shell8:
      shell8(actfield,actpart,actintra,actele,
             &estif_global,&emass_global,&intforce_global,
             kstep,0,NULL,action);
   break;
   case el_brick1:
      brick1(actpart,actintra,actele,
             &estif_global,&emass_global,
             action);
   break;
   case el_wall1:
      wall1(actpart,actintra,actele,
            &estif_global,&emass_global,&intforce_global,0,NULL,
            action);
   break;
   case el_fluid2:
      fluid2(actpart,actintra,actele,
             &estif_global,&emass_global,
             &etforce_global,&eiforce_global,&edforce_global,
	     action,&hasdirich); 
   break;
   case el_fluid3:
      fluid3(actpart,actintra,actele,
             &estif_global,&emass_global,
             &etforce_global,&eiforce_global,&edforce_global,
	     action,&hasdirich);    
   break;
   case el_ale3:
   break;
   case el_none:
      dserror("Typ of element unknown");
   break;
   default:
      dserror("Typ of element unknown");
   }/* end of calling elements */


   switch(*action)/*=== call assembly dependent on calculation-flag */
   {
   case calc_struct_linstiff     : assemble_action = assemble_one_matrix; break;
   case calc_struct_nlnstiff     : assemble_action = assemble_one_matrix; break;
   case calc_struct_nlnstiffmass : assemble_action = assemble_two_matrix; break;
   case calc_struct_eleload      : assemble_action = assemble_do_nothing; break;
   case calc_struct_stress       : assemble_action = assemble_do_nothing; break;
   case calc_struct_update_istep : assemble_action = assemble_do_nothing; break;
   case calc_fluid               : assemble_action = assemble_one_matrix; break;
   default: dserror("Unknown type of assembly"); break;
   }
   /*------------------------------ assemble one or two system matrices */
   assemble(sysarray1,
            &estif_global,
            sysarray2,
            NULL,
            actpart,
            actsolv,
            actintra,
            actele,
            assemble_action);
   /*------------------- assemble the vector etforce_global to time-rhs */
   if (nif!=0)
      assemble_intforce(actele,ftimerhs,global_numeq,&etforce_global);
   /*-------------- assemble the vector eiforce_global to iteration rhs */
   if (nii!=0)
      assemble_intforce(actele,fiterhs,global_numeq,&eiforce_global); 
   /*-------------- assemble the vector edforce_global to iteration rhs */
   if (hasdirich!=0)
      assemble_intforce(actele,fiterhs,global_numeq,&edforce_global);
}/* end of loop over elements */
/*----------------------------------------------------------------------*/
/*                    in parallel coupled dofs have to be exchanged now */
/*             (if there are any inter-proc couplings, which is tested) */
/*----------------------------------------------------------------------*/
#ifdef PARALLEL 
switch(*action)
{
case calc_struct_linstiff      : assemble_action = assemble_one_exchange; break;
case calc_struct_nlnstiff      : assemble_action = assemble_one_exchange; break;
case calc_struct_nlnstiffmass  : assemble_action = assemble_two_exchange; break;
case calc_struct_eleload       : assemble_action = assemble_do_nothing; break;
case calc_struct_stress        : assemble_action = assemble_do_nothing; break;
case calc_struct_update_istep  : assemble_action = assemble_do_nothing; break;
case calc_fluid                : assemble_action = assemble_one_exchange; break;
default: dserror("Unknown type of assembly"); break;
}
/*------------------------------ exchange coupled dofs, if there are any */
assemble(sysarray1,
         NULL,
         sysarray2,
         NULL,
         actpart,
         actsolv,
         actintra,
         actele,
         assemble_action);
#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of calelm_fluid */

/*!---------------------------------------------------------------------                                         
\brief call fluid elements and initialise them

<pre>                                                         genk 04/02       
</pre>
\param *actfield      FIELD	  (i)	 active field		 
\param *actpart       PARTITION   (i)	 my part. of this field  
\param *action	      CALC_ACTION (i)	 calculation option passed to element routines  
\return void                                                                       

------------------------------------------------------------------------*/
void calinit_fluid(FIELD       *actfield,    
                   PARTITION   *actpart,    
                   CALC_ACTION *action)
{
int i;                        /* a counter */
int is_shell8=0;              /* flags to check for presents of certain element types */
int is_brick1=0;
int is_wall1 =0;
int is_fluid2=0;
int is_fluid3=0;
int is_ale3=0;
ELEMENT *actele;              /* active element */

#ifdef DEBUG 
dstrc_enter("calinit");
#endif

/*-------------------------- define dense element matrices for assembly */
amdef("estif",&estif_global,(MAXNOD*MAXDOFPERNODE),(MAXNOD*MAXDOFPERNODE),"DA");
amdef("emass",&emass_global,(MAXNOD*MAXDOFPERNODE),(MAXNOD*MAXDOFPERNODE),"DA");
amdef("etforce",&etforce_global,(MAXNOD*MAXDOFPERNODE),1,"DV");
amdef("eiforce",&eiforce_global,(MAXNOD*MAXDOFPERNODE),1,"DV");
amdef("edforce",&edforce_global,(MAXNOD*MAXDOFPERNODE),1,"DV");
amdef("inforce",&intforce_global,(MAXNOD*MAXDOFPERNODE),1,"DV");

/*--------------------what kind of elements are there in this example ? */
for (i=0; i<actfield->dis[0].numele; i++)
{
   actele = &(actfield->dis[0].element[i]);
   switch(actele->eltyp)
   {
   case el_shell8:
      is_shell8=1;
   break;
   case el_brick1:
      is_brick1=1;
   break;
   case el_wall1:
      is_wall1=1;
   break;
   case el_fluid2:
      is_fluid2=1;
   break;
   case el_fluid3:
      is_fluid3=1;
   break;
   case el_ale3:
      is_ale3=1;
   break;
   default:
      dserror("Unknown typ of element");
   break;   
   }
}/* end of loop over all elements */
/*--------------------- init the element routines for all present types */
/*------------------------------- init all kind of routines for shell8  */
if (is_shell8==1)
{
   shell8(actfield,actpart,NULL,NULL,&estif_global,&emass_global,&intforce_global,0,0,NULL,action);
}
/*-------------------------------- init all kind of routines for brick1 */
if (is_brick1==1)
{
   brick1(actpart,NULL,NULL,&estif_global,&emass_global,action);
}
/*-------------------------------- init all kind of routines for wall1 */ 
if (is_wall1==1)
{
   wall1(actpart,NULL,NULL,&estif_global,&emass_global,&intforce_global,0,NULL,action);
}
/*-------------------------------- init all kind of routines for fluid2 */
if (is_fluid2==1)
{
   fluid2(actpart,NULL,NULL,
          &estif_global,&emass_global,
          &etforce_global,&eiforce_global,&edforce_global,
          action,NULL);
}
/*-------------------------------- init all kind of routines for fluid3 */
if (is_fluid3==1)
{
   fluid3(actpart,NULL,NULL,
          &estif_global,&emass_global,
          &etforce_global,&eiforce_global,&edforce_global,
          action,NULL);
}
/*----------------------------------- init all kind of routines for ale */
if (is_ale3==1)
{
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of calinit_fluid */



