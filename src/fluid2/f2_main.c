/*!----------------------------------------------------------------------
\file
\brief main routine fluid2 element

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
#include "../headers/standardtypes.h"
#include "fluid2.h"
#include "fluid2_prototypes.h"
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
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

/*!---------------------------------------------------------------------
\brief main fluid2 control routine

<pre>                                                         genk 03/02
</pre>
\param  *actpart	 PARTITION    (i)	    
\param	*actintra	 INTRA        (i)
\param	*ele		 ELEMENT      (i)    actual element
\param	*estif_global	 ARRAY        (o)    element stiffness matrix
\param  *emass_global	 ARRAY        (o)    element mass matrix
\param	*etforce_global  ARRAY        (o)    element time force vector
\param  *eiforce_global  ARRAY        (o)    element iter force vecotr
\param	*edforce_global  ARRAY        (o)    ele dirichl. force vector
\param	*action	         CALC_ACTION  (i)
\param	*hasdirich	 INT          (o)    flag
\param  *hasext          INT          (o)    flag
\param  *container       CONTAINER    (i)    container
\return void

------------------------------------------------------------------------*/
void fluid2(PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *ele,             
            ELEMENT     *eleke,             
            ARRAY       *estif_global,   
            ARRAY       *emass_global,   
	    ARRAY       *etforce_global, 
	    ARRAY       *eiforce_global, 
	    ARRAY       *edforce_global, 
            CALC_ACTION *action,
	    INT         *hasdirich,
	    INT         *hasext,
            CONTAINER   *container       
	   )
{
/*----------------------------------------------------------------------*/
#ifdef D_FLUID2 
/*----------------------------------------------------------------------*/

static INT              numff;      /* actual number of fluid field     */
static INT              viscstr;
static FLUID_DATA      *data;      
FLUID_DYNAMIC          *fdyn;
static FLUID_DYN_CALC  *dynvar;
static FLUID_DYN_ML    *mlvar;
static FLUID_ML_SMESH  *submesh;
static FLUID_ML_SMESH  *ssmesh;
static INT              ndum;       /* dummy variable                   */
static INT              xele,yele,zele;/* numb. of subm. ele. in x,y,z  */
FIELD                  *actfield;   /* actual field                     */
int smisal;

#ifdef DEBUG 
dstrc_enter("fluid2");
#endif

fdyn = alldyn[0].fdyn;
/*------------------------------------------------- switch to do option */
switch (*action)
{
/*------------------------------------------------------ initialization */
case calc_fluid_init:
/* ----------------------------------------- find number of fluid field */
   dynvar = &(alldyn[genprob.numff].fdyn->dynvar);
   data   = &(alldyn[genprob.numff].fdyn->dynvar.data);
   viscstr= alldyn[genprob.numff].fdyn->viscstr;
/*------------------------------------------- init the element routines */   
   f2_intg(data,0);
   f2_calele(data,dynvar,NULL,NULL,
             estif_global,emass_global,
	     etforce_global,eiforce_global,edforce_global,
	     NULL,NULL,0,0,1);
   f2_iedg(NULL,ele,-1,1);
/*---------------------------------------------------- multi-level FEM? */   
  if (fdyn->mlfem==1) 
  {  
    mlvar   = &(alldyn[numff].fdyn->mlvar);
    submesh = &(alldyn[numff].fdyn->mlvar.submesh);
    ssmesh  = &(alldyn[numff].fdyn->mlvar.ssmesh);
/*------- determine number of submesh elements in coordinate directions */   
    math_intextract(mlvar->smelenum,&ndum,&xele,&yele,&zele);
/*------------------------------------- create submesh on parent domain */   
    f2_pdsubmesh(submesh,xele,yele,mlvar->smorder,0);
/*-------------------- three-level FEM, i.e. dynamic subgrid viscosity? */   
    if (mlvar->smsgvi>2) 
    {
/*--- determine number of sub-submesh elements in coordinate directions */   
      math_intextract(mlvar->ssmelenum,&ndum,&xele,&yele,&zele);
/*--------------------------------- create sub-submesh on parent domain */   
      f2_pdsubmesh(ssmesh,xele,yele,mlvar->ssmorder,1);
    }
/*----------------------- init the element routines for multi-level FEM */   
    f2_lsele(data,dynvar,mlvar,submesh,ssmesh,ele,estif_global,emass_global,
             etforce_global,eiforce_global,edforce_global,hasdirich,hasext,1); 
  }     	    
break;

case calc_fluid_initvort:
/* ----------------------------------------- find number of fluid field */
   dynvar = &(alldyn[genprob.numff].fdyn->dynvar);
   data   = &(alldyn[genprob.numff].fdyn->dynvar.data);
/*------------------------------------------- init the element routines */
   f2_intg(data,0); 
   f2_calvort(data,dynvar,ele,1);  
break;

/*------------------------------------------- call the element routines */
case calc_fluid:
/*---------------------------------------------------- multi-level FEM? */   
if (fdyn->mlfem==1) 
{
  smisal = ele->e.f2->smisal;
  if (smisal!=1) 
  {
/*------------------------------ create element submesh if not yet done */   
    f2_elesubmesh(ele,submesh,0);
/*-------------------------- create element sub-submesh if not yet done */   
    if (mlvar->smsgvi>2) f2_elesubmesh(ele,ssmesh,1);
  }  
  f2_lsele(data,dynvar,mlvar,submesh,ssmesh,ele,estif_global,emass_global,
           etforce_global,eiforce_global,edforce_global,hasdirich,hasext,0);
}	      
else  
   f2_calele(data,dynvar,ele,eleke,
             estif_global,emass_global,
	     etforce_global,eiforce_global,edforce_global,
	     hasdirich,hasext,actintra->intra_rank,container->is_relax,0);
break;

/*------------------------------------------- calculate fluid vorticity */
case calc_fluid_vort:
   f2_calvort(data,dynvar,ele,0);
break;

/*-------------------------------------------- calculate fluid stresses */
case calc_fluid_stress:
    f2_stress(container->str,viscstr,data,ele,container->is_relax);
break;

/*------------------- calculate fluid stresses for lift&drag calculation */
case calc_fluid_liftdrag:
   if (ele->proc == actintra->intra_rank)
   {
      f2_stress(container->str,viscstr,data,ele,container->is_relax);
      f2_calliftdrag(ele,data,container);
   }
break;

/*--------------------------------- calculate curvature at free surface */
case calc_fluid_curvature:
   f2_curvature(data,dynvar,ele,actintra->intra_rank);
break;

/*---------------------------------------- calculate the shear stresses */
case calc_fluid_shearvelo:
   f2_shearstress(ele,dynvar);
break;

/*------------------------------------------------------- write restart */
case write_restart:
   f2_write_restart(ele,container->handsize,container->handles);
break;   
/*-------------------------------------------------------- read restart */
case read_restart:
   f2_read_restart(ele,container->handsize,container->handles);
break;   
/*----------------------------------------------------------------------*/
default:
   dserror("action unknown\n");
break;
} /* end swtich (*action) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
#endif
/*----------------------------------------------------------------------*/
return; 
} /* end of fluid2 */
/*! @} (documentation module close)*/
