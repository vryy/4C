/*!----------------------------------------------------------------------
\file
\brief main routine fluid3 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "../fluid3ml/fluid3ml_prototypes.h"
#include "fluid3.h"

#ifdef FLUID3_ML
#include "../fluid3ml/fluid3ml_prototypes.h"
#endif

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;   
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
 | global dense matrices for element routines             genk 04/02    |
 *----------------------------------------------------------------------*/
static FLUID_DATA      *data;

/*!---------------------------------------------------------------------                                         
\brief main fluid3 control routine

<pre>                                                         genk 05/02
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
\return void

------------------------------------------------------------------------*/
void fluid3(PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *ele,
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
#ifdef D_FLUID3 
static INT              viscstr;
static FLUID_DYNAMIC   *fdyn;
#ifdef FLUID3_ML
static FLUID_DYN_ML    *mlvar;
static FLUID_ML_SMESH  *submesh;
static FLUID_ML_SMESH  *ssmesh;
static INT              ndum;       /* dummy variable                   */
static INT              xele,yele,zele;/* numb. of subm. ele. in x,y,z  */
INT smisal;
#endif

INT       i;        /* simply a counter */
INT       ldflag;
GVOL     *actgvol;
GSURF    *actgsurf;
DSURF    *actdsurf;

#ifdef DEBUG 
dstrc_enter("fluid3");
#endif

/*------------------------------------------------- switch to do option */
switch (*action)
{
/*------------------------------------------------------ initialization */
case calc_fluid_init:
/* ----------------------------------------- find number of fluid field */
   fdyn   = alldyn[genprob.numff].fdyn;
   data   = alldyn[genprob.numff].fdyn->data;
   viscstr= alldyn[genprob.numff].fdyn->viscstr;

/*------------------------------------------- init the element routines */   
  f3_intg(data,0);
  f3_calele(data,NULL,estif_global,emass_global,etforce_global,
            eiforce_global,edforce_global,NULL,NULL,1);
/*---------------------------------------------------- multi-level FEM? */   
#ifdef FLUID3_ML
  if (fdyn->mlfem==1) 
  {  
    mlvar   = alldyn[genprob.numff].fdyn->mlvar;
    submesh = &(alldyn[genprob.numff].fdyn->mlvar->submesh);
    ssmesh  = &(alldyn[genprob.numff].fdyn->mlvar->ssmesh);
/*------- determine number of submesh elements in coordinate directions */   
    math_intextract(mlvar->smelenum,&ndum,&xele,&yele,&zele);
/*------------------------------------- create submesh on parent domain */   
    f3_pdsubmesh(submesh,xele,yele,zele,mlvar->smorder,0);
/*-------------------- three-level FEM, i.e. dynamic subgrid viscosity? */   
    if (mlvar->smsgvi>2) 
    {
/*--- determine number of sub-submesh elements in coordinate directions */   
      math_intextract(mlvar->ssmelenum,&ndum,&xele,&yele,&zele);
/*--------------------------------- create sub-submesh on parent domain */   
      f3_pdsubmesh(ssmesh,xele,yele,zele,mlvar->ssmorder,1);
    }
/*----------------------- init the element routines for multi-level FEM */   
    f3_lsele(data,mlvar,submesh,ssmesh,ele,estif_global,emass_global,
          etforce_global,eiforce_global,edforce_global,hasdirich,hasext,1); 
  }     	    
#endif
break;


/* call the element routines */
case calc_fluid:

/* multi-level FEM? */   
#ifdef FLUID3_ML
if (fdyn->mlfem==1) 
{
  smisal = ele->e.f3->smisal;
  if (smisal!=1) 
  {
/* create element submesh if not yet done */   
    f3_elesubmesh(ele,submesh,0);
/* create element sub-submesh if not yet done */   
    if (mlvar->smsgvi>2) f3_elesubmesh(ele,ssmesh,1);
  }  
  f3_lsele(data,mlvar,submesh,ssmesh,ele,estif_global,emass_global,
           etforce_global,eiforce_global,edforce_global,hasdirich,hasext,0);
}	      
else  
#endif
  f3_calele(data,ele,estif_global,emass_global,etforce_global,
                 eiforce_global,edforce_global,hasdirich,hasext,0);
break;

/*-------------------------------------------- calculate fluid stresses */
case calc_fluid_stress:
   /*------ calculate stresses only for elements belonging to this proc */
   if (par.myrank==ele->proc)
      f3_stress(container->str,viscstr,data,ele,container->is_relax);
break;

/* calculate fluid stresses for lift&drag calculation */
case calc_fluid_liftdrag:
  if (ele->proc == actintra->intra_rank)
  {
    /* check if element is on liftdrag-dline */
    actgvol=ele->g.gvol;
    ldflag=0;
    for (i=0;i<actgvol->ngsurf;i++)
    {
      actgsurf=actgvol->gsurf[i];
      actdsurf=actgsurf->dsurf;
      if (actdsurf==NULL) continue;
      if (actdsurf->liftdrag==NULL) continue;
      ldflag++;
      break;
    }
    if (ldflag>0)
    {
      f3_stress(container->str,viscstr,data,ele,container->is_relax);
      f3_liftdrag(ele,data,container);
    }
  }
break;

/*--------------------------------- calculate height function matrices */
case calc_fluid_heightfunc:
   f3_heightfunc(data,ele,estif_global,
                 eiforce_global,container);
break;

/*--------------------------- calculate element stabilisation parameter */
case calc_fluid_stab:
   f3_calstab(ele,data);
break;

/*----------------------------------------------------------------------*/
default:
   dserror("action unknown\n");
break;
} /* end switch (*action) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
#endif
/*----------------------------------------------------------------------*/
return; 
} /* end of fluid3 */

