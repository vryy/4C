/*!----------------------------------------------------------------------
\file
\brief main routine fluid2 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID2
*//*! @{ (documentation module open)*/
#include "../headers/standardtypes.h"
#include "fluid2.h"
#include "fluid2_prototypes.h"
#include "../fluid2ml/fluid2ml_prototypes.h"
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

/*!---------------------------------------------------------------------
\brief main fluid2 control routine

<pre>                                                         genk 03/02
</pre>
\param  *actpart	 PARTITION    (i)
\param	*actintra	 INTRA        (i)
\param	*ele		 ELEMENT      (i)    actual element
\param	*estif_global	 ARRAY        (o)    element stiffness matrix
\param  *emass_global	 ARRAY        (o)    element mass matrix
\param  *eforce_global   ARRAY        (o)    element internal force vector
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
	    ARRAY       *eforce_global,
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

static INT              viscstr;
static FLUID_DYNAMIC   *fdyn;
#ifdef FLUID2_ML
static FLUID_DYN_ML    *mlvar;
static FLUID_ML_SMESH  *submesh;
static FLUID_ML_SMESH  *ssmesh;
static INT              ndum;       /* dummy variable                   */
static INT              xele,yele,zele;/* numb. of subm. ele. in x,y,z  */
INT smisal;
#endif
ARRAY_POSITION* ipos;

#ifdef DEBUG
dstrc_enter("fluid2");
#endif

if (container!=NULL)
  ipos = &(field[genprob.numff].dis[container->actndis].ipos);
else
  ipos = NULL;

/*------------------------------------------------- switch to do option */
switch (*action)
{
/*------------------------------------------------------ initialization */
case calc_fluid_init:
   fdyn   = alldyn[genprob.numff].fdyn;
   viscstr= alldyn[genprob.numff].fdyn->viscstr;
/*------------------------------------------- init the element routines */
   f2_intg(0);
   f2_calele(NULL,NULL,
             estif_global,emass_global,
	     eforce_global,edforce_global,
	     NULL,NULL,0,ipos,0,1);
   f2_iedg(NULL,ele,-1,1);
/*---------------------------------------------------- multi-level FEM? */
#ifdef FLUID2_ML
  if (fdyn->mlfem==1)
  {
    mlvar   = alldyn[genprob.numff].fdyn->mlvar;
    submesh = &(alldyn[genprob.numff].fdyn->mlvar->submesh);
    ssmesh  = &(alldyn[genprob.numff].fdyn->mlvar->ssmesh);
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
    f2_lsele(fdyn->data,mlvar,submesh,ssmesh,ele,estif_global,emass_global,
             eforce_global,edforce_global,ipos,hasdirich,hasext,1);
  }
#endif
break;

case calc_fluid_initvort:
/*------------------------------------------- init the element routines */
   f2_intg(0);
   f2_calvort(ele,1);
break;

/*------------------------------------------- call the element routines */
case calc_fluid:
/*---------------------------------------------------- multi-level FEM? */
#ifdef FLUID2_ML
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
  f2_lsele(fdyn->data,mlvar,submesh,ssmesh,ele,estif_global,emass_global,
           eforce_global,edforce_global,ipos,hasdirich,hasext,0);
}
else
#endif
   f2_calele(ele,eleke,
             estif_global,emass_global,
	     eforce_global,edforce_global,
	     hasdirich,hasext,actintra->intra_rank,ipos,container->is_relax,0);
break;

/*------------------------------------------- calculate fluid vorticity */
case calc_fluid_vort:
   f2_calvort(ele,0);
break;

/*-------------------------------------------- calculate fluid stresses */
case calc_fluid_stress:
   /*------ calculate stresses only for elements belonging to this proc */
   if (par.myrank==ele->proc)
      f2_stress(container->str,viscstr,ele,ipos,container->is_relax);
break;

/*------------------- calculate fluid stresses for lift&drag calculation */
case calc_fluid_liftdrag:
   if (ele->proc == actintra->intra_rank)
   {
      f2_stress(container->str,viscstr,ele,ipos,container->is_relax);
      f2_calliftdrag(ele,container);
   }
break;

/*-------------------- do stress projection for lower order elements ---*/
case calc_fluid_stressprojection:
   f2_calstresspro(ele,hasext,estif_global,eforce_global, ipos);
break;

/*--------------------------------- calculate curvature at free surface */
case calc_fluid_curvature:
   f2_curvature(ele,actintra->intra_rank);
break;

/*--------------------------------- calculate height function matrices */
case calc_fluid_heightfunc:
   f2_heightfunc(ele,estif_global,
                 eforce_global,ipos,container);
break;

/*---------------------------------------- calculate the shear stresses */
case calc_fluid_shearvelo:
   f2_shearstress(ele, ipos);
break;

case calc_fluid_normal:
   f2_calnormal(ele);
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
