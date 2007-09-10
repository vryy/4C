/*!----------------------------------------------------------------------
\file
\brief Prepare calculation of stabilisation parameter

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6572
</pre>

------------------------------------------------------------------------*/
#ifndef CCADISCRET
/*!
\addtogroup FLUID2
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
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
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*!---------------------------------------------------------------------
\brief routine to prepare the determination of the stabilisation parameter

<pre>                                                        chfoe 10/04
This routine prepares the evaluation of the stabilisation parameter for
USFEM stabilised fluid elements in 2D if tau is calculated once per
element which currently is the default implementation.

</pre>
\param  *ele     ELEMENT	       (i)   actual element
\param **xzye    DOUBLE                (i)   nodal coordinates
\param  *funct   DOUBLE 	       (-)   shape functions
\param **deriv   DOUBLE 	       (-)   deriv. of shape funcs
\param **xjm     DOUBLE 	       (-)   jacobian matrix
\param **evelng  DOUBLE                (i)   (newest) element velocities
\param  *visc    DOUBLE 	       (-)   viscosity
\return void

------------------------------------------------------------------------*/
void f2_caltau(
	       ELEMENT         *ele,
	       DOUBLE         **xyze,
	       DOUBLE          *funct,
	       DOUBLE         **deriv,
	       DOUBLE         **xjm,
	       DOUBLE         **evelng,
               DOUBLE           visc
              )
{

INT            whichtau;
INT            which_hk;
INT            iel;

DOUBLE         velint[2];
DOUBLE         det;

FLUID_DYNAMIC *fdyn;
DIS_TYP        typ;

#ifdef DEBUG
dstrc_enter("f2_caltau");
#endif

/* The following includes the decision which stabilisation parameter
   will be used. It should rather be moved to the input. Current
   observations indicate however that the chosen selection behaves
   slightly better in most cases. It is therefore used preferabely while
   the other possibilities will be kept for research purpose
   for detailed explainations of the different tau concepts see
   preamble to routine 'f2_get_tau'                                     */

    whichtau = 0;   /* Franca  */
  /*whichtau = 1;*/ /* Whiting */


which_hk = 0; /* area square root */
/*which_hk = 1;*/ /* length in flow direction (Wall) */
/*which_hk = 2;*/ /* length in flow direction approximative (Codina) */
/*which_hk = 3;*/ /* length for anisotropic meshes (Codina) */

/*------------------------------------------------------- initialise ---*/
fdyn    = alldyn[genprob.numff].fdyn;

/*----------------------------------- set parameter for this element ---*/
iel    = ele->numnp;
typ    = ele->distyp;

/*------------ get shape functions and derivatives at element center ---*/
switch (ele->distyp)
{
  case quad4: case quad8: case quad9:
    f2_rec(funct,deriv,NULL,0.0,0.0,typ,2);
    break;
  case tri3: case tri6:
    f2_tri(funct,deriv,NULL,0.0,0.0,typ,2);
    break;
  default:
    dserror("typ unknown!");
} /* end switch(typ) */


/*--------------------------------- get velocities at element center ---*/
f2_veci(velint,funct,evelng,iel);

/*----------------------------- get Jacobian matrix and determinante ---*/
f2_jaco(xyze,deriv,xjm,&det,iel,ele);

/*--------------------------------------------- evaluate tau finally ---*/
f2_get_tau(ele,xjm,xyze,funct,det,velint,visc,whichtau,which_hk);



/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_caltau */

#endif
#endif
