/*!----------------------------------------------------------------------
\file
\brief integration loop for one fluid2 element

<pre>
Maintainer: Thomas Hettich
            hettich@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hettich/
            0771 - 685-6575
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID2TU
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2_tu.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
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

static FLUID_DYNAMIC   *fdyn;
/*----------------------------------------------------------------------*
 | integration loop for one fluid element                               |
 |                                                                      |
 |                                                          he    12/02 |
 *----------------------------------------------------------------------*/
/*!---------------------------------------------------------------------
\brief integration loop for one fluid2_tu element

<pre>                                                       he   12/02

In this routine the element stiffness matrix, iteration-RHS and
time-RHS for one fluid2 element is calculated

</pre>
\param  *ele	 ELEMENT	   (i)    actual element
\param  *elev	 ELEMENT	   (i)    actual rans-element
\param **estif     DOUBLE	   (o)    element stiffness matrix
\param **emass     DOUBLE	   (o)    element mass matrix
\param  *etforce   DOUBLE	   (o)    element time force vector
\param  *eiforce   DOUBLE	   (o)    element iter force vector
\param  *eproforce DOUBLE	   (o)    element production force vector
\param  *funct     DOUBLE	   (-)    natural shape functions
\param **deriv     DOUBLE	   (-)	  deriv. of nat. shape funcs
\param **deriv2    DOUBLE	   (-)    2nd deriv. of nat. shape f.
\param **xjm	 DOUBLE	   (-)    jacobian matrix
\param **xyze      DOUBLE        (-)    nodal coordinates
\param **derxy     DOUBLE	   (-)	  global derivatives
\param **derxy2    DOUBLE	   (-)    2nd global derivatives
\param  *kapepsn   DOUBLE	   (i)    kapeps at time n
\param  *kapepsg   DOUBLE	   (i)    kapeps at time n+g
\param  *eddyg     DOUBLE	   (i)    eddy-viscosity at time n+1
\param  *eddypro   DOUBLE	   (i)    eddy-viscosity for prod. term
\param  *kappa     DOUBLE	   (i)    kappa
\param  *kappan    DOUBLE	   (i)    kappa at time n
\param  *epsilon   DOUBLE	   (-)    epsilon at nodes
\param  *kapepspro    DOUBLE    (i)     kapeps for production term
\param  *kapepsderxy  DOUBLE	   (-)    kapeps derivates
\param  *kapepsderxy2 DOUBLE	   (-)    second kapeps derivatives
\param  *velint    DOUBLE	   (-)    vel at integration point
\param  *velint_dc    DOUBLE	   (-)    vel at integration point for DISC. CAPT.
\param **evel	 DOUBLE	   (i)    velocity at nodes
\param **vderxy	 DOUBLE	   (-)    velocity derivates
\param **vderxy2	 DOUBLE	   (-)    velocity 2nd derivates
\param **wa1	 DOUBLE	   (-)    working array
\param **wa2	 DOUBLE	   (-)    working array
\return void

------------------------------------------------------------------------*/
void f2_calint_tu(
	         ELEMENT         *ele,
	         ELEMENT         *elev,
               DOUBLE         **estif,
	         DOUBLE         **emass,
	         DOUBLE          *etforce,
	         DOUBLE          *eiforce,
	         DOUBLE          *eproforce,
	         DOUBLE          *funct,
	         DOUBLE         **deriv,
	         DOUBLE         **deriv2,
	         DOUBLE         **xjm,
	         DOUBLE         **xyze,
	         DOUBLE         **derxy,
	         DOUBLE         **derxy2,
	         DOUBLE          *kapepsn,
	         DOUBLE          *kapepsg,
               DOUBLE          *eddyg,
               DOUBLE          *eddypro,
               DOUBLE          *kappa,
               DOUBLE          *kappan,
               DOUBLE          *epsilon,
               DOUBLE          *kapepspro,
               DOUBLE          *kapepsderxy,
               DOUBLE          *kapepsderxy2,
	         DOUBLE          *velint,
	         DOUBLE          *velint_dc,
               DOUBLE         **evel,
               DOUBLE         **vderxy,
               DOUBLE         **vderxy2,
               DOUBLE         **wa1,
	         DOUBLE         **wa2
	        )
{
INT       i,j;        /* simply a counter                               */
INT       iel;        /* number of nodes                                */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       nir,nis;    /* number of integration nodesin r,s direction    */
INT       actmat;     /* material number of the element                 */
INT       icode=2;    /* flag for eveluation of shape functions         */
INT       ihoel=0;    /* flag for higher order elements                 */
INT       lr, ls;     /* counter for integration                        */
DOUBLE    dens;       /* density                                        */
DOUBLE    visc;       /* viscosity                                      */
DOUBLE    C_u,C_2;
DOUBLE    vderxy_12;
DOUBLE    fac;
DOUBLE    factor,factor1;
DOUBLE    factor2,sig,production;
DOUBLE    facr, facs; /* integration weights                            */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
DOUBLE    eddyint;
DOUBLE    eddynint;
DOUBLE    kappaint;
DOUBLE    kappanint;
DOUBLE    kapepsint;
DOUBLE    kapepsnint;
DOUBLE    eps_proint;
DIS_TYP   typ;	    /* element type                                   */
FLUID_DATA  *data;

#ifdef DEBUG
dstrc_enter("f2_calint_tu");
#endif

/*----------------------------------------------------- initialisation */
iel=ele->numnp;
actmat=ele->mat-1;
dens = mat[actmat].m.fluid->density;
visc = mat[actmat].m.fluid->viscosity;
typ  = ele->distyp;

fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;

/*------- get integraton data and check if elements are "higher order" */
switch (typ)
{
case quad4: case quad8: case quad9:  /* --> quad - element */
   icode   = 3;
   ihoel   = 1;
   /* initialise integration */
   nir = ele->e.f2_tu->nGP[0];
   nis = ele->e.f2_tu->nGP[1];
break;
case tri3: case tri6: /* --> tri - element */
   if (iel>3)
   {
      icode   = 3;
      ihoel   = 1;
   }
   /* initialise integration */
   nir  = ele->e.f2_tu->nGP[0];
   nis  = 1;
   intc = ele->e.f2_tu->nGP[1];
break;
default:
   dserror("typ unknown!");
} /* end switch(typ) */

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
for (lr=0;lr<nir;lr++)
{
   for (ls=0;ls<nis;ls++)
   {
/*--------------- get values of  shape functions and their derivatives */
      switch(typ)
      {
      case quad4: case quad8: case quad9:   /* --> quad - element */
	 e1   = data->qxg[lr][nir-1];
       facr = data->qwgt[lr][nir-1];
	 e2   = data->qxg[ls][nis-1];
	 facs = data->qwgt[ls][nis-1];
       f2_rec(funct,deriv,deriv2,e1,e2,typ,icode);
      break;
      case tri3: case tri6:   /* --> tri - element */
	 e1   = data->txgr[lr][intc];
	 facr = data->twgt[lr][intc];
	 e2   = data->txgs[lr][intc];
	 facs = ONE;
	 f2_tri(funct,deriv,deriv2,e1,e2,typ,icode);
      break;
      default:
         dserror("typ unknown!");
      } /* end switch(typ) */
/*-------------------------------------------- compute Jacobian matrix */
      f2_jaco(xyze,funct,deriv,xjm,&det,iel,ele);
      fac = facr*facs*det;

/*------------------------------------------- compute global derivates */
      f2_gder(derxy,deriv,xjm,det,iel);
      f2_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);

/*--------------------- get element vel. and vel. at integration point */
      f2_veci(velint,funct,evel,iel);

/*-------------------------------- get eddy-visc. at integration point */
      f2_eddyi(&eddyint,funct,eddyg,iel);

/*---------------------------------- compute global derivates for vel. */
      f2_vder(vderxy,derxy,evel,iel);

/*---------------------------- get kapeps (n+1,i)  at integraton point */
      f2_kapepsi(&kapepsint,funct,kapepsg,iel);

/*------------------ calculate stabilisation parameter for DISC. CAPT. */
      f2_kapepsder(kapepsderxy,derxy,kapepsn,iel);
      f2_vel_dc(velint,velint_dc,kapepsderxy);

/*---------------- get kapeps (n+1,i) derivatives at integration point */
      f2_kapepsder(kapepsderxy,derxy,kapepsg,iel);

/*------------- get factors for LOW-REYNOLD's Model for kappa-equation */
      if(fdyn->kapeps_flag==0)
      {
/*--------- calculate C_u with R_t and factors for LOW-REYNOLD's Model */
       f2_C_kappa(kapepsint,epsilon,funct,visc,&C_u,iel);
       f2_fac_kappa(C_u,eddyint,&factor,&factor1,&factor2,&sig);
/*--------------------------------get the Residuum at integraton point */
      }

/*----------- get factors for LOW-REYNOLD's Model for epsilon-equation */
      if(fdyn->kapeps_flag==1)
      {
/*------- calculate grad(grad u): grad(grad u) for LOW-REYNOLD's Model */
       f2_vder2(vderxy2,derxy2,evel,iel);
       f2_v(vderxy2,&vderxy_12);
/*--------- calculate C_2 with R_t and factors for LOW-REYNOLD's Model */
       f2_kappai_tu(&kappaint,&kappanint,&eps_proint,funct,kappa,kappan,
                    kapepspro,iel);
       f2_C_eps(kapepsint,kappaint,visc,&C_2,iel);
       f2_fac_eps(C_2,eps_proint,kappaint,kappanint,&factor,&factor1,
                  &factor2,&sig);
      }

/*----------------------------------------------------------------------*
 |         compute "Standard Galerkin" matrices                         |
 | NOTE:                                                                |
 |  Standard Galerkin matrices are all stored in one matrix "estif"     |
 |  Standard Galerkin mass matrix is stored in "emass"                  |
 *----------------------------------------------------------------------*/

 /*--------------------------------------------- compute matrix Kkapeps */
       f2_calkkapeps(estif,kapepsint,velint,eddyint,kapepsderxy,
                     funct,derxy,fac,visc,factor,sig,iel);

/*---------------------------------------------- compute matrix Mkapeps */
	 if (fdyn->nis==0)
          f2_calmkapeps(emass,funct,fac,iel);

/*----------------------------------------------------------------------*
 |         compute Stabilisation matrices                               |
 | NOTE:                                                                |
 |  Stabilisation matrices are all stored in one matrix "estif"         |
 |  Stabilisation mass matrices are all stored in one matrix "emass"    |
 *----------------------------------------------------------------------*/
/*------------------------------------ stabilisation for matrix Kkapeps */
       f2_calstabkkapeps(ele,elev,estif,kapepsint,velint,velint_dc,
                         eddyint,kapepsderxy,funct,derxy,derxy2,fac,visc,
                         factor,sig,iel);
/*------------------------------------ stabilisation for matrix Mkapeps */

       if (fdyn->nis==0)
          f2_calstabmkapeps(ele,emass,velint,velint_dc,funct,derxy,fac,iel);
/*----------------------------------------------------------------------*
 |         compute "Iteration" Force Vectors                            |
 *----------------------------------------------------------------------*/

/*------------- calculate galerkin part of "Iter-RHS" (kapeps dof) */
         f2_calgalifkapeps(eiforce,eddyint,kapepsint,funct,fac,
                           factor2,vderxy_12,visc,iel);

/*---------------- calculate stabilisation for "Iter-RHS" (kapeps dof) */
         f2_calstabifkapeps(ele,eiforce,kapepsint,velint,velint_dc,
                            eddyint,funct,derxy,fac,factor2,vderxy_12,visc,iel);

/*----------------------------------------------------------------------*
 |         compute Production "Time" Force Vector                        |
 *----------------------------------------------------------------------*/
      if (fdyn->niturbu_pro!=0)
      {
/*-------------------------------- get eddy-visc. at integration point */
        f2_eddyi(&eddynint,funct,eddypro,iel);

/*------------ calculate production-term for "PROTime-RHS" (kapeps-dofs)*/
        f2_production(vderxy,&production);

/*--------------- calculate galerkin part of "PROTime-RHS" (kapeps-dofs)*/
        f2_calgalprofkapeps(eproforce,eddynint,funct,visc,fac,factor1,
                            production,iel);

/*------------- calculate stabilisation for "PROTime-RHS" (kapeps-dofs) */
        f2_calstabprofkapeps(ele,eproforce,eddynint,funct,visc,fac,
                             factor1,production,velint,velint_dc,derxy,iel);
      }
/*----------------------------------------------------------------------*
 |         compute "Time" Force Vector                                  |
 *----------------------------------------------------------------------*/
      if (fdyn->niturbu_n!=0)
      {
/*-------------------------------- get kapeps (n) at integration point */
	    f2_kapepsi(&kapepsint,funct,kapepsn,iel);

/*--------------------- get kapepsderivatives (n) at integration point */
          f2_kapepsder(kapepsderxy,derxy,kapepsn,iel);

/*-------------------- get kapepsderivatives2 (n) at integration point */
          f2_kapepsder2(kapepsderxy2,derxy2,kapepsn,iel);

/*----------------- calculate galerkin part of "Time-RHS" (kapeps-dofs)*/
          f2_calgaltfkapeps(etforce,kapepsint,velint,eddynint,funct,
                            derxy,vderxy,kapepsderxy,visc,fac,factor,factor1,
                            factor2,sig,vderxy_12,production,iel);

/*--------------- calculate stabilisation for "Time-RHS" (kapeps-dofs) */
          f2_calstabtfkapeps(ele,etforce,kapepsint,velint,velint_dc,
                             eddynint,derxy,kapepsderxy2,vderxy,kapepsderxy,
                             fac,visc,factor,factor1,factor2,sig,vderxy_12,
                             production,iel);
     } /* endif fdyn->niturbu_n */
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */

  /*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calint */



#endif
