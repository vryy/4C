/*!----------------------------------------------------------------------
\file
\brief integration loop for one fluid2 element using USFEM

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/gammi/
            +49-(0)89-289-15235
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID2
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid2/fluid2_prototypes.h"
#ifdef D_FLUID2_TDS
#include "fluid2_TDS_prototypes.h"

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

/*!---------------------------------------------------------------------
\brief integration loop for one fluid2 element using time dependent
subscales

<pre>                                                         gammi 11/06

In this routine the element 'stiffness' matrix and RHS for one
fluid2 element is calculated.
Stabilisation is performed using time dependent subscales, see
Codina, Principe, Guasch, Badia:
"Time dependent subscales in the stabilized finite element approximation
 of incompressible flow problems"
 

</pre>
\param  *ele	      ELEMENT	   (i)    actual element
\param  *hasext       INT          (i)    element flag
\param **estif        DOUBLE	   (o)    element stiffness matrix
\param  *eforce       DOUBLE	   (o)    element iter force vector
\param **xyze         DOUBLE       (-)    nodal coordinates
\param  *funct        DOUBLE	   (-)    natural shape functions
\param **deriv        DOUBLE	   (-)	  deriv. of nat. shape funcs
\param **deriv2       DOUBLE	   (-)    2nd deriv. of nat. shape f.
\param **xjm	      DOUBLE	   (-)    jacobian matrix
\param **derxy        DOUBLE	   (-)	  global derivatives
\param **derxy2       DOUBLE	   (-)    2nd global derivatives
\param **evelng       DOUBLE	   (i)    ele vel. at time n+g
\param **evhist       DOUBLE	   (i)    lin. combination of recent vel and acc
\param **egridv       DOUBLE	   (i)    grid velocity of element
\param  *epreng       DOUBLE	   (-)    ele pres. at time n+1
\param  *epren        DOUBLE	   (-)    ele pres. at time n
\param  *edeadng      DOUBLE	   (-)    ele dead load (selfweight) at n+1
\param  *edeadn       DOUBLE	   (-)    ele dead load (selfweight) at n
\param  *velint       DOUBLE	   (-)    vel at integration point
\param  *vel2int      DOUBLE	   (-)    vel at integration point
\param  *covint       DOUBLE	   (-)    conv. vel. at integr. point
\param **vderxy       DOUBLE	   (-)    global vel. derivatives
\param **vderxy2      DOUBLE	   (-)    2nd global vel. deriv.
\param **vderxy_old   DOUBLE	   (-)    global vel. derivatives
\param **vderxy2_old  DOUBLE	   (-)    2nd global vel. deriv.
\param **wa1	      DOUBLE	   (-)    working array
\param **wa2	      DOUBLE	   (-)    working array
\return void

------------------------------------------------------------------------*/
void f2_int_tds(
	              ELEMENT         *ele,
                      INT             *hasext,
                      DOUBLE         **estif,
	              DOUBLE          *eforce,
	              DOUBLE         **xyze,
	              DOUBLE          *funct,
	              DOUBLE         **deriv,
	              DOUBLE         **deriv2,
	              DOUBLE         **xjm,
	              DOUBLE         **derxy,
	              DOUBLE         **derxy2,
	              DOUBLE         **evelng,
	              DOUBLE         **eveln,
	              DOUBLE         **evhist,
	              DOUBLE         **egridv,
	              DOUBLE          *epreng,
	              DOUBLE          *epren,
	              DOUBLE          *edeadng,
	              DOUBLE          *edeadn,
	              DOUBLE         **vderxy,
                      DOUBLE         **vderxy2,
	              DOUBLE         **vderxy_old,
                      DOUBLE         **vderxy2_old,
	              DOUBLE         **eacc,
		      DOUBLE           visc,
	              DOUBLE         **wa1,
	              DOUBLE         **wa2,
                      DOUBLE           estress[3][MAXNOD_F2],
                      INT              is_relax
	             )
{
INT       i;          /* a counter                                       */
INT       dim;        /* a counter                                       */
INT       iel;        /* number of nodes                                */
INT       intc=0;     /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       nir=0,nis=0;/* number of integration nodesin r,s direction    */
INT       ihoel;      /* flag for higher order elements                 */
INT       icode;      /* flag for eveluation of shape functions         */
INT       lr, ls;     /* counter for integration                        */
INT       is_ale;
DOUBLE    fac;        /* total integration vactor                       */
DOUBLE    facr=0;     /* integration weights                            */
DOUBLE    facs=0;
DOUBLE    det;        /* determinant of jacobian matrix at time (n+1)   */
DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
DOUBLE    gradp[2];   /* pressure gradient at integration point         */
DOUBLE    gradp_old[2];/* pressure gradient at integration point        */
DOUBLE    hot    [2];   /* old higher order terms at integration point  */
DOUBLE    hot_old[2];   /* old higher order terms at integration point  */
DOUBLE    velint[2];  /* velocity vector at integration point           */
DOUBLE    sub_pres;   /* subscale pressure at integration point         */
DOUBLE    histvec[2]; /* history data at integration point              */
DOUBLE    gridvelint[2]; /* grid velocity                               */
DOUBLE    divuold;
DOUBLE    press;

DOUBLE    velint_old[2];
DOUBLE    sub_vel[2];


/* time increment of velocities */
double  time_der[2];

/* trial value for subscale velocity without factor tau/(tau+theta*dt)*/
DOUBLE    sub_vel_trial_wo_facMtau[2];

DOUBLE    acc_old[2];

DOUBLE    res    [2];
DOUBLE    res_old[2];

DIS_TYP   typ;	      /* element type                                   */

FLUID_DYNAMIC   *fdyn;
FLUID_DATA      *data;


/* time algorithm */
double   theta;
double   dt;

#ifdef DEBUG
    dstrc_enter("f2_int_tds");
#endif

/*--------------------------------------------------- initialisation ---*/
iel    = ele->numnp;
typ    = ele->distyp;
fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;
is_ale = ele->e.f2->is_ale;

dt     =fdyn->dt;
theta  =fdyn->theta;


/*------- get integraton data and check if elements are "higher order" */
switch (typ)
{
case quad4: case quad8: case quad9:  /* --> quad - element */
   icode   = 3; /* flag for higher order elements                 */
   ihoel   = 1;	/* flag for eveluation of shape functions         */
   /* initialise integration */
   nir = ele->e.f2->nGP[0];
   nis = ele->e.f2->nGP[1];
break;
case tri6: /* --> tri - element */
   icode   = 3; /* flag for higher order elements                 */
   ihoel   = 1;	/* flag for eveluation of shape functions         */
   /* initialise integration */
   nir  = ele->e.f2->nGP[0];
   nis  = 1;
   intc = ele->e.f2->nGP[1];
   break;
case tri3:
    ihoel  =0;  /* flag for higher order elements                 */
    icode  =2;  /* flag for eveluation of shape functions         */
   /* initialise integration */
   nir  = ele->e.f2->nGP[0];
   nis  = 1;
   intc = ele->e.f2->nGP[1];
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
      /*------- get values of  shape functions and their derivatives ---*/
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

      /*------------------------ compute Jacobian matrix at time n+1 ---*/
      f2_jaco(xyze,deriv,xjm,&det,iel,ele);
      fac = facr * facs * det;

      /*----------------------------------- compute global derivates ---*/
      f2_gder(derxy,deriv,xjm,det,iel);

      /*---------------- get velocities (n+1,i) at integration point ---*/
      f2_veci(velint,funct,evelng,iel);

      /*---------------- get velocities (n) at integration point ---*/
      f2_veci(velint_old,funct,eveln,iel);
      
      /*--------------------------- calculate the velocity increment ---*/
      for(dim=0;dim<2;dim++)
      {
	  time_der[dim]=velint[dim]-velint_old[dim];
      }
      
      /*---------------- get history data (n,i) at integration point ---*/
      f2_veci(histvec,funct,evhist,iel);

      /*----------------- get accelerations (n) at integration point ---*/
      f2_veci(acc_old,funct,eacc,iel);

      
      /*--------------------- get grid velocity at integration point ---*/
      if (is_ale)
      {
	f2_veci(gridvelint,funct,egridv,iel);
      }
      else
      {
        gridvelint[0] = 0;
        gridvelint[1] = 0;
      }

      /*-------- get velocity (n,i) derivatives at integration point ---*/
      f2_vder(vderxy_old,derxy,eveln,iel);

      /*---------- set old divergence */
      divuold = vderxy_old[0][0] + vderxy_old[1][1];

      /*- get velocity derivatives (n+1,i) and (n) at integration point */
      f2_vder(vderxy,derxy,evelng,iel);

      if (ihoel!=0)
      {
	  f2_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);

	  /*- get second velocity derivatives (n) at integration point */
	  f2_vder2(vderxy2_old,derxy2,eveln ,iel);
	  /*- get second velocity derivatives (n+1,i) at integration point */
	  f2_vder2(vderxy2    ,derxy2,evelng,iel);
      }

      /*------------------------------------- get pressure gradients ---*/
      gradp    [0] = gradp    [1] = 0.0;
      gradp_old[0] = gradp_old[1] = 0.0;

      for (i=0; i<iel; i++)
      {
         gradp    [0] += derxy[0][i] * epreng[i];
         gradp    [1] += derxy[1][i] * epreng[i];

         gradp_old[0] += derxy[0][i] * epren [i];
         gradp_old[1] += derxy[1][i] * epren [i];
      }

      press = 0;
      for (i=0;i<iel;i++)
      {
        press += funct[i]*epren[i];
      }

      /*---------- set subscale pressure */
      sub_pres=ele->e.f2->sub_pres.a.dv[lr*nis+ls];


      /*---------- set old subscale velocity */
      for(dim=0;dim<2;dim++)
      {
	  sub_vel    [dim]=ele->e.f2->sub_vel.a.da[dim][lr*nis+ls];
      }
	      
      if(ihoel!=0)
      {
	  hot_old[0]=0.5 * (2.0*vderxy2_old[0][0]
			    +
			    (vderxy2_old[0][1] + vderxy2_old[1][2]));
	  hot_old[1]=0.5 * (2.0*vderxy2_old[1][1]
			    +
			    (vderxy2_old[1][0] + vderxy2_old[0][2]));

	  hot[0]=0.5 * (2.0*vderxy2[0][0]
			    +
			    (vderxy2[0][1] + vderxy2[1][2]));
	  hot[1]=0.5 * (2.0*vderxy2[1][1]
			    +
			    (vderxy2[1][0] + vderxy2[0][2]));
	  
      }	
      else
      {
	  hot_old[0]=0;
	  hot_old[1]=0;

	  hot    [0]=0;
	  hot    [1]=0;
      }
	    
      /* calculate old and new residual without time derivative       */
      for(dim=0;dim<2;dim++)
      {
	  res_old[dim] =0;
	  res_old[dim]+=(velint_old[0]*vderxy_old[dim][0]
			 +
			 velint_old[1]*vderxy_old[dim][1]);
  
	  res_old[dim]-=2*visc*hot_old[dim];
	  
	  res_old[dim]+=gradp_old[dim];
	  
	  res_old[dim]-=edeadn [dim];

	  res    [dim] =0;
	  res    [dim]+=(velint[0]*vderxy[dim][0]
			 +
			 velint[1]*vderxy[dim][1]);
  
	  res    [dim]-=2*visc*hot[dim];
	  
	  res    [dim]+=gradp    [dim];
	  
	  res    [dim]-=edeadng[dim];
 
      }

      /*---------------------- get new estimate for subscale velocities */
      for(dim=0;dim<2;dim++)
      {
	  sub_vel_trial_wo_facMtau[dim]=
	      sub_vel[dim]
	      -time_der[dim]
	      +theta    *dt* res    [dim]
	      +(1-theta)*dt* res_old[dim]
	      -1./fdyn->tau_old[0] * (1-theta)*dt* sub_vel[dim];
      }
      /*-------------- perform integration for entire matrix and rhs ---*/
      f2_calmat_tds(estif,eforce,velint,histvec,gridvelint,press,vderxy,
		    vderxy2,gradp,funct,derxy,derxy2,edeadng,fac,
		    visc,iel,hasext,is_ale, is_relax, sub_pres, divuold,
		    sub_vel,sub_vel_trial_wo_facMtau,velint_old,acc_old,res_old);

   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */



/*------------------------------------------- assure assembly of rhs ---*/
if (!is_relax)
*hasext = 1;

#ifdef DEBUG
    dstrc_exit();
#endif
    return;
}


#endif /*D_FLUID2_TDS*/
#endif /*D_FLUID2    */


/*! @} (documentation module close)*/
