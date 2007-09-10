/*!----------------------------------------------------------------------
\file
\brief evaluate 2D fluid coefficient matrix

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/gammi/
            +49-(0)89-289-15235
</pre>

------------------------------------------------------------------------*/
#ifndef CCADISCRET
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

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

static FLUID_DYNAMIC *fdyn;


/*!---------------------------------------------------------------------
\brief update of time dependent subscales

<pre>                                                        gammi 11/06

The time dependent pressure subscales are updated according to a
one step theta timestepping scheme.

                  dp_sub       1
                  ------ = - ----- * p_sub + res_C(u)
                    dt       tau_C

Here, res_C(u)=div(u) is the residual of the continuity equation.

</pre>
\param  *actpart       PARTITION        (i)
\param  *actintra      INTRA            (i)
\param  *actfield      FIELD            (i)
\param  *ipos          ARRAY_POSITION   (i)
\param   disnum_calc   INT              (i)
\return void

------------------------------------------------------------------------*/


void f2_update_subscale_pres(
    PARTITION      *actpart,
    INTRA          *actintra,
    FIELD          *actfield,
    ARRAY_POSITION *ipos,
    INT             disnum_calc)
{
/* multi purpose counter */
int      i;

/* counter for elements */
int      nele;

/* element related data */
int      ls,lr;
int      nis=0,nir=0;
int      icode,intc=0;
int      ihoel;
double   det;
DOUBLE   e1,e2;             /* natural coordinates of integr. point */
FLUID_DATA      *data;      /* Information on Gausspoints etc.      */

ELEMENT *ele;
NODE    *actnode;

/* material properties */
double   visc;

/* time algorithm */
double   theta;
double   dt;

/* constants containing dt and the stabilisation parameter */
double   facC,facCtau;
double   old_facC,old_facCtau;

/* new and old residual of the coninuity equation */
double   divu,divu_old;

/* the old subscale pressure */
double   sp_old;

/* a constant, */
double   CRHS;

static ARRAY     xyze_a;   /* actual element coordinates                */
static DOUBLE  **xyze;
static ARRAY     xjm_a;    /* jocobian matrix                           */
static DOUBLE  **xjm;
static ARRAY     funct_a;  /* shape functions                           */
static DOUBLE   *funct;
static ARRAY     deriv_a;  /* first natural derivatives                 */
static DOUBLE  **deriv;
static ARRAY     eveln_a;  /* element velocities at (n)                 */
static DOUBLE  **eveln;
static ARRAY     evelng_a; /* element velocities at (n+gamma)           */
static DOUBLE  **evelng;
static ARRAY     derxy_a;  /* coordinate - derivatives                  */
static DOUBLE  **derxy;
static ARRAY     vderxy_a; /* vel - derivatives                         */
static DOUBLE  **vderxy;
static ARRAY     deriv2_a; /* second natural derivatives                */
static DOUBLE  **deriv2;


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("f2_update_subscale_pres");
#endif

/*========================== initialisation ============================*/
fdyn = alldyn[genprob.numff].fdyn;
data   = fdyn->data;

dt   =fdyn->dt;
theta=fdyn->theta;

xyze      = amdef("xyze"     ,&xyze_a     ,2,MAXNOD_F2,"DA");
xjm       = amdef("xjm"      ,&xjm_a      ,2,2        ,"DA");
funct     = amdef("funct"    ,&funct_a    ,MAXNOD_F2,1,"DV");
deriv     = amdef("deriv"    ,&deriv_a    ,2,MAXNOD_F2,"DA");
eveln     = amdef("eveln"    ,&eveln_a    ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
evelng    = amdef("evelng"   ,&evelng_a   ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
vderxy    = amdef("vderxy"   ,&vderxy_a   ,2,2,"DA");
derxy     = amdef("derxy"    ,&derxy_a    ,2,MAXNOD_F2,"DA");
deriv2     = amdef("deriv2"   ,&deriv2_a   ,3,MAXNOD_F2,"DA");



for(nele=0;nele<actpart->pdis[disnum_calc].numele;nele++)
{
    ele=actpart->pdis[disnum_calc].element[nele];

/*------- get integraton data and check if elements are "higher order" */
    switch (ele->distyp)
    {
	case quad4: case quad8: case quad9:  /* --> quad - element */
	    icode   = 3; /* flag for higher order elements                 */
	    ihoel   = 1; /* flag for eveluation of shape functions         */
	    nir = ele->e.f2->nGP[0];
	    nis = ele->e.f2->nGP[1];
	    break;
	case tri6: /* --> tri - element */
	    icode   = 3; /* flag for higher order elements                 */
	    ihoel   = 1; /* flag for eveluation of shape functions         */
	    nir  = ele->e.f2->nGP[0];
	    nis  = 1;
	    break;
	case tri3:
	    ihoel  =0;  /* flag for higher order elements                 */
	    icode  =2;  /* flag for eveluation of shape functions         */
	    nir  = ele->e.f2->nGP[0];
	    nis  = 1;
	    break;
	default:
	    dserror("typ unknown!");
    } /* end switch(typ) */


    /*------------------------------------------ set element coordinates -*/
    for(i=0;i<ele->numnp;i++)
    {
	xyze[0][i]=ele->node[i]->x[0];
	xyze[1][i]=ele->node[i]->x[1];
    }

    /* -> implicit time integration method ---------*/
    for(i=0;i<ele->numnp;i++) /* loop nodes of element */
    {
	actnode=ele->node[i];
        /*----------------------------- set recent element velocities */
	evelng[0][i]=actnode->sol_increment.a.da[ipos->velnp][0];
	evelng[1][i]=actnode->sol_increment.a.da[ipos->velnp][1];
	eveln [0][i]=actnode->sol_increment.a.da[ipos->veln ][0];
	eveln [1][i]=actnode->sol_increment.a.da[ipos->veln ][1];

    } /* end of loop over nodes of element */


    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    f2_get_time_dependent_sub_tau(ele,xyze,funct,deriv,evelng,eveln,visc);

    facC   =1./(fdyn->tau[2]+theta*dt);
    facCtau=fdyn->tau[2]*facC;

    old_facC   = 1./(fdyn->tau_old[2]+theta*dt);
    old_facCtau= fdyn->tau_old[2]*facC;


/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
    for (lr=0;lr<nir;lr++)
    {
	for (ls=0;ls<nis;ls++)
	{
	    /*------- get values of  shape functions and their derivatives ---*/
	    switch(ele->distyp)
	    {
		case quad4: case quad8: case quad9:   /* --> quad - element */
		    e1   = data->qxg[lr][nir-1];
		    e2   = data->qxg[ls][nis-1];
		    f2_rec(funct,deriv,deriv2,e1,e2,ele->distyp,icode);
		    break;
		case tri3: case tri6:   /* --> tri - element */
		    e1   = data->txgr[lr][intc];
		    e2   = data->txgs[lr][intc];
		    f2_tri(funct,deriv,deriv2,e1,e2,ele->distyp,icode);
		    break;
		default:
		    dserror("typ unknown!");
	    } /* end switch(typ) */


	    /*------------------ compute Jacobian matrix at time n+1 ---*/
	    f2_jaco(xyze,deriv,xjm,&det,ele->numnp,ele);

	    /*----------------------------- compute global derivates ---*/
	    f2_gder(derxy,deriv,xjm,det,ele->numnp);

	    /*--- get velocity (n+1,i) derivatives at integration point */
	    f2_vder(vderxy,derxy,evelng,ele->numnp);

	    divu     = vderxy[0][0] + vderxy[1][1];

	    /*--- get velocity (n+1,i) derivatives at integration point */
	    f2_vder(vderxy,derxy,eveln,ele->numnp);

	    divu_old = vderxy[0][0] + vderxy[1][1];

	    sp_old=ele->e.f2->sub_pres.a.dv[lr*nis+ls];

	    CRHS=sp_old
		-
		(1.-theta)*dt*(sp_old/fdyn->tau_old[2]+divu_old);

	    ele->e.f2->sub_pres.a.dv[lr*nis+ls]=
		facCtau*(CRHS-theta*dt*divu);


	} /* end of loop over integration points ls*/
    } /* end of loop over integration points lr */
}


/* clean up ------------------------------------------------------------*/

amdel(&xyze_a);
amdel(&xjm_a);
amdel(&funct_a);
amdel(&deriv_a);
amdel(&eveln_a);
amdel(&evelng_a);
amdel(&derxy_a);
amdel(&vderxy_a);
amdel(&deriv2_a);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
}


/*!---------------------------------------------------------------------
\brief update of time dependent subscale velocities

<pre>                                                        gammi 11/06

The time dependent velocity subscales are updated according to a
one step theta time integration scheme of the equation

                  du_sub       1
                  ------ = - ----- * u_sub + res_M(u)
                    dt       tau_M

Here, res_M(u) is the residual of the momentum equation and contains a
time derivative, a convective term, a diffusion term, the pressure
gradient and the volume force.

</pre>
\param  *actpart       PARTITION        (i)
\param  *actintra      INTRA            (i)
\param  *actfield      FIELD            (i)
\param  *ipos          ARRAY_POSITION   (i)
\param   disnum_calc   INT              (i)
\return void

------------------------------------------------------------------------*/


void f2_update_subscale_vel(
    PARTITION      *actpart,
    INTRA          *actintra,
    FIELD          *actfield,
    ARRAY_POSITION *ipos,
    INT             disnum_calc)
{
/* multi purpose counter */
int      i;

/* counter for dimensions (x,y) */
int     dim;

/* counter for elements */
int      nele;

/* element related data */
DOUBLE   e1,e2;             /* natural coordinates of integr. point */
int      ls,lr;             /* counters for Gausspoints             */
int      nis=0,nir=0;       /* total number of integration points   */
int      icode=0,intc=0;    /* required for tri elements            */
int      ihoel=0;           /* higher order element flag            */
double   det;               /* Jacobideterminant                    */

ELEMENT *ele;
NODE    *actnode;

FLUID_DATA      *data;      /* Information on Gausspoints etc.      */

/* material properties */
double   visc;

/* time algorithm */
double  theta;
double  dt;

/* constants containing dt and the stabilisation parameter */
double  facM,facMtau;
double  old_facM,old_facMtau;

/* vector of subscale velocities */
double  sv_old[2];
/* residuals of both timesteps (without time derivative) */
double  res_old[2];
double  res_new[2];

/* time increment of velocities */
double  time_der[2];

/* factors from curves for timedependent deadload */
double  acttimefacn,acttimefac;

/* velocity vectors at integration point      */
DOUBLE  velint[2],velint_old[2];
/* pressure and pressure gradients */
DOUBLE  press    ,press_old    ;
DOUBLE  gradp [2],gradp_old [2];
/* higher order terms */
DOUBLE  hot   [2],hot_old   [2];

static ARRAY     xyze_a;   /* actual element coordinates                */
static DOUBLE  **xyze;
static ARRAY     xjm_a;    /* jocobian matrix                           */
static DOUBLE  **xjm;
static ARRAY     funct_a;  /* shape functions                           */
static DOUBLE   *funct;
static ARRAY     deriv_a;  /* first natural derivatives                 */
static DOUBLE  **deriv;
static ARRAY     eveln_a;  /* element velocities at (n)                 */
static DOUBLE  **eveln;
static ARRAY     evelng_a; /* element velocities at (n+gamma)           */
static DOUBLE  **evelng;
static ARRAY     derxy_a;  /* coordinate - derivatives                  */
static DOUBLE  **derxy;
static ARRAY     vderxy_a; /* vel - derivatives                         */
static DOUBLE  **vderxy;
static ARRAY     vderxy_old_a; /* vel - derivatives                     */
static DOUBLE  **vderxy_old;
static ARRAY     epren_a;  /* element pressures at (n)	                */
static DOUBLE   *epren;
static ARRAY     epreng_a; /* element pressures at (n)	                */
static DOUBLE   *epreng;
static ARRAY     edeadn_a; /* element dead load (selfweight)            */
static DOUBLE   *edeadng;
static ARRAY     edeadng_a;/* element dead load (selfweight)            */
static DOUBLE   *edeadn;
static ARRAY     w1_a;     /* working array of arbitrary chosen size    */
static DOUBLE  **wa1;      /* used in different element routines        */
static ARRAY     w2_a;     /* working array of arbitrary chosen size    */
static DOUBLE  **wa2;      /* used in different element routines        */
static ARRAY     vderxy2_a;/* vel - 2nd derivatives                     */
static DOUBLE  **vderxy2;
static ARRAY     vderxy2_old_a;/* vel - 2nd derivatives                     */
static DOUBLE  **vderxy2_old;
static ARRAY     derxy2_a; /* 2nd coordinate - derivatives              */
static DOUBLE  **derxy2;
static ARRAY     deriv2_a; /* second natural derivatives                */
static DOUBLE  **deriv2;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("f2_update_subscale_vel");
#endif

/*========================== initialisation ============================*/
fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;

dt   =fdyn->dt;
theta=fdyn->theta;

xyze       = amdef("xyze"      ,&xyze_a      ,2,MAXNOD_F2,"DA");
xjm        = amdef("xjm"       ,&xjm_a       ,2,2        ,"DA");
funct      = amdef("funct"     ,&funct_a     ,MAXNOD_F2,1,"DV");
deriv      = amdef("deriv"     ,&deriv_a     ,2,MAXNOD_F2,"DA");
deriv2     = amdef("deriv2"   ,&deriv2_a   ,3,MAXNOD_F2,"DA");
eveln      = amdef("eveln"     ,&eveln_a     ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
evelng     = amdef("evelng"     ,&evelng_a   ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
vderxy     = amdef("vderxy"    ,&vderxy_a    ,2,2,"DA");
vderxy_old = amdef("vderxy_old",&vderxy_old_a,2,2,"DA");
vderxy2    = amdef("vderxy2"  ,&vderxy2_a  ,2,3,"DA");
vderxy2_old = amdef("vderxy2_old"  ,&vderxy2_old_a  ,2,3,"DA");
derxy      = amdef("derxy"     ,&derxy_a     ,2,MAXNOD_F2,"DA");
derxy2     = amdef("derxy2"   ,&derxy2_a   ,3,MAXNOD_F2,"DA");
epren      = amdef("epren"     ,&epren_a     ,MAXNOD_F2,1,"DV");
epreng     = amdef("epreng"    ,&epreng_a    ,MAXNOD_F2,1,"DV");
edeadn     = amdef("edeadn"    ,&edeadn_a    ,2,1,"DV");
edeadng    = amdef("edeadng"   ,&edeadng_a   ,2,1,"DV");
wa1        = amdef("wa1"      ,&w1_a       ,MAXDOFPERELE,MAXDOFPERELE,"DA");
wa2        = amdef("wa2"      ,&w2_a       ,MAXDOFPERELE,MAXDOFPERELE,"DA");



for(nele=0;nele<actpart->pdis[disnum_calc].numele;nele++)
{
    ele=actpart->pdis[disnum_calc].element[nele];

/*------- get integraton data and check if elements are "higher order" */
    switch (ele->distyp)
    {
	case quad4: case quad8: case quad9:  /* --> quad - element */
	    icode   = 3; /* flag for higher order elements                 */
	    ihoel   = 1; /* flag for eveluation of shape functions         */
	    nir  = ele->e.f2->nGP[0];
	    nis  = ele->e.f2->nGP[1];
	    intc = ele->e.f2->nGP[1];
	    break;
	case tri6: /* --> tri - element */
	    icode   = 3; /* flag for higher order elements                 */
	    ihoel   = 1; /* flag for eveluation of shape functions         */
	    nir  = ele->e.f2->nGP[0];
	    nis  = 1;
	    intc = ele->e.f2->nGP[1];
	    break;
	case tri3:
	    ihoel  =0;  /* flag for higher order elements                 */
	    icode  =2;  /* flag for eveluation of shape functions         */
	    nir  = ele->e.f2->nGP[0];
	    nis  = 1;
	    intc = ele->e.f2->nGP[1];
	    break;
	default:
	    dserror("typ unknown!");
    } /* end switch(typ) */


    /*------------------------------------------ set element coordinates -*/
    for(i=0;i<ele->numnp;i++)
    {
	xyze[0][i]=ele->node[i]->x[0];
	xyze[1][i]=ele->node[i]->x[1];
    }

    /* -> implicit time integration method ---------*/
    for(i=0;i<ele->numnp;i++) /* loop nodes of element */
    {
	actnode=ele->node[i];
        /*----------------------------- set recent element velocities */
	evelng[0][i]=actnode->sol_increment.a.da[ipos->velnp][0];
	evelng[1][i]=actnode->sol_increment.a.da[ipos->velnp][1];
	eveln [0][i]=actnode->sol_increment.a.da[ipos->veln ][0];
	eveln [1][i]=actnode->sol_increment.a.da[ipos->veln ][1];

        /*--------------------------------------------- and pressures */
	epreng   [i]=actnode->sol_increment.a.da[ipos->velnp][2];
	epren    [i]=actnode->sol_increment.a.da[ipos->veln ][2];
    } /* end of loop over nodes of element */

    if(ele->g.gsurf->neum!=NULL)
    {

	if (ele->g.gsurf->neum->curve<1)
	{
	    acttimefac =ONE;
	    acttimefacn=ONE;
	}
	else
	{
	    dyn_facfromcurve(ele->g.gsurf->neum->curve-1,
			     fdyn->acttime,
			     &acttimefac) ;

	    dyn_facfromcurve(ele->g.gsurf->neum->curve-1,
			     fdyn->acttime-fdyn->dta,
			     &acttimefacn) ;
	}

	for (i=0;i<2;i++)
	{
	    if (ele->g.gsurf->neum->neum_onoff.a.iv[i]==0)
	    {
		edeadn[i]  = 0.0;
		edeadng[i] = 0.0;
	    }
	    if (ele->g.gsurf->neum->neum_type==neum_dead  &&
		ele->g.gsurf->neum->neum_onoff.a.iv[i]!=0)
	    {
		edeadn [i] = ele->g.gsurf->neum->neum_val.a.dv[i]*acttimefacn;
		edeadng[i] = ele->g.gsurf->neum->neum_val.a.dv[i]*acttimefac ;
	    }
	}
    }
    else
    {
	for (i=0;i<2;i++)
	{
	    edeadn [i] = 0.0;
	    edeadng[i] = 0.0;
	}
    }

    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    f2_get_time_dependent_sub_tau(ele,xyze,funct,deriv,evelng,eveln,visc);

    facM   =         1.0/(fdyn->tau[0]+theta*dt);
    facMtau=fdyn->tau[0]/(fdyn->tau[0]+theta*dt);


    old_facM   = 1./(fdyn->tau_old[0]+theta*dt);
    old_facMtau= fdyn->tau_old[0]*old_facM;


/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
    for (lr=0;lr<nir;lr++)
    {
	for (ls=0;ls<nis;ls++)
	{
	    /*------- get values of  shape functions and their derivatives ---*/
	    switch(ele->distyp)
	    {
		case quad4: case quad8: case quad9:   /* --> quad - element */
		    e1   = data->qxg[lr][nir-1];
		    e2   = data->qxg[ls][nis-1];
		    f2_rec(funct,deriv,deriv2,e1,e2,ele->distyp,icode);
		    break;
		case tri3: case tri6:   /* --> tri - element */
		    e1   = data->txgr[lr][intc];
		    e2   = data->txgs[lr][intc];
		    f2_tri(funct,deriv,deriv2,e1,e2,ele->distyp,icode);
		    break;
		default:
		    dserror("typ unknown!");
	    } /* end switch(typ) */

	    /*-------------- get velocities (n) at integration point ---*/
	    f2_veci(velint_old,funct,eveln ,ele->numnp);
	    /*------------ get velocities (n+1) at integration point ---*/
	    f2_veci(velint    ,funct,evelng,ele->numnp);
	    /*-------------------- get pressure at integration point ---*/
	    press_old = 0;
	    press = 0;
	    for (i=0;i<ele->numnp;i++)
	    {
		actnode=ele->node[i];

		/*--------------------------------- get pressure (n) ---*/
		press_old +=funct[i]*epren[i];
		/*------------------------------- get pressure (n+1) ---*/
		press     +=funct[i]*epreng[i];
	    }

	    /*------------------ compute Jacobian matrix at time n+1 ---*/
	    f2_jaco(xyze,deriv,xjm,&det,ele->numnp,ele);

	    /*----------------------------- compute global derivates ---*/
	    f2_gder(derxy,deriv,xjm,det,ele->numnp);

	    /*--- get velocity (n+1,i) derivatives at integration point */
	    f2_vder(vderxy,derxy,evelng,ele->numnp);

	    if (ihoel!=0)
	    {
		f2_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,ele->numnp);
		f2_vder2(vderxy2,derxy2,evelng,ele->numnp);
	    }

	    /*------- get velocity (n) derivatives at integration point */
	    f2_vder(vderxy_old,derxy,eveln,ele->numnp);

	    if (ihoel!=0)
	    {
		f2_vder2(vderxy2_old,derxy2,eveln,ele->numnp);
	    }

	    /*------------------------------- get pressure gradients ---*/
	    gradp[0] = gradp[1] = 0.0;

	    for (i=0; i<ele->numnp; i++)
	    {
		gradp[0] += derxy[0][i] * epreng[i];
		gradp[1] += derxy[1][i] * epreng[i];
	    }

	    gradp_old[0] = gradp_old[1] = 0.0;

	    for (i=0; i<ele->numnp; i++)
	    {
		gradp_old[0] += derxy[0][i] * epren[i];
		gradp_old[1] += derxy[1][i] * epren[i];
	    }
	    /*------------------------------ get higher order terms ---*/
	    if(ihoel!=0)
	    {
		hot    [0]=0.5 * (2.0*vderxy2[0][0]
				  +
				  (vderxy2[0][1] + vderxy2[1][2]));
		hot    [1]=0.5 * (2.0*vderxy2[1][1]
				  +
				  (vderxy2[1][0] + vderxy2[0][2]));

		hot_old[0]=0.5 * (2.0*vderxy2_old[0][0]
			   	  +
			   	  (vderxy2_old[0][1] + vderxy2_old[1][2]));
		hot_old[1]=0.5 * (2.0*vderxy2_old[1][1]
			   	  +
			   	  (vderxy2_old[1][0] + vderxy2_old[0][2]));
	    }
	    else
	    {
		hot    [0]=0;
		hot    [1]=0;
		hot_old[0]=0;
		hot_old[1]=0;
	    }

	    for (dim=0;dim<2;dim++)
	    {
		/* the current subscales velocities become the most
		 * recent subscale velocities for the next timestep    */
		/* sv_old is just an abbreviation                      */
		sv_old[dim]=ele->e.f2->sub_vel.a.da[dim][lr*nis+ls];

		/* calculate new residual without time derivative       */

		res_new[dim]=0;
		res_new[dim]+=(velint[0]*vderxy[dim][0]
			       +
			       velint[1]*vderxy[dim][1]);
		res_new[dim]-=2*visc*hot[dim];
		res_new[dim]+=gradp[dim];
		res_new[dim]-=edeadng[dim];

		/* calculate old residual without time derivative       */
		res_old[dim]=0;
		res_old[dim]+=(velint_old[0]*vderxy_old[dim][0]
			       +
			       velint_old[1]*vderxy_old[dim][1]);

		res_old[dim]-=2*visc*hot_old[dim];
		res_old[dim]+=gradp_old[dim];
		res_old[dim]-=edeadn[dim];

		/* calculate the time derivative                        */
		time_der[dim]=velint[dim]-velint_old[dim];


		/* set new subscale velocities                          */
		ele->e.f2->sub_vel.a.da[dim][lr*nis+ls]=
		     facMtau*                sv_old[dim]
		    -facMtau*(               time_der[dim]
			      +theta    *dt* res_new[dim]
			      +(1-theta)*dt* res_old[dim])
		    -facMtau/fdyn->tau_old[0] * (1-theta)*dt* sv_old[dim];

	    }
	} /* end of loop over integration points ls*/
    } /* end of loop over integration points lr */
}

/* clean up ------------------------------------------------------------*/
amdel(&xyze_a);
amdel(&xjm_a);
amdel(&funct_a);
amdel(&deriv_a);
amdel(&eveln_a);
amdel(&evelng_a);
amdel(&derxy_a);
amdel(&vderxy_a);
amdel(&vderxy_old_a);
amdel(&epren_a);
amdel(&epreng_a);
amdel(&edeadn_a);
amdel(&edeadng_a);
amdel(&w1_a);
amdel(&w2_a);
amdel(&derxy2_a);
amdel(&deriv2_a);
amdel(&vderxy2_a);
amdel(&vderxy2_old_a);


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
}


/*!---------------------------------------------------------------------
\brief update of time dependent subscales

<pre>                                                        gammi 11/06

The time dependent pressure subscales are updated according to a
generalzed alpha timestepping scheme.

                  dp_sub       1
                  ------ = - ----- * p_sub + res_C(u)
                    dt       tau_C

Here, res_C(u)=div(u) is the residual of the continuity equation.
The time discrete expression is



     p_sub^{n+1} =

               /dp_sub\ ^{n}
      = C1 *  |------ |     + C2 * p_sub^{n} - C3 * div(u^{n+alpha_F})
       	       \  dt  /


where Ci are constants from the time integration algorithm

   	               (alpha_M-theta)*(dt*tau_C)
   	       C1 =  ------------------------------
   	             alpha_M*tau_C+alpha_F*dt*theta



   	             (alpha_F-1)*dt*theta-alpha_M*tau_C
   	       C2 =  ----------------------------------
   	               alpha_M*tau_C+alpha_F*dt*theta


   	                   theta*dt*tau_C
   	       C3 =  ------------------------------
   	             alpha_M*tau_C+alpha_F*dt*theta


The acceleration of the subscale pressure is updated by use of the value
p_sub^{n+1} according to the linear relation between accelerations and
function values in the generalised alpha scheme.


</pre>
\param  *actpart       PARTITION        (i)
\param  *actintra      INTRA            (i)
\param  *actfield      FIELD            (i)
\param  *ipos          ARRAY_POSITION   (i)
\param   disnum_calc   INT              (i)
\return void

------------------------------------------------------------------------*/


void f2_update_subscale_pres_for_inc_gen_alpha(
    PARTITION      *actpart,
    INTRA          *actintra,
    FIELD          *actfield,
    ARRAY_POSITION *ipos,
    INT             disnum_calc)
{
/* multi purpose counter */
int      i;

/* counter for elements */
int      nele;

/* element related data */
int      ls,lr;
int      nis=0,nir=0;
int      icode=0,intc=0;
int      ihoel;
double   det;

DOUBLE    e1,e2;      /* natural coordinates of integr. point           */


ELEMENT *ele;
NODE    *actnode;

/* material properties */
double   visc;

/* new and old residual of the coninuity equation */
double   divu;

/* the old subscale pressure */
double   sp_old;
double   sp_acc_old;

double tau_C;

double  theta,alpha_F,alpha_M;

double  dt;

double    aftdt;

double  C1,C2,C3;


FLUID_DYNAMIC   *fdyn;
FLUID_DATA      *data;

static ARRAY     xyze_a;   /* actual element coordinates                */
static DOUBLE  **xyze;
static ARRAY     xjm_a;    /* jocobian matrix                           */
static DOUBLE  **xjm;
static ARRAY     funct_a;  /* shape functions                           */
static DOUBLE   *funct;
static ARRAY     deriv_a;  /* first natural derivatives                 */
static DOUBLE  **deriv;
static ARRAY     eveln_a;  /* element velocities at (n)                 */
static DOUBLE  **eveln;
static ARRAY     evelng_a; /* element velocities at (n+gamma)           */
static DOUBLE  **evelng;
static ARRAY     derxy_a;  /* coordinate - derivatives                  */
static DOUBLE  **derxy;
static ARRAY     vderxy_a; /* vel - derivatives                         */
static DOUBLE  **vderxy;
static ARRAY     deriv2_a; /* second natural derivatives                */
static DOUBLE  **deriv2;



#ifdef DEBUG
dstrc_enter("f2_update_subscale_pres_for_inc_gen_alpha");
#endif


xyze      = amdef("xyze"     ,&xyze_a     ,2,MAXNOD_F2,"DA");
xjm       = amdef("xjm"      ,&xjm_a      ,2,2        ,"DA");
funct     = amdef("funct"    ,&funct_a    ,MAXNOD_F2,1,"DV");
deriv     = amdef("deriv"    ,&deriv_a    ,2,MAXNOD_F2,"DA");
eveln     = amdef("eveln"    ,&eveln_a    ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
evelng    = amdef("evelng"   ,&evelng_a   ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
vderxy    = amdef("vderxy"   ,&vderxy_a   ,2,2,"DA");
derxy     = amdef("derxy"    ,&derxy_a    ,2,MAXNOD_F2,"DA");
deriv2    = amdef("deriv2"   ,&deriv2_a   ,3,MAXNOD_F2,"DA");


fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;


theta  = fdyn->theta;
alpha_F= fdyn->alpha_f;
alpha_M= fdyn->alpha_m;

dt     = fdyn->dt;
aftdt  = alpha_F*dt*theta;

for(nele=0;nele<actpart->pdis[disnum_calc].numele;nele++)
{
    ele=actpart->pdis[disnum_calc].element[nele];

    /*------------------------------------------ set element coordinates -*/
    for(i=0;i<ele->numnp;i++)
    {
	xyze[0][i]=ele->node[i]->x[0];
	xyze[1][i]=ele->node[i]->x[1];
    }

    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*------- get integraton data and check if elements are "higher order" */
    switch (ele->distyp)
    {
	case quad4: case quad8: case quad9:  /* --> quad - element */
	    icode   = 3; /* flag for higher order elements                 */
	    ihoel   = 1; /* flag for eveluation of shape functions         */
	    nir = ele->e.f2->nGP[0];
	    nis = ele->e.f2->nGP[1];
	    intc = ele->e.f2->nGP[1];
	    break;
	case tri6: /* --> tri - element */
	    icode   = 3; /* flag for higher order elements                 */
	    ihoel   = 1; /* flag for eveluation of shape functions         */
	    nir  = ele->e.f2->nGP[0];
	    nis  = 1;
	    intc = ele->e.f2->nGP[1];
	    break;
	case tri3:
	    ihoel  =0;  /* flag for higher order elements                 */
	    icode  =2;  /* flag for eveluation of shape functions         */
	    nir  = ele->e.f2->nGP[0];
	    nis  = 1;
	    intc = ele->e.f2->nGP[1];
	    break;
	default:
	    dserror("typ unknown!");
    } /* end switch(typ) */

    /* -> implicit time integration method ---------*/
    for(i=0;i<ele->numnp;i++) /* loop nodes of element */
    {
	actnode=ele->node[i];
        /*--------------------------- set element velocities at n+alpha_F */
	evelng[0][i]=actnode->sol_increment.a.da[ipos->velnm][0];
	evelng[1][i]=actnode->sol_increment.a.da[ipos->velnm][1];

    } /* end of loop over nodes of element */


    /*--------------------------------------------- stab-parameter ---*/
    f2_get_time_dependent_sub_tau(ele,xyze,funct,deriv,evelng,NULL,visc);

    tau_C  = fdyn->tau[2];

    C1 = (alpha_M-theta)*dt*(tau_C/(alpha_M*tau_C+aftdt));
    C2 = ((alpha_F-1.0)*dt*theta+alpha_M*tau_C)/(alpha_M*tau_C+aftdt);
    C3 = (tau_C/(alpha_M*tau_C+aftdt))*theta*dt;

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
    for (lr=0;lr<nir;lr++)
    {
	for (ls=0;ls<nis;ls++)
	{
	    /*------- get integraton data and check if elements are
	                                                 "higher order" */
	    switch (ele->distyp)
	    {
		case quad4: case quad8: case quad9:  /* --> quad - element */
		    e1   = data->qxg[lr][nir-1];
		    e2   = data->qxg[ls][nis-1];
		    f2_rec(funct,deriv,deriv2,e1,e2,ele->distyp,icode);
		    break;
		case tri6: /* --> tri - element */
		    e1   = data->txgr[lr][intc];
		    e2   = data->txgs[lr][intc];
		    f2_tri(funct,deriv,deriv2,e1,e2,ele->distyp,icode);
		    break;
		case tri3:
		    e1   = data->txgr[lr][intc];
		    e2   = data->txgs[lr][intc];
		    f2_tri(funct,deriv,deriv2,e1,e2,ele->distyp,icode);
		    break;
		default:
		    dserror("typ unknown!");
	    } /* end switch(typ) */

	    /*------------------ compute Jacobian matrix at time n+1 ---*/
	    f2_jaco(xyze,deriv,xjm,&det,ele->numnp,ele);

	    /*----------------------------- compute global derivates ---*/
	    f2_gder(derxy,deriv,xjm,det,ele->numnp);

	    /*--- get velocity (n+alpha_F,i) derivatives at integration
	                                                          point */
	    f2_vder(vderxy,derxy,evelng,ele->numnp);

	    divu     = vderxy[0][0] + vderxy[1][1];

	    sp_old    =ele->e.f2->sub_pres.a.dv[lr*nis+ls];

	    sp_acc_old=ele->e.f2->sub_pres_acc.a.dv[lr*nis+ls];

	    ele->e.f2->sub_pres_trial.a.dv[lr*nis+ls]=
		C1*sp_acc_old +	C2*sp_old - C3*divu;

	    ele->e.f2->sub_pres_acc_trial.a.dv[lr*nis+ls]=
		(theta-1.)/theta*sp_acc_old
		+(ele->e.f2->sub_pres_trial.a.dv[lr*nis+ls]-sp_old)/(theta*dt);

	} /* end of loop over integration points ls*/
    } /* end of loop over integration points lr */
} /* end of loop over elements nele */

/* clean up ------------------------------------------------------------*/

amdel(&xyze_a);
amdel(&xjm_a);
amdel(&funct_a);
amdel(&deriv_a);
amdel(&eveln_a);
amdel(&evelng_a);
amdel(&derxy_a);
amdel(&vderxy_a);
amdel(&deriv2_a);


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;

}



/*!---------------------------------------------------------------------
\brief update of time dependent subscales

<pre>                                                        gammi 11/06

The time dependent pressure subscales are updated according to a
generalzed alpha timestepping scheme.

                  du_sub       1
                  ------ = - ----- * u_sub + res_M(u,p)
                    dt       tau_M

Here, res_M(u,p) is the residual of the momentum equation.
The time discrete expression is




The subscale velocities are updated by use of the value of the subscale
acceleration according to the linear relation between accelerations and
function values in the generalised alpha scheme.


</pre>
\param  *actpart       PARTITION        (i)
\param  *actintra      INTRA            (i)
\param  *actfield      FIELD            (i)
\param  *ipos          ARRAY_POSITION   (i)
\param   disnum_calc   INT              (i)
\return void

------------------------------------------------------------------------*/


void f2_update_subscale_vel_for_inc_gen_alpha(
    PARTITION      *actpart,
    INTRA          *actintra,
    FIELD          *actfield,
    ARRAY_POSITION *ipos,
    INT             disnum_calc)
{
/* multi purpose counter */
int      i;

/* counter for elements */
int      nele;

/* element related data */
int      ls,lr;
int      nis=0,nir=0;
int      icode=0,intc=0;
int      ihoel=0;
double   det;

DOUBLE    e1,e2;      /* natural coordinates of integr. point           */


ELEMENT *ele;

/* material properties */
double   visc;


/* the old subscale velocities and accelerations */
double   sv_old    [2];
double   sv_acc_old[2];

/* the new subscale velocities and accelerations */
double   sv_new        [2];
double   sv_acc_new    [2];
double   sv_acc_mod_new[2];

/* temporary variables */
double   res_mod[2]; /* the residual of the momentum equation without f */

/* higher order terms */
DOUBLE  hot   [2];

/* pressure gradient */
DOUBLE  gradp [2];


/* intermediate acceleration (n+alpha_M) and velocities (n+alpha_F)    */
double   accint[2],velint[2];

double   tau_M;

double   theta,alpha_F,alpha_M;

double   dt;

double   aftdt;
double   amtau_M;


/* constants containing dt and the stabilisation parameter */
double  facM;


FLUID_DYNAMIC   *fdyn;
FLUID_DATA      *data;

static ARRAY     xyze_a;   /* actual element coordinates                */
static DOUBLE  **xyze;
static ARRAY     xjm_a;    /* jocobian matrix                           */
static DOUBLE  **xjm;
static ARRAY     funct_a;  /* shape functions                           */
static DOUBLE   *funct;
static ARRAY     deriv_a;  /* first natural derivatives                 */
static DOUBLE  **deriv;
static ARRAY     eveln_a;  /* element velocities at (n)                 */
static DOUBLE  **eveln;
static ARRAY     evelng_a; /* element velocities at (n+gamma)           */
static DOUBLE  **evelng;
static ARRAY     eaccng_a;  /* element accelerations at (n+alpha_M)     */
static DOUBLE  **eaccng;
static ARRAY     derxy_a;  /* coordinate - derivatives                  */
static DOUBLE  **derxy;
static ARRAY     vderxy_a; /* vel - derivatives                         */
static DOUBLE  **vderxy;
static ARRAY     deriv2_a; /* second natural derivatives                */
static DOUBLE  **deriv2;
static ARRAY     w1_a;     /* working array of arbitrary chosen size    */
static DOUBLE  **wa1;      /* used in different element routines        */
static ARRAY     w2_a;     /* working array of arbitrary chosen size    */
static DOUBLE  **wa2;      /* used in different element routines        */
static ARRAY     derxy2_a; /* 2nd coordinate - derivatives              */
static DOUBLE  **derxy2;
static ARRAY     vderxy2_a;/* vel - 2nd derivatives                     */
static DOUBLE  **vderxy2;
static ARRAY     epreng_a; /* element pressures at (n)	                */
static DOUBLE   *epreng;
static DOUBLE   *edeadng;
static ARRAY     edeadng_a;/* element dead load (selfweight)            */



#ifdef DEBUG
dstrc_enter("f2_update_subscale_vel_for_inc_gen_alpha");
#endif


xyze      = amdef("xyze"     ,&xyze_a     ,2,MAXNOD_F2,"DA");
xjm       = amdef("xjm"      ,&xjm_a      ,2,2        ,"DA");
funct     = amdef("funct"    ,&funct_a    ,MAXNOD_F2,1,"DV");
deriv     = amdef("deriv"    ,&deriv_a    ,2,MAXNOD_F2,"DA");
eveln     = amdef("eveln"    ,&eveln_a    ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
evelng    = amdef("evelng"   ,&evelng_a   ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
eaccng    = amdef("eaccng"   ,&eaccng_a   ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
vderxy    = amdef("vderxy"   ,&vderxy_a   ,2,2,"DA");
derxy     = amdef("derxy"    ,&derxy_a    ,2,MAXNOD_F2,"DA");
deriv2    = amdef("deriv2"   ,&deriv2_a   ,3,MAXNOD_F2,"DA");
wa1       = amdef("wa1"      ,&w1_a       ,MAXDOFPERELE,MAXDOFPERELE,"DA");
wa2       = amdef("wa2"      ,&w2_a       ,MAXDOFPERELE,MAXDOFPERELE,"DA");
derxy2    = amdef("derxy2"   ,&derxy2_a   ,3,MAXNOD_F2,"DA");
vderxy2   = amdef("vderxy2"  ,&vderxy2_a  ,2,3,"DA");
epreng    = amdef("epreng"   ,&epreng_a   ,MAXNOD_F2,1,"DV");
edeadng   = amdef("edeadng"  ,&edeadng_a  ,2,1,"DV");


fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;


theta  = fdyn->theta;
alpha_F= fdyn->alpha_f;
alpha_M= fdyn->alpha_m;

dt     = fdyn->dt;
aftdt  = alpha_F*theta*dt;

for(nele=0;nele<actpart->pdis[disnum_calc].numele;nele++)
{
    ele=actpart->pdis[disnum_calc].element[nele];

    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

/*------- get integraton data and check if elements are "higher order" */
    switch (ele->distyp)
    {
	case quad4: case quad8: case quad9:  /* --> quad - element */
	    icode   = 3; /* flag for higher order elements                 */
	    ihoel   = 1; /* flag for eveluation of shape functions         */
	    nir = ele->e.f2->nGP[0];
	    nis = ele->e.f2->nGP[1];
	    intc = ele->e.f2->nGP[1];
	    break;
	case tri6: /* --> tri - element */
	    icode   = 3; /* flag for higher order elements                 */
	    ihoel   = 1; /* flag for eveluation of shape functions         */
	    nir  = ele->e.f2->nGP[0];
	    nis  = 1;
	    intc = ele->e.f2->nGP[1];
	    break;
	case tri3:
	    ihoel  =0;  /* flag for higher order elements                 */
	    icode  =2;  /* flag for eveluation of shape functions         */
	    nir  = ele->e.f2->nGP[0];
	    nis  = 1;
	    intc = ele->e.f2->nGP[1];
	    break;
	default:
	    dserror("typ unknown!");
    } /* end switch(typ) */



    /*------------------------------------------ set element coordinates -*/
    f2_inc_gen_alpha_calset(
	ele,
	xyze,
	eaccng,
	evelng,
	epreng,
	edeadng,
	ipos,
	&visc
	);


    /*--------------------------------------------- stab-parameter ---*/
    f2_get_time_dependent_sub_tau(ele,xyze,funct,deriv,evelng,NULL,visc);

    tau_M  = fdyn->tau[0];

    amtau_M= alpha_M*tau_M;
    facM   = 1./(amtau_M+aftdt);


/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
    for (lr=0;lr<nir;lr++)
    {
	for (ls=0;ls<nis;ls++)
	{
            /*------- get integraton data and check if elements are
	                                                 "higher order" */
	    switch (ele->distyp)
	    {
		case quad4: case quad8: case quad9:  /* --> quad - element */
		    e1   = data->qxg[lr][nir-1];
		    e2   = data->qxg[ls][nis-1];
		    f2_rec(funct,deriv,deriv2,e1,e2,ele->distyp,icode);
		    break;
		case tri6: /* --> tri - element */
		    e1   = data->txgr[lr][intc];
		    e2   = data->txgs[lr][intc];
		    f2_tri(funct,deriv,deriv2,e1,e2,ele->distyp,icode);
		    break;
		case tri3:
		    e1   = data->txgr[lr][intc];
		    e2   = data->txgs[lr][intc];
		    f2_tri(funct,deriv,deriv2,e1,e2,ele->distyp,icode);
		    break;
		default:
		    dserror("typ unknown!");
	    } /* end switch(typ) */

	    /*------------------ compute Jacobian matrix at time n+1 ---*/
	    f2_jaco(xyze,deriv,xjm,&det,ele->numnp,ele);

	    /*----------------------------- compute global derivates ---*/
	    f2_gder(derxy,deriv,xjm,det,ele->numnp);

	    /*------ get velocities (n+alpha_F) at integration point ---*/
	    f2_veci(velint,funct,evelng,ele->numnp);

	    /*--- get accelerations (n+alpha_M) at integration point ---*/
	    f2_veci(accint,funct,eaccng,ele->numnp);

	    /*--- get velocity (n+alpha_F,i) derivatives at integration
	                                                          point */
	    f2_vder(vderxy,derxy,evelng,ele->numnp);

	    if (ihoel!=0)
	    {
		f2_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,ele->numnp);
		f2_vder2(vderxy2,derxy2,evelng,ele->numnp);
	    }

	    for(i=0;i<2;i++)
	    {
		sv_old    [i]=ele->e.f2->sub_vel.a.da    [i][lr*nis+ls];

		sv_acc_old[i]=ele->e.f2->sub_vel_acc.a.da[i][lr*nis+ls]
		    +1./alpha_M*edeadng[i];
	    }

	    /*------------------------------ get higher order terms ---*/
	    if(ihoel!=0)
	    {
		hot    [0]=0.5 * (2.0*vderxy2[0][0]
				  +
				  (vderxy2[0][1] + vderxy2[1][2]));
		hot    [1]=0.5 * (2.0*vderxy2[1][1]
				  +
				  (vderxy2[1][0] + vderxy2[0][2]));
	    }
	    else
	    {
		hot    [0]=0;
		hot    [1]=0;
	    }

	    /*------------------------------- get pressure gradients ---*/
	    gradp[0] = gradp[1] = 0.0;

	    for (i=0; i<ele->numnp; i++)
	    {
		gradp[0] += derxy[0][i] * epreng[i];
		gradp[1] += derxy[1][i] * epreng[i];
	    }

	    for(i=0;i<2;i++)
	    {
		/* residual without body force */
		res_mod[i] = accint[i] + 0 - 2*visc*hot[i]
		    + gradp[i];


		sv_acc_mod_new[i] = - facM * sv_old[i]
		    - facM * (tau_M+alpha_F*dt-amtau_M-aftdt)*sv_acc_old[i]
		    - facM * (tau_M*res_mod[i]+aftdt/alpha_M*edeadng[i]);

		ele->e.f2->sub_vel_acc_trial.a.da[i][lr*nis+ls]
		    =sv_acc_mod_new[i];

		sv_acc_new[i] = sv_acc_mod_new[i]+ 1./alpha_M*edeadng[i];


		sv_new[i] = sv_old[i] + dt*sv_acc_old[i]
		    + dt*theta *(sv_acc_new[i]-sv_acc_old[i]);

		ele->e.f2->sub_vel_trial.a.da[i][lr*nis+ls]=sv_new [i];

	    }
	} /* end of loop over integration points ls*/
    } /* end of loop over integration points lr */
} /* end of loop over elements nele */

/* clean up ------------------------------------------------------------*/

amdel(&xyze_a);
amdel(&xjm_a);
amdel(&funct_a);
amdel(&deriv_a);
amdel(&eveln_a);
amdel(&evelng_a);
amdel(&eaccng_a);
amdel(&derxy_a);
amdel(&vderxy_a);
amdel(&deriv2_a);
amdel(&w1_a);
amdel(&w2_a);
amdel(&derxy2_a);
amdel(&vderxy2_a);
amdel(&epreng_a);
amdel(&edeadng_a);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;

}


/*!---------------------------------------------------------------------
\brief update of time dependent subscales

<pre>                                                        gammi 11/06

converged time dependent subscale trial values are copied to the new
subscale values
</pre>
\param  *actpart       PARTITION        (i)
\param  *actintra      INTRA            (i)
\param  *actfield      FIELD            (i)
\param  *ipos          ARRAY_POSITION   (i)
\param   disnum_calc   INT              (i)
\return void

------------------------------------------------------------------------*/


void f2_time_update_subscales_for_incr_gen_alpha (
    PARTITION      *actpart,
    INTRA          *actintra,
    FIELD          *actfield,
    ARRAY_POSITION *ipos,
    INT             disnum_calc)
{
/* dimension counter */
  int      dim;

/* counter for elements */
  int      nele;

/* element related data */
  ELEMENT *ele;
  int      ls,lr;
  int      nis=0,nir=0;


#ifdef DEBUG
dstrc_enter("f2_time_update_subscales_for_incr_gen_alpha");
#endif


for(nele=0;nele<actpart->pdis[disnum_calc].numele;nele++)
{
    ele=actpart->pdis[disnum_calc].element[nele];

/*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
    for (lr=0;lr<nir;lr++)
    {
      for (ls=0;ls<nis;ls++)
      {

        ele->e.f2->sub_pres.a.dv[lr*nis+ls]
          =ele->e.f2->sub_pres_trial.a.dv[lr*nis+ls];

        ele->e.f2->sub_pres_acc.a.dv[lr*nis+ls]
          =ele->e.f2->sub_pres_acc_trial.a.dv[lr*nis+ls];

        for(dim=0;dim<2;dim++)
        {
          ele->e.f2->sub_vel.a.da[dim][lr*nis+ls]
            =ele->e.f2->sub_vel_trial.a.da[dim][lr*nis+ls];

          ele->e.f2->sub_vel_acc.a.da[dim][lr*nis+ls]=
            ele->e.f2->sub_vel_acc_trial.a.da[dim][lr*nis+ls];
        }

      } /* end of loop over integration points ls*/
    } /* end of loop over integration points lr */
} /* end of loop over elements nele */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;

}

/*!---------------------------------------------------------------------
\brief update of time dependent subscales for inc acc genalpha

<pre>                                                        gammi 04/07

time dependent subscale trial values are calculated in one gausspoint of
one element. This function is called during the loop over the
gausspoints before calculating the element matrix and the element force
vector!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The time dependent pressure subscales are updated according to a
one step theta timestepping scheme.

                  dp_sub       1
                  ------ = - ----- * p_sub + res_C(u)
                    dt       tau_C

Here, res_C(u)=div(u) is the residual of the continuity equation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The time dependent velocity subscales are updated according to a
one step theta time integration scheme of the equation

                  du_sub       1
                  ------ = - ----- * u_sub + res_M(u)
                    dt       tau_M

Here, res_M(u) is the residual of the momentum equation and contains a
time derivative, a convective term, a diffusion term, the pressure
gradient and the volume force.


</pre>
\param  sp_trial              *double             (o)
\param  sp_acc_trial          *double             (o)
\param  sub_u_trial           *double             (o)
\param  sub_u_acc_trial       *double             (o)
\param  sp_old                 double[2]          (i)
\param  sp_acc_old             double[2]          (i)
\param  sub_u                  double[2]          (i)
\param  sub_u_acc              double[2]          (i)
\param  ihoel                  int                (i)
\param  alpha_M                double             (i)
\param  alpha_F                double             (i)
\param  theta                  double             (i)
\param  dt                     double             (i)
\param  tau_M                  double             (i)
\param  tau_C                  double             (i)
\param  visc                   double             (i)
\param  velint                 double[2]          (i)
\param  accint                 double[2]          (i)
\param  gradp                  double[2]          (i)
\param  vderxy                 double**           (i)
\param  vderxy2                double**           (i)
\param  edeadng                double*            (i)

\return void

------------------------------------------------------------------------*/
void f2_up_tds_at_gp_genalpha (
  double  *sp_trial       , /* trial value for subscale pressure        */
  double  *sp_acc_trial   , /* trial value for acc. of subscale p       */
  double  *su_trial       , /* trial value for subscale velocity        */
  double  *su_acc_mod     , /* subscale u acceleration wo bodyforce     */
  double   sp_old         , /* old subscale pressure                    */
  double   sp_acc_old     , /* old acceleration of subscale pressure    */
  double   su_old      [2], /* old subscale velocity                    */
  double   su_acc_old  [2], /* old acceleration of subscale velocity    */
  int      ihoel        ,/* flag for higher order elements              */
  double   alpha_M      ,/* momentum parameter of genalpha              */
  double   alpha_F      ,/* force parameter of genalpha                 */
  double   theta        ,/* parameter relating acc and vel of genalpha  */
  double   dt           ,/* timestepsize                                */
  double   tau_M        ,/* stabilisation parameter, velocity           */
  double   tau_C        ,/* stabilisation parameter, continuity         */
  double   visc         ,/* kinematic viscosity                         */
  double   velint [2]   ,/* velocity at t=t^{n+alpha_F}                 */
  double   accint [2]   ,/* acceleration at t=t^{n+alpha_M}             */
  double   gradp  [2]   ,/* pressure gradient at t=t^{n+1}              */
  double **vderxy       ,/* velocity derivatives at t=t^{n+alpha_F}     */
  double **vderxy2      ,/* second der. of velocity at t=t^{n+alpha_F}  */
  double  *edeadng       /* body force at time t=???                    */
  )
{

  int       dim;         /* counter for dimensions (x,y) */

  double    aftdt;
  double    amtauC;

  double    facM;
  double    amtauM;


  double    C1,C2,C3;

  double    hot         [2]; /* higher order terms                     */
  double    res_mod     [2]; /* residual of momentum equation without
                                                            body force */
#ifdef DEBUG
  dstrc_enter("f2_up_tds_at_gp_genalpha");
#endif

  aftdt   = alpha_F*dt*theta;


  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
   *                                                                   *
   *                  update of SUBSCALE PRESSURE                      *
   *                                                                   *
   *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

  amtauC  = alpha_M*tau_C;

  C1 = (alpha_M-theta)*dt*(tau_C/(amtauC+aftdt));
  C2 = (dt*theta)/(amtauC+aftdt);
  C3 = dt*tau_C*theta/(amtauC+aftdt);


/*  printf("C1 %22.15e C2 %22.15e  C3  %22.15e\n",C1,C2,C3);*/


  sp_trial[0]  = C1*sp_acc_old                   ;
  sp_trial[0] -= C2*sp_old                       ;
  sp_trial[0] -= C3*(vderxy[0][0] + vderxy[1][1]);
  sp_trial[0] += sp_old;

  sp_acc_trial[0] = (1.0-1./theta)*sp_acc_old;
  sp_acc_trial[0]+= (sp_trial[0]-sp_old)/(theta*dt);

  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
   *                                                                   *
   *                  update of SUBSCALE VELOCITY                      *
   *                                                                   *
   *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

  amtauM   = alpha_M*tau_M;
  facM     = 1./(amtauM+aftdt);

  /*------------------------------ get higher order terms ---*/
  if(ihoel!=0)
  {
    hot    [0]=0.5 * (2.0*vderxy2[0][0]
                      +
                      (vderxy2[0][1] + vderxy2[1][2]));
    hot    [1]=0.5 * (2.0*vderxy2[1][1]
                      +
                      (vderxy2[1][0] + vderxy2[0][2]));
  }
  else
  {
    hot    [0]=0;
    hot    [1]=0;
  }

  for(dim=0;dim<2;dim++)
  {
    /* residual without body force */
    res_mod[dim] =
      accint[dim]
      + 0
      - 2*visc*hot[dim]
      + gradp[dim];

    /* this is a modified subscale acceleration --- it doesn't contain
     * the full contribution of the body force!!                       */
    su_acc_mod[dim] =
      su_old[dim]
      + (tau_M*(1-alpha_M)+alpha_F*dt*(1-theta))*su_acc_old[dim]
      + (tau_M*res_mod[dim]+aftdt/alpha_M*edeadng[dim]);
    su_acc_mod[dim] *=-facM;

    if(0)
    {
      printf("res_mod              [%d]               %22.15e\n",
             dim,res_mod[dim]);
      printf("su_acc_old           [%d]               %22.15e\n",
             dim,su_acc_old [dim]);
      printf("aftdt/alpha_M*edeadng[%d]               %22.15e\n",
             dim,aftdt/alpha_M*edeadng[dim]);
      printf("tau_M*(1-alpha_M)+alpha_F*dt*(1-theta) %22.15e\n",
             tau_M*(1-alpha_M)+alpha_F*dt*(1-theta));
    }
#if 0
    if(alldyn->fdyn->step == 1 || 0)
    {
    su_trial[dim] = su_old[dim]
      + dt*(1.-theta)*(su_acc_old[dim])
      + dt*    theta *(su_acc_mod[dim])
      + dt*(theta)*(edeadng[dim]/alpha_M);
    }
    else
#endif
    {
    su_trial[dim] = su_old[dim]
      + dt*(1.-theta)*(su_acc_old[dim])
      + dt*    theta *(su_acc_mod[dim])
      + dt*(edeadng[dim]/alpha_M);
    }
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

}

#endif /*D_FLUID2_TDS*/
#endif /*D_FLUID2*/
/*! @} (documentation module close)*/
#endif
