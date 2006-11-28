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

</pre>
\param  *actpart       PARTITION        (i)   
\param  *actintra      INTRA            (i)   
\param  *actfield      FIELD            (i)   
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
int      i;    
int      nele;

int      ls,lr;
int      nis=0,nir=0;
int      icode;
int      ihoel;

double   det;
double   visc;
double   divu,divu_old;
ELEMENT *ele;
NODE    *actnode;

double  facC;
double  theta;
double  dt;
double  sp_old;
double  CRHS;


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
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("f2_update_subscale_pres");
#endif

/*========================== initialisation ============================*/
fdyn = alldyn[genprob.numff].fdyn;

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
    f2_caltau(ele,xyze,funct,deriv,xjm,evelng,visc);

    facC=fdyn->tau[2]/(fdyn->tau[2]+theta*dt);
    
/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
    for (lr=0;lr<nir;lr++)
    {
	for (ls=0;ls<nis;ls++)
	{
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

	    sp_old=ele->e.f2->sub_pres.a.dv[lr+nis*ls];
	    CRHS=sp_old-(1.-theta)*dt*(1./fdyn->tau[2]*sp_old+divu_old);

	    ele->e.f2->sub_pres.a.dv[lr+nis*ls]=
		facC*(CRHS-theta*dt*divu);

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

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
}

#endif /*D_FLUID2_TDS*/
#endif /*D_FLUID2*/
/*! @} (documentation module close)*/
