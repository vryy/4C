/*!----------------------------------------------------------------------
\file
\brief service routines for fluid time algorithms

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fluid_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"
#include "../fluid2_is/fluid2_is.h"
#include "../fluid3_is/fluid3_is.h"
#include "../io/io.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;
/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum _CALC_ACTION calc_action[MAXFIELD];


/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

static FLUID_DYNAMIC *fdyn;
/*!---------------------------------------------------------------------
\brief routine to check starting algorithm

<pre>                                                         genk 04/02

this routine conrols the starting algorithms schemes. For the first 'nums'
iteration steps a different time integration scheme may be used then
during the rest of the simulation.

</pre>
\param *nfrastep 	INT           (o)  number of fract. steps
\param  init            INT           (i)  init flag

\return void
\warning this routine is not completely tested yet!

------------------------------------------------------------------------*/
void fluid_startproc(
		          INT               *nfrastep,
                          INT                init
		    )
{
static INT     iop_s;
static INT     iops_s;
static INT     itemax_s;
static DOUBLE  theta_s ;
static DOUBLE  thetas_s;

#ifdef DEBUG
dstrc_enter("fluid_startproc");
#endif

fdyn = alldyn[genprob.numff].fdyn;

if (init==1)  /* save parameter from input */
{
   iop_s    = fdyn->iop;
   iops_s   = fdyn->iops;
   itemax_s = fdyn->itemax;
   theta_s  = fdyn->theta;
   thetas_s = fdyn->thetas;
   if (theta_s < EPS5)
   {
     theta_s = 1.0; /* hier bei Gelegenheit mal noch 'ne dswarning */
   }
   if (iop_s == 7) /* BDF2 */
     theta_s = 1.0;
   goto end;
}

if (fdyn->step<=fdyn->nums) /* set parameter for starting phase */
{
   fdyn->iop    = iops_s;
   fdyn->theta  = thetas_s;
   *nfrastep    = 1;
   if (iops_s==5)
      *nfrastep = 3;
   if (iops_s==2 || iops_s==3)
      fdyn->itemax = 1;
   if (fdyn->step==fdyn->nums && iop_s==3)
       dserror("starting algo for semi-impl. two-step method not implemented yet\n");
       /* set U(n-1) in last step of start-algo */
}
else if (fdyn->step==(fdyn->nums+1)) /* set original parameter */
{
   fdyn->iop    = iop_s;
   fdyn->theta  = theta_s;
   *nfrastep    = 1;
   if (iop_s==5)
      *nfrastep = 3;
   if (iop_s==2 || iop_s==3)
      fdyn->itemax = 1;
   else
      fdyn->itemax = itemax_s;
}

/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_startproc*/


/*!---------------------------------------------------------------------
\brief calculating time independent time integration constants

<pre>                                                       chfoe 09/03

in this routine the constants for the time integration algorithms are
calculated as far as they are independent of the time

</pre>

\return void
\warning only ONE-STEP-THETA implemented up to now!

------------------------------------------------------------------------*/
void fluid_cons( void )
{

#ifdef DEBUG
dstrc_enter("fluid_cons");
#endif

/*----------------------------------------------------- initialisation */
fdyn = alldyn[genprob.numff].fdyn;
fdyn->nir=0;
fdyn->nil=0;
fdyn->nii=0;
fdyn->nis=0;
fdyn->totarea=ZERO;

/*----------------------------------------------------------------------*/
fdyn->dtp = fdyn->dt; /* set starting 'previous' step to prescribed value */
fluid_startproc(NULL,1);


/*--- the following is to distinguish beween Newton and fixed ----------
      point like iteration it has been implemented for GLS- ------------
      Stabilisation (Wall/Genkinger) only! -----------------------------
      Within USFEM a full Newton Iteration is used always!!! -----------*/
switch (fdyn->ite)
{
case 0:              /* no iteration */
   dserror("results with no nonlin. iteration not checked yet!\n");
break;
case 1:              /* fixed point like iteration */
break;
case 2:              /* Newton iteration */
   fdyn->nir = 1;
   fdyn->nii = 1;
break;
case 3:              /* fixed point iteration */
   dserror("fixed point iteration removed!!!\n");
break;
default:
   dserror("Unknown nonlinear iteration");
}

/*----------------------------------------------------- check algorithm */
switch(fdyn->iop)
{

case 1:		/* gen alpha implementation */
   if (fabs(fdyn->theta*TWO*fdyn->alpha_f/fdyn->alpha_m - ONE) > EPS10)
   {
       fdyn->theta = fdyn->alpha_m / fdyn->alpha_f * 0.5;
       printf("\nWarning: Theta, Alpha_m and Alpha_f do not satisfy 2nd order condition.\n");
       printf("         Theta is recalculated.\n");
       printf("\n Theta = Alpha_m / (2 Alpha_f) = %6.4f \n\n", fdyn->theta);
   }
   fdyn->dta = 0.0;
break;

case 4:		/* one step theta */
   fdyn->dta = 0.0;
break;

case 7:		/* 2nd order backward differencing BDF2 */
   fdyn->dta = 0.0;
   /* set number of starting steps for restarts */
   if (genprob.restart == 0) /* restart using the pss-file */
     fdyn->nums += fdyn->step;

case 8:		/* incremental accelerations gen alpha implementation */
    
  if (fabs((1./2.+fdyn->alpha_m-fdyn->alpha_f) - fdyn->theta)> EPS10)
  {
      printf("\nWarning: Theta, Alpha_m and Alpha_f do not satisfy 2nd order condition.\n");
  }

  if(fdyn->alpha_m<fdyn->alpha_f || 0.5 > fdyn->alpha_m)
  {
      dserror("\n Unstable choice of Alpha_m and Alpha_f!");
  }  
   fdyn->dta = 0.0;
break;
break;


default:
   dserror("constants for time algorithm not implemented yet!\n");
} /* end switch */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_cons*/


/*!---------------------------------------------------------------------
\brief calculating time integration constants

<pre>                                                         genk 03/02

in this routine the constants for the time integration algorithms are
calculated, here time dependent values only are calculated.

NOTE: The time discretisation scheme itself is time dependent since
      starting steps using another scheme are possible!

</pre>

\return void
\warning adaptive time stepping does not work for ALE calculations!

------------------------------------------------------------------------*/
void fluid_tcons( void )
{

#ifdef DEBUG
dstrc_enter("fluid_tcons");
#endif

fdyn = alldyn[genprob.numff].fdyn;


/*----------------------------------------------------- check algorithm */
switch (fdyn->iop)
{
case 1:		/* generalised alpha */
   if (fdyn->adaptive)
   {
      if(fdyn->dta == 0.0) fdyn->dta = fdyn->dt;
   }
   else /* constant time step size */
   {
      fdyn->dta  = fdyn->dt;
   }
   fdyn->thsl = fdyn->dta * fdyn->theta * fdyn->alpha_f / fdyn->alpha_m;
   fdyn->thpl = fdyn->thsl;
   fdyn->thsr = ZERO;
   fdyn->thpr = fdyn->thsr;
break;
case 4:		/* one step theta */
   if (fdyn->adaptive)
   {
      if(fdyn->dta == 0.0) fdyn->dta = fdyn->dt;
   }
   else /* constant time step size */
   {
      fdyn->dta  = fdyn->dt;
   }
   fdyn->thsl = fdyn->dta*fdyn->theta;
   fdyn->thpl = fdyn->thsl;
   fdyn->thsr = (ONE - fdyn->theta)*fdyn->dta;
   fdyn->thpr = fdyn->thsr;
break;
case 7:		/* 2nd order backward differencing BDF2 */
   if (fdyn->adaptive)
   {
      if(fdyn->dta == 0.0) fdyn->dta = fdyn->dt;
      fdyn->thsl = (DSQR(fdyn->dta) + fdyn->dta*fdyn->dtp)
                     / (2.0*fdyn->dta + fdyn->dtp);
      fdyn->thpl = fdyn->thsl;
      fdyn->thsr = 0.0;
      fdyn->thpr = fdyn->thsr;

   }
   else /* constant time step size */
   {
      fdyn->dta  = fdyn->dt;
      fdyn->thsl = fdyn->dta*TWO/THREE;
      fdyn->thpl = fdyn->thsl;
      fdyn->thsr = 0.0;
      fdyn->thpr = fdyn->thsr;
   }
break;
#ifdef D_FLUID2_TDS
case 8:	  /* incremental acceleration generalised alpha */
   /* constant time step size !!!*/
   fdyn->dta  = fdyn->dt;
    
   fdyn->thsl = fdyn->dta * fdyn->theta * fdyn->alpha_f;
   fdyn->thpl = fdyn->thsl;
   fdyn->thsr = ZERO;
   fdyn->thpr = fdyn->thsr;
break;
#endif
default:
   dserror("constants for time algorithm not implemented yet!\n");
}

/*----------------------------------------------- treatment of pressure */
if (fdyn->iprerhs!=1)
    dserror("treatment of pressure not implemented yet!\n");

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_tcons*/


/*!---------------------------------------------------------------------
\brief setting flags for nonlinear iteration schemes for free surface

<pre>                                                         genk 03/02

in this routine the free surface flags for the different nonlinear
iteration schemes are set. Depending on the iteration schemes the
evaluation of some termes in the LHS or the RHS have to turned on or off.

NOTE: This routine is not called any more in the Euler case! It is required
      for ALE-free-surface problems only!                    chfoe 02/05
</pre>
\param  itnum      INT  	   (i)     actual number of iterations
\return void
\warning up to now, only fixed-point like iteration checked!!!

------------------------------------------------------------------------*/
void fluid_icons( INT itnum )
{

#ifdef DEBUG
dstrc_enter("fluid_icons");
#endif

fdyn = alldyn[genprob.numff].fdyn;

fdyn->totarea=ZERO;

/*------------------------------ flags for free surface tension effects */
switch (fdyn->freesurf)
{
case 0: /* do nothing */
break;
case 1: /* explicit local langrage: include only for first iteration
           step in "time rhs" */
   if (itnum>1)
   {
      fdyn->surftens= 0;
      fdyn->fsstnif = 0;
      fdyn->fsstnii = 0;
   }
   else
   {
      fdyn->fsstnif = fdyn->surftens;
      fdyn->fsstnii = fdyn->surftens;
   }
break;
case 2: case 6:/* partioned implicit local lagrange*/
   if (itnum>1)
   {
      fdyn->fsstnif = 0;
      fdyn->fsstnii = fdyn->surftens;
   }
   else
   {
      fdyn->fsstnif = fdyn->surftens;
      fdyn->fsstnii = fdyn->surftens;
   }
break;
case 3: case 5:/* vertical heightfunction */
   /*----------------------- evaluate "time rhs" in each iteration step */
   fdyn->fsstnif  = fdyn->surftens;
   fdyn->fsstnii  = fdyn->surftens;
break;
default:
      dserror("parameter fdyn->freesurf out of range!\n");
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_icons*/

/*!---------------------------------------------------------------------
\brief initialisation of solution history

<pre>                                                         genk 04/02

in this routine the solution history of the fluid field is initialised.
The time-integration schemes require the solutions at different time-
steps. They are stored in
node->sol_incement: solution history used for calculations
      sol_increment.a.da[ipos.velnm][i]:  solution at (n-1)
      sol_increment.a.da[ipos.veln][i]:   solution at (n)
      sol_increment.a.da[ipos.hist][i]:   sol. history
      sol_increment.a.da[ipos.velnp][i]:  solution at (n+1)
      sol_increment.a.da[ipos.gridv][i]:  grid velocity
      sol_increment.a.da[ipos.convn][i]:  convective velocity at (n)
      sol_increment.a.da[ipos.convnp][i]: convective velocity at (n+1)
In node->sol one findes the converged solutions of the time-steps, which
are written to the output-file and the pss-file.

If the initial fluid fuild comes from the input-file 'fluid-start-data',
these data are copied to sol_increment.

</pre>
\param *actpart    PARTITION       (i)     partition in this field
\param *actintra   INTRA           (i)
\param *actfield   FIELD           (i)     the fluid field
\param  disnum     INT             (i)     discretization number
\param *fdyn	   FLUID_DYNAMIC   (i)
\param  numr	   INT             (i/o)   number of rows in sol_incr
\param *ipos                       (i)     node array positions
\param  str        FLUID_STRESS    (i)     flag for stress calculation
\return void

------------------------------------------------------------------------*/
void fluid_init(
    PARTITION          *actpart,
    INTRA              *actintra,
    FIELD              *actfield,
    INT                 disnum_calc,
    INT                 disnum_io,
    CALC_ACTION        *action,
    CONTAINER          *container,
    INT                 numr,
    ARRAY_POSITION     *ipos,
    FLUID_STRESS        str
    )
{
  INT    i,j;             /* simply counters                              */
  INT    actmat;          /* number of actual material                    */
  INT    numdf;           /* number of dofs in this discretisation        */
  INT    numnp_total=0;   /* total number of nodes in this discretisation */
  INT    numele_total;    /* total number of elements in this discr.      */
  INT    numnp=0;         /* number of nodes at actual element            */
  INT    ldflag;
  NODE  *actnode;         /* the actual node                              */
  ELEMENT *actele;        /* the actual element                           */
  GSURF   *actgsurf;

#ifdef D_FSI
  INT    addcol=0;
  DOUBLE dens;            /* density                                      */
  DOUBLE phi;
  NODE  *actanode;        /* the actual ale ode                           */
  GNODE *actgnode;        /* the actual gnode                             */
#endif

#ifdef D_FLUID3
  GVOL    *actgvol;
  DSURF   *actdsurf;
#endif

#ifdef D_FLUID2
  GLINE   *actgline;
  DLINE   *actdline;
#endif

  DOUBLE d,t,p,u1,u2,u3;

  /*----------------------------------------- variables for solitary wave */
#ifdef D_FSI
  DOUBLE eta,c,g,H,x,y,fac,fac1,sech;
#endif

#ifdef D_FLUID3_F
  FAST_ELES        *act_fast_eles;
  INT               l;
#endif

#ifdef D_FLUID2_TDS
  INT               nis;
#endif
  
  /* variables for beltrami */
  DOUBLE    visc,a,x1,x2,x3;

#ifdef DEBUG
  dstrc_enter("fluid_init");
#endif

  /*----------------------- set control variables for element evaluation */
  fdyn = alldyn[genprob.numff].fdyn;

  fdyn->ishape = 1;

  numdf = fdyn->numdf;
  numele_total = actfield->dis[disnum_calc].numele;

  /*---------------------------------------------------------------------*
    INITIALISE NODAL ARRAYS
   *---------------------------------------------------------------------*/
  /*-------------------------------- allocate space for solution history */
  numnp_total=actfield->dis[disnum_calc].numnp;
  for (i=0;i<numnp_total;i++)
  {
    actnode=&(actfield->dis[disnum_calc].node[i]);
    amredef(&(actnode->sol_increment),numr,actnode->numdf,"DA");
    amzero(&(actnode->sol_increment));
  }


  /*----------------------------------- set initial free surface position */
#ifdef D_FSI
  if (fdyn->freesurf>0)
  {
    for (i=0;i<numnp_total;i++)
    {
      actnode=&(actfield->dis[disnum_calc].node[i]);
      actgnode = actnode->gnode;
      if (actgnode->freesurf!=NULL)
      {
        actnode->xfs = (DOUBLE*)CCACALLOC(3,sizeof(DOUBLE));
        if (genprob.restart==0)
        {
          actnode->xfs[0] = actnode->x[0];
          actnode->xfs[1] = actnode->x[1];
          actnode->xfs[2] = actnode->x[2];
        }
        else /* get free surf position from ale-displacements */
        {
          dserror("restart for free surface not checked up to now!\n");
          actanode=actgnode->mfcpnode[genprob.numaf];
          actnode->xfs[0] = actanode->x[0]+actanode->sol_mf.a.da[1][0];
          actnode->xfs[1] = actanode->x[0]+actanode->sol_mf.a.da[1][1];
          actnode->xfs[2] = actanode->x[0]+actanode->sol_mf.a.da[1][2];
        }
      }
      else
        actnode->xfs=NULL;
    }
  }

  if (fdyn->freesurf==3 || fdyn->freesurf==5)
  {
    if (fdyn->freesurf==3) addcol=1;
    for (i=0;i<numnp_total;i++)
    {
      actnode=&(actfield->dis[disnum_calc].node[i]);
      if (actnode->xfs!=NULL)
      {
        amredef(&(actnode->sol_increment),numr,actnode->numdf+addcol,"DA");
        amzero(&(actnode->sol_increment));
        amredef(&(actnode->sol),actnode->sol.fdim,actnode->numdf+addcol,"DA");
        amzero(&(actnode->sol));
        amredef(&(actnode->sol_mf),actnode->sol.fdim,actnode->numdf+addcol,"DA");
        amzero(&(actnode->sol_mf));
        phi = actnode->x[numdf-2];
        actnode->sol.a.da[0][numdf] = phi;
        actnode->sol_increment.a.da[ipos->velnm][numdf] = phi;
        actnode->sol_increment.a.da[ipos->veln][numdf] = phi;
        actnode->sol_increment.a.da[ipos->velnp][numdf] = phi;
        actnode->sol_mf.a.da[0][numdf] = phi;
      }
    }
  }
#endif

  /*---------------------------------------------------------------------*
    INITIALISE ELEMENT ARRAYS
   *---------------------------------------------------------------------*/
  /*------------- allocate array for stabilisation parameter (only pdis) */
  if (fdyn->dyntyp==dyntyp_nln_time_int)
    for (i=0;i<actpart->pdis[disnum_calc].numele;i++)
    {
      actele = actpart->pdis[disnum_calc].element[i];
      switch(actele->eltyp)
      {
#ifdef D_FLUID2
        case el_fluid2:
          amdef("Stabpar",&(actele->e.f2->tau_old),3,1,"DV");
          amzero(&(actele->e.f2->tau_old));
#ifdef D_FLUID2_TDS
	  if(actele->distyp==tri3 || actele->distyp==tri6)
	  {
	      nis=1;
	  }
	  else
	  {
	      nis=actele->e.f2->nGP[1];
	  }

	  
	  amdef("subscale_pressure",
		&(actele->e.f2->sub_pres),
		actele->e.f2->nGP[0]*nis,1,
		"DV");
	  amzero(&(actele->e.f2->sub_pres));


	  amdef("subscale_velocities",
		&(actele->e.f2->sub_vel),
		3/*numdf (xyz)*/,
		actele->e.f2->nGP[0]*nis,
		"DA");

	  amzero(&(actele->e.f2->sub_vel));

	  if(fdyn->iop == 8)
	  {
	      amdef("subscale_pressure_acceleration",
		    &(actele->e.f2->sub_pres_acc),
		    actele->e.f2->nGP[0]*nis,1,
		    "DV");
	      amzero(&(actele->e.f2->sub_pres));
	      
	      
	      amdef("subscale_velocity_accelerations",
		    &(actele->e.f2->sub_vel_acc),
		    3/*numdf (xyz)*/,
		    actele->e.f2->nGP[0]*nis,
		    "DA");
	      
	      amzero(&(actele->e.f2->sub_vel));
	  }

	  
#endif
          break;
#endif

#ifdef D_FLUID3
        case el_fluid3:
          amdef("Stabpar",&(actele->e.f3->tau_old),3,1,"DV");
          amzero(&(actele->e.f3->tau_old));
          break;
#endif

#ifdef D_FLUID3_F
        case el_fluid3_fast:
          amdef("Stabpar",&(actele->e.f3->tau_old),3,1,"DV");
          amzero(&(actele->e.f3->tau_old));
          break;
#endif

#ifdef D_FLUID2_IS
        case el_fluid2_is:
#if 0
          amdef("Stabpar",&(actele->e.f2is->tau_old),3,1,"DV");
          amzero(&(actele->e.f2is->tau_old));
#endif
          break;
#endif

#ifdef D_FLUID3_IS
        case el_fluid3_is:
#if 0
          amdef("Stabpar",&(actele->e.f3is->tau_old),3,1,"DV");
          amzero(&(actele->e.f3is->tau_old));
#endif
          break;
#endif

        default:
          dserror("wrong type of element for fluid init!");
          break;
      }
    }

  /* check for stress calculation and allocate arrays */
  switch (str)
  {

    case str_none: /* do nothing */
      break;


#ifdef D_FSI
    case str_fsicoupling: /* allocate stress field for elements with fsi-coupled nodes */
      for (i=0;i<numele_total;i++)
      {
        INT found;
        actele = &(actfield->dis[disnum_calc].element[i]);
        numnp=actele->numnp;
        found=0;
        for (j=0;j<numnp;j++)
        {
          actnode=actele->node[j];
          actgnode = actnode->gnode;
          /* this approach does not work with a nonconforming discretization
             of the interface, thus it is replaced by the second one */
          /*
             if (actgnode->mfcpnode[0]!=NULL) found=1;
             */
          if (actgnode->fsicouple != NULL) found = 1;
        }

        if (found==1)
        {
#ifdef D_FLUID2
          if (actele->eltyp == el_fluid2)
          {
            amdef("stress_ND",&(actele->e.f2->stress_ND),numnp,3,"DA");
            amzero(&(actele->e.f2->stress_ND));
          }
#endif

#ifdef D_FLUID3
          if (actele->eltyp == el_fluid3)
          {
            amdef("stress_ND",&(actele->e.f3->stress_ND),numnp,6,"DA");
            amzero(&(actele->e.f3->stress_ND));
          }
#endif

#ifdef D_FLUID3_F
          if (actele->eltyp == el_fluid3_fast)
          {
            amdef("stress_ND",&(actele->e.f3->stress_ND),numnp,6,"DA");
            amzero(&(actele->e.f3->stress_ND));
          }
#endif

        }
      }


#ifdef D_FLUID3_F
      for (i=0; i<actpart->pdis[disnum_calc].num_fele; i++)
      {
        INT found;

        act_fast_eles = &(actpart->pdis[disnum_calc].fast_eles[i]);

        switch(act_fast_eles->fast_ele_typ)
        {
          case fele_f3f_hex8_e:
          case fele_f3f_hex8_a:
          case fele_f3f_hex20_e:
          case fele_f3f_hex20_a:
          case fele_f3f_tet4_e:
          case fele_f3f_tet4_a:

            found=0;
            for (l=0;l<act_fast_eles->aloopl;l++)
            {
              actele = act_fast_eles->ele_vec[l];
              numnp=actele->numnp;
              for (j=0;j<numnp;j++)
              {
                actnode=actele->node[j];
                actgnode = actnode->gnode;
                if (actgnode->fsicouple != NULL)
                  found = 1;
              }
            }

            if (found == 1)
            {
              for (l=0;l<act_fast_eles->aloopl;l++)
              {
                actele = act_fast_eles->ele_vec[l];
                if (actele->e.f3->stress_ND.Typ == cca_XX)
                {
                  amdef("stress_ND",&(actele->e.f3->stress_ND),numnp,6,"DA");
                  amzero(&(actele->e.f3->stress_ND));
                }

              }
            }
            break;

          default:
            break;

        } /* switch(act_fast_eles->fast_ele_typ) */
      }
#endif

      break;
#endif


    case str_all: /* allocate stress field for all elements */
      for (i=0;i<numele_total;i++)
      {
        actele = &(actfield->dis[disnum_calc].element[i]);
        numnp=actele->numnp;

#ifdef D_FLUID2
        if (actele->eltyp == el_fluid2)
        {
          amdef("stress_ND",&(actele->e.f2->stress_ND),numnp,3,"DA");
          amzero(&(actele->e.f2->stress_ND));
        }
#endif

#ifdef D_FLUID3
        if (actele->eltyp == el_fluid3)
        {
          amdef("stress_ND",&(actele->e.f3->stress_ND),numnp,6,"DA");
          amzero(&(actele->e.f3->stress_ND));
        }
#endif
      }

#ifdef D_FLUID3_F
      for (i=0; i<actpart->pdis[disnum_calc].num_fele; i++)
      {

        act_fast_eles = &(actpart->pdis[disnum_calc].fast_eles[i]);

        switch(act_fast_eles->fast_ele_typ)
        {
          case fele_f3f_hex8_e:
          case fele_f3f_hex8_a:
          case fele_f3f_hex20_e:
          case fele_f3f_hex20_a:
          case fele_f3f_tet4_e:
          case fele_f3f_tet4_a:

            for (l=0;l<act_fast_eles->aloopl;l++)
            {
              actele = act_fast_eles->ele_vec[l];
              if (actele->e.f3->stress_ND.Typ == cca_XX)
              {
                amdef("stress_ND",&(actele->e.f3->stress_ND),numnp,6,"DA");
                amzero(&(actele->e.f3->stress_ND));
              }

            }
            break;

          default:
            break;

        } /* switch(act_fast_eles->fast_ele_typ) */
      }
#endif

      break;


    case str_liftdrag:
      for (i=0;i<numele_total;i++)
      {
        actele = &(actfield->dis[disnum_calc].element[i]);
#ifdef D_FLUID2
        if (actele->eltyp == el_fluid2)
        {
          actgsurf=actele->g.gsurf;
          ldflag=0;
          for (j=0;j<actgsurf->ngline;j++)
          {
            actgline=actgsurf->gline[j];
            actdline=actgline->dline;
            if (actdline==NULL) continue;
            if (actdline->liftdrag==NULL) continue;
            ldflag++;
            break;
          }
          if (ldflag>0)
          {
            amdef("stress_ND",&(actele->e.f2->stress_ND),actele->numnp,6,"DA");
            amzero(&(actele->e.f2->stress_ND));
          }
        }
#endif

#ifdef D_FLUID3
        if (actele->eltyp == el_fluid3)
        {
          actgvol=actele->g.gvol;
          ldflag=0;
          for (j=0;j<actgvol->ngsurf;j++)
          {
            actgsurf=actgvol->gsurf[j];
            actdsurf=actgsurf->dsurf;
            if (actdsurf==NULL) continue;
            if (actdsurf->liftdrag==NULL) continue;
            ldflag++;
            break;
          }
          if (ldflag>0)
          {
            amdef("stress_ND",&(actele->e.f3->stress_ND),actele->numnp,6,"DA");
            amzero(&(actele->e.f3->stress_ND));
          }
        }
#endif
      } /* for (i=0;i<numele_total;i++) */


#ifdef D_FLUID3_F
      for (i=0; i<actpart->pdis[disnum_calc].num_fele; i++)
      {

        act_fast_eles = &(actpart->pdis[disnum_calc].fast_eles[i]);

        switch(act_fast_eles->fast_ele_typ)
        {
          case fele_f3f_hex8_e:
          case fele_f3f_hex8_a:
          case fele_f3f_hex20_e:
          case fele_f3f_hex20_a:
          case fele_f3f_tet4_e:
          case fele_f3f_tet4_a:

            ldflag=0;
            for (l=0;l<act_fast_eles->aloopl;l++)
            {
              actele = act_fast_eles->ele_vec[l];
              actgvol=actele->g.gvol;
              for (j=0;j<actgvol->ngsurf;j++)
              {
                actgsurf=actgvol->gsurf[j];
                actdsurf=actgsurf->dsurf;
                if (actdsurf==NULL) continue;
                if (actdsurf->liftdrag==NULL) continue;
                ldflag++;
                break;
              }
            }

            if (ldflag > 0)
            {
              for (l=0;l<act_fast_eles->aloopl;l++)
              {
                actele = act_fast_eles->ele_vec[l];
                if (actele->e.f3->stress_ND.Typ == cca_XX)
                {
                  amdef("stress_ND",&(actele->e.f3->stress_ND),actele->numnp,6,"DA");
                  amzero(&(actele->e.f3->stress_ND));
                }

              }
            }
            break;

          default:
            break;

        } /* switch(act_fast_eles->fast_ele_typ) */
      }
#endif

      break;


    default:
      dserror("option for 'str' not implemented yet!\n");
  }

  /*--------------- allocate curvature field for elements at free surface */
  if (fdyn->surftens>0)
  {
    for (i=0;i<numele_total;i++)
    {
      actele = &(actfield->dis[disnum_calc].element[i]);
#ifdef D_FLUID2
      if (actele->eltyp == el_fluid2)
      {
        if (actele->e.f2->fs_on==0) continue;
        {
          amdef("kappa_ND",&(actele->e.f2->kappa_ND),numnp,2,"DA");
          amzero(&(actele->e.f2->kappa_ND));
        }
      }
#endif

#ifdef D_FLUID3
      if (actele->eltyp == el_fluid3)
        dserror("curvature not implemented for 3D fluid elements!\n");
#endif

#ifdef D_FLUID3_F
      if (actele->eltyp == el_fluid3_fast)
        dserror("curvature not implemented for 3D fluid fast elements!\n");
#endif
    }
  }

  /*-------- inherit the Neumann conditions from design to discretization */
  inherit_design_dis_neum(&(actfield->dis[disnum_calc]));

  /*-------------------------------------------- create the initial field */
  if (fdyn->init>=1)
  {
    if (fdyn->init==1) /*------------ initial data from fluid_start.data */
      inp_fluid_start_data(actfield, ipos, fdyn);

    if (fdyn->init==2) { /*------------------ initial data from pss file */
#if defined(BINIO)
      restart_read_bin_fluiddyn(fdyn, NULL, NULL, actfield, actpart, disnum_calc,
          actintra, fdyn->resstep);
#else
      restart_read_fluiddyn(fdyn->resstep,fdyn,actfield,actpart,actintra,
          action,container);
#endif
    }

#ifdef D_FSI
    if (fdyn->init==6) /*--------------------------------- solitary wave */
    {
      actele = &(actfield->dis[disnum_calc].element[0]);
      actmat = actele->mat-1;
      numnp_total=actfield->dis[disnum_calc].numnp;
      /*-------------------------------------------- set some constants */
      dens   = mat[actmat].m.fluid->density;
      g      = -actele->g.gsurf->neum->neum_val.a.dv[1];
      d      = TEN;
      H      = ONE/FIVE*d;
      fac1   = sqrt((THREE*H)/(FOUR*d*d*d));
      c      = sqrt(g*d*(1+H/d));
      t      = ZERO;
      for (i=0;i<numnp_total;i++) /* loop nodes */
      {
        actnode=&(actfield->dis[disnum_calc].node[i]);
        x   = actnode->x[0];
        y   = actnode->x[1];
        fac = fac1*(x-c*t);
        /* sech = sech(fac*xct)**2 = 1/cosh(fac*xct)**2 */
        sech = cosh(fac);
        sech = ONE/sech/sech;
        eta  = d + H*sech;
        /*---------------------------------------- modify initial values */
        p    = dens*g*(eta-y);
        u1   = sqrt(g*d)*H/d*sech;
        u2   = sqrt(THREE*g*d)*(H/d)*sqrt(H/d)*(y/d)*sech*tanh(fac);
        actgnode = actnode->gnode;
        /*---------------------------------------- write values to nodes */
        actnode->sol.a.da[0][0] = u1;
        actnode->sol.a.da[0][1] = u2;
        actnode->sol.a.da[0][2] = p ;
        actnode->sol_increment.a.da[ipos->veln][0] = u1;
        actnode->sol_increment.a.da[ipos->veln][1] = u2;
        actnode->sol_increment.a.da[ipos->veln][2] = p ;
        actnode->sol_increment.a.da[ipos->velnp][0] = u1;
        actnode->sol_increment.a.da[ipos->velnp][1] = u2;
        actnode->sol_increment.a.da[ipos->velnp][2] = p ;
        if (fdyn->freesurf==2 && actnode->xfs!=NULL)
        {
          actnode->sol_increment.a.da[ipos->veln][3] = u1;
          actnode->sol_increment.a.da[ipos->veln][4] = u2;
          actnode->sol_increment.a.da[ipos->velnp][3] = u1;
          actnode->sol_increment.a.da[ipos->velnp][4] = u2;
          actgnode = actnode->gnode;
          if (actgnode->dirich!=NULL)
          {
            if (actgnode->dirich->dirich_type!=dirich_none)
              dserror("dirich type not allowed for solwave init\n");
            for (j=3;j<actnode->numdf;j++) /* loop all grid dofs */
            {
              if (actgnode->dirich->dirich_onoff.a.iv[j]==1)
              {
                actnode->sol_increment.a.da[ipos->veln][j] = ZERO;
                actnode->sol_increment.a.da[ipos->velnp][j] = ZERO;
              }
            } /* end loop over dofs */
          }
        }
      }
    }
#endif

#ifdef D_FSI
    if (fdyn->init==7) /*-------------------------- wavebreaking problem */
    {
      actele = &(actfield->dis[disnum_calc].element[0]);
      actmat = actele->mat-1;
      numnp_total=actfield->dis[disnum_calc].numnp;
      /*-------------------------------------------- set some constants */
      dens   = mat[actmat].m.fluid->density;
      g      = -actele->g.gsurf->neum->neum_val.a.dv[1];
      d      = TEN;
      H      = ONE/FIVE*d;
      fac1   = sqrt((THREE*H)/(FOUR*d*d*d));
      c      = sqrt(g*d*(1+H/d));
      t      = ZERO;
      for (i=0;i<numnp_total;i++) /* loop nodes */
      {
        actnode=&(actfield->dis[disnum_calc].node[i]);
        x   = actnode->x[0];
        y   = actnode->x[1];
        fac = fac1*(x-c*t);
        /* sech = sech(fac*xct)**2 = 1/cosh(fac*xct)**2 */
        sech = cosh(fac);
        sech = ONE/sech/sech;
        eta  = d + H*sech;
        /*---------------------------------------- modify initial values */
        if (x>=80.0)
        {
          p    = dens*g*(d-y);
          u1 = ZERO;
          u2 = ZERO;
        }
        else
        {
          p    = dens*g*(eta-y);
          u1   = sqrt(g*d)*H/d*sech;
          u2   = sqrt(THREE*g*d)*(H/d)*sqrt(H/d)*(y/d)*sech*tanh(fac);
        }
        /*---------------------------------------- write values to nodes */
        actnode->sol.a.da[0][0] = u1;
        actnode->sol.a.da[0][1] = u2;
        actnode->sol.a.da[0][2] = p ;
        actnode->sol_increment.a.da[ipos->veln][0] = u1;
        actnode->sol_increment.a.da[ipos->veln][1] = u2;
        actnode->sol_increment.a.da[ipos->veln][2] = p ;
        actnode->sol_increment.a.da[ipos->velnp][0] = u1;
        actnode->sol_increment.a.da[ipos->velnp][1] = u2;
        actnode->sol_increment.a.da[ipos->velnp][2] = p ;
        if (fdyn->freesurf==2 && actnode->xfs!=NULL)
        {
          actnode->sol_increment.a.da[ipos->veln][3] = u1;
          actnode->sol_increment.a.da[ipos->veln][4] = u2;
          actnode->sol_increment.a.da[ipos->velnp][3] = u1;
          actnode->sol_increment.a.da[ipos->velnp][4] = u2;
          actgnode = actnode->gnode;
          if (actgnode->dirich!=NULL)
          {
            if (actgnode->dirich->dirich_type!=dirich_none)
              dserror("dirich type not allowed for solwave init\n");
            for (j=3;j<actnode->numdf;j++) /* loop all grid dofs */
            {
              if (actgnode->dirich->dirich_onoff.a.iv[j]==1)
              {
                actnode->sol_increment.a.da[ipos->veln][j] = ZERO;
                actnode->sol_increment.a.da[ipos->velnp][j] = ZERO;
              }
            } /* end loop over dofs */
          }
        }
      }
    }
#endif
    if (fdyn->init==8) /* Beltrami flow */
    {
      actele = &(actfield->dis[disnum_calc].element[0]);
      actmat = actele->mat-1;
      numnp_total=actfield->dis[disnum_calc].numnp;
      /* set some constants */
      visc   = mat[actmat].m.fluid->viscosity;
      a      = PI/4.0;
      d      = PI/2.0;
      t      = 0.0;
      for (i=0;i<numnp_total;i++) /* loop nodes */
      {
        actnode=&(actfield->dis[disnum_calc].node[i]);
        x1   = actnode->x[0];
        x2   = actnode->x[1];
        x3   = actnode->x[2];
        /* calculate initial values */
        p    = -a*a/2.0 * ( exp(2.0*a*x1) + exp(2.0*a*x2) + exp(2.0*a*x3)
            + 2.0 * sin(a*x1 + d*x2) * cos(a*x3 + d*x1) * exp(a*(x2+x3))
            + 2.0 * sin(a*x2 + d*x3) * cos(a*x1 + d*x2) * exp(a*(x3+x1))
            + 2.0 * sin(a*x3 + d*x1) * cos(a*x2 + d*x3) * exp(a*(x1+x2)) )
          * exp(-2.0*visc*d*d*t);
        u1   = -a * ( exp(a*x1) * sin(a*x2 + d*x3) +
            exp(a*x3) * cos(a*x1 + d*x2) ) * exp(-visc*d*d*t);
        u2   = -a * ( exp(a*x2) * sin(a*x3 + d*x1) +
            exp(a*x1) * cos(a*x2 + d*x3) ) * exp(-visc*d*d*t);
        u3   = -a * ( exp(a*x3) * sin(a*x1 + d*x2) +
            exp(a*x2) * cos(a*x3 + d*x1) ) * exp(-visc*d*d*t);
        /* write values to nodes */
        actnode->sol.a.da[0][0] = u1;
        actnode->sol.a.da[0][1] = u2;
        actnode->sol.a.da[0][2] = u3;
        actnode->sol.a.da[0][3] = p ;
        actnode->sol_increment.a.da[ipos->velnm][0] = u1;
        actnode->sol_increment.a.da[ipos->velnm][1] = u2;
        actnode->sol_increment.a.da[ipos->velnm][2] = u3;
        actnode->sol_increment.a.da[ipos->velnm][3] = p ;
        actnode->sol_increment.a.da[ipos->veln][0] = u1;
        actnode->sol_increment.a.da[ipos->veln][1] = u2;
        actnode->sol_increment.a.da[ipos->veln][2] = u3;
        actnode->sol_increment.a.da[ipos->veln][3] = p ;
        actnode->sol_increment.a.da[ipos->velnp][0] = u1;
        actnode->sol_increment.a.da[ipos->velnp][1] = u2;
        actnode->sol_increment.a.da[ipos->velnp][2] = u3;
        actnode->sol_increment.a.da[ipos->velnp][3] = p ;
      }
    }

    if (fdyn->init==9) /* Kim-Moin flow */
    {
      actele = &(actfield->dis[disnum_calc].element[0]);
      actmat = actele->mat-1;
      numnp_total=actfield->dis[disnum_calc].numnp;
      /* set some constants */
      visc   = mat[actmat].m.fluid->viscosity;
      a      = 2.0;
      t      = 0.0;
      for (i=0;i<numnp_total;i++) /* loop nodes */
      {
        actnode=&(actfield->dis[disnum_calc].node[i]);
        x1   = actnode->x[0];
        x2   = actnode->x[1];
        /* calculate initial values */
        p    = -1.0/4.0 * ( cos(2.0*a*PI*x1) + cos(2.0*a*PI*x2) ) * exp(-4.0*a*a*PI*PI*t*visc);
        u1   = - cos(a*PI*x1) * sin(a*PI*x2) * exp(-2.0*a*a*PI*PI*t*visc);
        u2   = + sin(a*PI*x1) * cos(a*PI*x2) * exp(-2.0*a*a*PI*PI*t*visc);
        /* write values to nodes */
        actnode->sol.a.da[0][0] = u1;
        actnode->sol.a.da[0][1] = u2;
        actnode->sol.a.da[0][2] = p ;
        actnode->sol_increment.a.da[ipos->velnm][0] = u1;
        actnode->sol_increment.a.da[ipos->velnm][1] = u2;
        actnode->sol_increment.a.da[ipos->velnm][2] = p ;
        actnode->sol_increment.a.da[ipos->veln][0] = u1;
        actnode->sol_increment.a.da[ipos->veln][1] = u2;
        actnode->sol_increment.a.da[ipos->veln][2] = p ;
        actnode->sol_increment.a.da[ipos->velnp][0] = u1;
        actnode->sol_increment.a.da[ipos->velnp][1] = u2;
        actnode->sol_increment.a.da[ipos->velnp][2] = p ;
      }
    }
  }

  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of fluid_init */


/*!---------------------------------------------------------------------
\brief storing results in solution history

<pre>                                                         genk 05/02

in this routine the results in the DIST_VECTOR are put to the nodes in
a certain place in ARRAY sol_increment.
Result has to be allreduced and is put to the whole field on each proc.
If necassary the norms for the iteration convergence check of the
nonlinear iteration scheme are calculated.

</pre>
\param **actfield      FIELD	      (i)    actual field
\param  *actintra      INTRA	      (i)    actual intra comm.
\param	*sol 	       DIST_VECTOR    (i)    solution vector
\param	 place         INT	      (i)    place in sol_incr.
\param	*sysarray      SPARSE_ARRAY   (i)
\param	*sysarray_typ  SPARSE_TYP     (i)
\param	*vrat          DOUBLE	      (o)    vel.  conv. ratio
\param	*prat          DOUBLE	      (o)    pre.  conv. ratio
\param  *grat          DOUBLE         (o)    grid  conv. ratio
\return void

------------------------------------------------------------------------*/
void fluid_result_incre(
                          FIELD             *actfield,
                          INT                disnum,
                          INTRA             *actintra,
			  DIST_VECTOR       *sol,
			  DIST_VECTOR       *rhs,
                          INT                place,
			  SPARSE_ARRAY      *sysarray,
			  SPARSE_TYP        *sysarray_typ,
			  DOUBLE            *vrat,
			  DOUBLE            *prat,
                          DOUBLE            *grat
		       )
{
INT      i,j;          /* simply some counters                         */
INT      max;
INT      diff;
INT      dof;          /* actual dof number                            */
INT      predof;       /* number of pressure dof (2=2D; 3=3D)          */
INT      numeq_total;  /* total number of equations                    */
DOUBLE   dvnorm=ZERO;
DOUBLE   dpnorm=ZERO;
DOUBLE   dgnorm=ZERO;
DOUBLE    vnorm=ZERO;
DOUBLE    pnorm=ZERO;  /* values for norm calculation                  */
DOUBLE    gnorm=ZERO;
NODE    *actnode;      /* actual node                                  */
static ARRAY    result_a;
static DOUBLE  *result;       /* redundandent result vector                   */

#ifdef DEBUG
dstrc_enter("fluid_result_incre");
#endif

fdyn = alldyn[genprob.numff].fdyn;
/*----------------------------------------------------------------------*/
numeq_total = sol->numeq_total;
predof      = fdyn->numdf-1;

/*------------------------- allocate space to allreduce the DIST_VECTOR */
if(result_a.Typ==cca_XX)
   result = amdef("result",&result_a,numeq_total,1,"DV");
/* Chimera has fluid discretizations of different size. */
if (result_a.fdim < numeq_total) {
  amdel(&result_a);
  result = amdef("result",&result_a,numeq_total,1,"DV");
}
amzero(&result_a);

/*------------------ copy distributed result to redundant result vector */
solserv_reddistvec(
                      sol,
                      sysarray,
                      sysarray_typ,
                      result,
                      sol->numeq_total,
                      actintra
                     );

/* values in sol_increment are in XYZ co-sys - result is in xyz* co-sys
   so we have to tranform sol_increment for the convergence check       */
locsys_trans_sol(actfield,disnum,1,place,0);

switch (fdyn->itnorm) /* switch to norm */
{
case fncc_no:
   dvnorm=ONE;
   dpnorm=ONE;
   dgnorm=ONE;
   /*--------  loop nodes and put the result back to the node structure */
   for (i=0; i<actfield->dis[disnum].numnp; i++)
   {
      actnode = &(actfield->dis[disnum].node[i]);
      /*--------------------------- enlarge sol_increment, if necessary */
      if (place >= actnode->sol_increment.fdim)
      {
         diff = place - actnode->sol_increment.fdim;
         max  = IMAX(diff,5);
         amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,
	           actnode->sol_increment.sdim,"DA");
      } /* endif enlargement */
      for (j=0; j<actnode->numdf; j++) /* loop dofs */
      {
         dof = actnode->dof[j];
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
         if (actnode->gnode->dirich!=NULL &&
             actnode->gnode->dirich->dirich_onoff.a.iv[j]!=0)
           continue;
#else
         if (dof>=numeq_total) continue;
#endif
#ifdef FLUID_INCREMENTAL
         actnode->sol_increment.a.da[place][j] += result[dof];
#else
         actnode->sol_increment.a.da[place][j] = result[dof];
#endif
      } /* end of loop over dofs */

   } /* end of loop over nodes */
break;
/*----------------------------------------------------------------------*/
case fncc_Linf: /* L_infinity norm */
   /*-----  loop nodes and put the result back to the node structure */
   for (i=0; i<actfield->dis[disnum].numnp; i++)
   {
      actnode = &(actfield->dis[disnum].node[i]);
      /*------------------------ enlarge sol_increment, if necessary */
      if (place >= actnode->sol_increment.fdim)
      {
         diff = place - actnode->sol_increment.fdim;
         max  = IMAX(diff,5);
         amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,
                   actnode->sol_increment.sdim,"DA");
      } /* endif enlargement */
      for (j=0; j<actnode->numdf; j++) /* loop dofs and calculate the norms */
      {
         dof = actnode->dof[j];
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
         if (actnode->gnode->dirich!=NULL &&
             actnode->gnode->dirich->dirich_onoff.a.iv[j]!=0)
           continue;
#else
         if (dof>=numeq_total) continue;
#endif
         if (j==predof) /* pressure dof */
         {
#ifdef FLUID_INCREMENTAL
            pnorm = DMAX(pnorm,
                FABS(result[dof]+actnode->sol_increment.a.da[place][j]));
            dpnorm  = DMAX(dpnorm, FABS(result[dof]));
            actnode->sol_increment.a.da[place][j] += result[dof];
#else
            dpnorm = DMAX(dpnorm,
                FABS(result[dof]-actnode->sol_increment.a.da[place][j]));
            pnorm  = DMAX(pnorm, FABS(result[dof]));
            actnode->sol_increment.a.da[place][j] = result[dof];
#endif
         } /* endif pressure dof */
         else if (j>predof) /* grid velocity or height function dof */
         {
#ifdef FLUID_INCREMENTAL
            gnorm = DMAX(gnorm,
                FABS(result[dof]+actnode->sol_increment.a.da[place][j]));
            dgnorm  = DMAX(dgnorm, FABS(result[dof]));
            actnode->sol_increment.a.da[place][j] += result[dof];
#else
            dgnorm = DMAX(dgnorm,
                FABS(result[dof]-actnode->sol_increment.a.da[place][j]));
            gnorm  = DMAX(gnorm,FABS(result[dof]));
            actnode->sol_increment.a.da[place][j] = result[dof];
#endif
         }
         else /* vel - dof */
         {
#ifdef FLUID_INCREMENTAL
            vnorm = DMAX(vnorm,
                FABS(result[dof]+actnode->sol_increment.a.da[place][j]));
            dvnorm  = DMAX(dvnorm, FABS(result[dof]));
            actnode->sol_increment.a.da[place][j] += result[dof];
#else
            dvnorm = DMAX(dvnorm,
                FABS(result[dof]-actnode->sol_increment.a.da[place][j]));
            vnorm  = DMAX(vnorm, FABS(result[dof]));
            actnode->sol_increment.a.da[place][j] = result[dof];
#endif
         } /* endif vel dof */
      } /* end of loop over dofs */
   } /* end of loop over nodes */
break;
/*-------------------------------------------------------------------------*/
case fncc_L1: /* L_1 norm */
   /*-----------  loop nodes and put the result back to the node structure */
   for (i=0; i<actfield->dis[disnum].numnp; i++)
   {
      actnode = &(actfield->dis[disnum].node[i]);
      /*------------------------------ enlarge sol_increment, if necessary */
      if (place >= actnode->sol_increment.fdim)
      {
         diff = place - actnode->sol_increment.fdim;
         max  = IMAX(diff,5);
         amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,
                   actnode->sol_increment.sdim,"DA");
      } /* endif enlargement */
      for (j=0; j<actnode->numdf; j++) /* loop dofs and calculate the norms */
      {
         dof = actnode->dof[j];
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
         if (actnode->gnode->dirich!=NULL &&
             actnode->gnode->dirich->dirich_onoff.a.iv[j]!=0)
           continue;
#else
         if (dof>=numeq_total) continue;
#endif
         if (j==predof) /* pressure dof */
         {
#ifdef FLUID_INCREMENTAL
            pnorm += FABS(result[dof]+actnode->sol_increment.a.da[place][j]);
            dpnorm  += FABS(result[dof]);
            actnode->sol_increment.a.da[place][j] += result[dof];
#else
            dpnorm += FABS(result[dof]-actnode->sol_increment.a.da[place][j]);
            pnorm  += FABS(result[dof]);
            actnode->sol_increment.a.da[place][j] = result[dof];
#endif
         } /* endif pressure dof */
         else if (j>predof) /* grid velocity or height function dof */
         {
#ifdef FLUID_INCREMENTAL
            gnorm += FABS(result[dof]+actnode->sol_increment.a.da[place][j]);
            dgnorm  += FABS(result[dof]);
            actnode->sol_increment.a.da[place][j] += result[dof];
#else
            dgnorm += FABS(result[dof]-actnode->sol_increment.a.da[place][j]);
            gnorm  += FABS(result[dof]);
            actnode->sol_increment.a.da[place][j] = result[dof];
#endif
         } /* endif grid velocity dofs */
         else /* vel - dof */
         {
#ifdef FLUID_INCREMENTAL
            vnorm += FABS(result[dof]+actnode->sol_increment.a.da[place][j]);
            dvnorm  += FABS(result[dof]);
            actnode->sol_increment.a.da[place][j] += result[dof];
#else
            dvnorm += FABS(result[dof]-actnode->sol_increment.a.da[place][j]);
            vnorm  += FABS(result[dof]);
            actnode->sol_increment.a.da[place][j] = result[dof];
#endif
         } /* endif vel dof */
      } /* end of loop over dofs */
   } /* end of loop over nodes */
break;
/*-------------------------------------------------------------------------*/
case fncc_L2: /* L_2 norm */
   /*-----------  loop nodes and put the result back to the node structure */
   for (i=0; i<actfield->dis[disnum].numnp; i++)
   {
      actnode = &(actfield->dis[disnum].node[i]);
      /*------------------------------ enlarge sol_increment, if necessary */
      if (place >= actnode->sol_increment.fdim)
      {
         diff = place - actnode->sol_increment.fdim;
         max  = IMAX(diff,5);
         amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,
                   actnode->sol_increment.sdim,"DA");
      }
      for (j=0; j<actnode->numdf; j++) /* loop dofs and calculate the norms */
      {
         dof = actnode->dof[j];
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
         if (actnode->gnode->dirich!=NULL &&
             actnode->gnode->dirich->dirich_onoff.a.iv[j]!=0)
           continue;
#else
         if (dof>=numeq_total) continue;
#endif
         if (j==predof) /* pressure dof */
         {
#ifdef FLUID_INCREMENTAL
            pnorm  += DSQR(result[dof]+actnode->sol_increment.a.da[place][j]);
            dpnorm += DSQR(result[dof]);
            actnode->sol_increment.a.da[place][j] += result[dof];
#else
            dpnorm += DSQR(result[dof]-actnode->sol_increment.a.da[place][j]);
            pnorm  += DSQR(result[dof]);
            actnode->sol_increment.a.da[place][j] = result[dof];
#endif
         } /* endif pressure dof */
         else if (j>predof) /* grid velocity or height function dof */
         {
#ifdef FLUID_INCREMENTAL
            gnorm  += DSQR(result[dof]+actnode->sol_increment.a.da[place][j]);
            dgnorm += DSQR(result[dof]);
            actnode->sol_increment.a.da[place][j] += result[dof];
#else
            dgnorm += DSQR(result[dof]-actnode->sol_increment.a.da[place][j]);
            gnorm  += DSQR(result[dof]);
            actnode->sol_increment.a.da[place][j] = result[dof];
#endif
         } /* endif grid velocity dofs */
         else /* vel - dof */
         {
#ifdef FLUID_INCREMENTAL
            vnorm  += DSQR(result[dof]+actnode->sol_increment.a.da[place][j]);
            dvnorm += DSQR(result[dof]);
            actnode->sol_increment.a.da[place][j] += result[dof];
#else
            dvnorm += DSQR(result[dof]-actnode->sol_increment.a.da[place][j]);
            vnorm  += DSQR(result[dof]);
            actnode->sol_increment.a.da[place][j] = result[dof];
#endif
         } /* endif vel dof */
      } /* end of loop over dofs */
   } /* end of loop over nodes */
   dvnorm = sqrt(dvnorm);
    vnorm = sqrt( vnorm);
   dpnorm = sqrt(dpnorm);
    pnorm = sqrt( pnorm);
   dgnorm = sqrt(dgnorm);
    gnorm = sqrt( gnorm);
break;
/*-------------------------------------------------------------------*/
default:
   dserror("unknown norm for convergence check!\n");
} /* end of switch(fdyn->itnorm) */

/*------------------------------------------- check for "ZERO-field" */
if (vnorm<EPS5)
{
   vnorm = ONE;
#ifdef DEBUG
    printf("ATTENTION: zero vel field - norm <= 1.0e-5 set to 1.0!! \n");
#endif
}
/* If prat==NULL we are not interessted in the pressure norm. This is
 * the case in the projection method where we use discontinous
 * pressures. */
if (prat!=NULL)
{
  if (pnorm<EPS5)
  {
    pnorm = ONE;
#ifdef DEBUG
    printf("ATTENTION: zero pre field - norm <= 1.0e-5 set to 1.0!! \n");
#endif
  }
}
if (gnorm<EPS5)
{
   gnorm = ONE;
}
/*------------------------------------- set final convergence ratios */
*vrat = dvnorm/vnorm;
if (prat!=NULL)
  *prat = dpnorm/pnorm;
if (fdyn->freesurf==2 || fdyn->freesurf==5 || fdyn->freesurf==6)
   *grat = dgnorm/gnorm;

/*------------------------------------------- local co-ordinate system:
  the values in sol_increment[3][] are given in the xyz* co-system,
  however everything else is in the XYZ co-system. So we have to
  tranform them from xyz* to XYZ                                        */
locsys_trans_sol(actfield,disnum,1,place,1);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of  fluid_result_incre */


/*!---------------------------------------------------------------------
\brief calculating norms for steady state check

<pre>                                                         genk 05/02

in this routine the velocity and pressure norms for the steady state
check are calculated:
   norm = ||U(n+1) - U(n)|| / ||U(n)||
      solution at (n+1): node->sol_increment[3][j]
      solution at (n)  : node->sol_increment[1][j]

</pre>
\param *fdyn 	      FLUID_DYNAMIC  (i)
\param *actfield      FIELD	     (i)   actual field
\param  numeq_total   INT	     (i)   total number of equations
\param *vrat	      DOUBLE	     (o)   vel.  conv. ratio
\param *prat	      DOUBLE	     (o)   pres. conv. ratio
\param *ipos                         (i)   node array positions
\return void

------------------------------------------------------------------------*/
void fluid_norm(
    FIELD             *actfield,
    INT                disnum,
    INT                numeq_total,
    DOUBLE            *vrat,
    DOUBLE            *prat,
    ARRAY_POSITION    *ipos
    )
{
INT         i,j;           /* simply some counters                      */
INT         numdf;         /* number of fluid dofs                      */
INT         numvel;        /* total number of vel-dofs                  */
INT         predof;        /* actual number of pres dof                 */
INT         numnp_total;   /* total number of fluid nodes               */
INT         actdof;        /* actual dof number                         */
DOUBLE      dvnorm=ZERO;   /* norms					*/
DOUBLE       vnorm=ZERO;   /* norms 					*/
DOUBLE      dpnorm=ZERO;   /* norms					*/
DOUBLE       pnorm=ZERO;   /* norms                                     */
NODE       *actnode;       /* actual node                               */

#ifdef DEBUG
dstrc_enter("fluid_norm");
#endif

/*---------------------------------------------------- set some values */
fdyn = alldyn[genprob.numff].fdyn;

numdf        = fdyn->numdf;
numnp_total  = actfield->dis[disnum].numnp;
predof       = numdf-1;
numvel       = numdf-1;

switch (fdyn->stnorm)
{
case fnst_no: /* do nothing */
break;
case fnst_Linf: /* L_infinity norm */
   /*-------------------------------------------------- loop all nodes */
   for (i=0;i<numnp_total;i++) /* loop nodes */
   {
      actnode=&(actfield->dis[disnum].node[i]);
      for (j=0;j<numvel;j++) /* loop vel-dofs */
      {
	 actdof = actnode->dof[j];
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
         if (actnode->gnode->dirich!=NULL &&
             actnode->gnode->dirich->dirich_onoff.a.iv[j]!=0)
           continue;
#else
         if (actdof>=numeq_total) continue;
#endif
         dvnorm = DMAX(dvnorm,FABS(actnode->sol_increment.a.da[ipos->velnp][j]  \
	                          -actnode->sol_increment.a.da[ipos->veln][j]));
          vnorm = DMAX( vnorm,FABS(actnode->sol_increment.a.da[ipos->velnp][j]));
      } /* end of loop over vel-dofs */
      actdof = actnode->dof[predof];
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
      if (actnode->gnode->dirich!=NULL &&
          actnode->gnode->dirich->dirich_onoff.a.iv[j]!=0)
        continue;
#else
      if (actdof>=numeq_total) continue;
#endif
      dpnorm = DMAX(dpnorm,FABS(actnode->sol_increment.a.da[ipos->velnp][predof]  \
	                       -actnode->sol_increment.a.da[ipos->veln][predof]));
       pnorm = DMAX( pnorm,FABS(actnode->sol_increment.a.da[ipos->velnp][predof]));
   } /* end of loop over nodes */
break;
/*----------------------------------------------------------------------*/
case fnst_L1: /* L_1 norm */
   /*--------------------------------------------------- loop all nodes */
   for (i=0;i<numnp_total;i++) /* loop nodes */
   {
      actnode=&(actfield->dis[disnum].node[i]);
      for (j=0;j<numvel;j++) /* loop vel-dofs */
      {
	 actdof = actnode->dof[j];
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
         if (actnode->gnode->dirich!=NULL &&
             actnode->gnode->dirich->dirich_onoff.a.iv[j]!=0)
           continue;
#else
         if (actdof>=numeq_total) continue;
#endif
         dvnorm += FABS(actnode->sol_increment.a.da[ipos->velnp][j]  \
	               -actnode->sol_increment.a.da[ipos->veln][j]);
          vnorm += FABS(actnode->sol_increment.a.da[ipos->velnp][j]);
      } /* end of loop over vel-dofs */
      actdof = actnode->dof[predof];
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
      if (actnode->gnode->dirich!=NULL &&
          actnode->gnode->dirich->dirich_onoff.a.iv[j]!=0)
        continue;
#else
      if (actdof>=numeq_total) continue;
#endif
      dpnorm += FABS(actnode->sol_increment.a.da[ipos->velnp][predof]  \
	            -actnode->sol_increment.a.da[ipos->veln][predof]);
       pnorm += FABS(actnode->sol_increment.a.da[ipos->velnp][predof]);
   } /* end of loop over nodes */
break;
/*----------------------------------------------------------------------*/
case fnst_L2: /* L_2 norm */
   /*--------------------------------------------------- loop all nodes */
   for (i=0;i<numnp_total;i++) /* loop nodes */
   {
      actnode=&(actfield->dis[disnum].node[i]);
      for (j=0;j<numvel;j++) /* loop vel-dofs */
      {
	 actdof = actnode->dof[j];
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
         if (actnode->gnode->dirich!=NULL &&
             actnode->gnode->dirich->dirich_onoff.a.iv[j]!=0)
           continue;
#else
         if (actdof>=numeq_total) continue;
#endif
         dvnorm += DSQR(actnode->sol_increment.a.da[ipos->velnp][j]  \
	               -actnode->sol_increment.a.da[ipos->veln][j]);
          vnorm += DSQR(actnode->sol_increment.a.da[ipos->velnp][j]);
      } /* end of loop over vel-dofs */
      if (actnode->numdf <= predof)
        continue;
      actdof = actnode->dof[predof];
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
      if (actnode->gnode->dirich!=NULL &&
          actnode->gnode->dirich->dirich_onoff.a.iv[j]!=0)
        continue;
#else
      if (actdof>=numeq_total) continue;
#endif
      dpnorm += DSQR(actnode->sol_increment.a.da[ipos->velnp][predof]  \
	            -actnode->sol_increment.a.da[ipos->veln][predof]);
       pnorm += DSQR(actnode->sol_increment.a.da[ipos->velnp][predof]);
   } /* end of loop over nodes */
   dvnorm = sqrt(dvnorm);
    vnorm = sqrt( vnorm);
   dpnorm = sqrt(dpnorm);
    pnorm = sqrt( pnorm);
break;
/*----------------------------------------------------------------------*/
default:
   dserror("unknown norm for steady state check!\n");
} /* end of switch(fdyn->stnorm) */

/*---------------------------------------------- check for "ZERO-field" */
if (vnorm<EPS5)
{
   vnorm = ONE;
   printf("ATTENTION: zero vel field - norm <= 1.0e-5 set to 1.0!! \n");
}
if (pnorm<EPS5)
{
   pnorm = ONE;
   printf("ATTENTION: zero pre field - norm <= 1.0e-5 set to 1.0!! \n");
}

/*---------------------------------------- set final convergence ratios */
*vrat = dvnorm/vnorm;
*prat = dpnorm/pnorm;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_norm*/

/*!---------------------------------------------------------------------
\brief copy solution history

<pre>                                                         genk 05/02

in this routine the solution at postion 'from' in the nodal solution
history is copied to the positon 'to'.
The pressure is transformed from kinematic to real pressure

</pre>
\param *actfield      FIELD	     (i)  actual field
\param  disnum        INT	     (i)  number of the discr.
\param  arrayfrom     INT	     (i)  index of the "from-array"
\param  arrayto       INT	     (i)  index of the "to-array"
\param  from          INT            (i)  place in "from-array"
\param  to            INT            (i)  place in "to-array"
\param  numdf         INT            (i)  number of dofs -> 2D or 3D probl.
\return void
\sa fluid_init()

------------------------------------------------------------------------*/
void fluid_sol_copy(
                          FIELD             *actfield,
			  INT                disnum,
			  INT                arrayfrom,
			  INT                arrayto,
                          INT                from,
		          INT                to,
		          INT                numdf
		  )
{
INT         i,j;	/* simply some counters				*/
INT         diff,max;	/* integers for amredef				*/
INT         predof;	/* pressure dof					*/
DOUBLE      dens;	/* density					*/
NODE       *actnode;	/* actual node					*/
ELEMENT    *actele;	/* actual element				*/
DISCRET    *actdis;
ARRAY      *arrayf,*arrayt;

#ifdef DEBUG
dstrc_enter("fluid_sol_copy");
#endif

predof=numdf-1;

/* since different materials are not allowed  one can work with the
   material parameters of any element ---------------------------------*/
actele = &(actfield->dis[disnum].element[0]);
dens  = mat[actele->mat-1].m.fluid->density;

actdis = &(actfield->dis[disnum]);
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   /*----------------------------------------- select correct arrayfrom */
   switch(arrayfrom)
   {
   case 0:
      arrayf = &(actnode->sol);
   break;
   case 1:
      arrayf = &(actnode->sol_increment);
   break;
   case 2:
      arrayf = &(actnode->sol_residual);
   break;
   case 3:
      arrayf = &(actnode->sol_mf);
   break;
   default:
      arrayf = NULL;
      dserror("Only 0,1,2,3 allowed for arrayfrom to select sol, sol_increment, sol_residual, sol_mf");
   }
   /*----------------------------------------- select correct arrayfrom */
   switch(arrayto)
   {
   case 0:
      arrayt = &(actnode->sol);
   break;
   case 1:
      arrayt = &(actnode->sol_increment);
   break;
   case 2:
      arrayt = &(actnode->sol_residual);
   break;
   case 3:
      arrayt = &(actnode->sol_mf);
   break;
   default:
      arrayt = NULL;
      dserror("Only 0,1,2,3 allowed for arrayfrom to select sol, sol_increment, sol_residual, sol_mf");
   }
   /* check the size of arrayf */
   if (from >= arrayf->fdim)
      dserror("Cannot copy from array, because place doesn't exist");
   /* check the size of arrayt */
   if (to >= arrayt->fdim)
   {
      diff = to - arrayt->fdim;
      /*max  = IMAX(diff,5);*/
      max = diff;
      amredef(arrayt,arrayt->fdim+max+1,arrayt->sdim,"DA");
   }
   for (j=0;j<arrayf->sdim;j++)
      arrayt->a.da[to][j] = arrayf->a.da[from][j];
   /*----------------------------------------------- transform pressure */
   arrayt->a.da[to][predof]*=dens;
} /* end of loop over nodes */



/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_sol_copy*/

/*!---------------------------------------------------------------------
\brief transform pressure

<pre>                                                         genk 05/02

in this routine the pressure is transformed in the solution history.

  index = 0: actnode->sol
  index = 1: actnode->sol_increment
  index = 2: actnode->sol_residual
  index = 3: actnode->sol_mf

  option = 0: kinematic pressure -> real pressure
  option = 1: real pressure -> kinematic pressure

  (kintematic pressure) = (real pressure)/density

</pre>
\param *actfield      FIELD	     (i)  actual field
\param  disnum        INT	     (i)  number of the discr.
\param  index         INT	     (i)  index of the array
\param  actpos        INT            (i)  position where to transform
\param  predof        INT            (i)  postion of pressure dof
\return void
\sa

------------------------------------------------------------------------*/
void fluid_transpres(
    FIELD             *actfield,
    INT                disnum,
    INT                index,
    INT                actpos,
    INT                predof,
    INT                option
    )
{
INT         i;		/* simply some counters				*/
DOUBLE      dens;	/* density					*/
NODE       *actnode;	/* actual node					*/
ELEMENT    *actele;	/* actual element				*/
DISCRET    *actdis;	/* actual discretisation			*/
ARRAY      *array;	/* pointer to solution array			*/

#ifdef DEBUG
dstrc_enter("fluid_transpres");
#endif

/* since different materials are not allowed  one can work with the
   material parameters of any element ---------------------------------*/
actele = &(actfield->dis[disnum].element[0]);
dens  = mat[actele->mat-1].m.fluid->density;

actdis = &(actfield->dis[disnum]);
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   /*----------------------------------------- select correct arrayfrom */
   switch(index)
   {
   case 0:
      array = &(actnode->sol);
   break;
   case 1:
      array = &(actnode->sol_increment);
   break;
   case 2:
      array = &(actnode->sol_residual);
   break;
   case 3:
      array = &(actnode->sol_mf);
   break;
   default:
      array = NULL;
      dserror("Only 0,1,2,3 allowed for index to select sol, sol_increment, sol_residual, sol_mf");
   }
   dsassert(actpos<array->fdim,"cannot transform pressure\n");
   /*dsassert(predof<array->sdim,"cannot transform pressure\n");*/
   if (predof<array->sdim)
   {
     /*------------------------------------------ pressure transformation */
     switch (option)
     {
     case 0: /*kinematic pressure -> real pressure  */
       array->a.da[actpos][predof]*=dens;
       break;
     case 1: /*real pressure -> kinematic pressure  */
       array->a.da[actpos][predof]/=dens;
       break;
     default:
       dserror("option out of range: don't know what to do!\n");
     }
   }
} /* end of loop over nodes */


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_transpres*/


/*!---------------------------------------------------------------------
\brief steady state check

<pre>                                                         genk 05/02

in this routine the convergence ratios for the steady state check are
calculated and the result is printed to the screen.

</pre>
\param *actfield      FIELD	     (i)  actual field
\param *ipos                         (i)  node array positions
\param  numeq_total   INT	     (i)  total number of equations
\return INT steady

------------------------------------------------------------------------*/
INT fluid_steadycheck(
    FIELD              *actfield,
    INT                 disnum,
    ARRAY_POSITION     *ipos,
    INT                 numeq_total
    )

{

  INT         steady=0;   /* flag for steady state                        */
  DOUBLE      vrat,prat;  /* vel. & pres. ratios                          */
  FILE       *out = allfiles.out_out;


#ifdef DEBUG
  dstrc_enter("fluid_steadycheck");
#endif


  fdyn = alldyn[genprob.numff].fdyn;
  if (fdyn->stnorm==fnst_no) goto end;


  /* determine the conv. ratios */
  fluid_norm(actfield, disnum, numeq_total,&vrat,&prat,ipos);


  /* output to the screen and to .out */
  if (par.myrank==0)
  {
    switch (fdyn->stnorm)
    {
      case fnst_Linf:
        printf("   --> steady state check   (tolerance[norm]):  %10.3E [L_in] \n",
            fdyn->sttol);
        break;
      case fnst_L1:
        printf("   --> steady state check   (tolerance[norm]):  %10.3E [L_1 ] \n",
            fdyn->sttol);
        break;
      case fnst_L2:
        printf("   --> steady state check   (tolerance[norm]):  %10.3E [L_2 ] \n",
            fdyn->sttol);
        break;
      case fnst_no: /* do nothing */
        break;
      default:
        dserror("Norm for steady state check unknwon!\n");
    } /* end switch (fdyn->stnorm)  */

    printf("         velocities: %10.3E	   pressures:   %10.3E  \n",vrat,prat);
    fprintf(out," %10.3E | %10.3E |",vrat,prat);

  } /* endif (par.myrank) */


  /* check if the ratios are smaller than the given tolerance and set flag */
  if (vrat<fdyn->sttol && prat<fdyn->sttol)
  {
    steady=1;
    if (par.myrank==0)
    {
      printf("\n");
      printf("    >>>>>> STEADY STATE REACHED <<<<<< \n");
      printf("\n");
      fprintf(out,">>> STEADY STATE REACHED <<<");
    }
  }

end:


#ifdef DEBUG
  dstrc_exit();
#endif

  return (steady);

} /* end of fluid_steadycheck*/

/*!---------------------------------------------------------------------
\brief iteration convergence check

<pre>                                                         genk 05/02

in this routine the iteration convergence ratios are compared with
the given tolerance. The result is printed out to the screen.

</pre>
\param  vrat          DOUBLE  	     (i)  vel. conv. ratio
\param  prat          DOUBLE         (i)  pres. conv. ratio
\param  grat          DOUBLE         (i)  grid. conv. ratio
\param	itnum 	      INT	     (i)  actual numb. of iter steps
\param	te 	      DOUBLE	     (i)  time for element calcul.
\param	ts	      DOUBLE	     (i)  solver time
\return INT converged

------------------------------------------------------------------------*/
INT fluid_convcheck(
                          DOUBLE             vrat,
		          DOUBLE             prat,
			  DOUBLE             grat,
                          INT                itnum,
		          DOUBLE             te,
		          DOUBLE             ts
		   )
{
FILE       *out = allfiles.out_out;
INT         converged=0;  /* flag for convergence check                  */

#ifdef DEBUG
dstrc_enter("fluid_convcheck");
#endif

fdyn = alldyn[genprob.numff].fdyn;

if (fdyn->itnorm!=fncc_no)
{
   if (par.myrank==0) /* output to the screen */
   {
      switch(fdyn->itnorm)
      {
      case fncc_Linf: /* infinity norm */
        if (fdyn->dyntyp==dyntyp_nln_time_int)
          printf("|  %3d/%3d   | %10.3E[L_in]  | %10.3E   | %10.3E   | {te: %10.3E} {ts:%10.3E} \n",
                 itnum,fdyn->itemax,fdyn->ittol,vrat,prat,te,ts);
        else
          printf("|  %3d/%3d   | %10.3E[L_in]  | %10.3E   | {te: %10.3E} {ts:%10.3E} \n",
                 itnum,fdyn->itemax,fdyn->ittol,vrat,te,ts);
      break;
      case fncc_L1: /* L_1 norm */
        if (fdyn->dyntyp==dyntyp_nln_time_int)
          printf("|  %3d/%3d   | %10.3E[L_1 ]  | %10.3E   | %10.3E   | {te: %10.3E} {ts:%10.3E} \n",
                 itnum,fdyn->itemax,fdyn->ittol,vrat,prat,te,ts);
        else
          printf("|  %3d/%3d   | %10.3E[L_1 ]  | %10.3E   | {te: %10.3E} {ts:%10.3E} \n",
                 itnum,fdyn->itemax,fdyn->ittol,vrat,te,ts);
      break;
      case fncc_L2: /* L_2 norm */
         if (fdyn->freesurf>1)
         printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | {te: %10.3E} {ts:%10.3E} \n",
                 itnum,fdyn->itemax,fdyn->ittol,vrat,prat,grat,te,ts);
         else if (fdyn->dyntyp==dyntyp_nln_time_int)
           printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | {te: %10.3E} {ts:%10.3E} \n",
                  itnum,fdyn->itemax,fdyn->ittol,vrat,prat,te,ts);
         else
           printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | {te: %10.3E} {ts:%10.3E} \n",
                  itnum,fdyn->itemax,fdyn->ittol,vrat,te,ts);
      break;
      default:
         dserror("Norm for nonlin. convergence check unknown!!\n");
      } /*end of switch(fdyn->itnorm) */
   } /* endif (par.myrank) */
   /*------------------------------------------------ convergence check */
   if (vrat<fdyn->ittol && prat<fdyn->ittol && grat<fdyn->ittol)
      converged=2;
   if (itnum==fdyn->itemax)
      converged++;
   if (converged==1 && par.myrank==0)
   {
      if (fdyn->freesurf>1)
      printf("------------------------------------------------------------------------------- \n");
      else
      printf("---------------------------------------------------------------- \n");
      printf("|          >>>>>> not converged in itemax steps!               | \n");
      printf("---------------------------------------------------------------- \n");
      fprintf(out," %3d | %10.3E | %10.3E |",
          itnum,vrat,prat);
      fprintf(out," not converged in itemax steps\n");
   }
   else if (converged>0 && par.myrank==0)
   {
      if (fdyn->freesurf>1)
      printf("------------------------------------------------------------------------------- \n");
      else if (fdyn->dyntyp==dyntyp_nln_time_int)
        printf("----------------------------------------------------------------\n");
      else
        printf("-------------------------------------------------\n");
      printf("\n");
      fprintf(out," %3d | %10.3E | %10.3E |",
          itnum,vrat,prat);
   }
} /* endif (fdyn->itchk) */
else
{
   if (par.myrank==0)
   {
      printf("      iteration step: %3d / %3d \n",
                              itnum, fdyn->itemax) ;
      if(itnum==fdyn->itemax)
      {
	 printf("----------------------------------------------------------------- \n");
         printf("\n");
      }
   }
} /* endif (par.myrank) */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return (converged);
} /* end of fluid_convcheck*/

/*!---------------------------------------------------------------------
\brief print out to the screen

<pre>                                                         genk 05/02

time-integration parameters are printed out to the screen

</pre>
\return void

------------------------------------------------------------------------*/
void fluid_algoout( void )
{

#ifdef DEBUG
dstrc_enter("fluid_algoout");
#endif

fdyn = alldyn[genprob.numff].fdyn;

printf("\n");

switch(fdyn->iop)
{
case 1:
   printf("TIME: %11.4E/%11.4E  DT = %11.4E  Generalised Alpha  STEP = %4d/%4d \n",
          fdyn->acttime,fdyn->maxtime,fdyn->dta,fdyn->step,fdyn->nstep);
break;
case 4:
   printf("TIME: %11.4E/%11.4E  DT = %11.4E  One-Step-Theta  STEP = %4d/%4d \n",
          fdyn->acttime,fdyn->maxtime,fdyn->dta,fdyn->step,fdyn->nstep);
break;
case 7:
   printf("TIME: %11.4E/%11.4E  DT = %11.4E     BDF2         STEP = %4d/%4d \n",
          fdyn->acttime,fdyn->maxtime,fdyn->dta,fdyn->step,fdyn->nstep);
break;
#ifdef D_FLUID2_TDS
case 8:
   printf("TIME: %11.4E/%11.4E  DT = %11.4E incr. Gen-Alpha  STEP = %4d/%4d \n",
          fdyn->acttime,fdyn->maxtime,fdyn->dta,fdyn->step,fdyn->nstep);
break;
#endif /* D_FLUID2_TDS */
default:
   dserror("parameter out of range: IOP\n");
} /* end of swtich(fdyn->iop) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_algoout*/
/*!---------------------------------------------------------------------
\brief reduce shearstresses in parallel case

<pre>                                                       he 03/03

shearstress is calcuted within a parellel loop over the elements.
So the result is not available on all procs. Here they are reduced by
MPI to all procs!

</pre>
\param *actintra          INTRA         (i)
\param *actfield          FIELD         (i)   the actual field
\return void

------------------------------------------------------------------------*/
void fluid_reduceshstr(INTRA             *actintra,
                         FIELD             *actfield)
{
#ifdef PARALLEL
INT      numnp_total;
INT      i;
NODE    *actnode;


#ifdef DEBUG
dstrc_enter("fluid_reduceshstr");
#endif

fdyn = alldyn[genprob.numff].fdyn;

/*------------------------------------------- get total number of nodes */
numnp_total = actfield->dis[0].numnp;
for (i=0; i<numnp_total; i++) /* loop nodes */
{
   actnode = &(actfield->dis[0].node[i]);
   MPI_Bcast(&(actnode->fluid_varia->c_f_shear),1,MPI_DOUBLE,actnode->proc,
             actintra->MPI_INTRA_COMM);
/*------------------- compute shearvelocity for the scaned coordinates */
   if (FABS(actnode->x[0]-fdyn->coord_scale[0])<EPS7 &&
       FABS(actnode->x[1]-fdyn->coord_scale[1])<EPS15)
   MPI_Bcast(&(fdyn->washvel),1,MPI_DOUBLE,actnode->proc,
             actintra->MPI_INTRA_COMM);
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

#endif
return;
} /* end of fluid_reduceshstr*/
/*!---------------------------------------------------------------------
\brief reduce shearstresses in parallel case

<pre>                                                       he 03/03

shearstress is calcuted within a parellel loop over the elements.
So the result is not available on all procs. Here they are reduced by
MPI to all procs!

</pre>
\param *actintra         INTRA           (i)
\param *actfield	       FIELD           (i)    the actual field
\return void

------------------------------------------------------------------------*/
void fluid_nullshstr(INTRA             *actintra,
                        PARTITION         *actpart,
                        FIELD             *actfield)
{
#ifdef PARALLEL
INT      numnp_total;
INT      i;
INT      myrank,imyrank;
NODE    *actnode;


#ifdef DEBUG
dstrc_enter("fluid_nullshstr");
#endif

fdyn = alldyn[genprob.numff].fdyn;

/*------------------------------------------- get total number of nodes */
myrank = par.myrank;
imyrank= actintra->intra_rank;
numnp_total = actfield->dis[0].numnp;

if (imyrank!=0 || myrank!=0)
for (i=0; i<actfield->dis[0].numnp; i++)
{
 actnode = &(actfield->dis[0].node[i]);
 actnode->fluid_varia->c_f_shear = ZERO;
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

#endif
return;
} /* end of fluid_nullshstr*/
/*!---------------------------------------------------------------------
\brief reduce fluid stresses

<pre>                                                         genk 01/03

fluid stresses are calcuted within a parellel loop over the elements.
So the results are not available on all procs. Here they are reduced by
MPI to all procs!

</pre>
\param *actintra         INTRA           (i)
\param *actfield	 FIELD           (i)    the actual field
\param *actsolv          SOLVAR          (i)    the actual solver variables
\param  numdf            INT             (i)    number of fluid dofs
\param  flag             INT             (i)    a flag
\return void

------------------------------------------------------------------------*/
void fluid_reducestress(
    INTRA              *actintra,
    PARTITION          *actpart,
    FIELD              *actfield,
    INT                 disnum,
    INT                 numdf,
    FLUID_STRESS        str
    )

{

#ifdef PARALLEL
  INT      numele_total,numnp;
  INT      i,j, coupled;
  static INT dim;
  INT      start;
  INT      mone=-1;
  static ARRAY    index_a;
  static INT     *index;
  static DOUBLE  *buffer;
  static ARRAY    buffer_a;
  static DOUBLE  *recvbuffer;
  static ARRAY    recvbuffer_a;
  ELEMENT *actele;
  GNODE   *actgnode;

#ifdef DEBUG
  dstrc_enter("fluid_reducestress");
#endif

  numele_total = actfield->dis[disnum].numele;

  if (index_a.Typ!= cca_IV)
  {
    index=amdef("index",&index_a,numele_total,1,"IV");
    aminit(&index_a,&mone);
    switch(str)
    {
#ifdef D_FSI
      case str_fsicoupling: /*- sent only if element has fsi-coupled nodes */
        coupled=0;
        for(i=0;i<numele_total;i++) /* loop elements */
        {
          actele=&(actfield->dis[disnum].element[i]);
          numnp = actele->numnp;
          /*----- loop nodes and check if element has a fsi-coupled node */
          for (j=0;j<numnp;j++)
          {
            actgnode = actele->node[j]->gnode;
            /* check if there is a coupled struct node */
            if (actgnode->mfcpnode[genprob.numsf]==NULL) continue;
            index[actele->Id_loc]=coupled;
            coupled++;
            break;
          } /* end of loop over nodes */
        }
        break;
#endif
      default:
        dserror("str not valid for reduction\n");
    } /* end switch */

    /*------------------------------------------ allocate stress buffers */
    if (numdf==3)
    {
      dim=MAXNOD*coupled*3;
      buffer     = amdef("buffer",&buffer_a,MAXNOD*coupled*3,1,"DV");
      recvbuffer = amdef("recvbuffer",&recvbuffer_a,dim,1,"DV");
    }
    else if (numdf==4)
    {
      dim=MAXNOD*coupled*6;
      buffer     = amdef("buffer",&buffer_a,MAXNOD*coupled*6,1,"DV");
      recvbuffer = amdef("recvbuffer",&recvbuffer_a,dim,1,"DV");
    }
    else
      dserror("numdf out of range!\n");
  } /* endif */

  amzero(&buffer_a);
  amzero(&recvbuffer_a);

  for (i=0; i<actpart->pdis[disnum].numele; i++)
  {
    actele = actpart->pdis[disnum].element[i];
    coupled=index[actele->Id_loc];
    /* continue if element belongs to other proc or is not coupled */
    if (par.myrank!=actele->proc || coupled==-1) continue;
    numnp=actele->numnp;

    /*------------------------------------------- set buffer and dim */
#ifdef D_FLUID2
    if (actele->eltyp == el_fluid2)
    {
      start=coupled*MAXNOD*3;
      for (j=0;j<numnp;j++)
      {
        buffer[start]  =actele->e.f2->stress_ND.a.da[j][0];
        buffer[start+1]=actele->e.f2->stress_ND.a.da[j][1];
        buffer[start+2]=actele->e.f2->stress_ND.a.da[j][2];
        start+=3;
      }
    }
#endif

#ifdef D_FLUID3
    if (actele->eltyp == el_fluid3)
    {
      start=coupled*MAXNOD*6;
      for (j=0;j<numnp;j++)
      {
        buffer[start]  =actele->e.f3->stress_ND.a.da[j][0];
        buffer[start+1]=actele->e.f3->stress_ND.a.da[j][1];
        buffer[start+2]=actele->e.f3->stress_ND.a.da[j][2];
        buffer[start+3]=actele->e.f3->stress_ND.a.da[j][3];
        buffer[start+4]=actele->e.f3->stress_ND.a.da[j][4];
        buffer[start+5]=actele->e.f3->stress_ND.a.da[j][5];
        start+=6;
      }
    }
#endif

#ifdef D_FLUID3_F
    if (actele->eltyp == el_fluid3_fast)
    {
      start=coupled*MAXNOD*6;
      for (j=0;j<numnp;j++)
      {
        buffer[start]  =actele->e.f3->stress_ND.a.da[j][0];
        buffer[start+1]=actele->e.f3->stress_ND.a.da[j][1];
        buffer[start+2]=actele->e.f3->stress_ND.a.da[j][2];
        buffer[start+3]=actele->e.f3->stress_ND.a.da[j][3];
        buffer[start+4]=actele->e.f3->stress_ND.a.da[j][4];
        buffer[start+5]=actele->e.f3->stress_ND.a.da[j][5];
        start+=6;
      }
    }
#endif

  } /* end of loop over all elements */

  /*------------------------------------------------------ allreduce the vector */
  MPI_Allreduce(buffer,recvbuffer,dim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);

  /* now recbuffer contains all stresses, which have to be distrubuted to the
     elements */
  for(i=0;i<numele_total;i++) /* loop elements */
  {
    actele=&(actfield->dis[disnum].element[i]);
    coupled=index[actele->Id_loc];
    if (coupled==-1) continue; /* no stress reduction */
    numnp=actele->numnp;
#ifdef D_FLUID2
    if (actele->eltyp == el_fluid2)
    {
      start=coupled*MAXNOD*3;
      for (j=0;j<numnp;j++)
      {
        actele->e.f2->stress_ND.a.da[j][0]=recvbuffer[start]  ;
        actele->e.f2->stress_ND.a.da[j][1]=recvbuffer[start+1];
        actele->e.f2->stress_ND.a.da[j][2]=recvbuffer[start+2];
        start+=3;
      }
    }
#endif

#ifdef D_FLUID3
    if (actele->eltyp == el_fluid3)
    {
      start=coupled*MAXNOD*6;
      for (j=0;j<numnp;j++)
      {
        actele->e.f3->stress_ND.a.da[j][0]=recvbuffer[start]  ;
        actele->e.f3->stress_ND.a.da[j][1]=recvbuffer[start+1];
        actele->e.f3->stress_ND.a.da[j][2]=recvbuffer[start+2];
        actele->e.f3->stress_ND.a.da[j][3]=recvbuffer[start+3];
        actele->e.f3->stress_ND.a.da[j][4]=recvbuffer[start+4];
        actele->e.f3->stress_ND.a.da[j][5]=recvbuffer[start+5];
        start+=6;
      }
    }
#endif

#ifdef D_FLUID3_F
    if (actele->eltyp == el_fluid3_fast)
    {
      start=coupled*MAXNOD*6;
      for (j=0;j<numnp;j++)
      {
        actele->e.f3->stress_ND.a.da[j][0]=recvbuffer[start]  ;
        actele->e.f3->stress_ND.a.da[j][1]=recvbuffer[start+1];
        actele->e.f3->stress_ND.a.da[j][2]=recvbuffer[start+2];
        actele->e.f3->stress_ND.a.da[j][3]=recvbuffer[start+3];
        actele->e.f3->stress_ND.a.da[j][4]=recvbuffer[start+4];
        actele->e.f3->stress_ND.a.da[j][5]=recvbuffer[start+5];
        start+=6;
      }
    }
#endif

  } /* end of loop over all elements */



#ifdef DEBUG
  dstrc_exit();
#endif


#endif  /* ifdef PARALLEL */


  return;

} /* end of fluid_reducestress*/









/*!---------------------------------------------------------------------
\brief calculating errors for beltrami and kim-moin

<pre>                                                         mn 02/04

In this routine different velocity and pressure error norms for the
beltramii and kim-moin flow are calculated.

</pre>
\param *actfield      FIELD          (i)   actual field
\param *ipos                         (i)   node array positions
\param  index         INT            (i)   index for type of flow
                                           index = 1: beltrami
                                           index = 2: kim-moin
\return void

------------------------------------------------------------------------*/
void fluid_cal_error(
    FIELD             *actfield,
    INT                disnum,
    ARRAY_POSITION    *ipos,
    INT                index
    )
{
  INT         i;            /* simply some counters */
  INT         numdf;          /* number of fluid dofs */
  INT         numvel;         /* total number of vel-dofs */
  INT         predof;         /* actual number of pres dof */
  INT         numnp_total;    /* total number of fluid nodes */

  DOUBLE      dvinorm=ZERO;   /* infinity norms */
  DOUBLE       vinorm=ZERO;   /* infinity norms */
  DOUBLE      dpinorm=ZERO;   /* infinity norms */
  DOUBLE       pinorm=ZERO;   /* infinity norms */


  DOUBLE      dv1norm=ZERO;   /* L1 norms */
  DOUBLE       v1norm=ZERO;   /* L1 norms */
  DOUBLE      dp1norm=ZERO;   /* L1 norms */
  DOUBLE       p1norm=ZERO;   /* L1 norms */

  DOUBLE      dv2norm=ZERO;   /* L2 norms */
  DOUBLE       v2norm=ZERO;   /* L2 norms */
  DOUBLE      dp2norm=ZERO;   /* L2 norms */
  DOUBLE       p2norm=ZERO;   /* L2 norms */

  SOLVAR      *actsolv;       /* pointer to active sol. structure */
  PARTITION   *actpart;       /* pointer to active partition      */
  INTRA       *actintra;      /* pointer to active intra-communic.*/

  CONTAINER    container;     /* contains variables defined in container.h */
  CALC_ACTION *action;        /* pointer to the cal_action enum   */

  INT            numeq_total;
  FLUID_DYNAMIC *fdyn;                /* pointer to fluid dyn. inp.data   */



#ifdef DEBUG
  dstrc_enter("fluid_cal_error");
#endif

  /*---------------------------------------------------- set some values */
#ifdef PARALLEL
  actintra    = &(par.intra[0]);
#else
  actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
  actintra->intra_fieldtyp = fluid;
  actintra->intra_rank     = 0;
  actintra->intra_nprocs   = 1;
#endif

  action      = &(calc_action[0]);
  *action     = calc_struct_none;

  actsolv     = &(solv[0]);
  actpart     = &(partition[0]);

  /* Zero everything so we can rely on the values we find. */
  memset(&container, 0, sizeof(CONTAINER));

  container.disnum   = disnum;
  container.fieldtyp = fluid;

  fdyn = alldyn[genprob.numff].fdyn;
  numdf        = fdyn->numdf;
  numnp_total  = actfield->dis[disnum].numnp;
  numeq_total  = actfield->dis[disnum].numeq;
  predof       = numdf-1;
  numvel       = numdf-1;


  /* new error calculation, or better integration */
  for (i=0;i<=2;i++)
  {
    container.vel_error  = 0.0;
    container.pre_error  = 0.0;
    container.vel_norm   = 0.0;
    container.pre_norm   = 0.0;
    container.error_norm = i;

    *action = calc_fluid_error;
    calelm(actfield,actsolv,actpart,actintra,-1,-1,
        &container,action);


    if (i==2)
    {
      container.vel_error = sqrt(container.vel_error);
      container.pre_error = sqrt(container.pre_error);
    }

    switch (i)
    {
      case 0:
        dvinorm = container.vel_error;
        dpinorm = container.pre_error;
        vinorm  = container.vel_norm;
        pinorm  = container.pre_norm;
        break;

      case 1:
        dv1norm = container.vel_error;
        dp1norm = container.pre_error;
        v1norm  = container.vel_norm;
        p1norm  = container.pre_norm;
        break;

      case 2:
        dv2norm = container.vel_error;
        dp2norm = container.pre_error;
        v2norm  = container.vel_norm;
        p2norm  = container.pre_norm;
        break;

      default:
        dserror("Error norm %2i not possible!!", i);
        break;

    }  /* switch (i) */

  }  /* for (i=0;i<=2;i++) */


  switch (index)
  {
    case 1: /* BELTRAMI */
      printf("\nErrors for Beltrami flow:\n");
      break;
      /* end of BELTRAMI */



    case 2: /* KIM-MOIN */
      printf("\nErrors for Kim-Moin flow:\n");
      break;
      /* end of KIM-MOIN */

    default:
      dserror("Unknown type of flow for error calculation");
      break;
  }

  printf("norm | abs. vel.  | rel. vel.  | abs. pre.  | rel. pre.  |\n");
  printf("----------------------------------------------------------\n");
  printf("inf  |%11.4E |%11.4E |%11.4E |%11.4E |\n",dvinorm,dvinorm/vinorm,dpinorm,dpinorm/pinorm);
  printf("L1   |%11.4E |%11.4E |%11.4E |%11.4E |\n",dv1norm,dv1norm/v1norm,dp1norm,dp1norm/p1norm);
  printf("L2   |%11.4E |%11.4E |%11.4E |%11.4E |\n\n",dv2norm,dv2norm/v2norm,dp2norm,dp2norm/p2norm);


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of fluid_cal_error */



/*!---------------------------------------------------------------------
\brief init positions in sol_increment

<pre>                                                        chfoe 11/04

This routine inits the positions in sol_increment in the case of a pure
Eulerian problem.

\param  *ipos                           (i)   node array positions
</pre>
\return void

------------------------------------------------------------------------*/
void fluid_init_pos_euler(ARRAY_POSITION    *ipos)
{
#ifdef DEBUG
  dstrc_enter("fluid_init_pos_euler");
#endif

fdyn = alldyn[genprob.numff].fdyn;

if (fdyn->adaptive)
{
   /*---------------------------------------- adaptive time stepping ---*/
   ipos->velnm = 0;
   ipos->veln  = 1;
   ipos->hist  = 2;
   ipos->velnp = 3;
   ipos->accnm = 4;
   ipos->accn  = 5;
   ipos->pred  = 6;
   ipos->terr  = 7;

   ipos->numincr = 8;
}
else
{
   /*------------------------------------------ fixed time step size ---*/
   ipos->velnp = 0;  /* most recent solution */
   ipos->veln  = 1;  /* previous soulution (at time n) */
   ipos->velnm = 2;  /* last historical solution (at time n-1) */
   ipos->accnm = 3;  /* acceleration at time n-1 */
   ipos->accn  = 4;  /* acceleration at time n */
   ipos->hist  = 5;  /* linear combination of old solutions */

   ipos->numincr = 6;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif


/*! @} (documentation module close)*/
