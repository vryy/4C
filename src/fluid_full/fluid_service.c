#include "../headers/standardtypes.h"
#include "../headers/solution.h"
#include "fluid_prototypes.h"
#define ONE (1.0) 
#define ZERO (0.0)
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
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
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | ranks and communicators                                              |
 | This structure struct _PAR par; is defined in main_ccarat.c
 *----------------------------------------------------------------------*/
 extern struct _PAR   par; 
/*----------------------------------------------------------------------*
 | routine to get the split (vel - pre) dof numbers                     |
 *----------------------------------------------------------------------*/
/* void fluid_splitdof(FLUID_DYNAMIC *fdyn)
{
int     numdf;              /* number of dofs in this discretisation */
/*int     numvel;             /* number of vel-dofs */
/* int     predof;             /* number of pressure dof in dof / spldof field */
/*int     j,k;
int     counter;
int     numnp_total;
NODE   *actnode;           /* the actual node */
/*FIELD  *actfield;          /* the actual field */
/*
#ifdef DEBUG 
dstrc_enter("fluid_splitdof");
#endif

/*------------------------------------------------ set actfield = fluid */
/*actfield = &(field[0]);

/*--------------------------- for fluid problems numdf is either 3 or 4 */
/*numdf=fdyn->numdf;
numvel=numdf-1;
predof=numdf-1;
numnp_total=actfield->dis[0].numnp;
fdyn->dynvar.npeqmax=numnp_total;
fdyn->dynvar.nveqmax=(numdf-1)*numnp_total;

/*---------------------------------- loop all nodes and allocate spldof */
/*for (j=0; j<actfield->dis[0].numnp; j++)
{
   actfield->dis[0].node[j].spldof  = (int*)CALLOC(numdf,sizeof(int));
   if (!(actfield->dis[0].node[j].spldof)) 
         dserror("Allocation of dof in NODE failed");
}
/*---------------------------- loop all nodes and assign split pre-dofs */
/*counter=0;
for (j=0; j<actfield->dis[0].numnp; j++)
{
   actnode = &(actfield->dis[0].node[j]);
   actnode->spldof[predof] = counter;
   counter++;
}
/*---------------------------- loop all nodes and assign split vel-dofs */
/*counter=0;
for (j=0; j<actfield->dis[0].numnp; j++)
{
   actnode = &(actfield->dis[0].node[j]);
   for (k=0; k<numvel; k++)
   {   
      actnode->spldof[k] = counter;
      counter++;
   }
}

/*----------------------------------------------------------------------*/
/*#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_splitodf*/ 

/*----------------------------------------------------------------------*
 | routine to check starting algorithm                      genk 04/02  |
 *----------------------------------------------------------------------*/
void fluid_startproc(
                     FLUID_DYNAMIC  *fdyn,
		     int            *nfrastep 
		    )
{
static int     iopfsi_s;
static int     iopfss_s;
static int     itemax_s;
static double  theta_s ;
static double  thetas_s;

 #ifdef DEBUG 
dstrc_enter("fluid_startproc");
#endif 

if (fdyn->step==1)
{
   iopfsi_s = fdyn->iopfsi;
   iopfss_s = fdyn->iopfss;
   itemax_s = fdyn->itemax;
   theta_s  = fdyn->theta;
   thetas_s = fdyn->thetas;
}
if (fdyn->step<=fdyn->numfss)
{
   fdyn->iopfsi = iopfss_s;
   fdyn->theta  = theta_s;
   *nfrastep    = 1;
   if (iopfss_s==5)
      *nfrastep = 3;
   if (iopfss_s==2 || iopfss_s==3)
      fdyn->itemax = 1;
   if (fdyn->step==fdyn->numfss && iopfsi_s==3)
       dserror ("starting algo for semi-impl. two-step method not implemented yet");
       /* set U(n-1) in last step of start-algo */
}
else if (fdyn->step==(fdyn->numfss+1))
{
   fdyn->iopfsi = iopfsi_s;
   fdyn->theta  = theta_s;
   *nfrastep    = 1;
   if (iopfsi_s==5)
      *nfrastep = 3;
   if (iopfsi_s==2 || iopfsi_s==3)
      fdyn->itemax = 1;
   else
      fdyn->itemax = itemax_s;
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_startproc*/

/*----------------------------------------------------------------------*
 | routine to calculate time integration constants                      |
 | for fluid time algorithms                                genk  03/02 |
 *----------------------------------------------------------------------*/
void fluid_tcons(FLUID_DYNAMIC *fdyn,
                 FLUID_DYN_CALC *dynvar)
{
 #ifdef DEBUG 
dstrc_enter("fluid_tcons");
#endif
/*----------------------------------------------------- check algorithm */
if (fdyn->iopfsi==4)   /* one step theta */
{
    dynvar->dta  = fdyn->dt;
    dynvar->thsl = fdyn->dt*fdyn->theta;
    dynvar->thpl = dynvar->thsl;
    dynvar->thsr = (ONE - fdyn->theta)*fdyn->dt;
    dynvar->thpr = dynvar->thsr;
}
else
   dserror ("constants for time algorithm not implemented yet!");
/*----------------------------------------------- treatment of pressure */ 
if (fdyn->iprerhs!=1)
    dserror ("treatment of pressure not implemented yet!");

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_tcons*/ 

/*----------------------------------------------------------------------*
 | routine to calculate constants for nonlinear iteration   genk  03/02 |
 |									|
 |   nik <-> 'K' --> EVALUATION OF LHS-MATRICES (w/o NONLINEAR TERM)    |
 |   nic <-> 'C' --> EVALUATION OF NONLINEAR LHS N-CONVECTIVE	        |
 |   nir <-> 'R' --> EVALUATION OF NONLINEAR LHS N-REACTION	        |
 |   nie <-> 'E' --> EVALUATE ONLY LHS-TERMS FOR EXPLICIT VELOCITY      |
 |   nil <-> 'L' --> EVALUATION OF LUMPED MASS MATRIX (Mvv-lumped)      |
 |   nif <-> 'F' --> EVALUATION OF "TIME - RHS" (F-hat)  	        |
 |   nii <-> 'I' --> EVALUATION OF "ITERATION - RHS"		        |
 |   nis <-> 'S' --> STATIONARY CASE (NO TIMEDEPENDENT TERMS)	        |
 |									|
 *----------------------------------------------------------------------*/
void fluid_icons(FLUID_DYNAMIC *fdyn,
                FLUID_DYN_CALC *dynvar,
		int itnum)
{
int i;

 #ifdef DEBUG 
dstrc_enter("fluid_icons");
#endif

/*----------------------------------------------------- initialisation */
dynvar->nik=0;
dynvar->nic=0;
dynvar->nir=0;
dynvar->nie=0;
dynvar->nil=0;
dynvar->nif=0;
dynvar->nii=0;
dynvar->nis=0;

if(fdyn->ite==0)           /* no iteration */
{
   dynvar->sigma=ZERO;
   dynvar->nik=1;       
   dynvar->nic=2; 
   dynvar->nif=3;      /* KCF */
}
else if (fdyn->ite==1)    /* fixed point like iteration */
{
   dynvar->sigma=ZERO;
   if (itnum>1)
   {
      dynvar->nik=1;  
      dynvar->nic=2;   /* KC */  
   }
   else
   {
      dynvar->nik=1;
      dynvar->nic=2;
      dynvar->nif=3;   /* KCF */
   }
}
else if (fdyn->ite==2)    /* Newton iteration */
{
   dynvar->sigma=ONE;
   if (itnum>1)
   {
      dynvar->nik=1;
      dynvar->nic=2;
      dynvar->nir=3;
      dynvar->nii=4;  /* KCRI */
   }
   else
   {
      dynvar->nik=1;
      dynvar->nic=2;
      dynvar->nir=3;
      dynvar->nif=4;
      dynvar->nii=5;  /* KCRFI */
   }
}
else if (fdyn->ite==3)    /* fixed point iteration */
{
   dynvar->sigma=-ONE;
   if (itnum>1)
   {
      dynvar->nik=1;
      dynvar->nii=2;  /* KI */     
   }
   else
   {
      dynvar->nik=1;
      dynvar->nif=2;      
      dynvar->nii=3;  /* KFI */      
   }   
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_icons*/

/*----------------------------------------------------------------------*
 | routine to initialise solution history                    genk 04/02 |
 *----------------------------------------------------------------------*/
void fluid_init(
		FIELD  *actfield,  
                FLUID_DYNAMIC *fdyn	
	       )      
{
int    i,j;
int    actmat;
int    numdf;        /* number of dofs in this discretisation */
int    numnp_total;  /* total number of nodes in this discretisation */
int    numele_total; /* total number of elements in this discr. */  
int    predof;       /* pressure dof number */
double dens;
NODE  *actnode;      /* the actual node */
ELEMENT * actele;    /* the actual element */

#ifdef DEBUG 
dstrc_enter("fluid_init");
#endif

/*----------------------- set control variables for element evaluation */
fdyn->dynvar.itwost = 0;
fdyn->dynvar.isemim = 0;
fdyn->dynvar.ishape = 1;
fdyn->dynvar.iprerhs= fdyn->iprerhs;

if(fdyn->iopfsi==3) 
   fdyn->dynvar.itwost = 1;
if(fdyn->iopfsi==2 || fdyn->iopfsi==3)   
  fdyn->dynvar.isemim = 1;  
  
/*---------------------------------------------------- set some values */
numdf        = fdyn->numdf;
numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;
predof       = numdf-1;

/*-------------------------------- allocate space for solution history *
 | node->sol : solutions for pss-file and output but not used during   |
 |             calculations!!!                                         |
 | node->sol_incement: solution history used for calculations          |
 |       sol_increment[0][i]: solution at (n-1)                        |
 |	 sol_increment[1][i]: solution at (n)                          |
 |	 sol_increment[2][i]: solution at (n+g)                        |
 |	 sol_increment[3][i]: solution at (n+1)                        |
 *---------------------------------------------------------------------*/
for (i=0;i<numnp_total;i++)
{
   actnode=&(actfield->dis[0].node[i]);
   amredef(&(actnode->sol_increment),4,numdf,"DA");
   amzero(&(actnode->sol_increment));
}


if (fdyn->init>=1)
{
/*------------------------------ copy initial data to solution history */
   for (i=0;i<numnp_total;i++)
   {
      actnode=&(actfield->dis[0].node[i]);
      for (j=0;j<numdf;j++)
      {      
         actnode->sol.a.da[0][j]          =fdyn->start.a.da[actnode->Id_loc][j];
	 actnode->sol_increment.a.da[1][j]=fdyn->start.a.da[actnode->Id_loc][j];
	 actnode->sol_increment.a.da[3][j]=fdyn->start.a.da[actnode->Id_loc][j];
      }
   }
/*--------------------------- transform pressure values of initial data 
             from real pressure to kinematic pressure ------------------*/
    for (i=0;i<numele_total;i++)
    {
       actele = &(actfield->dis[0].element[i]);
       actmat = actele->mat-1;
       dens   = mat[actmat].m.fluid->density;
       for(j=0;j<actele->numnp;j++)
       {
          actnode=actele->node[j];
	  actnode->sol_increment.a.da[1][predof] /= dens;
	  actnode->sol_increment.a.da[3][predof] /= dens;
       }
    }
                               

/*--------------------- free starting field and initial solution vector */
   amdel(&(fdyn->start));
}


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_init*/ 


/*----------------------------------------------------------------------*
 | routine to determine the absolute maximum velocity for use as global |
 | scaling velocity in some stability parameter definitions (Tau_mp)    |
 |                                                          genk  03/02 |
 | at the moment not used and only sequentiell version                  |
 *----------------------------------------------------------------------*/
void fluid_maxvel(
                  FLUID_DYNAMIC *fdyn, 
                  FIELD  *actfield,               		 
		  int actpos,
		  double *res
		 )      
{
int    i,j;
int    numvel;       /* number of vel dofs */
int    numdf;        /* number of dofs in this discretisation */
int    numnp_total;  /* total number of nodes in this discretisation */
int    actdof;       /* actual dofnumber */     
double dm=0.0;
NODE  *actnode;      /* the actual node */

#ifdef DEBUG 
dstrc_enter("fluid_maxvel");
#endif

numdf=fdyn->numdf;
numvel=numdf-1;
numnp_total=actfield->dis[0].numnp;
/*------------------------------------ loop over total number of nodes */
/* --------------------- its possible to parallise this loop --------- */
for (i=0;i<numnp_total;i++)
{
   actnode=&(actfield->dis[0].node[i]);
   for (j=0;j<numvel;j++)
   {
      actdof =actnode->dof[j];             
      dm=DMAX(dm,FABS(actnode->sol.a.da[actpos][j]));
   }
}

*res=dm;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_maxvel*/ 



/*----------------------------------------------------------------------*
 | routine to extract digits from integer number                        |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void intextract(
                int num,    /* integer number */
                int *it,    /* integer on position "thousand" */
		int *ih,    /* integer on position "hundred" */
		int *id,    /* integer on position "ten" */
		int *io     /* integer on position "one" */
	       )
{
int nit, nih, nid, nio;
#ifdef DEBUG 
dstrc_enter("intextract");
#endif

nit = num/1000;
nih = (num-nit*1000)/100;
nid = (num-nit*1000-nih*100)/10;
nio = num -nit*1000-nih*100-nid*10;

*it=nit;
*ih=nih;
*id=nid;
*io=nio;

 /*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of intextract*/ 
   
/*----------------------------------------------------------------------*
 |  Put the results of a DIST_VECTOR to the nodes in a       m.gee 11/01|
 |  certain place in ARRAY sol_increment                                |
 |  Result have to bee allreduced and are put to the whole              |
 |  field on each proc                                                  |
 |  Calculates the norms for the iteration check if necessary           |
 |                                                          genk 05/02  |
 *----------------------------------------------------------------------*/
void fluid_result_incre(FIELD *actfield,INTRA *actintra,DIST_VECTOR *sol,
                          int place,SPARSE_ARRAY *sysarray,SPARSE_TYP *sysarray_typ,
			  double *vrat, double *prat, FLUID_DYNAMIC *fdyn )
{
int      i,j;
int      max;
int      diff;
int      dof;
int      predof;

int      numeq_total;
NODE    *actnode;
ARRAY    result_a;
double  *result;

double   dvnorm=ZERO;
double   dpnorm=ZERO;
double    vnorm=ZERO;
double    pnorm=ZERO;

#ifdef DEBUG 
dstrc_enter("fluid_result_incre");
#endif
/*----------------------------------------------------------------------*/
numeq_total = sol->numeq_total;
predof      = fdyn->numdf-1;
/*------------------------- allocate space to allreduce the DIST_VECTOR */
result = amdef("result",&result_a,numeq_total,1,"DV");
         amzero(&result_a);
solserv_reddistvec(
                      sol,
                      sysarray,
                      sysarray_typ,
                      result,
                      sol->numeq_total,
                      actintra
                     );
switch (fdyn->itchk)
{
case 0: /* no convergence check; */
   /*-----------  loop nodes and put the result back to the node structure */
   for (i=0; i<actfield->dis[0].numnp; i++)
   {
      actnode = &(actfield->dis[0].node[i]);
      /*------------------------------ enlarge sol_increment, if necessary */
      if (place >= actnode->sol_increment.fdim)
      {
         diff = place - actnode->sol_increment.fdim;
         max  = IMAX(diff,5);
         amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,actnode->sol_increment.sdim,"DA");
      }
      for (j=0; j<actnode->numdf; j++)
      {
         dof = actnode->dof[j];
         if (dof>=numeq_total) continue;
         actnode->sol_increment.a.da[place][j] = result[dof];
      }   
     
   }
   break;
case 1: /* convergence check  */
   switch (fdyn->itnorm)
   {
   case 0: /* L_infinity norm */
      /*-----------  loop nodes and put the result back to the node structure */
      for (i=0; i<actfield->dis[0].numnp; i++)
      {
         actnode = &(actfield->dis[0].node[i]);
         /*------------------------------ enlarge sol_increment, if necessary */
         if (place >= actnode->sol_increment.fdim)
         {
            diff = place - actnode->sol_increment.fdim;
            max  = IMAX(diff,5);
            amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,actnode->sol_increment.sdim,"DA");
         }
         for (j=0; j<actnode->numdf; j++) /* loop dofs and calculate the norms */
         {
	    dof = actnode->dof[j];
            if (dof>=numeq_total) continue;
            if (j==predof) /* pressure dof */
	    {
               dpnorm = DMAX(dpnorm, FABS(result[dof]-actnode->sol_increment.a.da[place][j]));
	       pnorm  = DMAX(pnorm, FABS(result[dof]));
	       actnode->sol_increment.a.da[place][j] = result[dof];
	    }
	    else /* vel - dof */
	    {	       
	       dvnorm = DMAX(dvnorm, FABS(result[dof]-actnode->sol_increment.a.da[place][j]));
	       vnorm  = DMAX(vnorm, FABS(result[dof]));
               actnode->sol_increment.a.da[place][j] = result[dof];
	    }
         }        
      }      
      break;
   case 1: /* L_1 norm */
      /*-----------  loop nodes and put the result back to the node structure */
      for (i=0; i<actfield->dis[0].numnp; i++)
      {
         actnode = &(actfield->dis[0].node[i]);
         /*------------------------------ enlarge sol_increment, if necessary */
         if (place >= actnode->sol_increment.fdim)
         {
            diff = place - actnode->sol_increment.fdim;
            max  = IMAX(diff,5);
            amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,actnode->sol_increment.sdim,"DA");
         }
         for (j=0; j<actnode->numdf; j++) /* loop dofs and calculate the norms */
         {
	    dof = actnode->dof[j];
            if (dof>=numeq_total) continue;
            if (j==predof) /* pressure dof */
	    {
               dpnorm += FABS(result[dof]-actnode->sol_increment.a.da[place][j]);
	       pnorm  += FABS(result[dof]);
	       actnode->sol_increment.a.da[place][j] = result[dof];
	    }
	    else /* vel - dof */
	    {	       
	       dvnorm += FABS(result[dof]-actnode->sol_increment.a.da[place][j]);
	       vnorm  += FABS(result[dof]);
               actnode->sol_increment.a.da[place][j] = result[dof];
	    }
         }        
      }   
      break;
   case 2: /* L_2 norm */
      /*-----------  loop nodes and put the result back to the node structure */
      for (i=0; i<actfield->dis[0].numnp; i++)
      {
         actnode = &(actfield->dis[0].node[i]);
         /*------------------------------ enlarge sol_increment, if necessary */
         if (place >= actnode->sol_increment.fdim)
         {
            diff = place - actnode->sol_increment.fdim;
            max  = IMAX(diff,5);
            amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,actnode->sol_increment.sdim,"DA");
         }
         for (j=0; j<actnode->numdf; j++) /* loop dofs and calculate the norms */
         {
	    dof = actnode->dof[j];
            if (dof>=numeq_total) continue;
            if (j==predof) /* pressure dof */
	    {
               dpnorm += pow(result[dof]-actnode->sol_increment.a.da[place][j],2);
	       pnorm  += pow(result[dof],2);
	       actnode->sol_increment.a.da[place][j] = result[dof];
	    }
	    else /* vel - dof */
	    {	       
	       dvnorm += pow(result[dof]-actnode->sol_increment.a.da[place][j],2);
	       vnorm  += pow(result[dof],2);
               actnode->sol_increment.a.da[place][j] = result[dof];
	    }
         }        
      }
      dvnorm = sqrt(dvnorm);
       vnorm = sqrt( vnorm);
      dpnorm = sqrt(dpnorm);
       pnorm = sqrt( pnorm);  
      break;
   default:
      dserror("unknown norm for convergence check!");
   }   
   /*------------------------------------------- check for "ZERO-field" */
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
   *vrat = dvnorm/vnorm;
   *prat = dpnorm/pnorm;
   break;
default:
   dserror("parameter itchk out of range!");
}
/*----------------------------------------------------------------------*/
amdel(&result_a);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of  fluid_result_incre */


/*----------------------------------------------------------------------*
 | calculate vel- and pre-norms for steady state check                  |
 | norm = ||U(n+1) - U(n)|| / ||U(n)||                      genk 05/02  |
 |    solution at (n+1): node->sol_increment[3][j]                      |
 |    solution at (n)  : node->sol_increment[1][j]                      |
 *----------------------------------------------------------------------*/
void fluid_norm( FLUID_DYNAMIC *fdyn, FIELD  *actfield, int numeq_total,
                 double *vrat, double *prat )
{
int         i,j;
int         numdf, numvel, predof;
int         numnp_total;
int         actdof;
double      dvnorm=ZERO;
double       vnorm=ZERO;
double      dpnorm=ZERO;
double       pnorm=ZERO;
NODE       *actnode;

#ifdef DEBUG 
dstrc_enter("fluid_norm");
#endif
/*---------------------------------------------------- set some values */
numdf        = fdyn->numdf;
numnp_total  = actfield->dis[0].numnp;
predof       = numdf-1;
numvel       = numdf-1;

switch (fdyn->stnorm)
{
case 0: /* L_infinity norm */
   /*-------------------------------------------------- loop all nodes */
   for (i=0;i<numnp_total;i++)
   {
      actnode=&(actfield->dis[0].node[i]);
      for (j=0;j<numvel;j++)
      {
	 actdof = actnode->dof[j];
         if (actdof>=numeq_total) continue;
         dvnorm = DMAX(dvnorm,FABS(actnode->sol_increment.a.da[3][j]  \
	                          -actnode->sol_increment.a.da[1][j]));
          vnorm = DMAX( vnorm,FABS(actnode->sol_increment.a.da[3][j]));
      }
      actdof = actnode->dof[predof];
      if (actdof>=numeq_total) continue;
      dpnorm = DMAX(dpnorm,FABS(actnode->sol_increment.a.da[3][predof]  \
	                       -actnode->sol_increment.a.da[1][predof]));
       pnorm = DMAX( pnorm,FABS(actnode->sol_increment.a.da[3][predof]));                
   }
   break;
case 1: /* L_1 norm */   
   /*-------------------------------------------------- loop all nodes */
   for (i=0;i<numnp_total;i++)
   {
      actnode=&(actfield->dis[0].node[i]);
      for (j=0;j<numvel;j++)
      {
	 actdof = actnode->dof[j];
         if (actdof>=numeq_total) continue;
         dvnorm += FABS(actnode->sol_increment.a.da[3][j]  \
	               -actnode->sol_increment.a.da[1][j]);
          vnorm += FABS(actnode->sol_increment.a.da[3][j]);
      }
      actdof = actnode->dof[predof];
      if (actdof>=numeq_total) continue;
      dpnorm += FABS(actnode->sol_increment.a.da[3][predof]  \
	            -actnode->sol_increment.a.da[1][predof]);
       pnorm += FABS(actnode->sol_increment.a.da[3][predof]);                
   }     
   break;
case 2: /* L_2 norm */   
   /*-------------------------------------------------- loop all nodes */
   for (i=0;i<numnp_total;i++)
   {
      actnode=&(actfield->dis[0].node[i]);
      for (j=0;j<numvel;j++)
      {
	 actdof = actnode->dof[j];
         if (actdof>=numeq_total) continue;
         dvnorm += pow(actnode->sol_increment.a.da[3][j]  \
	               -actnode->sol_increment.a.da[1][j],2);
          vnorm += pow(actnode->sol_increment.a.da[3][j],2);
      }
      actdof = actnode->dof[predof];
      if (actdof>=numeq_total) continue;
      dpnorm += pow(actnode->sol_increment.a.da[3][predof]  \
	            -actnode->sol_increment.a.da[1][predof],2);
       pnorm += pow(actnode->sol_increment.a.da[3][predof],2);                
   }
   dvnorm = sqrt(dvnorm);
    vnorm = sqrt( vnorm);
   dpnorm = sqrt(dpnorm);
    pnorm = sqrt( pnorm);   
   break;
default:
   dserror("unknown norm for steady state check!");
}   
/*------------------------------------------- check for "ZERO-field" */
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
*vrat = dvnorm/vnorm;
*prat = dpnorm/pnorm;


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_norm*/

/*----------------------------------------------------------------------*
 | routine to copy solution history                         genk 05/02  |
 |    solution at (n+1): node->sol_increment[3][j]                      |
 |    solution at (n)  : node->sol_increment[1][j]                      |
 *----------------------------------------------------------------------*/
void fluid_copysol(FLUID_DYNAMIC *fdyn, FIELD  *actfield, 
                   int from, int to, int flag)
{
int         i,j;
int         numnp_total;
int         diff,max;
int         numdf;
NODE       *actnode;

#ifdef DEBUG 
dstrc_enter("fluid_copysol");
#endif

numnp_total  = actfield->dis[0].numnp;
numdf        = fdyn->numdf;

switch(flag)
{
case 0: /* copy from sol_increment[from] to sol_increment[to] */
   for (i=0;i<numnp_total;i++)
   {
    /*------------------------------ enlarge sol_increment, if necessary */
      actnode=&(actfield->dis[0].node[i]);
      if (to >= actnode->sol.fdim)
      {
         diff = to - actnode->sol_increment.fdim;
         max  = IMAX(diff,5);
         amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,actnode->sol_increment.sdim,"DA");
      }      
      for (j=0;j<numdf;j++)
        actnode->sol_increment.a.da[to][j]=actnode->sol_increment.a.da[from][j];

   }
case 1: /* copy from sol_increment[from] to sol[to]*/
   for (i=0;i<numnp_total;i++)
   {
    /*------------------------------ enlarge sol, if necessary */
      actnode=&(actfield->dis[0].node[i]);
      if (to >= actnode->sol.fdim)
      {
         diff = to - actnode->sol_increment.fdim;
         max  = IMAX(diff,5);
         amredef(&(actnode->sol),actnode->sol.fdim+max,actnode->sol.sdim,"DA");
      }        
      for (j=0;j<numdf;j++)
         actnode->sol.a.da[to][j]=actnode->sol_increment.a.da[from][j];
   }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_copysol*/

/*----------------------------------------------------------------------*
 | routine to do the steady state check                     genk 05/02  |
 *----------------------------------------------------------------------*/
int fluid_steadycheck(FLUID_DYNAMIC *fdyn, FIELD  *actfield, int numeq_total)
{
int         steady=0;
double      vrat,prat;

#ifdef DEBUG 
dstrc_enter("fluid_steadycheck");
#endif

fluid_norm(fdyn,actfield,numeq_total,&vrat,&prat);
/*--------------------------------------------- output to the screen */
if (par.myrank==0)
{
   switch (fdyn->stnorm) 
   {  
   case 0:
      printf("   --> steady state check   (tolerance[norm]):  %10.3#E [L_in] \n",  
	        fdyn->sttol);
      break;
   case 1:
      printf("   --> steady state check   (tolerance[norm]):  %10.3#E [L_1 ] \n",
	        fdyn->sttol);
   case 2:
      printf("   --> steady state check   (tolerance[norm]):  %10.3#E [L_2 ] \n",
	        fdyn->sttol);
   }
   printf("         velocities: %10.3#E	   pressures:   %10.3#E  \n", 
          vrat,prat);
}	  
if (vrat<fdyn->sttol && prat<fdyn->sttol)    
{
   steady=1;
   if (par.myrank==0)
   {
      printf("\n");
      printf("    >>>>>> STEADY STATE REACHED <<<<<< \n");
      printf("\n");      
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return ((int)(steady));
} /* end of fluid_steadycheck*/

/*----------------------------------------------------------------------*
 | routine to do the iteration convergence check            genk 05/02  |
 *----------------------------------------------------------------------*/
int fluid_convcheck(FLUID_DYNAMIC *fdyn, double vrat, double prat, 
                    int itnum, double te, double ts)
{
int         converged=0;

#ifdef DEBUG 
dstrc_enter("fluid_convcheck");
#endif
if (fdyn->itchk!=0)
{
   if (par.myrank==0)
   { 
      switch(fdyn->itnorm)
      {
      case 0:
         printf("|  %3d/%3d   | %10.3#E[L_in]  | %10.3#E   | %10.3#E  | {te: %10.3#E} {ts:%10.3#E} \n", 
                 itnum,fdyn->itemax,fdyn->ittol,vrat,prat,te,ts);
         printf("|            |                   |              |             | \n");
         break;
      case 1:
         printf("|  %3d/%3d   | %10.3#E[L_1 ]  | %10.3#E   | %10.3#E  | {te: %10.3#E} {ts:%10.3#E} \n", 
              itnum,fdyn->itemax,fdyn->ittol,vrat,prat,te,ts);
         printf("|            |                   |              |             | \n");
         break;
      case 2:
         printf("|  %3d/%3d   | %10.3#E[L_2 ]  | %10.3#E   | %10.3#E  | {te: %10.3#E} {ts:%10.3#E} \n", 
                 itnum,fdyn->itemax,fdyn->ittol,vrat,prat,te,ts);
         printf("|            |                   |              |             | \n");
         break;
      }
   }         
   if (vrat<fdyn->ittol && prat<fdyn->ittol)
      converged=1;
   if (itnum==fdyn->itemax)
      converged+=1;
   if (converged==1 && par.myrank==0)
   {
      printf("|            |                   |              |             | \n");
      printf("|          >>>>>> not converged in itemax steps!              |\n");        
      printf("|            |                   |              |             | \n");  
   }
}
else if (par.myrank==0)
{
    printf("      iteration step: %3d / 3d \n",  
            itnum, fdyn->itemax) ;          
}   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return ((int)(converged));
} /* end of fluid_steadycheck*/

/*----------------------------------------------------------------------*
 | routine to print out time and algo to the screen         genk 05/02  |
 *----------------------------------------------------------------------*/
void fluid_algoout(FLUID_DYNAMIC *fdyn, FLUID_DYN_CALC *dynvar)
{

#ifdef DEBUG 
dstrc_enter("fluid_algoout");
#endif

printf("\n");

switch(fdyn->iopfsi)
{
case 2:
   printf("TIME: %11.4#E/%11.4#E  DT = %11.4#E  Semi-Impl-One-Step  STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->step,fdyn->nstep);
   break;
case 3:
   printf("TIME: %11.4#E/%11.4#E  DT = %11.4#E  Semi-Impl-Two-Step  STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->step,fdyn->nstep);
   break;
case 4:
   printf("TIME: %11.4#E/%11.4#E  DT = %11.4#E  One-Step-Theta  STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->step,fdyn->nstep);
   break;
case 5:
   printf("TIME: %11.4#E/%11.4#E  DT = %11.4#E  Fract-Step-Theta  STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->step,fdyn->nstep);
   break;         
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_algoout*/
