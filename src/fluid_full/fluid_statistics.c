/*!----------------------------------------------------------------------
\file
\brief statistical evaluation of fluid solution

<pre>
Maintainer:  name
             e-mail
             homepage
             telephone number


</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../headers/prototypes.h"
#include "fluid_prototypes.h"
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
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;                      

/*----------------------------------------------------------------------*
 | arrays and variables for statistical evaluation    v.gravem 02/03    |
 *----------------------------------------------------------------------*/
static ARRAY     statistic_a;    /* array for solution statistics      */
static DOUBLE  **statistic;
static ARRAY     indicator_a;    /* indicator vector          	       */
static INT      *indicator;
static INT       numnhd;         /* number of non-homogeneous direction*/
static INT       numsnhd=0;         
static INT       numstat=31;     /* number of entries in statistic arr.*/

/*!---------------------------------------------------------------------                                         
\brief initialization for gathering solution values and statistics

<pre>                                                       gravem 01/03

First of all, in this routine homogeneous coordinate directions are 
identified and the respective nodes are assigned.	   

Afterwards, the initialization of the array for gathering single values 
and correlations as well as derivatives at nodes after homogenization 
for statistical evaluation is prepared. Actual initialization in 
fluid_isi. 
The array looks as follows:	   
gathering[i][j]: i=ID of homogenized node (node-ID for no homogenization)
                 j=single values of velocity and pressure (2-D:3/3-D:4)		    
                  +double-correlations of velocity and pressure (3/4)			    
                  +triple-correlations of velocity (2/3)			    
                  +quadruple-correlations of velocity (2/3)		    
                  +uxuy 			    
                  +velocity derivative dux/dy 		    
                  +double-correlation of velocity derivative dux/dy		    
j := 13 for 2-D and 17 for 3-D 			    

Furthermore, the array for statistical evaluation of different measures
is initialized. They are stored in 	   
statistic[i][j]: i=ID of homogenized node (node-ID for no homogenization)
                 j=mean values of velocity and pressure (2-D:3/3-D:4)		    
                  +mv double-correlations of velocity and pressure (3/4)			    
                  +mv triple-correlations of velocity (2/3)			    
                  +mv quadruple-correlations of velocity (2/3)		    
                  +mv uxuy 			    
                  +mv velocity derivative dux/dy 		    
                  +mv double-correlation of velocity derivative dux/dy
		  +mv shear stress 12
		  +mv Reynolds stress 12
		  +mv total stress 12
		  +rms velocities and pressures (3/4)
		  +rms vorticity
		  +skewness velocity (2/3)
		  +flatness velocity (2/3) 		    
j := 24 for 2-D and 31 for 3-D
Both arrays will be allocated for the 3D-case with the respective entries
missing in the 2D-case. At first, all correlations of the velocity, then, 
all correlations of the pressure and, afterwards, uxuy and the 
correlations of the velocity derivatives are stored in the gathering 
array. In the statistic-array the additional entries follow accordingly.  			    

Finally, the nodal-based velocity derivative dux/dy is set to zero.			    
			     
</pre>   
\param *actfield   FIELD	   (i)  actual field
\param *fdyn	   FLUID_DYNAMIC   (i)  
\param  numdco	   INT             (o)  numb. of diff. coordinates   
\return void 

------------------------------------------------------------------------*/
void fluid_initstat(FIELD         *actfield,  
                    FLUID_DYNAMIC *fdyn,
		    INT           *numdco)      
{
INT    i;            /* simply a counter                                */
INT    numdf;        /* number of dofs in this discretisation           */
INT    numnp;        /* total number of nodes in this discretisation    */
INT    ndum;         /* dummy value                                     */
INT    dir[3],homdir,dirsave;
INT    flag;
NODE  *actnode;      /* the actual node                                 */

#ifdef DEBUG 
dstrc_enter("fluid_initstat");
#endif

/*-------------------------------------- set and initialize some values */
numdf = fdyn->numdf;
numnp = actfield->dis[0].numnp;
fdyn->dynvar.avwastr =ZERO;
fdyn->dynvar.avflrate=ZERO;
/*-------- initialize timestep counter for statistics by stress profile */
if (fdyn->turbstat==2) fdyn->statimst = 0;

if (fdyn->spavhom!=0)
/*------------------------------------- determine homogeneous direction */
{
  homdir = fdyn->spavhom;
  math_intextract(homdir,&ndum,&dir[0],&dir[1],&dir[2]);
  for (i=0; i<3; i++) 
  {
    if (dir[i]==0) numnhd = i;
  }  

/*------- define indicator array for nodes in homogeneous lines/planes */
  indicator = amdef("indicator",&indicator_a,numnp,1,"DV");

  if (fdyn->turbstat==4) flag=1; /* turbulent 3-D driven cavity */
  else                   flag=0;
/*------- indicator array identifying nodes in homogeneous directions  */
  fluid_homogen(actfield,indicator,numnp,numdco,flag);
}  
/*------------------------------------------- no homogeneous direction */
else *numdco = numnp;

/*----------------- allocate space for gathering and statistics arrays */
statistic = amdef("statistic",&statistic_a,*numdco,numstat,"DA");
amzero(&statistic_a);

/*------------------- set velocity derivative dux/dy at nodes to zero  */
for (i=0;i<numnp;i++)
{
   actnode=&(actfield->dis[0].node[i]);
   actnode->deruxy=ZERO;
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_initstat */ 

/*!---------------------------------------------------------------------                                         
\brief create indicator array to distinguish homogeneous directions

<pre>                                                       gravem 01/03

In this routine, an indicator array to distinguish homogeneous directions
is created.			    
			     
</pre>   
\param *actfield      FIELD	     (i)  actual field
\param *indicator     INT	     (o)  indicator array
\param  numnp         INT	     (i)  number of nodal points
\param *numdco        INT	     (o)  number of diff. coordinates
\param  flag          INT	     (i)  flag for num. of non-hom. dir.
\return void 

------------------------------------------------------------------------*/
void fluid_homogen(FIELD *actfield, 
                   INT   *indicator,
		   INT    numnp,
		   INT   *numdco,
		   INT    flag)
{
INT        i,j,k;
NODE      *firnode,*secnode;

#ifdef DEBUG 
dstrc_enter("fluid_homogen");
#endif

/*------------------------ initialize entries of indicator-array to -1  */
for (i=0;i<numnp;i++)
{
  indicator[i]=-1;
}
/*-- initialize counter for number of different coordinates of homogen. */
*numdco=0;

for (j=0; j<numnp; j++)
{
  if (indicator[j]>=0) continue;
/*---------------------------- set first node if not already registered */
  firnode = &(actfield->dis[0].node[j]);
  if (flag!=0)
  {
/*- identify nodes in y-dir. (x=0.5 and z=0.5) for turb. 3-D driv. cav. */
    if (FABS(firnode->x[0]-0.5)<EPS5 && FABS(firnode->x[2]-0.5)<EPS5) 
    {
      indicator[j] = *numdco;
      ++(*numdco);
    }  
  }
  else
  {
    indicator[j] = *numdco;
    for (k=j+1; k<numnp; k++)
    {  
      if (indicator[k]>=0) continue;
/* identify second node with identical coordinate if not already regist. */
      secnode = &(actfield->dis[0].node[k]);
      if (FABS(firnode->x[numnhd]-secnode->x[numnhd])<EPS5) 
      indicator[k] = *numdco;
    }
    ++(*numdco);
  }
}

/* there is a 2nd non-homogen. direction for turbulent 3-D driven cavity */
if (flag!=0)
{
  numsnhd=*numdco;
  for (j=0; j<numnp; j++)
  {
/*---------------------------- set first node if not already registered */
    firnode = &(actfield->dis[0].node[j]);
/*- identify nodes in x-dir. (y=0.5 and z=0.5) for turb. 3-D driv. cav. */
    if (FABS(firnode->x[1]-0.5)<EPS5 && FABS(firnode->x[2]-0.5)<EPS5)
    { 
      indicator[j] = *numdco;
      ++(*numdco);
    }  
  }  
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_homogen*/

/*!---------------------------------------------------------------------                                         
\brief calculate vorticity thickness for plane mixing layer

<pre>                                                       gravem 07/03

In this routine, the vorticity thickness for the plane mixing layer is
calculated. -> see dissertation for definition of vorticity thickness		    
			     
</pre>   
\param  *fdyn 	      FLUID_DYNAMIC  (i)   
\param  *dynvar       FLUID_DYN_CALC (o)   
\param  *actfield     FIELD	     (i)  actual field
\param  *actpart      PARTITION	     (i)  
\param  *actintra     INTRA	     (i)  
\param   numdco       INT	     (i)  number of diff. coordinates
\return void 

------------------------------------------------------------------------*/
void fluid_vorthick(FLUID_DYNAMIC  *fdyn,
                    FLUID_DYN_CALC *dynvar, 
                    FIELD	   *actfield,  
                    PARTITION	   *actpart, 
		    INTRA	   *actintra,
		    INT 	    numdco)
{
INT        j,k;               /* simply some counters                   */
INT        myrank,imyrank;    /* processor rank identification          */
INT	   numnp;             /* number of nodal points                 */
INT	   numnod;            /* node counter                           */
DOUBLE     imvortx;           /* integral mean of vorticity in x-dir.   */
DOUBLE     maxvort;           /* maximum value of vorticity             */
DOUBLE     uinf=ONE;          /* maximum velocity in x-direction        */
DOUBLE     delta=ONE/(FOUR*SEVEN); /* initial vorticity thickness       */ 
NODE      *actnode;           /* actual node                            */

#ifdef DEBUG 
dstrc_enter("fluid_vorthick");
#endif

/*----------------------------------------------------- set some values */
numnp   = actfield->dis[0].numnp;
myrank  = par.myrank;
imyrank = actintra->intra_rank;
maxvort = ZERO;
/*----------------------------------------------------------------------*/
if (imyrank==0 && myrank==0)
{
  for (j=0; j<numdco; j++)
  {
/* calculate integral mean in x-direction as discrete sum of nod. values */
    imvortx=ZERO;
    numnod=-1;
    for (k=0; k<numnp; k++)
    {
      if (indicator[k]==j) 
      { 
        actnode = &(actfield->dis[0].node[k]);
        if ((actnode->x[0]-ZERO)<EPS5 || (actnode->x[0]-ONE)<EPS5)
/*----------------------------------------------------- boundary values */
          imvortx+=FABS((ONE/TWO)*actnode->deruxy);
/*------------------------------------------------- non-boundary values */
        else imvortx+=FABS(actnode->deruxy);
        numnod++; 
      }  
    }
/*-------------- determine maximum value of vorticity along y-direction */
    maxvort=DMAX(maxvort,imvortx);
  }
/*--------------------------------------- calculate vorticity thickness */
  dynvar->vorthick=TWO*uinf*numnod/(maxvort*delta);
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_vorthick */

/*!---------------------------------------------------------------------                                         
\brief calculate wall shear stress tauw

<pre>                                                       gravem 02/03

In this routine, the wall shear stress tauw is calculated as the mean 
value of the wall shear at the upper and at the lower wall. 			    
			     
</pre>   
\param *dynvar 	      FLUID_DYN_CALC (o)   
\param *actfield      FIELD	     (i)  actual field
\param *actpart       PARTITION	     (i)  
\param *actintra      INTRA	     (i)  
\return void 

------------------------------------------------------------------------*/
void fluid_washstr(FLUID_DYN_CALC *dynvar, 
                   FIELD          *actfield,
                   PARTITION      *actpart,    
                   INTRA          *actintra)    
{
INT         i,j,k,l;      /* simply some counters                       */
INT         numnp;        /* total number of nodes                      */
INT         numnlw=0;     /* node counter at lower wall                 */
INT         numnuw=0;     /* node counter at upper wall                 */
INT         numnlwrec;    /* node counter at lower wall                 */
INT         numnuwrec;    /* node counter at upper wall                 */
INT         actmat;       /* material number of active element          */
DOUBLE      visc;         /* viscosity                                  */
DOUBLE      dens;         /* density                                    */
DOUBLE      wsslw=0.0;    /* wall shear stress at lower wall            */
DOUBLE      wssuw=0.0;    /* wall shear stress at upper wall            */
DOUBLE      wsslwrec;     /* wall shear stress at lower wall            */
DOUBLE      wssuwrec;     /* wall shear stress at upper wall            */
NODE       *actnode;      /* actual node                                */

#ifdef DEBUG 
dstrc_enter("fluid_washstr");
#endif

/*------------------------------- assume constant viscosity and density */
actmat=actfield->dis[0].element[0].mat-1;
dens = mat[actmat].m.fluid->density;
visc = mat[actmat].m.fluid->viscosity*dens;/* here we need dynamic visc.! */

for (i=0; i<actpart->pdis[0].numnp; i++)
{
   /*------------------------------------ set pointer to active node */
  actnode = actpart->pdis[0].node[i];
  if (FABS(actnode->x[numnhd]+ONE)<EPS5) 
  {
    numnlw++;
/*--------------------------- shear stress 12 of mean flow at lower wall */
    wsslw += FABS(actnode->deruxy);
  }
  if (FABS(actnode->x[numnhd]-ONE)<EPS5) 
  {
    numnuw++;
/*-------------------------- shear stress 12 of mean flow at upper wall */
    wssuw += FABS(actnode->deruxy);
  }
/*------------------------------------- set velocity derivative to zero */
  actnode->deruxy=ZERO;
} /* end of loop over nodes */

/*------------------------------------- wall shear stress as mean value */
#ifdef PARALLEL 
MPI_Allreduce(&wsslw,&wsslwrec,1,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
MPI_Allreduce(&wssuw,&wssuwrec,1,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
MPI_Allreduce(&numnlw,&numnlwrec,1,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
MPI_Allreduce(&numnuw,&numnuwrec,1,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
dynvar->washstr = (visc/TWO)*((wsslwrec/numnlwrec)+(wssuwrec/numnuwrec));
#else
dynvar->washstr = (visc/TWO)*((wsslw/numnlw)+(wssuw/numnuw));
#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_washstr*/

/*!---------------------------------------------------------------------                                         
\brief gather solution values at nodes every time step

<pre>                                                       gravem 01/03

In this routine, single values and correlations at nodes are gathered for
statistical evaluation. They are stored in 	   
gathering[i][j]: i=ID of homogenized node (node-ID for no homogenization)
                 j=single values of velocity and pressure (2-D:3/3-D:4)		    
                  +double-correlations of velocity and pressure (3/4)			    
                  +triple-correlations of velocity (2/3)			    
                  +quadruple-correlations of velocity (2/3)		    
                  +uxuy			    
                  +velocity derivative dux/dy 	    
                  +double-correlation of velocity derivative dux/dy		    

</pre>   
\param  *fdyn 	      FLUID_DYNAMIC  (i)   
\param  *actfield     FIELD	     (i)  actual field
\param **gathering    DOUBLE	     (o)  gathering array
\param   numdco       INT	     (i)  number of diff. coordinates
\return void 

------------------------------------------------------------------------*/
void fluid_solgath(FLUID_DYNAMIC *fdyn, 
                   FIELD         *actfield,
		   DOUBLE       **gathering,
		   INT            numdco)
{
INT         i,j,k;        /* simply some counters                       */
INT         nent;         /* simply a counter for the entries           */
INT         numco;        /* number of coordinate in homogeneous plane  */
INT         numnp;        /* total number of nodes                      */
INT         numnop;       /* number of nodes in one homogeneous plane   */
INT         numdf;        /* number of fluid dofs                       */
NODE       *actnode;      /* actual node                                */

#ifdef DEBUG 
dstrc_enter("fluid_solgath");
#endif

/*----------------------------------------------------- set some values */
numnp = actfield->dis[0].numnp;
numdf = fdyn->numdf;
/*---------------------------- number of nodes in one homogeneous plane */
if (fdyn->turbstat==4) numnop = 1;
else                   numnop = numnp/numdco;

for (i=0;i<numnp;i++) /* loop nodes */
{
  actnode = &(actfield->dis[0].node[i]);
  numco=indicator[i];
  if (indicator[i]==-1) continue;
  nent=0;
/* gather single, double, triple and quadruple correlations of velocity */
  for (j=0; j<(numdf-1); j++) /* loop dofs */
  {
    gathering[numco][nent]  +=actnode->sol_increment.a.da[3][j]/numnop;
    gathering[numco][nent+1]+=(actnode->sol_increment.a.da[3][j]\
                              *actnode->sol_increment.a.da[3][j])/numnop;
    gathering[numco][nent+2]+=(actnode->sol_increment.a.da[3][j]\
                              *actnode->sol_increment.a.da[3][j]\
                              *actnode->sol_increment.a.da[3][j])/numnop;
    gathering[numco][nent+3]+=(actnode->sol_increment.a.da[3][j]\
                              *actnode->sol_increment.a.da[3][j]\
                              *actnode->sol_increment.a.da[3][j]\
                              *actnode->sol_increment.a.da[3][j])/numnop;
    nent=nent+4;
  } /* end of loop over dofs */
/*---------------------------------------- increase counter for 2D-case */
  if (numdf<4) nent=nent+4;
  
/*------------------- gather single and double correlations of pressure */
  gathering[numco][nent]  +=actnode->sol_increment.a.da[3][numdf-1]/numnop;
  gathering[numco][nent+1]+=(actnode->sol_increment.a.da[3][numdf-1]\
			    *actnode->sol_increment.a.da[3][numdf-1])/numnop;
  nent=nent+2;

/*--------------------------------------------------------- gather uxuy */
  gathering[numco][nent]+=(actnode->sol_increment.a.da[3][0]\
			  *actnode->sol_increment.a.da[3][1])/numnop;
  nent++;

/*- gather single and double correlations of velocity derivative dux/dy */
  gathering[numco][nent]  +=actnode->deruxy/numnop;
  gathering[numco][nent+1]+=(actnode->deruxy*actnode->deruxy)/numnop;
  actnode->deruxy=ZERO;
} /* end of loop over nodes */
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_solgath */

/*!---------------------------------------------------------------------                                         
\brief calculate and check stress profile

<pre>                                                       gravem 02/03

In this routine, the stress profile (shear stresses + Reynolds stresses)
is calculated and checked subject to its linearity. Therefore, the
statistic array is partially used. 			    
statistic[i][j]: i=ID of homogenized node (node-ID for no homogenization)
                 j=mean values of velocity and pressure (2-D:3/3-D:4)		    
                  +mv double-correlations of velocity and pressure (3/4)			    
                  +mv triple-correlations of velocity (2/3)			    
                  +mv quadruple-correlations of velocity (2/3)		    
                  +mv uxuy			    
                  +mv velocity derivative dux/dy 		    
                  +mv double-correlation of velocity derivative dux/dy
		  +mv shear stress 12
		  +mv Reynolds stress 12
		  +mv total stress 12
		  +rms velocities and pressures (3/4)
		  +rms vorticity
		  +skewness velocity (2/3)
		  +flatness velocity (2/3) 		    
j := 24 for 2-D and 31 for 3-D 
			     
</pre>   
\param  *fdyn 	      FLUID_DYNAMIC  (i)   
\param  *actfield     FIELD	     (i)  actual field
\param **gathering    DOUBLE	     (i)  gathering array
\param   numdco       INT	     (i)  number of diff. coordinates
\return INT statstart 
\sa fluid_profcheck()

------------------------------------------------------------------------*/
INT fluid_profcheck(FLUID_DYNAMIC *fdyn, 
                    FIELD         *actfield,
		    DOUBLE       **gathering,
		    INT            numdco)
{
INT         i;            /* simply a counter                           */
INT         numts;        /* number of time steps for statistics        */
INT         numnod;       /* node counter                               */
INT         numco;        /* number of coordinate of homogeneous plane  */
INT         numdf;        /* number of fluid dofs                       */
INT         actmat;       /* material number of active element          */
INT         statstart=0;  /* flag for start of gathering for statistics */
DOUBLE      visc;         /* viscosity                                  */
DOUBLE      dens;         /* density                                    */
DOUBLE      stprof;       /* value sum of stress profile                */
DOUBLE      stmax;        /* maximum averaged value of stress profile   */
DOUBLE      sttol;        /* tolerance for sum of stress profile        */
NODE       *actnode;      /* actual node                                */

#ifdef DEBUG 
dstrc_enter("fluid_profcheck");
#endif

/*----------------------------------------------------- set some values */
numts = fdyn->step;
stprof=ZERO;
stmax=ZERO;
/*----------- assuming constant material parameters all over the domain */
actmat=actfield->dis[0].element[0].mat-1;
dens = mat[actmat].m.fluid->density;
visc = mat[actmat].m.fluid->viscosity*dens;/* here we need dynamic visc.! */

for (i=0;i<numdco;i++) /* loop over different coordinates */
{
/*------------------------------------- mean values of velocity 1 and 2 */
  statistic[i][0]=gathering[i][0]/numts;
  statistic[i][4]=gathering[i][4]/numts;

/*----------------------- mean value of correlation of velocity 1 and 2 */
  statistic[i][14]=gathering[i][14]/numts;

/*----------------------------------------- shear stress 12 of mean flow */
  statistic[i][17]=(visc/numts)*gathering[i][15];

/*-------------------------------------------------- Reynolds stress 12 */
  statistic[i][18]=statistic[i][14]-(statistic[i][0]*statistic[i][4]);

/*------------------------------------------------------ total stress 12 */
  statistic[i][19]=statistic[i][17]-statistic[i][18];

/*--------------- summation to get a representative value of the profile */
  stprof += statistic[i][19];
  stmax = DMAX(stmax,FABS(statistic[i][19]));
} /* end of loop over different coordinates */

/*--------- check if sum of profile within 0.1%-deviation of max. stress */
sttol = stmax/(TEN*TEN*TEN);
if (stprof<=sttol) statstart = 1;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return ((INT)(statstart));
} /* end of fluid_profcheck*/

/*!---------------------------------------------------------------------                                         
\brief calculate statistical measures

<pre>                                                       gravem 01/03

In this routine statistical measures are calculated. They are stored in 	   
statistic[i][j]: i=ID of homogenized node (node-ID for no homogenization)
                 j=mean values of velocity and pressure (2-D:3/3-D:4)		    
                  +mv double-correlations of velocity and pressure (3/4)			    
                  +mv triple-correlations of velocity (2/3)			    
                  +mv quadruple-correlations of velocity (2/3)		    
                  +mv uxuy			    
                  +mv velocity derivative dux/dy 		    
                  +mv double-correlation of velocity derivative dux/dy
		  +mv shear stress 12
		  +mv Reynolds stress 12
		  +mv total stress 12
		  +rms velocities and pressures (3/4)
		  +rms vorticity
		  +skewness velocity (2/3)
		  +flatness velocity (2/3) 		    
			     
</pre>   
\param  *fdyn 	      FLUID_DYNAMIC  (i)   
\param  *actfield     FIELD	     (i)  actual field
\param  *actpart      PARTITION	     (i)  
\param  *actintra     INTRA	     (i)  
\param **gathering    DOUBLE	     (i)  gathering array
\param   numgath      INT	     (i)  num. of entries in gath. array
\param   numdco       INT	     (i)  number of diff. coordinates
\return void 

------------------------------------------------------------------------*/
void fluid_statistics(FLUID_DYNAMIC *fdyn, 
                      FIELD         *actfield,  
                      PARTITION     *actpart, 
		      INTRA         *actintra,
		      DOUBLE       **gathering,
		      INT            numgath,
		      INT            numdco)
{
INT         i,j;          /* simply some counters                       */
INT         numdf;        /* number of fluid dofs                       */
INT         numts;        /* number of time steps for statistics        */
INT         actmat;       /* material number of active element          */
DOUBLE      visc;         /* viscosity                                  */
DOUBLE      dens;         /* density                                    */
NODE       *actnode;      /* actual node                                */

#ifdef DEBUG 
dstrc_enter("fluid_statistics");
#endif

/*----------------------------------------------------- set some values */
numdf = fdyn->numdf;
numts = (fdyn->step)-(fdyn->statimst)+1;
/*----------- assuming constant material parameters all over the domain */
actmat=actfield->dis[0].element[0].mat-1;
dens = mat[actmat].m.fluid->density;
visc = mat[actmat].m.fluid->viscosity*dens;/* here we need dynamic visc.! */

for (i=0;i<numdco;i++) /* loop nodes */
{
  for (j=0; j<numgath; j++) /* loop over entries of gathering array */
  {
/*- m. v. of single, double, triple and quadruple correlations velocity */
/*----------- mean values of single and double correlations of pressure */
/*-------------------------------------------------- mean value of uxuy */
/* m.v. of single and double correlations of velocity derivative dux/dy */
      statistic[i][j]=gathering[i][j]/numts;
  } /* end of loop over entries of gathering array */
  
/*--------------------------------------------------- additional entries */
/*----------------------------------------- shear stress 12 of mean flow */
  statistic[i][17]=(visc/numts)*gathering[i][15];

/*-------------------------------------------------- Reynolds stress 12 */
  statistic[i][18]=statistic[i][14]-(statistic[i][0]*statistic[i][4]);

/*----------------------------------------------------- total stress 12 */
  statistic[i][19]=statistic[i][17]-statistic[i][18];

/*------------- rms of velocities, pressure and vorticity omega/omega3  */
  statistic[i][20]=sqrt(statistic[i][1]-(statistic[i][0]*statistic[i][0]));
  statistic[i][21]=sqrt(statistic[i][5]-(statistic[i][4]*statistic[i][4]));
  if (numdf==4) 
    statistic[i][22]=sqrt(statistic[i][9]-(statistic[i][8]*statistic[i][8]));
  statistic[i][23]=sqrt(statistic[i][13]-(statistic[i][12]*statistic[i][12]));
  statistic[i][24]=sqrt(statistic[i][16]-(statistic[i][15]*statistic[i][15]));

/*---------------------------------------------- skewness of velocities */
  statistic[i][25]=(statistic[i][2]-(statistic[i][0]*statistic[i][20]\
                   *statistic[i][20])-\
		   (statistic[i][0]*statistic[i][0]*statistic[i][0]))/\
		   (statistic[i][20]*statistic[i][20]*statistic[i][20]);
  statistic[i][26]=(statistic[i][6]-(statistic[i][4]*statistic[i][21]\
                   *statistic[i][21])-\
                   (statistic[i][4]*statistic[i][4]*statistic[i][4]))/\
		   (statistic[i][21]*statistic[i][21]*statistic[i][21]);
  if (numdf==4) 
    statistic[i][27]=(statistic[i][10]-(statistic[i][8]*statistic[i][22]\
                     *statistic[i][22])-\
                     (statistic[i][8]*statistic[i][8]*statistic[i][8]))/\
		     (statistic[i][22]*statistic[i][22]*statistic[i][22]);

/*---------------------------------------------- flatness of velocities */
  statistic[i][28]=(statistic[i][3]-(statistic[i][0]*statistic[i][25])-\
                   (statistic[i][0]*statistic[i][0]*statistic[i][20]\
		   *statistic[i][20])-(statistic[i][0]*statistic[i][0]\
                   *statistic[i][0]*statistic[i][0]))/(statistic[i][20]\
		   *statistic[i][20]*statistic[i][20]*statistic[i][20]);
  statistic[i][29]=(statistic[i][7]-(statistic[i][4]*statistic[i][26])-\
                   (statistic[i][4]*statistic[i][4]*statistic[i][21]\
		   *statistic[i][21])-(statistic[i][4]*statistic[i][4]\
                   *statistic[i][4]*statistic[i][4]))/(statistic[i][21]\
		   *statistic[i][21]*statistic[i][21]*statistic[i][21]);
  if (numdf==4) 
    statistic[i][30]=(statistic[i][11]-(statistic[i][8]*statistic[i][27])-\
                     (statistic[i][8]*statistic[i][8]*statistic[i][22]\
		     *statistic[i][22])-(statistic[i][8]*statistic[i][8]\
                     *statistic[i][8]*statistic[i][8]))/(statistic[i][22]\
		     *statistic[i][22]*statistic[i][22]*statistic[i][22]);
} /* end of loop over nodes */

/*---------------------------------------- print out statistics to .out */
out_turbstat(actfield,actpart,actintra,statistic,indicator,numdco,numnhd,
             numsnhd,fdyn->statimst,numts);

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_statistics*/


#endif
