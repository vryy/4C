/*!----------------------------------------------------------------------
\file
\brief Call and control VISUAL2 

------------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#ifdef VISUAL2_PACKAGE
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD         *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 *----------------------------------------------------------------------*/
extern struct _GENPROB       genprob;
/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION    *partition;
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum _CALC_ACTION calc_action[MAXFIELD];
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn; 
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                        genk 07/02    |
 | all variables needed by VISUAL2 are defined extern                   |
 *----------------------------------------------------------------------*/
static int      MPTRI,MPPTRI,MFACE;   /* see VISUAL2 manual             */
static int      MPFACE,MEDGE,MPEDGE;  /* see VISUAL2 manual             */
static int      KCELL;                /* number of visualised cells - 
                                 different from number of elements!     */
static int      IOPT;                 /* program mode                   */
static int      IVORT;                /* flag for vorticity calculation */
static int      DSTEP;                /* increment of visualised steps  */
static int      NUMA;                 /* number of ALE-FIELD            */
static int      CMNCOL=300;           /* see VISUAL2 manual             */
static int      CMUNIT=37;            /* see VISUAL2 manual             */
static int      MNODE;                /* maximum number of nodes        */
static int      NKEYS=11;              /* see VISUAL2 manual             */
static int      XYPIX[2];             /* size of VISUAL2 window         */
static int      IKEYS[]={112,115,118,102,120,121,97,116,102,92,109};
static int      FKEYS[]={1,1,1,3,1,1,1,1,3,1,1};             
static float    FLIMS[11][2];          /* data limits                    */
static float    XYMIN[2], XYMAX[2];   /* min. and max. coordinates      */
static struct  _ARRAY PCELL_A;
static struct  _ARRAY WCELL_A;
static struct  _ARRAY WNODE_A;
static struct  _ARRAY WFACE_A;
static struct  _ARRAY CEDGE_A;        /* arrays needed in qat2v2        */
static int    **PCELL;
static int     *WCELL;
static int     *WNODE;
static int     *WFACE; 
static int     *CEDGE;                /* pointers to arrays             */
/*------------------------------- variables needed for own calculations */
static int     numnp;	        /* number of nodes of actual field	*/
static int     numele;	        /* number of elements of actual field	*/
static int     ncols=1;         /* number of sol steps stored in sol	*/
static int     nacols=1;        /* number of sol steps of ALE field     */
static int     isstr=0;         /* flag for streamline calculation	*/
static int     icol=-1;         /* act. num. of sol step to be visual.  */
static int     hsize, vsize;    /* size for VISUAL2 window		*/
static int     true=-1;         /* flag for v2_cursor			*/
static int     false=0;         /* flag for v2_cursor			*/
static int     STSTEP=-2;       /* stopping step                        */
static int     IMOVIE=0;        /* counter to slow down for movie creat.*/
static int     bgcolour;        /* background colour                    */
static float   sstrval=0.0;     /* stationary streamline value		*/
static float   dx,dy,dxy;       /* for determination of coord. limits*/
static float   ratio;	        /*					*/
static double  yscale;	        /* scaling factors for y-direction	*/
static double  velx;	        /*					*/
static double  vely;	        /*					*/
static double  pres;	        /*					*/
static double  absv;	        /*					*/
static double  vort;            /*                                      */
static double  xy;
static double  minxy,maxxy;     /*					*/
static double  minpre,maxpre;   /*					*/
static double  minvx,maxvx;     /*					*/
static double  minvy,maxvy;     /*					*/
static double  minabsv,maxabsv; /*					*/
static double  minvort,maxvort; /*					*/
static double  mingv,maxgv;     /* for determination of data limits	*/
static double  minvort,maxvort; /*					*/
ELEMENT       *actele;          /* actual element			*/
NODE          *actnode;         /* actual node				*/
NODE          *actanode;        /* actual ale-node                      */
GNODE         *actgnode;        /* actual gnode  			*/
DIS_TYP        distyp;          /* element type  			*/
FIELDTYP       actfieldtyp;     /* actual fieldtyp			*/
FIELD         *actfield;        /* actual field  			*/ 
FIELD         *alefield;
FIELD         *structfield;
ARRAY          time_a ;         /* time array				*/
FLUID_DYNAMIC *fdyn;                /* pointer to fluid dyn. inp.data   */
FLUID_DYN_CALC *dynvar;             /* pointer to fluid_dyn_calc        */ 
/*!---------------------------------------------------------------------                                         
\brief call of visual2 for fluid 

<pre>                                                         genk 07/02      
</pre>  
\param  numf   int      (i)       actual number of fluid field
\return void                                                                       

------------------------------------------------------------------------*/
void vis2caf(int numff, int numaf, int numsf) 
{

int iscan, i,j,k;
CALC_ACTION    *action;             /* pointer to the cal_action enum   */
INTRA          *actintra;           /* pointer to active intra-communic.*/
PARTITION      *actpart;            /* pointer to active partition      */
CONTAINER       container;          /* contains variables defined in container.h */

#ifdef DEBUG 
dstrc_enter("vis2caf");
#endif

actfield=&(field[numff]);
actfieldtyp=fluid;
actpart= &(partition[numff]);
action = &(calc_action[numff]);
NUMA=numaf;
/*------------------------- read solution data from pss-file of proc 0 */
visual_readpss(actfield,&ncols,&time_a);

/*---------------------------------------------------------------------*
 * input from the screen                                               *
 *---------------------------------------------------------------------*/

/*----------------------------------------------- data structure sizes */
MPTRI  = 100000;
MPPTRI = 1500;
MFACE  = 175000;
MPFACE = 10000;
MEDGE  = 10000;
MPEDGE = 1000;
MNODE  = actfield->dis[0].numnp;

input1:
printf("\n");
printf("   Actual size of data structures:\n");
printf("      MNODE  = %d\n",MNODE);
printf("      MPTRI  = %d\n",MPTRI);
printf("      MPPTRI = %d\n",MPPTRI);
printf("      MFACE  = %d\n",MFACE);
printf("      MPFACE = %d\n",MPFACE);
printf("      MEDGE  = %d\n",MEDGE);
printf("      MPEDGE = %d\n",MPEDGE);
printf("\n");
printf("   New size of data structures? (0: no; 1: yes)\n");
scanf("%d",&iscan);
if (iscan!=0)
{
   printf("MPTRI  =\n");  
   scanf("%d",&MPTRI);
   printf("MPPTRI =\n"); 
   scanf("%d",&MPPTRI);
   printf("MFACE  =\n");
   scanf("%d",&MFACE);
   printf("MPFACE =\n");
   scanf("%d",&MPFACE);
   printf("MEDGE  =\n");
   scanf("%d",&MEDGE);
   printf("MPEDGE =\n");
   scanf("%d",&MPEDGE);
   goto input1;
} /* endif (iscan!=0) */

/*----------------------------------------------- calculate vorticity */
/* READ in decision from screen !

/* idea for storing of vorticity:
   amredef of node->sol and then store it there! *gg*                 */

printf("\n");
printf("   Do you want to calculate the vorticity? (0: no; 1: yes)\n");
scanf("%d",&IVORT);

if (IVORT==1)
{
   printf("\n");
   printf("   Calculating the vorticity...\n");

   fdyn = alldyn[numff].fdyn;
   dynvar      = &(fdyn->dynvar);
   dynvar->ncols=ncols;
   *action = calc_fluid_initvort;
   calinit(actfield,actpart,action,&container);
   *action = calc_fluid_vort; 
   container.nif=0;
   container.nii=0;
   container.fieldtyp=fluid;
   container.dvec=NULL;
   calelm(actfield,NULL,actpart,NULL,0,0,
          &container,action); 
} /* endif (IVORT==1) */

input2:
printf("\n");
printf("   Please give the mode in which you want to run VISUAL2:\n");
printf("     0: steady data structure, grid and variables\n");
printf("     1: steady data structure and grid, unsteady variables\n");
printf("     2: steady data structure, unsteady grid and variables\n");
printf("     3: unsteady data structure, grid and variables\n");
scanf("%d",&IOPT);
if (IOPT==0) icol++;
if (IOPT==3)
{
   printf("   sorry, mode not implemented yet - new input!\n");
   goto input2;
}
if (IOPT==2) /*-------------------------- read ALE field from pss-file */
{
   if (genprob.probtyp==prb_fsi)
   {
      if (numaf<0) dserror("ALE-field does not exist!");
      if (numsf<0) dserror("STRUCTURE-field does not exist!");
      alefield=&(field[numaf]);
      structfield=&(field[numsf]);
      fsi_initcoupling(structfield,actfield,alefield);
      visual_readpss(alefield,&nacols,NULL);
      if (ncols!=nacols)
      {
         printf("\n");
	 printf("WARNING:number of columns different in ALE and FLUID field \n");
         printf("\n");
	 ncols=IMIN(ncols,nacols);
      }
   }
   if (genprob.probtyp==prb_fluid)
   {
      if (numaf<0) dserror("ALE-field does not exist!");
      alefield=&(field[numaf]);
      fluid_initmfcoupling(actfield,alefield);
      visual_readpss(alefield,&nacols,NULL);
      if (ncols!=nacols)
      {
         printf("\n");
	 printf("WARNING:number of columns different in ALE and FLUID field \n");
         printf("\n");
         ncols=IMIN(ncols,nacols);
      }
   }
}
/*--------------------------------------------------------- y-scaling */
printf("\n");
printf("   Scaling parameter for y-coordinates?\n");
scanf("%lf",&yscale);

/*--------------------------------------------------------- increment */
input3:
printf("\n");
printf("   Increment step number by ...? (<%d)\n",ncols);
scanf("%d",&DSTEP);
if (DSTEP>ncols-1)
{
   printf("   Increment too large --> Try again!\n");
   goto input3;
}

/*------------------------------------------------- background colour */
printf("\n");
printf("   Background colour? (0: black   1: white)\n",bgcolour);
scanf("%d",&bgcolour);


/*--------------------------------------------------- set some values */
numele = actfield->dis[0].numele;
numnp  = actfield->dis[0].numnp;
actele = &(actfield->dis[0].element[0]);
distyp = actele->distyp;

/*---------------- all element types are cut to 3- or 4-node elements */
switch (distyp)
{
case tri3:
   KCELL = numele;
break;
case tri6:
   KCELL = 4*numele;
break;
case quad4:
   KCELL = numele;
break;
case quad8:
   KCELL = 5*numele;
break;
case quad9:
   KCELL = 4*numele;
break;
default:
   dserror("distyp not implemented yet!\n");         
} /* end switch(distyp) */

/*----------------------------------------------------- define arrays */
PCELL = amdef("PCELL",&PCELL_A,KCELL,4,"IA");
WCELL = amdef("WCELL",&WCELL_A,KCELL,1,"IV");
WNODE = amdef("WNODE",&WNODE_A,numnp,1,"IV");
WFACE = amdef("WFACE",&WFACE_A,6*MFACE,1,"IV");
CEDGE = amdef("CEDGE",&CEDGE_A,MEDGE,1,"IV");
amzero(&PCELL_A);
amzero(&WCELL_A);
amzero(&WNODE_A);
amzero(&WFACE_A);
amzero(&CEDGE_A); 

/*------------------------- check if all elements are the same distyp */
for (i=1;i<numele;i++) /* loop all elements */
{
   actele=&(actfield->dis[0].element[i]);
   if (actele->distyp!=distyp)
      dserror ("up to now, all elements have to be the same distyp!\n");
} /* end loop over all elements */

/*---------------------------------- set up modified element topology */
v2cell(actfield);

/*------------------------------------------------ get the x/y limits */
/* find maximum x-coordinate */
if (IOPT>=2)
{
   actnode=&(actfield->dis[0].node[0]);
   minxy=actnode->x[0];
   maxxy=actnode->x[0];
   for (i=1;i<numnp;i++)
   {
      actnode=&(actfield->dis[0].node[i]);
      actgnode=actnode->gnode;
      actanode=actgnode->mfcpnode[genprob.numaf];
      if (actanode==NULL)
      {
         minxy=DMIN(minxy,actnode->x[0]);
         maxxy=DMAX(maxxy,actnode->x[0]);
      }
      else
      {
         xy=actanode->x[0];
         for (j=0;j<ncols;j++)
         {
            dxy=xy+actanode->sol.a.da[j][0];
	    minxy=DMIN(minxy,dxy);
	    maxxy=DMAX(maxxy,dxy);
         }
      }
   }
}
else
{
   actnode=&(actfield->dis[0].node[0]);
   minxy=actnode->x[0];
   maxxy=actnode->x[0];
   for (i=1;i<numnp;i++)
   {
      actnode=&(actfield->dis[0].node[i]);
      minxy=DMIN(minxy,actnode->x[0]);
      maxxy=DMAX(maxxy,actnode->x[0]);
   }
   
}
XYMIN[0] = minxy - (maxxy-minxy)*0.1;
XYMAX[0] = maxxy + (maxxy-minxy)*0.1;

/* find maximum y-coordinate */
if (IOPT>=2)
{
   actnode=&(actfield->dis[0].node[0]);
   minxy=actnode->x[1];
   maxxy=actnode->x[1];
   for (i=1;i<numnp;i++)
   {
      actnode=&(actfield->dis[0].node[i]);
      actgnode=actnode->gnode;
      actanode=actgnode->mfcpnode[numaf];
      if (actanode==NULL)
      {
         minxy=DMIN(minxy,actnode->x[1]);
         maxxy=DMAX(maxxy,actnode->x[1]);      
      }
      else
      {
         xy=actanode->x[1];
         for (j=0;j<ncols;j++)
         {
            dxy=xy+actanode->sol.a.da[j][1];
	    minxy=DMIN(minxy,dxy);
	    maxxy=DMAX(maxxy,dxy);
         }
      }
   }
}
else
{
   actnode=&(actfield->dis[0].node[0]);
   minxy=actnode->x[1];
   maxxy=actnode->x[1];
   for (i=1;i<numnp;i++)
   {
      actnode=&(actfield->dis[0].node[i]);
      minxy=DMIN(minxy,actnode->x[1]);
      maxxy=DMAX(maxxy,actnode->x[1]);
   }
   
}
XYMIN[1] = minxy - (maxxy-minxy)*0.1;
XYMAX[1] = maxxy + (maxxy-minxy)*0.1;

dx = XYMAX[0] - XYMIN[0];
dy = XYMAX[1] - XYMIN[1];
dxy= sqrt(dx*dx+dy*dy);

/*------------------------------------------------------ window size */
hsize = 700;
vsize = 525;
ratio = 700.0/525.0;

/*------------------------ modify x/y - limits to window aspect ratio */
if (dx>=dy)
{
   XYPIX[0] = hsize;
   XYPIX[1] = vsize;
   if (dx>ratio*dy)
   {
      dy = dx/ratio;
      XYMIN[1] = 0.5*(XYMIN[1]+XYMAX[1]) - 0.5*dy;
      XYMAX[1] = XYMIN[1] + dx;
   }
   else
   {
      dx = dy*ratio;
      XYMIN[0] = 0.5*(XYMIN[0]+XYMAX[0]) - 0.5*dx;
      XYMAX[0] = XYMIN[0] + dx;
   }         
} /*endif (dx>=dy) */
else 
{
   XYPIX[0] = vsize;
   XYPIX[1] = hsize;
   if (dy>ratio*dx)
   {
      dx = dy/ratio;
      XYMIN[0] = 0.5*(XYMIN[0]+XYMAX[0]) - 0.5*dx;
      XYMAX[0] = XYMIN[0] + dy;
   }
   else
   {
      dy = dx*ratio;
      XYMIN[1] = 0.5*(XYMIN[1]+XYMAX[1]) - 0.5*dy;
      XYMAX[1] = XYMIN[1] + dy;
   }   
} /* endif else */
XYMIN[1] = XYMIN[1]/yscale;
XYMAX[1] = XYMAX[1]/yscale;

/*----------------------------------------------- get the data limits */
actnode=&(actfield->dis[0].node[0]);
velx    = actnode->sol.a.da[0][0];
vely    = actnode->sol.a.da[0][1];
absv    = sqrt(velx*velx+vely*vely);
pres    = actnode->sol.a.da[0][2];
if (IVORT==1) 
{
   vort = actnode->sol.a.da[0][3];
}

minvx   = velx;
maxvx   = velx;
minvy   = vely;
maxvy   = vely; 
minpre  = pres;
maxpre  = pres;
minabsv = sqrt(velx*velx+vely*vely);
maxabsv = sqrt(velx*velx+vely*vely);
if (IVORT==1)
{
   minvort = vort;
   maxvort = vort;
}

for (i=0;i<numnp;i++) /* loop nodes */
{
   actnode=&(actfield->dis[0].node[i]);
   for (j=0;j<ncols;j++) /* loop columns in sol-history */
   {
      velx    = actnode->sol.a.da[j][0];
      vely    = actnode->sol.a.da[j][1];
      pres    = actnode->sol.a.da[j][2];
      absv    = sqrt(velx*velx+vely*vely);  
      if(IVORT==1)
      {
         vort = actnode->sol.a.da[j][3];
      }          
      minvx   = DMIN(minvx  ,velx);
      maxvx   = DMAX(maxvx  ,velx);
      minvy   = DMIN(minvy  ,vely);
      maxvy   = DMAX(maxvy  ,vely); 
      minpre  = DMIN(minpre ,pres);
      maxpre  = DMAX(maxpre ,pres);
      minabsv = DMIN(minabsv,absv);
      maxabsv = DMAX(maxabsv,absv);
      if (IVORT==1)
      { 
         minvort = DMIN(minvort,vort);
	 maxvort = DMAX(maxvort,vort);
      }
   } /* end of loop over columns */
} /* end of loop over nodes */

/* grid velocity not available yet !!!! */
mingv = 0.0;
maxgv = 1.0;

/*--------------------------------------- store max/min data in FLIMS */
if (maxpre==minpre) maxpre=minpre+1.0;
FLIMS[0][0]=minpre;
FLIMS[0][1]=maxpre;

FLIMS[1][0]=0.0;
FLIMS[1][1]=1.0;

if (minvort==maxvort) maxvort=minvort+1.0;
FLIMS[2][0]=minvort;
FLIMS[2][1]=maxvort;

if (minabsv==maxabsv) maxabsv=minabsv+1.0;
FLIMS[3][0]=dxy/maxabsv;
FLIMS[3][1]=0.0;
FLIMS[6][0]=minabsv;
FLIMS[6][1]=maxabsv;

if (minvx==maxvx) maxvx=minvx+1.0;
FLIMS[4][0]=minvx;
FLIMS[4][1]=maxvx;

if (minvx==maxvy) maxvy=minvy+1.0;
FLIMS[5][0]=minvy;
FLIMS[5][1]=maxvy;

FLIMS[7][0]=FLIMS[0][0];
FLIMS[7][1]=FLIMS[0][1];

FLIMS[8][0]=mingv;
FLIMS[8][1]=maxgv;

FLIMS[9][0]=ZERO;
FLIMS[9][1]=ZERO;

FLIMS[10][0]=ZERO;
FLIMS[10][1]=ZERO;

/*-------------------------- call Fortran routine which calls VISUAL2 */
v2call(&IOPT,&CMNCOL,&CMUNIT,
        XYPIX,XYMIN,XYMAX,
        &NKEYS,IKEYS,FKEYS,FLIMS,
	&MNODE,&MPTRI,&MPPTRI,
	&MFACE,&MPFACE,&MEDGE,&MPEDGE,&bgcolour);

#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of vis2caf */

/*!---------------------------------------------------------------------                                         
\brief set VISUAL2's constants and structures

<pre>                                                         genk 07/02      

use QAT2V2 to make VISUAL2 structures
this routine is called by VISUAL2 during the visualisation

</pre>  
\param  *KNODE     int   (o)    Number of nodes
\param	*PTRI	   int	 (o)	PolyTri pointers to nodes
\param	*CTRI	   int	 (o)	PolyTri pointers to neighbouring cells
\param	*ITRI	   int	 (o)	PolyTri pointers to original cells
\param	*KPTRI	   int	 (o)	Number of PolyTri strips
\param	*PPTRI	   int	 (o)	Pointers to ends of strips
\param	*PFACE	   int	 (o)	Face PolyLine pointers to nodes
\param	*KPFACE    int	 (o)	Number of face PolyLines
\param	*PPFACE    int	 (o)	Pointers to ends of face PolyLines
\param	*PEDGE	   int	 (o)	Edge PolyLine pointers to nodes
\param	*KPEDGE    int	 (o)	Number of Edge PolyLines
\param	*PPEDGE    int	 (o)	Pointers to ends of edge PolyLines
\return void                                                                       
------------------------------------------------------------------------*/
void v2data(int *KNODE, int *PTRI, int *CTRI, int *ITRI, int *KPPTRI, 
            int *PPTRI, int *PFACE, int *KPFACE, int *PPFACE, 
	    int *PEDGE, int *KPEDGE, int *PPEDGE) 
{
int KPTRI, KEDGE, KFACE;

#ifdef DEBUG 
dstrc_enter("V2DATA");
#endif

*KNODE = numnp;

/*------------------------------------------------ generate polystrips */
qat2v2(&PCELL[0][0],WCELL,&KCELL,WNODE,KNODE,WFACE,
       &MPTRI,&MPPTRI,&MEDGE,&MPEDGE,&MFACE,&MPFACE,
       PTRI,CTRI,ITRI,&KPTRI,PPTRI,KPPTRI,
       PEDGE,CEDGE,&KEDGE,PPEDGE,KPEDGE,
       PFACE,&KFACE,PPFACE,KPFACE); 

/*-------------------------------------------------------------- output */
printf("\n");
printf("   Number of Nodes           = %d \n",MNODE);
printf("   Number of PolyTri Strips  = %d \n",*KPPTRI);
printf("   Length of PolyTri Strips  = %d \n",KPTRI);
if (*KPPTRI!=0)
printf("   Ave PolyTri Strip Length  = %d \n",KPTRI/(*KPPTRI));
printf("   Number of PolyLine Strips = %d \n",*KPFACE);
printf("   Length of PolyLine Strips = %d \n",KFACE);
if (*KPFACE!=0)
printf("   Ave PolyLine Strip Length = %d \n",KFACE/(*KPFACE));
printf("   Number of edge PLn Strips = %d \n",*KPEDGE);
printf("   Length of edge PLn Strips = %d \n",KEDGE);
if (*KPEDGE!=0)
printf("   Ave Edge PLn Strip Length = %d \n",KEDGE/(*KPEDGE));
printf("\n");

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of V2DATA */

/*!---------------------------------------------------------------------                                         
\brief set VISUAL2's grid coordinates 

<pre>                                                         genk 07/02      

this routine is called by VISUAL2 during the visualisation

</pre>  
\param  *XY        float   (o)  xy coordinates of grid nodes		    
\param	*XYDELT    float   (o)  Maximum width (dx) and heigth dy of any cell
\param	*XYMINP    float   (o)  Minimum values of x,y in each PolyTri strip 
\param	*XYMAXP    float   (o)  Maximum values of x,y in each PolyTri strip 
\param	*XYMINF    float   (o)  Minimum values of x,y in each face PolyTri  
\param	*XYMAXF    float   (o)  Maximum values of x,y in each face PolyTri  
\param	*XYMINE    float   (o)  Minimum values of x,y in each edge PolyTri  
\param	*XYMAXE    float   (o)  Maximum values of x,y in each edge PolyTri  
\return void                                                                       
------------------------------------------------------------------------*/
void v2grid(float XY[][2], float *XYDELT, float *XYMINP, float *XYMAXP,
            float *XYMINF, float *XYMAXF, float *XYMINE, float *XYMAXE) 
{

int i;

#ifdef DEBUG 
dstrc_enter("V2GRID");
#endif

switch(actfieldtyp)
{
case fluid:
   if (NUMA>=0 && IOPT==0)
   {
      input4:
      v2_cursor(&true);  
      printf("Give the column number (min. 0; max. %d): ?\n",ncols-1);
      scanf("%d",&icol);      
      v2_cursor(&false);  
      if (icol<0 || icol>ncols-1)
      {
         printf("Column number out of range. New input!\n");
	 printf("\n");
	 goto input4;
      }
   }
   if (IOPT>=2)
   {
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);
         actgnode=actnode->gnode;
         actanode=actgnode->mfcpnode[NUMA];
	 if (actanode==NULL)
	 {
            XY[i][0] = actnode->x[0];
            XY[i][1] = actnode->x[1];	 
         }
	 else
	 {
	    XY[i][0] = actanode->x[0] + actanode->sol.a.da[icol][0];
	    XY[i][1] = actanode->x[1] + actanode->sol.a.da[icol][1];	   
         }
      }
   }
   else
   {
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);
         XY[i][0] = actnode->x[0];
         XY[i][1] = actnode->x[1];
      }
   }
break;
case structure:
   dserror("fieldtyp not implemented yet!");
case ale:
   dserror("fieldtyp not implemented yet!");
} /* end switch(actfieldtyp) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of V2GRID */

/*!---------------------------------------------------------------------                                         
\brief set VISUAL2's scalar function values

<pre>                                                         genk 07/02      

This routine is called by VISUAL2 during the visualisation.
It filters the scalar values out of the nodal sol-field at the actual
position icol.
icol is set every timestep in v2update


</pre>  
\param  *JKEY      int    (i)  Index of key			 
\param	*S         float  (o)  Scalar function at nodes or cells 
\return void                                                                       
\sa v2update

------------------------------------------------------------------------*/
void v2scal(int *JKEY, float *S) 
{

int i;
float vx,vy;

#ifdef DEBUG 
dstrc_enter("V2SCAL");
#endif

switch(actfieldtyp)
{
case fluid:
   if (IOPT==0) /* read in the column number */
   {
      input3:
      v2_cursor(&true);
      printf("Give the column number (min. 0; max. %d): ?\n",ncols-1);
      scanf("%d",&icol);      
      v2_cursor(&false);      
      if (icol<0 || icol>ncols-1)
      {
         printf("Column number out of range. New input!\n");
	 printf("\n");
	 goto input3;
      }
   }
   /*-------------------------------------------------- get scalar data */
   switch(*JKEY)
   {
   case 1: /* pressure */
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);
	 S[i]=actnode->sol.a.da[icol][2];
      }
   break;
   /*-------------------------------------------------------------------*/   
   case 2: /* ???? */
      printf("NOT IMPLEMTED SO FAR!!!!!!!!!!!!!!!!!!!!\n");
   break;
   /*-------------------------------------------------------------------*/   
   case 3: /* vorticity */
      if (IVORT==1)
      {
         for (i=0;i<numnp;i++)
         {
            actnode=&(actfield->dis[0].node[i]);
	    S[i]=actnode->sol.a.da[icol][3];
         }
      }
      else
      {
         printf("vorticity not calculated!!!! --> ZERO-field\n");
	 for (i=0;i<numnp;i++)
	 S[i]=ZERO;
      }        
   break;
   case 4: 
      printf("VECTOR DATA IN SCALAR ROUTINE ???????\n");
   break;   
   /*-------------------------------------------------------------------*/   
   case 5: /* horizontal velocity */
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);
	 S[i]=actnode->sol.a.da[icol][0];
      }      
   break;
   /*-------------------------------------------------------------------*/   
   case 6: /* vertical velocity */
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);
	 S[i]=actnode->sol.a.da[icol][1];
      }        
   break;
   /*-------------------------------------------------------------------*/  
   case 7: /* absolute velocity */
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);
	 vx=actnode->sol.a.da[icol][0];
	 vy=actnode->sol.a.da[icol][1];
	 S[i] = sqrt(vx*vx+vy*vy);
      }       
   break;
   /*-------------------------------------------------------------------*/   
   case 8: /* switch streamline - stationary streamline */
      if(isstr==0)
      {
         isstr=1;
	 v2_cursor(&true);
	 printf("Subtraction value for stationary streamlines: ?\n");
	 scanf("%f",&sstrval);
	 v2_cursor(&false);	 
      }
      else
         isstr=0;
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);
	 S[i] = actnode->sol.a.da[icol][2];
      }          
   break;
   /*-------------------------------------------------------------------*/    
   case 9:
      printf("VECTOR DATA IN SCALAR ROUTINE ???????\n");
   break;
   /*-------------------------------------------------------------------*/   
   case 10: /* stopping at step */
      printf("   Which step do you want to stop?\n");
      scanf("%d",&STSTEP); 
   break;
   /*-------------------------------------------------------------------*/
   case 11:
      if (IMOVIE==0)
      {
         printf("\n");
         printf("   Starting movie creation at time 0.0\n");
	 printf("\n");
         IMOVIE=1;
      }      
   /*-------------------------------------------------------------------*/                 
   } /* end switch(*JKEY) */
break;
case structure:
   dserror("fieldtyp not implemented yet!\n");
case ale:
   dserror("fieldtyp not implemented yet!\n");
} /* end switch(actfieldtyp) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of V2SCAL */

/*!---------------------------------------------------------------------                                         
\brief set VISUAL2's vector function values

<pre>                                                         genk 07/02      

This routine is called by VISUAL2 during the visualisation.
It filters the vector values out of the nodal sol-field at the actual
position icol.
icol is set every timestep in v2update

</pre>  
\param  *JKEY      int    (i)  Index of key			 
\param	*V         float  (o)  Vector function at nodes or cells 
\return void                                                                       
\sa v2update

------------------------------------------------------------------------*/
void v2vect(int *JKEY, float V[][2]) 
{

int i;
float vx,vy;

#ifdef DEBUG 
dstrc_enter("V2VECT");
#endif

switch(actfieldtyp)
{
case fluid:
   /*-------------------------------------------------- get vector data */
   if (*JKEY==4)
   {
      switch(isstr)
      {
      case 0: /* streamlines */
         for (i=0;i<numnp;i++)
         {
            actnode=&(actfield->dis[0].node[i]);
	    V[i][0] = actnode->sol.a.da[icol][0];
	    V[i][1] = actnode->sol.a.da[icol][1];
         }
         break;
   /*-------------------------------------------------------------------*/   
      default: /* stationary streamlines */
         printf("stationary streamlines not checked yet!!!!");
         for (i=0;i<numnp;i++)
         {
            actnode=&(actfield->dis[0].node[i]);
	    actgnode=actnode->gnode;
	    if (actgnode->dirich!=NULL) 
	    {        
	       if (actgnode->dirich->dirich_val.a.dv[0]==0.0)
	       {
	          V[i][0] = actnode->sol.a.da[icol][0];
	       }
	       else
	       {
	          V[i][0] = actnode->sol.a.da[icol][0]-sstrval;
	       }
	       V[i][1] = actnode->sol.a.da[icol][1];
	    }
	    else
	    {
	       V[i][0] = actnode->sol.a.da[icol][0]-sstrval;
	       V[i][1] = actnode->sol.a.da[icol][1];	 
	    }	 
         }      
         break;                 
      }
   }
break;
case structure:
   dserror("fieldtyp not implemented yet!\n");
case ale:
   dserror("fieldtyp not implemented yet!\n");
} /* end switch(actfieldtyp) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of V2VECT */

/*!---------------------------------------------------------------------                                         
\brief update time and icol 

<pre>                                                         genk 07/02      

This routine is called by VISUAL2 during the visualisation.
Beside the actual time the actual position in the nodal sol-history 
icol is also updated

</pre>  
\param  *TIME      float    (o)  actual time		 
\return void                                                                       

------------------------------------------------------------------------*/
void v2update(float *TIME) 
{
int i;

#ifdef DEBUG 
dstrc_enter("V2UPDATE");
#endif


/*-------------------------------------------------- check for stopping */
if (icol==STSTEP)
{
   printf("   Next Step to Stop (-1 to terminate)\n");
   scanf("%d",&STSTEP);        
}

/*------------------------------------------- check for movie creation */
if (icol==0 && IMOVIE==2) 
{   
   IMOVIE=0;
   printf("\n");
   printf("Movie creation finished\n");
}   
if (icol==0 && IMOVIE==1) IMOVIE=2;
if (IMOVIE==2) v2movie();


switch(actfieldtyp)
{
case fluid: 
   if (icol==-1 || icol+DSTEP > ncols-1) icol=0;
   else icol+=DSTEP;
   *TIME = time_a.a.dv[icol];
break;   
case structure:
   dserror("fieldtyp not implemented yet!\n");
case ale:
   dserror("fieldtyp not implemented yet!\n");
} /* end switch(actfieldtyp) */


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of V2UPDATE */

/*!---------------------------------------------------------------------                                         
\brief  additional label for  VISUAL2 plots 

<pre>                                                         genk 07/02      

This routine is called by VISUAL2 during the visualisation.

</pre>  
\param  STRING     char    (o)  plot on VISUAL2 window		 
\return void                                                                       

------------------------------------------------------------------------*/
void v2string(char STRING[81]) 
{

char *charpointer;
float t;

#ifdef DEBUG 
dstrc_enter("VSTRING");
#endif

switch(actfieldtyp)
{
case fluid:
   t=time_a.a.dv[icol];  
   strcpy(STRING,"Time: ");
   charpointer=STRING+strlen(STRING);
   sprintf(charpointer,"%8.4f",t);
   charpointer=STRING+strlen(STRING);
   strcpy(charpointer,"    Step: ");
   charpointer=STRING+strlen(STRING);
   sprintf(charpointer,"%5d",icol);  
   charpointer=STRING+strlen(STRING); 
   strcpy(charpointer,"/");
   charpointer=STRING+strlen(STRING); 
   sprintf(charpointer,"%-5d",ncols-1); 
   charpointer=STRING+strlen(STRING);      
   strcpy(charpointer,"                                             ");      
break;   
case structure:
   dserror("fieldtyp not implemented yet!\n");
case ale:
   dserror("fieldtyp not implemented yet!\n");
} /* end switch(actfieldtyp) */


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of V2STRING */

/*!---------------------------------------------------------------------                                         
\brief  compute the array WCELL for use in QAT2V2 

<pre>                                                         genk 07/02      
</pre>  
\param  *actfield     FIELD    (i)  actual field		 
\return void                                                                       
\warning QAT2V2 requires FORTRAN numbering, so increase node-numbers by one!!

------------------------------------------------------------------------*/
void v2cell(FIELD *actfield) 
{
int i,j,k;
int inel;

#ifdef DEBUG 
dstrc_enter("v2cell");
#endif

inel=-1;

switch (distyp)
{
case tri3: /* 3 node triangle */
   for (i=0;i<numele;i++)
   {
      actele=&(actfield->dis[0].element[i]); 
      inel++;
      for(j=0;j<3;j++) 
      {
         PCELL[inel][j] = actele->node[j]->Id_loc+1;
      }
      PCELL[inel][3] = PCELL[inel][0];
   }   
   break;
case tri6:
   for (i=0;i<numele;i++)
   {
      actele=&(actfield->dis[0].element[i]); 
      inel++;
      /*----------------------------------------- sub-element 1 */
      PCELL[inel][0] = actele->node[0]->Id_loc+1;    
      PCELL[inel][1] = actele->node[3]->Id_loc+1;    
      PCELL[inel][2] = actele->node[5]->Id_loc+1;    
      PCELL[inel][3] = actele->node[0]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 2 */
      PCELL[inel][0] = actele->node[3]->Id_loc+1;    
      PCELL[inel][1] = actele->node[1]->Id_loc+1;    
      PCELL[inel][2] = actele->node[4]->Id_loc+1;    
      PCELL[inel][3] = actele->node[3]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 3 */
      PCELL[inel][0] = actele->node[3]->Id_loc+1;    
      PCELL[inel][1] = actele->node[4]->Id_loc+1;    
      PCELL[inel][2] = actele->node[5]->Id_loc+1;    
      PCELL[inel][3] = actele->node[3]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 4 */
      PCELL[inel][0] = actele->node[4]->Id_loc+1;    
      PCELL[inel][1] = actele->node[3]->Id_loc+1;    
      PCELL[inel][2] = actele->node[5]->Id_loc+1;    
      PCELL[inel][3] = actele->node[4]->Id_loc+1;
    } 
   break;
case quad4: /* 4 node rectangle */
   for (i=0;i<numele;i++)
   {
      actele=&(actfield->dis[0].element[i]); 
      inel++;
      for(j=0;j<4;j++) 
      {
         PCELL[inel][j] = actele->node[j]->Id_loc+1;
      }
   }
   break;
case quad8:
   for (i=0;i<numele;i++)
   {
      actele=&(actfield->dis[0].element[i]); 
      inel++;
      /*----------------------------------------- sub-element 1 */
      PCELL[inel][0] = actele->node[0]->Id_loc+1;    
      PCELL[inel][1] = actele->node[4]->Id_loc+1;    
      PCELL[inel][2] = actele->node[7]->Id_loc+1;    
      PCELL[inel][3] = actele->node[0]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 2 */
      PCELL[inel][0] = actele->node[4]->Id_loc+1;    
      PCELL[inel][1] = actele->node[1]->Id_loc+1;    
      PCELL[inel][2] = actele->node[5]->Id_loc+1;    
      PCELL[inel][3] = actele->node[4]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 3 */
      PCELL[inel][0] = actele->node[5]->Id_loc+1;    
      PCELL[inel][1] = actele->node[2]->Id_loc+1;    
      PCELL[inel][2] = actele->node[6]->Id_loc+1;    
      PCELL[inel][3] = actele->node[5]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 4 */
      PCELL[inel][0] = actele->node[6]->Id_loc+1;    
      PCELL[inel][1] = actele->node[3]->Id_loc+1;    
      PCELL[inel][2] = actele->node[7]->Id_loc+1;    
      PCELL[inel][3] = actele->node[6]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 5 */
      PCELL[inel][0] = actele->node[6]->Id_loc+1;    
      PCELL[inel][1] = actele->node[7]->Id_loc+1;    
      PCELL[inel][2] = actele->node[4]->Id_loc+1;    
      PCELL[inel][3] = actele->node[5]->Id_loc+1;
   }
   break;
case quad9:
   for (i=0;i<numele;i++)
   {
      actele=&(actfield->dis[0].element[i]); 
      inel++;
      /*----------------------------------------- sub-element 1 */
      PCELL[inel][0] = actele->node[0]->Id_loc+1;    
      PCELL[inel][1] = actele->node[4]->Id_loc+1;    
      PCELL[inel][2] = actele->node[8]->Id_loc+1;    
      PCELL[inel][3] = actele->node[7]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 2 */
      PCELL[inel][0] = actele->node[4]->Id_loc+1;    
      PCELL[inel][1] = actele->node[1]->Id_loc+1;    
      PCELL[inel][2] = actele->node[5]->Id_loc+1;    
      PCELL[inel][3] = actele->node[8]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 3 */
      PCELL[inel][0] = actele->node[8]->Id_loc+1;    
      PCELL[inel][1] = actele->node[5]->Id_loc+1;    
      PCELL[inel][2] = actele->node[2]->Id_loc+1;    
      PCELL[inel][3] = actele->node[6]->Id_loc+1;
      inel++;
      /*----------------------------------------- sub-element 4 */
      PCELL[inel][0] = actele->node[7]->Id_loc+1;    
      PCELL[inel][1] = actele->node[8]->Id_loc+1;    
      PCELL[inel][2] = actele->node[6]->Id_loc+1;    
      PCELL[inel][3] = actele->node[3]->Id_loc+1;
   }
   break;           
default:
   dserror("disytp not implemented yet!\n");
}


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of v2cell*/

/*!---------------------------------------------------------------------                                         
\brief  movie creation 

<pre>                                                         genk 03/03      

this function							   
creates the file vis<icol>.xwd via system command (screen-shot     
xwd -out ...)							   
converts the xwd-file into a .gif file via system command (convert)
deletes the xwd-file via system command (rm)			   

</pre>  

\return void                                                                       
------------------------------------------------------------------------*/
void v2movie()
{
char string[100]=("xwd -name Visual2X -out vis          ");
char *charpointer;

char convert[100]=("convert vis                         ");
char remove[100] =("rm vis                               ");
char screen[100] =("Storing vis                          ");

#ifdef DEBUG 
dstrc_enter("v2movie");
#endif

/*--------------------------------------- add file counter to filenames */
sprintf(&string[27] ,"%-d" ,icol);
sprintf(&convert[11],"%-d" ,icol);
sprintf(&remove[6]  ,"%-d" ,icol);
sprintf(&screen[11 ],"%-d" ,icol);

/*---------------------------------------- add extensions to file names */
charpointer = string+strlen(&string[0]);
strncpy(charpointer,".xwd",4);

charpointer = remove+strlen(&remove[0]);
strncpy(charpointer,".xwd",4);

charpointer = screen+strlen(&screen[0]);
strncpy(charpointer,".gif\n",6);

charpointer = convert+strlen(&convert[0]);
strncpy(charpointer,".xwd vis",8);
charpointer +=8;
sprintf(charpointer,"%-d" ,icol);
charpointer = convert+strlen(&convert[0]);
strncpy(charpointer,".gif",4);

printf(screen);

/*---------------------------------------------------- call UNIX-system */
system(&string[0]);
system(&convert[0]);
system(&remove[0]);

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of v2movie*/

#else
/*!---------------------------------------------------------------------                                         
\brief dummy routine 

<pre>                                                         genk 07/02      

since VISUAL2 is called by a fortran routine, this one is necessary, if
VIS2-routines are not compiled into the program

</pre>  		 
\return void                                                                       

------------------------------------------------------------------------*/
void v2_init(
             char *titl, int *iopt, int *cmncol, char *cmfile, int *cmunit,
	     int *xypix, float *xymin, float *xymax,
	     int *nkeys, int *ikeys, int *tkeys, int *fkeys, float **flims,
	     int *mnode, int *mptri, int *mpptri,
	     int *mface, int *mpface, int *medge, int *mpedge
	    )   
{
dserror("VISUAL2 PACKAGE not compiled into programm\n");
return;
}
#endif
