/*!----------------------------------------------------------------------
\file
\brief Call and control VISUAL2 

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

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
/*----------------------------------------------------------------------*
 |                                                        genk 07/02    |
 | all variables needed by VISUAL2 are defined extern                   |
 *----------------------------------------------------------------------*/
static INT      MPTRI,MPPTRI,MFACE;   /* see VISUAL2 manual             */
static INT      MPFACE,MEDGE,MPEDGE;  /* see VISUAL2 manual             */
static INT      KCELL;                /* number of visualised cells - 
                                 different from number of elements!     */
static INT      IOPT;                 /* program mode                   */
static INT      IVORT;                /* flag for vorticity calculation */
static INT      ITURBU;               /* flag for turbulence calc.      */
static INT      SCAL;
static INT      DSTEP;                /* increment of visualised steps  */
static INT      INCRE;                /* increment of visualised steps  */
static INT      NUMA;                 /* number of ALE-FIELD            */
static INT      DATAFILE;             /* file-flag for results          */
static INT      CMNCOL=300;           /* see VISUAL2 manual             */
static INT      CMUNIT=37;            /* see VISUAL2 manual             */
static INT      MNODE;                /* maximum number of nodes        */
static INT      NKEYS=14;              /* see VISUAL2 manual             */
static INT      XYPIX[2];             /* size of VISUAL2 window         */
static INT      IKEYS[]={112,115,118,102,120,121,97,116,102,92,109,107,100,110};
static INT      FKEYS[]={1,1,1,3,1,1,1,1,3,1,1,1,1,1};             
static float    FLIMS[14][2];          /* data limits                    */
static float    XYMIN[2], XYMAX[2];   /* min. and max. coordinates      */
static struct  _ARRAY PCELL_A;
static struct  _ARRAY WCELL_A;
static struct  _ARRAY WNODE_A;
static struct  _ARRAY WFACE_A;
static struct  _ARRAY CEDGE_A;        /* arrays needed in qat2v2        */
static INT    **PCELL;
static INT     *WCELL;
static INT     *WNODE;
static INT     *WFACE; 
static INT     *CEDGE;                /* pointers to arrays             */
/*------------------------------- variables needed for own calculations */
static INT     numnp;	        /* number of nodes of actual field	*/
static INT     numele;	        /* number of elements of actual field	*/
static INT     ncols=1;         /* number of sol steps stored in sol	*/
static INT     nacols=1;        /* number of sol steps of ALE field     */
static INT     isstr=0;         /* flag for streamline calculation	*/
static INT     icol=-1;         /* act. num. of sol step to be visual.  */
static INT     FIRSTSTEP;
static INT     LASTSTEP;
static INT     ACTSTEP;
static INT     hsize, vsize;    /* size for VISUAL2 window		*/
static INT     true=-1;         /* flag for v2_cursor			*/
static INT     false=0;         /* flag for v2_cursor			*/
static INT     STSTEP=-2;       /* stopping step                        */
static INT     IMOVIE=0;        /* counter to slow down for movie creat.*/
static INT     bgcolour;        /* background colour                    */
static float   sstrval=0.0;     /* stationary streamline value		*/
static float   dx,dy,dxy;       /* for determination of coord. limits*/
static float   ratio;	        /*					*/
static DOUBLE  yscale;	        /* scaling factors for y-direction	*/
static DOUBLE  velx;	        /*					*/
static DOUBLE  vely;	        /*					*/
static DOUBLE  pres;	        /*					*/
static DOUBLE  absv;	        /*					*/
static DOUBLE  vort;            /*                                      */
static DOUBLE  kappa;	        /*					*/
static DOUBLE  dissi;	        /*					*/
static DOUBLE  visco;	        /*					*/
static DOUBLE  xy;              /*					*/
static DOUBLE  minxy,maxxy;     /*					*/
static DOUBLE  centerx,centery; /*                                      */
static DOUBLE  minpre,maxpre;   /*					*/
static DOUBLE  minvx,maxvx;     /*					*/
static DOUBLE  minvy,maxvy;     /*					*/
static DOUBLE  minabsv,maxabsv; /*					*/
static DOUBLE  minvort,maxvort; /*					*/
static DOUBLE  minkappa,maxkappa; /*					*/
static DOUBLE  mindissi,maxdissi; /*					*/
static DOUBLE  minvisco,maxvisco; /*					*/
static DOUBLE  mingv,maxgv;     /* for determination of data limits	*/
static DOUBLE  minvort,maxvort; /*					*/
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
ARRAY          step_a ;         /* time array				*/
FLUID_DYNAMIC *fdyn;                /* pointer to fluid dyn. inp.data   */
FLUID_DYN_CALC *dynvar;             /* pointer to fluid_dyn_calc        */ 
/*!---------------------------------------------------------------------                                         
\brief call of visual2 for fluid 

<pre>                                                         genk 07/02      
</pre>  
\param  numf   INT      (i)       actual number of fluid field
\return void                                                                       

------------------------------------------------------------------------*/
void vis2caf(INT numff, INT numaf, INT numsf) 
{

INT iscan, i,j,k;
INT screen, dummy;
CALC_ACTION    *action;             /* pointer to the cal_action enum   */
INTRA          *actintra;           /* pointer to active intra-communic.*/
PARTITION      *actpart;            /* pointer to active partition      */
CONTAINER       container;          /* contains variables defined in container.h */
FLUID_DYNAMIC *fdyn;                /* pointer to fluid dyn. inp.data   */

#ifdef DEBUG 
dstrc_enter("vis2caf");
#endif

actfield=&(field[numff]);
actfieldtyp=fluid;
fdyn = alldyn[genprob.numff].fdyn;

actpart= &(partition[numff]);
action = &(calc_action[numff]);
NUMA=numaf;

/*------ read solution data from pss-file or flavia.res-file of proc 0 */
input0:
printf("\n");
printf("     Where are the visualisation data stored?\n");
printf("      0 :   pss-file\n");
printf("     [1]:   flavia.res-file\n");
screen=getchar();
switch(screen)
{
case 10: DATAFILE=1; break;
case 48: DATAFILE=0; dummy=getchar(); break;
case 49: DATAFILE=1; dummy=getchar(); break;
default: 
   printf("\nTry again!\n");
   goto input0;
}

if (DATAFILE==0)
{
   visual_readpss(actfield,&ncols,&time_a);
   FIRSTSTEP=0;
   LASTSTEP=ncols;
}
else if (DATAFILE==1)
{
   visual_readflaviares(actfield,&ncols,&time_a,&step_a,
                        &FIRSTSTEP,&LASTSTEP,&DSTEP);
   INCRE=1;
}
else 
   goto input0;

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
printf("     Actual size of data structures:\n");
printf("      MNODE  = %d\n",MNODE);
printf("      MPTRI  = %d\n",MPTRI);
printf("      MPPTRI = %d\n",MPPTRI);
printf("      MFACE  = %d\n",MFACE);
printf("      MPFACE = %d\n",MPFACE);
printf("      MEDGE  = %d\n",MEDGE);
printf("      MPEDGE = %d\n",MPEDGE);
printf("\n");
printf("     New size of data structures?\n");
printf("     [0]:   no\n");
printf("      1 :   yes\n");
screen=getchar();
switch(screen)
{
case 10: iscan=0; break;
case 48: iscan=0; dummy=getchar(); break;
case 49: iscan=1; dummy=getchar(); break;
default: 
   printf("\nTry again!\n");
   goto input1;
}

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
input2:
printf("\n");
printf("     Do you want to calculate the vorticity?\n");
printf("     [0]:   no\n");
printf("      1 :   yes\n");
screen=getchar();
switch(screen)
{
case 10: IVORT=0; break;
case 48: IVORT=0; dummy=getchar(); break;
case 49: IVORT=1; dummy=getchar(); break;
default: 
   printf("\nTry again!\n");
   goto input2;
}

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
   container.actndis=0;
   container.nif=0;
   container.nii=0;
   container.nim=0;
   container.gen_alpha=0;
   container.fieldtyp=fluid;
   container.dvec=NULL;
   container.is_relax     = 0;
   container.dvec         = NULL;
   container.actndis  = 0;   
   calelm(actfield,NULL,actpart,NULL,0,0,
          &container,action); 
} /* endif (IVORT==1) */

input3:
printf("\n");
printf("     Please give the mode in which you want to run VISUAL2:\n");
printf("      0 : steady data structure, grid and variables\n");
printf("      1 : steady data structure and grid, unsteady variables\n");
printf("     [2]: steady data structure, unsteady grid and variables\n");
printf("      3 : unsteady data structure, grid and variables\n");
screen=getchar();
switch(screen)
{
case 10: IOPT=2; break;
case 48: IOPT=0; dummy=getchar(); break;
case 49: IOPT=1; dummy=getchar(); break;
case 50: IOPT=2; dummy=getchar(); break;
case 51: IOPT=3; dummy=getchar(); break;
default: 
   printf("\nTry again!\n");
   goto input3;
}

if (IOPT==0) icol++;
if (IOPT==3)
{
   printf("   sorry, mode not implemented yet - new input!\n");
   goto input3;
}
if (IOPT==2) /*-------------------------- read ALE field from pss-file */
{
#ifdef D_FSI
   if (genprob.probtyp==prb_fsi)
   {
      if (numaf<0) dserror("ALE-field does not exist!");
      if (numsf<0) dserror("STRUCTURE-field does not exist!");
      alefield=&(field[numaf]);
      structfield=&(field[numsf]);
      fsi_initcoupling(structfield,actfield,alefield);
      if (DATAFILE==0)
         visual_readpss(alefield,&nacols,NULL);
      else if (DATAFILE==1)
         visual_readflaviares(alefield,&nacols,NULL,NULL,&FIRSTSTEP,&LASTSTEP,&DSTEP);
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
      if (DATAFILE==0)
         visual_readpss(alefield,&nacols,NULL);
      else if (DATAFILE==1)
         visual_readflaviares(alefield,&nacols,NULL,NULL,&FIRSTSTEP,&LASTSTEP,&DSTEP);
      if (ncols!=nacols)
      {
         printf("\n");
	 printf("WARNING:number of columns different in ALE and FLUID field \n");
         printf("\n");
         ncols=IMIN(ncols,nacols);
      }
   }
#else
dserror("FSI functions not compiled in!\n");
#endif
}
/*--------------------------------------------------------- y-scaling */
input4:
printf("\n");
printf("     Scaling Geometry?\n");
printf("     [0]: no \n");
printf("      1 : yes\n");

screen=getchar();
switch(screen)
{
case 10: SCAL=0; break;
case 48: SCAL=0; dummy=getchar(); break;
case 49: SCAL=1; dummy=getchar(); break;
default: 
   printf("\nTry again!\n");
   goto input4;
}
if (SCAL==1)
{
   printf("\n");
   printf("     Scaling parameter for y-coordinates?\n");
   scanf("%lf",&yscale);
}
else
   yscale=ONE;

/*--------------------------------------------------------- increment */
input5:
if (DATAFILE==0)
{
   printf("\n");
   printf("   Increment step number by ...? (<%d)\n",ncols);
   scanf("%d",&INCRE);
   DSTEP=INCRE;
   if (INCRE>ncols-1)
   {
      printf("   Increment too large --> Try again!\n");
      goto input5;
   }
}
/*--------------------------------------------------- modiy last step */
else
   LASTSTEP=step_a.a.iv[ncols-1];

/*------------------------------------------------- background colour */
input6:
printf("\n");
printf("     Colours? \n");
printf("     [0]: colours    - black background\n");
printf("      1 : colours    - white background\n");
printf("      2 : grey scale - white background\n");
screen=getchar();
switch(screen)
{
case 10: bgcolour=0; break;
case 48: bgcolour=0; dummy=getchar(); break;
case 49: bgcolour=1; dummy=getchar(); break;
case 50: bgcolour=2; dummy=getchar(); break;
default: 
   printf("\nTry again!\n");
   goto input6;
}
#ifdef SUSE73
bgcolour+=10;
#endif
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
#ifdef D_FSI
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
#else
   dserror("FSI-functions not compiled in!\n");
#endif   
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
#ifdef D_FSI
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
#else
   dserror("FSI-functions not compiled in!\n");
#endif   
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
centerx =XYMIN[0]+dx/TWO;
centery =XYMIN[1]+dy/TWO;
if (dx>=dy)
{
   XYPIX[0] = hsize;
   XYPIX[1] = vsize;
   if (dx/dy<ratio)
   {
      dx = dy*ratio;
      XYMIN[0] = centerx - dx/TWO;
      XYMAX[0] = centerx + dx/TWO;           
   }
   else
   {
      dy = dx/ratio;
      XYMIN[1] = centery - dy/TWO;
      XYMAX[1] = centery + dy/TWO;
   }
}
else
{
   XYPIX[0]=vsize;
   XYPIX[1]=hsize;
   if (dy/dx<ratio)
   {
      dy = dx*ratio;
      XYMIN[1] = centery - dy/TWO;
      XYMAX[1] = centery + dy/TWO;
   }
   else
   {
      dx = dy/ratio;
      XYMIN[0] = centerx - dx/TWO;
      XYMAX[0] = centerx + dx/TWO;
   }
}
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

if (fdyn->turbu == 2 || fdyn->turbu == 3)
{
 ITURBU = 1;
 for (i=0;i<numnp;i++) /* loop nodes */
 {
   actnode=&(actfield->dis[1].node[i]);
   for (j=0;j<ncols;j++) /* loop columns in sol-history */
   {
      kappa      = actnode->sol.a.da[j][0];
      dissi      = actnode->sol.a.da[j][2];
      visco      = actnode->sol.a.da[j][1];

      minkappa   = DMIN(minkappa   ,kappa);
      maxkappa   = DMAX(maxkappa   ,kappa);
      mindissi   = DMIN(mindissi   ,dissi);
      maxdissi   = DMAX(maxdissi   ,dissi); 
      minvisco   = DMIN(minvisco   ,visco);
      maxvisco   = DMAX(maxvisco   ,visco); 
   } /* end of loop over columns */
 } /* end of loop over nodes */
} /* end of container.turbu==2 || container.turbu==3 */
else
{
 ITURBU = 0;
      minkappa   = 0.0;
      maxkappa   = 1.0;
      mindissi   = 0.0;
      maxdissi   = 1.0; 
      minvisco   = 0.0;
      maxvisco   = 1.0; 
}
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

if (minvy==maxvy) maxvy=minvy+1.0;
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

if (minkappa==maxkappa) maxkappa=minkappa+1.0;
FLIMS[11][0]=minkappa;
FLIMS[11][1]=maxkappa;

if (mindissi==maxdissi) maxdissi=mindissi+1.0;
FLIMS[12][0]=mindissi;
FLIMS[12][1]=maxdissi;

if (minvisco==maxvisco) maxvisco=minvisco+1.0;
FLIMS[13][0]=minvisco;
FLIMS[13][1]=maxvisco;

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
\param  *KNODE     INT   (o)    Number of nodes
\param	*PTRI	   INT	 (o)	PolyTri pointers to nodes
\param	*CTRI	   INT	 (o)	PolyTri pointers to neighbouring cells
\param	*ITRI	   INT	 (o)	PolyTri pointers to original cells
\param	*KPTRI	   INT	 (o)	Number of PolyTri strips
\param	*PPTRI	   INT	 (o)	Pointers to ends of strips
\param	*PFACE	   INT	 (o)	Face PolyLine pointers to nodes
\param	*KPFACE    INT	 (o)	Number of face PolyLines
\param	*PPFACE    INT	 (o)	Pointers to ends of face PolyLines
\param	*PEDGE	   INT	 (o)	Edge PolyLine pointers to nodes
\param	*KPEDGE    INT	 (o)	Number of Edge PolyLines
\param	*PPEDGE    INT	 (o)	Pointers to ends of edge PolyLines
\return void                                                                       
------------------------------------------------------------------------*/
void v2data(INT *KNODE, INT *PTRI, INT *CTRI, INT *ITRI, INT *KPPTRI, 
            INT *PPTRI, INT *PFACE, INT *KPFACE, INT *PPFACE, 
	    INT *PEDGE, INT *KPEDGE, INT *PPEDGE) 
{
INT KPTRI, KEDGE, KFACE;

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

INT i;

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
#ifdef D_FSI
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
#else
   dserror("FSI-functions not compiled in!\n");
#endif   
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
\param  *JKEY      INT    (i)  Index of key			 
\param	*S         float  (o)  Scalar function at nodes or cells 
\return void                                                                       
\sa v2update

------------------------------------------------------------------------*/
void v2scal(INT *JKEY, float *S) 
{

INT i;
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
   break;
   /*-------------------------------------------------------------------*/   
   case 12: /* kappa */
      if (ITURBU==1)
      {
       for (i=0;i<numnp;i++)
       {
        actnode=&(actfield->dis[1].node[i]);
	  S[i]=actnode->sol.a.da[icol][0];
       }
      }        
      else
      {
       printf("kappa not calculated!!!! --> ZERO-field\n");
	 for (i=0;i<numnp;i++)
	 S[i]=ZERO;
      }        
   break;
   /*-------------------------------------------------------------------*/  
   case 13: /* dissi */
      if (ITURBU==1)
      {
       for (i=0;i<numnp;i++)
       {
        actnode=&(actfield->dis[1].node[i]);
	  S[i] = actnode->sol.a.da[icol][2];
       }       
      }        
      else
      {
       printf("dissipation not calculated!!!! --> ZERO-field\n");
	 for (i=0;i<numnp;i++)
	 S[i]=ZERO;
      }        
   break;
   /*-------------------------------------------------------------------*/  
   case 14: /* visco */
      if (ITURBU==1)
      {
       for (i=0;i<numnp;i++)
       {
        actnode=&(actfield->dis[1].node[i]);
	  S[i] = actnode->sol.a.da[icol][1];
       }       
      }        
      else
      {
       printf("eddy-viscosity with ke- or kw-model not calculated!!!! --> ZERO-field\n");
	 for (i=0;i<numnp;i++)
	 S[i]=ZERO;
      }        
   break;
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
\param  *JKEY      INT    (i)  Index of key			 
\param	*V         float  (o)  Vector function at nodes or cells 
\return void                                                                       
\sa v2update

------------------------------------------------------------------------*/
void v2vect(INT *JKEY, float V[][2]) 
{

INT i;
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
INT i;

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
   if (icol==-1 || icol+INCRE > ncols-1) 
   {
      icol=0;
      if (DATAFILE==0)
         ACTSTEP=FIRSTSTEP;
      else         
         ACTSTEP=step_a.a.iv[icol];
   }
   else
   { 
      icol+=INCRE;
      if (DATAFILE==0)
         ACTSTEP+=DSTEP;
      else
         ACTSTEP=step_a.a.iv[icol];
   }
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
dstrc_enter("V2STRING");
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
   sprintf(charpointer,"%5d",ACTSTEP);  
   charpointer=STRING+strlen(STRING); 
   strcpy(charpointer,"/");
   charpointer=STRING+strlen(STRING); 
   sprintf(charpointer,"%-5d",LASTSTEP); 
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
INT i,j,k;
INT inel;

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
creates the file vis<ACTSTEP>.xwd via system command (screen-shot     
xwd -out ...)							   
converts the xwd-file into a .gif file via system command (convert)
deletes the xwd-file via system command (rm)			   

see also /bau/statik/info/README.visual2_ccarat_movie

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
#ifdef HPUX11
printf("movie creation ony possible under HPUX 10.20\n");
goto end;
#endif

sprintf(&string[27] ,"%-d" ,ACTSTEP);
sprintf(&convert[11],"%-d" ,ACTSTEP);
sprintf(&remove[6]  ,"%-d" ,ACTSTEP);
sprintf(&screen[11 ],"%-d" ,ACTSTEP);

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

/*--------------------------------------------------- call UNIX-system */
system(&string[0]);
system(&convert[0]);
system(&remove[0]);

/*----------------------------------------------------------------------*/
end:
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
             char *titl, INT *iopt, INT *cmncol, char *cmfile, INT *cmunit,
	     INT *xypix, float *xymin, float *xymax,
	     INT *nkeys, INT *ikeys, INT *tkeys, INT *fkeys, float **flims,
	     INT *mnode, INT *mptri, INT *mpptri,
	     INT *mface, INT *mpface, INT *medge, INT *mpedge
	    )   
{
dserror("VISUAL2 PACKAGE not compiled into programm\n");
return;
}
void v2_init_(
             char *titl, INT *iopt, INT *cmncol, char *cmfile, INT *cmunit,
	     INT *xypix, float *xymin, float *xymax,
	     INT *nkeys, INT *ikeys, INT *tkeys, INT *fkeys, float **flims,
	     INT *mnode, INT *mptri, INT *mpptri,
	     INT *mface, INT *mpface, INT *medge, INT *mpedge
	    )   
{
dserror("VISUAL2 PACKAGE not compiled into programm\n");
return;
}
void v2_init__(
             char *titl, INT *iopt, INT *cmncol, char *cmfile, INT *cmunit,
	     INT *xypix, float *xymin, float *xymax,
	     INT *nkeys, INT *ikeys, INT *tkeys, INT *fkeys, float **flims,
	     INT *mnode, INT *mptri, INT *mpptri,
	     INT *mface, INT *mpface, INT *medge, INT *mpedge
	    )   
{
dserror("VISUAL2 PACKAGE not compiled into programm\n");
return;
}
#endif
