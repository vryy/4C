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
#ifdef VISUAL3_PACKAGE
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
static INT      KSURF;              
static INT      KCEL1,KCEL2,KCEL3,KCEL4;
static INT      WIN3D,HWIN3D;
static INT      MIRROR;
static INT      LIMITS;
static INT      IOPT;                 /* program mode                   */
static INT      DSTEP;                /* increment of visualised steps  */
static INT      INCRE;                /* increment of visualised steps  */
static INT      NUMA;                 /* number of ALE-FIELD            */
static INT      DATAFILE;             /* file-flag for results          */
static INT      CMNCOL=300;           /* see VISUAL2 manual             */
static INT      CMUNIT=37;            /* see VISUAL2 manual             */
static INT      MNODE;                /* maximum number of nodes        */
static INT      NKEYS=5;              /* see VISUAL2 manual             */
static INT      IKEYS[]={120,121,122,112,102,97,109};
static INT      FKEYS[]={1,1,1,1,2,1,1};             
static float    FLIMS[7][2];          /* data limits                    */
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
static DOUBLE  velx;	        /*					*/
static DOUBLE  vely;	        /*					*/
static DOUBLE  velz;	        /*					*/
static DOUBLE  pres;	        /*					*/
static DOUBLE  absv;	        /*					*/
static DOUBLE  vort;            /*                                      */
static DOUBLE  minpre,maxpre;   /*					*/
static DOUBLE  minvx,maxvx;     /*					*/
static DOUBLE  minvy,maxvy;     /*					*/
static DOUBLE  minvz,maxvz;     /*					*/
static DOUBLE  minabsv,maxabsv; /*					*/
static DOUBLE  FACX,FACY,FACZ;
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
\brief call of visual3 for fluid 

<pre>                                                         genk 12/03      
</pre>  
\param  numf   INT      (i)       actual number of fluid field
\return void                                                                       

------------------------------------------------------------------------*/
void vis3caf(INT numff, INT numaf, INT numsf) 
{

INT iscan, i,j,k;
CALC_ACTION    *action;             /* pointer to the cal_action enum   */
INTRA          *actintra;           /* pointer to active intra-communic.*/
PARTITION      *actpart;            /* pointer to active partition      */
CONTAINER       container;          /* contains variables defined in container.h */
FLUID_DYNAMIC *fdyn;                /* pointer to fluid dyn. inp.data   */

#ifdef DEBUG 
dstrc_enter("vis3caf");
#endif

actfield=&(field[numff]);
actfieldtyp=fluid;
fdyn = alldyn[genprob.numff].fdyn;

actpart= &(partition[numff]);
action = &(calc_action[numff]);
NUMA=numaf;

/*--------------------------------------------------- set some values */
numele = actfield->dis[0].numele;
numnp  = actfield->dis[0].numnp;
actele = &(actfield->dis[0].element[0]);
distyp = actele->distyp;

/*------------------------- check if all elements are the same distyp */
for (i=1;i<numele;i++) /* loop all elements */
{
   actele=&(actfield->dis[0].element[i]);
   if (actele->distyp!=distyp)
      dserror ("up to now, all elements have to be the same distyp!\n");
} /* end loop over all elements */

/*---------------- all element types are cut to 3- or 4-node elements */
switch (distyp)
{
case hex8:
   KCEL1=0;
   KCEL2=0;
   KCEL3=0;
   KCEL4=numele;
   KSURF=numele*8;
break;
default:
   dserror("distyp not implemented yet!\n");         
} /* end switch(distyp) */

/*-------------------- read solution data from flavia.res-file of proc 0 */
visual_readflaviares(actfield,&ncols,&time_a,&FIRSTSTEP,&LASTSTEP,&DSTEP);
INCRE=1;

input2:
printf("\n");
printf("     Please give the mode in which you want to run VISUAL2:\n");
printf("     0: steady data structure, grid and variables\n");
printf("     1: steady data structure and grid, unsteady variables\n");
printf("     2: steady data structure, unsteady grid and variables\n");
scanf("%d",&IOPT);
if (IOPT==0) icol++;

/*----------------------------------------------------- get window size */
WIN3D = true;
printf("\n");
printf("     Size of Windows:\n");
printf("\n");
printf("     0: 3D window bigger than 2D window\n");
printf("     1: 2D window bigger than 3D window\n");
scanf("%d",&HWIN3D);
if (HWIN3D==0) WIN3D = false;

/*----------------------------------------------------- get mirror flag */
input3:
printf("\n");
printf("     Mirroring:\n");
printf("\n");
printf("     0: no mirroring\n");
printf("     1: mirroring about the plane x=0.0\n");
printf("     2: mirroring about the plane y=0.0\n");
printf("     3: mirroring about the plane z=0.0\n");
scanf("%d",&MIRROR);
if (MIRROR<0 || MIRROR>3)
{
   printf("Input out of range! - Try again!\n");
   goto input3;
}

/* ----------------------------------------------------- get limits flag */
LIMITS=1;
printf("\n");
printf("     Limits:\n");
printf("\n");
printf("     0: no calculation of limits\n");
printf("     1: calculation of limits\n");
scanf("%d",&LIMITS);

/* --------------------------------------------------- get scaling facs */
printf("\n");
printf("     Scaling Factors:\n");
printf("\n");
printf("     factor in x-direction:\n");
scanf("%lf",&FACX);
printf("     factor in y-direction:\n");
scanf("%lf",&FACY);
printf("     factor in z-direction:\n");
scanf("%lf",&FACZ);

/*------------------------------------------------- background colour */
bgcolour=0;
printf("\n");
printf("   Background colour? (0: black   1: white)\n");
scanf("%d",&bgcolour);
#ifdef SUSE73
bgcolour+=10;
#endif

/*---------------------------------------------------- get data limits */
/*----------------------------------------------- get the data limits */
actnode=&(actfield->dis[0].node[0]);
velx    = actnode->sol.a.da[0][0];
vely    = actnode->sol.a.da[0][1];
velz    = actnode->sol.a.da[0][2];
absv    = sqrt(velx*velx+vely*vely+velz*velz);
pres    = actnode->sol.a.da[0][3];

minvx   = velx;
maxvx   = velx;
minvy   = vely;
maxvy   = vely; 
minvz   = velz;
maxvz   = velz; 
minpre  = pres;
maxpre  = pres;
minabsv = absv;
maxabsv = absv;

for (i=0;i<numnp;i++) /* loop nodes */
{
   actnode=&(actfield->dis[0].node[i]);
   for (j=0;j<ncols;j++) /* loop columns in sol-history */
   {
      velx    = actnode->sol.a.da[j][0];
      vely    = actnode->sol.a.da[j][1];
      velz    = actnode->sol.a.da[j][2];
      pres    = actnode->sol.a.da[j][3];
      absv    = sqrt(velx*velx+vely*vely+velz*velz);  
      minvx   = DMIN(minvx  ,velx);
      maxvx   = DMAX(maxvx  ,velx);
      minvy   = DMIN(minvy  ,vely);
      maxvy   = DMAX(maxvy  ,vely); 
      minvz   = DMIN(minvz  ,velz);
      maxvz   = DMAX(maxvz  ,velz); 
      minpre  = DMIN(minpre ,pres);
      maxpre  = DMAX(maxpre ,pres);
      minabsv = DMIN(minabsv,absv);
      maxabsv = DMAX(maxabsv,absv);
   } /* end of loop over columns */
} /* end of loop over nodes */

/*--------------------------------------- store max/min data in FLIMS */
if (minvx==maxvx) maxvx=minvx+ONE;
FLIMS[0][0]=minvx;
FLIMS[0][1]=maxvx;

if (minvy==maxvy) maxvy=minvy+ONE;
FLIMS[1][0]=minvy;
FLIMS[1][1]=maxvy;

if (minvz==maxvz) maxvz=minvz+ONE;
FLIMS[2][0]=minvz;
FLIMS[2][1]=maxvz;

if (maxpre==minpre) maxpre=minpre+ONE;
FLIMS[3][0]=minpre;
FLIMS[3][1]=maxpre;

if (minabsv==maxabsv) maxabsv=minabsv+ONE;
FLIMS[5][0]=maxabsv;
FLIMS[5][1]=ZERO;

FLIMS[4][0]=maxabsv;
FLIMS[4][1]=ZERO;


/*-------------------------- call Fortran routine which calls VISUAL2 */
v3call(&IOPT,&WIN3D,
       &NKEYS,IKEYS,FKEYS,FLIMS,
       &MIRROR,&numnp,&KCEL1,&KCEL2,&KCEL3,&KCEL4,
       &KSURF,&bgcolour);


#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of vis2caf */

/*!---------------------------------------------------------------------                                         
\brief return nodes to VISUAL3

<pre>                                                         genk 01/04      


</pre>  
\return void                                                                       
------------------------------------------------------------------------*/
void v3cell(INT **cel1, INT **cel2, INT **cel3, INT cel4[][8], 
            INT **nptet, INT *ptet) 
{
INT i,j;

#ifdef DEBUG 
dstrc_enter("v3cell");
#endif

for (i=0;i<numele;i++) /* loop all elements */
{
   actele=&(actfield->dis[0].element[i]);
   switch (actele->distyp)
   {
   case hex8:
      for (j=0;j<8;j++)
      {
         actnode = actele->node[j];
	 cel4[i][j]=actnode->Id_loc+1;
      }
   break;      
   default:
      dserror("distyp not implemented yet!\n");         
   }
} /* end loop over all elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of v3cell */

/*!---------------------------------------------------------------------                                         
\brief return nodes to VISUAL3

<pre>                                                         genk 01/04      


</pre>  
\return void                                                                       
------------------------------------------------------------------------*/
void v3surface(INT **nsurf, INT *surf, INT **scel, char tsurf[][100]) 
{

#ifdef DEBUG 
dstrc_enter("v3surface");
#endif

/* do nothing at the moment! */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of v3surface */


/*!---------------------------------------------------------------------                                         
\brief set VISUAL3's grid coordinates 

<pre>                                                         genk 01/04      


</pre>  
\return void                                                                       
------------------------------------------------------------------------*/
void v3grid(float XYZ[][3]) 
{
INT i;

#ifdef DEBUG 
dstrc_enter("v3grid");
#endif

switch(actfieldtyp)
{
case fluid:
#if 0
   if (NUMA>=0 && IOPT==0)
   {
      input4:
      v3_cursor(&true);  
      printf("Give the column number (min. 0; max. %d): ?\n",ncols-1);
      scanf("%d",&icol);      
      v3_cursor(&false);  
      if (icol<0 || icol>ncols-1)
      {
         printf("Column number out of range. New input!\n");
	 printf("\n");
	 goto input4;
      }
   }
#endif   
   if (IOPT==2)
   {
#ifdef D_FSI
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);
         actgnode=actnode->gnode;
         actanode=actgnode->mfcpnode[NUMA];
	 if (actanode==NULL)
	 {
            XYZ[i][0] = actnode->x[0];
            XYZ[i][1] = actnode->x[1];	 
            XYZ[i][2] = actnode->x[2];	 
         }
	 else
	 {
	    XYZ[i][0] = actanode->x[0] + actanode->sol.a.da[icol][0];
	    XYZ[i][1] = actanode->x[1] + actanode->sol.a.da[icol][1];	   
	    XYZ[i][2] = actanode->x[2] + actanode->sol.a.da[icol][2];	   
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
         XYZ[i][0] = actnode->x[0];
         XYZ[i][1] = actnode->x[1];
         XYZ[i][2] = actnode->x[2];
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
} /* end of V3GRID */


/*!---------------------------------------------------------------------                                         
\brief set VISUAL3's scalar function values

<pre>                                                         genk 01/04      

This routine is called by VISUAL3 during the visualisation.
It filters the scalar values out of the nodal sol-field at the actual
position icol.
icol is set every timestep in v3update


</pre>  
\param  *JKEY      INT    (i)  Index of key			 
\param	*S         float  (o)  Scalar function at nodes or cells 
\return void                                                                       
\sa v2update

------------------------------------------------------------------------*/
void v3scal(INT *JKEY, float *S) 
{

INT i;
float vx,vy,vz;

#ifdef DEBUG 
dstrc_enter("v3scal");
#endif

switch(actfieldtyp)
{
case fluid:
#if 0
   if (IOPT==0) /* read in the column number */
   {
      input3:
      v3_cursor(&true);
      printf("Give the column number (min. 0; max. %d): ?\n",ncols-1);
      scanf("%d",&icol);      
      v3_cursor(&false);      
      if (icol<0 || icol>ncols-1)
      {
         printf("Column number out of range. New input!\n");
	 printf("\n");
	 goto input3;
      }
   }
#endif   
   /*-------------------------------------------------- get scalar data */
   switch(*JKEY)
   {
   /*-------------------------------------------------------------------*/   
   case 1: /* Ux */
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);
	 S[i]=actnode->sol.a.da[icol][0];
      }      
   break;
   /*-------------------------------------------------------------------*/   
   case 2: /* Uy */
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);
	 S[i]=actnode->sol.a.da[icol][1];
      }        
   break;
   /*-------------------------------------------------------------------*/   
   case 3: /* Uz */
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);
	 S[i]=actnode->sol.a.da[icol][2];
      }        
   break;
   case 4: /* pressure */
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);
	 S[i]=actnode->sol.a.da[icol][3];
      }
   break;
   /*-------------------------------------------------------------------*/  
   case 5: 
      printf("VECTOR DATA IN SCALAR ROUTINE ???????\n");
   break;
   /*-------------------------------------------------------------------*/  
   case 6: /* absolute velocity */
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);
	 vx=actnode->sol.a.da[icol][0];
	 vy=actnode->sol.a.da[icol][1];
	 vz=actnode->sol.a.da[icol][2];
	 S[i] = sqrt(vx*vx+vy*vy+vz*vz);
      }       
   break;
   /*-------------------------------------------------------------------*/
   case 7:
      if (IMOVIE==0)
      {
         printf("\n");
         printf("   Starting movie creation at time 0.0\n");
	 printf("\n");
         IMOVIE=1;
      }      
   break;
   /*-------------------------------------------------------------------*/                 
   } /* end switch(*JKEY) */
break;
case structure:
   dserror("fieldtyp not implemented yet!\n");
case ale:
   dserror("fieldtyp not implemented yet!\n");
default:
   dserror("fieldtyp unknown!\n");
} /* end switch(actfieldtyp) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of V3SCAL */

/*!---------------------------------------------------------------------                                         
\brief set VISUAL3's vector function values

<pre>                                                         genk 01/04      

This routine is called by VISUAL3 during the visualisation.
It filters the vector values out of the nodal sol-field at the actual
position icol.
icol is set every timestep in v2update

</pre>  
\param  *JKEY      INT    (i)  Index of key			 
\param	*V         float  (o)  Vector function at nodes or cells 
\return void                                                                       
\sa v2update

------------------------------------------------------------------------*/
void v3vect(INT *JKEY, float V[][3]) 
{

INT i;

#ifdef DEBUG 
dstrc_enter("V3VECT");
#endif

switch(actfieldtyp)
{
case fluid:
   /*-------------------------------------------------- get vector data */
   if (*JKEY==5)  /* streamlines */
   {
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);
     	 V[i][0] = actnode->sol.a.da[icol][0];
     	 V[i][1] = actnode->sol.a.da[icol][1];
     	 V[i][2] = actnode->sol.a.da[icol][2];
      }
   }   
   else
     printf("SCALAR DATA IN VECTOR ROUTINE ???????\n");
   
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
} /* end of V3VECT */

/*!---------------------------------------------------------------------                                         
\brief update time and icol 

<pre>                                                         genk 01/04      

This routine is called by VISUAL3 during the visualisation.
Beside the actual time the actual position in the nodal sol-history 
icol is also updated

</pre>  
\param  *TIME      float    (o)  actual time		 
\return void                                                                       

------------------------------------------------------------------------*/
void v3update(float *TIME) 
{

#ifdef DEBUG 
dstrc_enter("V3UPDATE");
#endif

#if 0
/*-------------------------------------------------- check for stopping */
if (icol==STSTEP)
{
   printf("   Next Step to Stop (-1 to terminate)\n");
   scanf("%d",&STSTEP);        
}
#endif

#if 0
/*------------------------------------------- check for movie creation */
if (icol==0 && IMOVIE==2) 
{   
   IMOVIE=0;
   printf("\n");
   printf("Movie creation finished\n");
}   
if (icol==0 && IMOVIE==1) IMOVIE=2;
if (IMOVIE==2) v3movie();

#endif

switch(actfieldtyp)
{
case fluid: 
   if (icol==-1 || icol+INCRE > ncols-1) 
   {
      icol=0;
      ACTSTEP=FIRSTSTEP;
   }
   else
   { 
      icol+=INCRE;
      ACTSTEP+=DSTEP;
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
} /* end of V3UPDATE */


/*!---------------------------------------------------------------------                                         
\brief  additional label for  VISUAL3 plots 

<pre>                                                         genk 01/04      

This routine is called by VISUAL2 during the visualisation.

</pre>  
\param  STRING     char    (o)  plot on VISUAL3 window		 
\return void                                                                       

------------------------------------------------------------------------*/
void v3string(char STRING[81]) 
{

char *charpointer;
float t;

#ifdef DEBUG 
dstrc_enter("V3STRING");
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
} /* end of V3STRING */
#else
/*!---------------------------------------------------------------------                                         
\brief dummy routine 

<pre>                                                         genk 07/02      

since VISUAL3 is called by a fortran routine, this one is necessary, if
VIS2-routines are not compiled into the program

</pre>  		 
\return void                                                                       

------------------------------------------------------------------------*/
void v3_init()   
{
dserror("VISUAL3 PACKAGE not compiled into programm\n");
return;
}

#endif
