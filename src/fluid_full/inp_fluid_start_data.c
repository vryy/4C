/*!----------------------------------------------------------------------
\file
\brief reading fluid data from fluid_start_data

------------------------------------------------------------------------*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
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
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;
/*----------------------------------------------------------------------*
 | variables needed for parallel comp.                    m.gee 8/00    |
 *----------------------------------------------------------------------*/
extern struct _PAR   par;
/*!---------------------------------------------------------------------                                         
\brief input from file 'fluid_start-data'

<pre>                                                         genk 07/02

in this routine the inital fluid data are read form 'fluid_start_data'
and stored in the array 'start'
			     
</pre>   
\return void  
\warning ONLY SEQUENTIEL VERSION !!!!!
------------------------------------------------------------------------*/
void inp_fluid_start_data()
{
int irc=1;                        /* flag for file opening              */
int k;                            /* counters                           */
int ierr1=0;                      /* flag for reading form file         */
int ierr2=0;                      /* flag for reading form file         */
int ierr3=0;                      /* flag for reading form file         */
int numnp, numdf;                 /* number of points, dof per node     */
int num;                          /* actual number of node during input */
int counter=0;                    
int actnf;                        /* actual field number                */
double time;                      
char line [500];
char *foundit  = NULL;
char *foundit1 = NULL;
char *foundit2 = NULL;
char *foundit3 = NULL;		  /* pointers for input                 */
char *end;
double **start;
FILE *file;
FLUID_DYNAMIC *fdyn;              /* pointer to fluid dynamic input data*/
FIELD *actfield;                  /* pointer to actual field            */

#ifdef DEBUG 
dstrc_enter("inp_fluid_start_data");
#endif
/* -------------------------------------------find number of fluid field */
for (actnf=0;actnf<genprob.numfld;actnf++)
{
   actfield=&(field[actnf]);
   if (actfield->fieldtyp==fluid)
   goto fluid;
}
goto end;
fluid:    
/* --------------------------------------------------------- set pointer */
fdyn = alldyn[actnf].fdyn;
/*------------------------------------------------ check initialiazation */
if(fdyn->init==1)                     /* read from file fluid_start.data */
{
/* ----------------------------------------- open file fluid_start.data */
   if (par.myrank==0)
   {
      file = fopen("fluid_start.data","r");
      if (file == NULL) irc=0;
   }
#ifdef PARALLEL
   if (par.myrank==0)
   dserror (" parallel version of inp_fluid_start_data not implemented yet!");
   printf("WARNING: parallel version of inp_fluid_start_data not checekd yet!\n");
/*   MPI_Bcas(&irc,1,MPI_INT,0,MPI_COMM_WORLD); */
#endif 	         
   if (irc==0)
   {
      if (par.myrank==0)
      printf("opening of file fluid_start.data failed\n");
#ifdef PARALLEL
      MPI_Finalize();
#else               
      exit(1);
#endif	 
   }      
   if (par.myrank==0)
/*------------------------------------------- read file fluid_start.data */
   {
      if (fgets(line,499,file)==NULL)
      dserror("An error occured reading a line from file fluid_start.data");
      while(strncmp(line,"------",6)!=0)
      {
	 foundit1 = strstr(line,"NUMNP");
	 if (foundit1 != NULL) ierr1 = 1;
	 foundit2 = strstr(line,"TIME");
	 if (foundit2 != NULL) ierr2 = 1;
	 foundit3 = strstr(line,"NUMDF");
	 if (foundit3 != NULL) ierr3 = 1;
	 if (ierr1==1)
	 {
	    foundit1 += strlen("NUMNP");
	    foundit1=strpbrk(foundit1,"-1234567890");
	    sscanf(foundit1," %d ",&numnp);
	    ierr1=0;
	    counter++;
	 }
	 if (ierr2==1)  	 
	 {
	    foundit2 += strlen("TIME");
	    foundit2=strpbrk(foundit2,"-1234567890");
	    time = strtod(foundit2,&end);
	    ierr2=0;
	    counter++;
	 }
	 if (ierr3==1)  	
	 {
	    foundit3 += strlen("NUMDF");
	    foundit3=strpbrk(foundit3,"-1234567890");
	    sscanf(foundit3," %d ",&numdf);
	    ierr3=0;
	    counter++;  		       
	 }
      if (fgets(line,499,file)==NULL)
	 dserror("An error occured reading a line from file fluid_start.data");	    
      }
      if (counter!=3)
      dserror("An error occured reading NUMDF, NUMNP, TIME from file fluid_start.data"); 	 	 
   }
#ifdef PARALLEL
/*   MPI_Bcas(&numnp,1,MPI_INT   ,0,MPI_COMM_WORLD);
   MPI_Bcas(&time ,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcas(&numdf,1,MPI_INT   ,0,MPI_COMM_WORLD); */
#endif
/*------------------------------------------------- store / check values */      
   if (numnp != field->dis[0].numnp)
      dserror("NUMNP not the same in input file and file fluid_start.data");
   if (numdf != fdyn->numdf)
      dserror("NUMDF not the same in input file and file fluid_start.data");
   fdyn->time=time;
/*------------------------------------- allocate vector for storing data */
    start=amdef("start",&(fdyn->start),numnp,numdf,"DA");
    amzero(&(fdyn->start));
/*------------------------------------------ read initial data from file */ 
    if (par.myrank==0)
    {
       if (fgets(line,499,file)==NULL)
          dserror("An error occured reading a line from file fluid_start.data");
       counter=0;
       while(strncmp(line,"------",6)!=0)
       {
	  foundit=strstr(line," ");
	  foundit += strlen(" ");
	  foundit =strpbrk(foundit,"-1234567890");
	  num = strtol(foundit,&foundit,10);
	  num--;
	  counter++;
	  for(k=0;k<numdf;k++)
	  {
	     foundit=strpbrk(foundit,"-.1234567890");		    
	     fdyn->start.a.da[num][k] = strtod(foundit,&foundit); 
	  }
       if (fgets(line,499,file)==NULL)
          dserror("An error occured reading a line from file fluid_start.data");	    
       }
    if(counter != numnp)
       dserror("NUMNP wrong in file fluid_start.data");
/*----------------------------------- close file fluid_start.data */
    fclose(file);
    printf("initial field read from    fluid_start.data\n"); 
   }            
}    

/*----------------------------------------------- broadcast array start */
#ifdef PARALLEL
/* MPI_Bcas(start,numnp*numdf,MPI_DOUBLE,0,MPI_COMM_WORLD); */
#endif
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} 
/* end of inp_fluid_start_data */ 


#endif
