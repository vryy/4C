/*!----------------------------------------------------------------------
\file
\brief read data from flavia.res and make them available for VISUAL2

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#if defined VISUAL2_PACKAGE || defined VISUAL3_PACKAGE
/*----------------------------------------------------------- prototype */
void vis_frfind(char string[]);

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 *----------------------------------------------------------------------*/
extern struct _GENPROB       genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD         *field;
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

/*--------------------------------------------------- global variables */
FILE  *in;                        /* pointer to flavia.res file        */
char  line [500];                 /* actual line                       */
INT   eof;                        /* flag for end of flavia.res file   */
/*!---------------------------------------------------------------------                                         
\brief read flavia.res file

<pre>                                                         genk 07/02      

in this function the result from the .flavia.res file is read in order
to visualise it with VISUAL2/VISUAL3

</pre>  
\param  *actfield   FIELD      (i)       actual field
\param  *ntsteps    INT        (o)       number of steps read
\param  *time_a     ARRAY      (o)       time values
\param  *FIRSTSTEP  INT        (o)       first step to visualise
\param  *LASTSTEP   INT        (o)       last step to visualise
\param  *DSTEP      INT        (o)       increment
\return void                                                                       

------------------------------------------------------------------------*/
void visual_readflaviares(FIELD   *actfield, 
                          INT     *ntsteps,  
		          ARRAY   *time_a,
			  ARRAY   *step_a,
			  INT     *FIRSTSTEP,
			  INT     *LASTSTEP,
			  INT     *DSTEP	 
		         )  
{
INT             i,k;              /* simply some counters	        */
INT             screen;
INT             dummy;
INT             stepin=0;         /* number of result steps read        */
INT             num;              /* nodenumber                         */
INT             numnp;            /* number of nodes in actfield        */
INT             numdf;            /* number of fluid dofs               */
INT             mone=-1;
INT             ALL;
INT             var;
#ifdef D_FSI
INT             numnp_struct;
INT             numnp_ale;        
INT             numsf,numff;
#endif
static INT      firststep;
static INT      laststep;
static INT      ncols;
static INT      incre;
INT             counter=0;
ARRAY           val_a;            /* temporary resutls                  */
double        **val;
ARRAY           globloc;          /* global - local node Ids            */
char           *foundit = NULL;   
char           *end;
FLUID_DYNAMIC  *fdyn;             /* pointer to fluid dynamic input data*/
NODE           *actnode;          /* actual node                        */
#ifdef D_FSI
NODE           *actfnode;
GNODE          *actgnode;
FIELD          *structfield;
FIELD          *fluidfield;
#endif

#ifdef DEBUG 
dstrc_enter("visual_readflaviares");
#endif

/*------------------------------------------------------ initialisation */
eof=0;
in = allfiles.gidres;
rewind(in);
if (fgets(line,499,in)==NULL)
   dserror("An error occured reading a line from .flavia.res");

/*------------------------------------- determine global local node Ids */
amdef("globloc",&globloc,genprob.nnode,1,"IV");
aminit(&globloc,&mone);

switch (actfield->fieldtyp)
{
case fluid:
   numnp = actfield->dis[0].numnp;   
   for (i=0;i<numnp;i++)
   {
      actnode = &(actfield->dis[0].node[i]);  
      globloc.a.iv[actnode->Id] = actnode->Id_loc;
   }
   /* ---------------------------------------------------initialisation */
   fdyn = alldyn[genprob.numff].fdyn;
   numdf = fdyn->numdf;
   /*------------------------------------- count results in .flavia.res */
   while (eof==0)
   {
      vis_frfind("# RESULT velocity on FIELD fluid");
      stepin++;
      if (fgets(line,499,in)==NULL && eof==0)
        dserror("An error occured reading a line from fluid_start.data");
   }
   stepin--;
   input1:
   printf("\n");
   printf("     There are %d set of results in .flavia.res\n",stepin);
   printf("     Do you want to see all of them?\n");
   printf("      0 : no \n");
   printf("     [1]: yes\n");

   screen=getchar();
   switch(screen)
   {
   case 10: ALL=1; break;
   case 48: ALL=0; dummy=getchar(); break;
   case 49: ALL=1; dummy=getchar(); break;
   default: 
      printf("\nTry again!\n");
      goto input1;
   }
   if (ALL==0)
   {
      input2:
      printf("     At which step shall the visualisation start? (min = 0) \n");
      scanf("%d",&firststep);
      printf("     At which step shall the visualisation end? (max = %d)\n",
                stepin-1);   
      scanf("%d",&laststep);
      if (laststep<=firststep || firststep<0 || laststep >= stepin)
      {
         printf("   Input out of range --> try again!\n");
         goto input2;
      }
      ncols = laststep - firststep + 1;
      /*------------------------------------------------------ increment */
      input3:
      printf("\n");
      printf("     Increment step number by ...? (<%d)\n",ncols-1);
      scanf("%d",&incre);
      if (incre>ncols-1 || incre==0)
      {
         printf("   Increment out of range --> Try again!\n");
         goto input3;
      }
      dummy=getchar();
   }
   else
   {
      firststep  = 0;
      laststep = stepin-1;
      incre    = 1;
   }
   /*-------------------------------------------------- determine ncols */
   counter = firststep;
   stepin = 0;
   while (counter<=laststep)
   {
      stepin++;
      counter+=incre;
   }
   ncols = stepin;
   laststep = counter-incre;
    
   /*---------------------------------------- allocate temporary arrays */
   val = amdef("val",&val_a,numnp,numdf,"DA");
   /*------------------------------------------------ create time array */
   amdef("time",time_a,ncols,1,"DV");
   /*------------------------------------------------ create step array */
   amdef("step",step_a,ncols,1,"IV");   
   /*----------------------------------------- allocate nodal sol field */
   for (i=0;i<numnp;i++)
   {
      actnode=&(actfield->dis[0].node[i]);      
      amredef(&(actnode->sol),ncols,numdf,"DA");
   }
   /*------------------------------------ find firststep in .flavia.res */
   rewind(in);
   stepin=-1;
   eof=0;
   while (stepin<firststep)
   {
      vis_frfind("# RESULT velocity on FIELD fluid");
      stepin++;
      if (eof==1)
         dserror("Cannot read from fluid_start.data: stepin non-existent!\n");
      if (fgets(line,499,in)==NULL)
         dserror("An error occured reading a line from fluid_start.data");
   }
   stepin = firststep;
   counter = 0;
   eof = 0;
   printf("\n     Reading FLUID results ... \n");
   while (stepin<=laststep && eof==0)
   {
      printf("        STEP %d\n",stepin);
      /*---------------------------------------------- find & read time */
      vis_frfind("# TIME");
      foundit=strpbrk(line,"-.1234567890");
      time_a->a.dv[counter] = strtod(foundit,&end);
      /*---------------------------------------------- find & read step */
      vis_frfind("# STEP");
      foundit=strpbrk(line,"1234567890");
      sscanf(foundit," %d ",&var);
      step_a->a.iv[counter]= var;
      /*---------------------------------- find & read velocity results */
      vis_frfind("VALUES");
      for (i=0;i<numnp;i++)
      {
         if (fgets(line,499,in)==NULL)
            dserror("An error occured reading a line from .flavia.res");
         foundit = strstr(line," ");
         foundit=strpbrk(foundit,"-.1234567890");
         /*---------------------------------------- read global node Id */
	 num = strtod(foundit,&end)-1;   
	 /*------------------------------------ determine local node Id */
	 num = globloc.a.iv[num];
         if (num<0) dserror("node number not valid!\n");
         for (k=0;k<numdf-1;k++)
         {
            foundit = strstr(foundit," ");
            foundit=strpbrk(foundit,"-.1234567890");
            val[num][k] = strtod(foundit,&end); 
         }
      }
      /*------------------------------------------ plausibility checkes */
      if (fgets(line,499,in)==NULL)
         dserror("An error occured reading a line from .flavia.res\n");
      if (strstr(line,"END VALUES") == NULL)
         dserror("Number of Fluid nodes not correct in .flavia.res\n");
      /*--------------------------------- find  & read pressure results */
      vis_frfind("# RESULT pressure on FIELD fluid");
      vis_frfind("VALUES");
      for (i=0;i<numnp;i++)
      {
         if (fgets(line,499,in)==NULL)
            dserror("An error occured reading a line from .flavia.res");
         foundit = strstr(line," ");
         foundit=strpbrk(foundit,"-.1234567890");
         /*---------------------------------------- read global node Id */
	 num = strtod(foundit,&end)-1;   
	 /*------------------------------------ determine local node Id */
	 num = globloc.a.iv[num];
         if (num<0) dserror("node number not valid!\n");
	 foundit = strstr(foundit," ");
         foundit=strpbrk(foundit,"-.1234567890");
         val[num][numdf-1] = strtod(foundit,&end); 
      }
      /*------------------------------------------ plausibility checkes */
      if (fgets(line,499,in)==NULL)
         dserror("An error occured reading a line from .flavia.res\n");
      if (strstr(line,"END VALUES") == NULL)
         dserror("Number of Fluid nodes not correct in .flavia.res\n");
      /*------------------------------------- copy results to the nodes */
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);      
         num = actnode->Id_loc;
         for (k=0;k<numdf;k++)
            actnode->sol.a.da[counter][k]=val[num][k];
      }
      counter++;
      stepin+=incre;
      /*----------------------- check if there's another set of results */
      if (stepin<=laststep)
      for (i=0;i<incre;i++)
      {
         vis_frfind("# RESULT velocity on FIELD fluid");      
         if (fgets(line,499,in)==NULL)
            dserror("An error occured reading a line from .flavia.res\n");      
      }
   }
   if (counter!=ncols)
      dserror("an error occured reading .flavie.res!\n");
   *ntsteps=ncols;
   if (*ntsteps==0) dserror("no fluid results in .flavia.res\n");
   *FIRSTSTEP = firststep;
   *LASTSTEP  = laststep;
   *DSTEP     = incre;   
   /*------------------------------------------- delete temporary arrays */
   amdel(&globloc);
   amdel(&val_a);
break;
case structure:
   /*------------------------------------- count results in .flavia.res */
   while (eof==0)
   {
      vis_frfind("# RESULT DISPLACEMENTS on FIELD FSI");
      stepin++;
      if (fgets(line,499,in)==NULL && eof==0)
        dserror("An error occured reading a line from fluid_start.data");
   }
   stepin--;
   input4:
   printf("\n");
   printf("     There are %d set of results in .flavia.res\n",stepin);
   printf("     Do you want to see all of them?\n");
   printf("      0 : no \n");
   printf("     [1]: yes\n");

   screen=getchar();
   switch(screen)
   {
   case 10: ALL=1; break;
   case 48: ALL=0; dummy=getchar(); break;
   case 49: ALL=1; dummy=getchar(); break;
   default: 
      printf("\nTry again!\n");
      goto input4;
   }
   if (ALL==0)
   {
      input5:
      printf("     At which step shall the visualisation start? (min = 0) \n");
      scanf("%d",&firststep);
      printf("     At which step shall the visualisation end? (max = %d)\n",
                stepin-1);   
      scanf("%d",&laststep);
      if (laststep<=firststep || firststep<0 || laststep >= stepin)
      {
         printf("   Input out of range --> try again!\n");
         goto input5;
      }
      ncols = laststep - firststep + 1;
      /*------------------------------------------------------ increment */
      input6:
      printf("\n");
      printf("     Increment step number by ...? (<%d)\n",ncols-1);
      scanf("%d",&incre);
      if (incre>ncols-1 || incre==0)
      {
         printf("   Increment out of range --> Try again!\n");
         goto input6;
      }
      dummy=getchar();
   }
   else
   {
      firststep  = 0;
      laststep = stepin-1;
      incre    = 1;
   }
   /*-------------------------------------------------- determine ncols */
   counter = firststep;
   stepin = 0;
   while (counter<=laststep)
   {
      stepin++;
      counter+=incre;
   }
   ncols = stepin;
   laststep = counter-incre;

   numnp = actfield->dis[0].numnp;
   /*----------------------------------------------- create index array */
   for (i=0;i<numnp;i++)
   {
      actnode = &(actfield->dis[0].node[i]);  
      globloc.a.iv[actnode->Id] = actnode->Id_loc;
   }
   actnode=&(actfield->dis[0].node[0]);      
   numdf=actnode->numdf;
   if (numdf==6) numdf=3;
   /*---------------------------------------- allocate temporary arrays */
   val = amdef("val",&val_a,numnp,numdf,"DA");
   /*------------------------------------------------ create time array */
   amdef("time",time_a,ncols,1,"DV");
   /*------------------------------------------------ create step array */
   amdef("step",step_a,ncols,1,"IV");   
   /*----------------------------------------- allocate nodal sol field */
   for (i=0;i<numnp;i++)
   {
      actnode=&(actfield->dis[0].node[i]);      
      amredef(&(actnode->sol),ncols,actnode->numdf,"DA");
   }
   /*------------------------------------ find firststep in .flavia.res */
   rewind(in);
   stepin=-1;
   eof=0;
   while (stepin<firststep)
   {
      vis_frfind("# RESULT DISPLACEMENTS on FIELD FSI");
      stepin++;
      if (eof==1)
         dserror("Cannot read from fluid_start.data: stepin non-existent!\n");
      if (fgets(line,499,in)==NULL)
         dserror("An error occured reading a line from fluid_start.data");
   }
   stepin = firststep;
   counter = 0;
   eof = 0;
   printf("\n     Reading STRUCTURE results ... \n");
   while (stepin<=laststep && eof==0)
   {
      printf("        STEP %d\n",stepin);
      /*---------------------------------------------- find & read time */
      vis_frfind("# TIME");
      foundit=strpbrk(line,"-.1234567890");
      time_a->a.dv[counter] = strtod(foundit,&end);
      /*---------------------------------------------- find & read step */
      vis_frfind("# STEP");
      foundit=strpbrk(line,"1234567890");
      sscanf(foundit," %d ",&var);
      step_a->a.iv[counter]= var;
      /*------------------------------ find & read displacement results */
      vis_frfind("VALUES");
      /*---------------------------- read displacements of struct nodes */
      for (i=0;i<numnp;i++)
      {
         if (fgets(line,499,in)==NULL)
            dserror("An error occured reading a line from .flavia.res");
         foundit = strstr(line," ");
         foundit=strpbrk(foundit,"-.1234567890");
         /*---------------------------------------- read global node Id */
	 num = strtod(foundit,&end)-1;   
	 /*------------------------------------ determine local node Id */
	 num = globloc.a.iv[num];
         if (num<0) dserror("node number not valid!\n");
         for (k=0;k<numdf;k++)
         {
            foundit = strstr(foundit," ");
            foundit=strpbrk(foundit,"-.1234567890");
            val[num][k] = strtod(foundit,&end); 
         }
      }
      /*------------------------------------------ plausibility checkes */
      if (fgets(line,499,in)==NULL)
         dserror("An error occured reading a line from .flavia.res\n");
     /*------------------------------------- copy results to the nodes */
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);      
	 num = actnode->Id_loc;
         for (k=0;k<numdf;k++)
            actnode->sol.a.da[counter][k]=val[num][k];
      }
      counter++;
      stepin+=incre;
      /*----------------------- check if there's another set of results */
      if (stepin<=laststep)
      for (i=0;i<incre;i++)
      {
         vis_frfind("# RESULT DISPLACEMENTS on FIELD FSI");      
         if (fgets(line,499,in)==NULL)
            dserror("An error occured reading a line from .flavia.res\n");      
      }
   }
   if (counter!=ncols)
      dserror("an error occured reading .flavie.res!\n");
   *ntsteps=ncols;
   if (*ntsteps==0) dserror("no displ. results in .flavia.res\n");


#if 0
/* THIS DOES NOT WORK UP TO NOW */
   /*------------------------------------ find firststep in .flavia.res */
   rewind(in);
   stepin=-1;
   eof=0;
   while (stepin<firststep)
   {
      vis_frfind("# RESULT FSI LOADS on FIELD STRUCTURE");
      stepin++;
      if (eof==1)
         dserror("Cannot read from fluid_start.data: stepin non-existent!\n");
      if (fgets(line,499,in)==NULL)
         dserror("An error occured reading a line from fluid_start.data");
   }
   stepin = firststep;
   counter = 0;
   eof = 0;
   printf("\n     Reading FSI-loads ... \n");
   while (stepin<=laststep && eof==0)
   {
      printf("        STEP %d\n",stepin);
      /*--------------------------------- find & read fsi-loads results */
      vis_frfind("VALUES");
      /*-------------------------------- read fsi loads of struct nodes */
      for (i=0;i<numnp;i++)
      {
         if (fgets(line,499,in)==NULL)
            dserror("An error occured reading a line from .flavia.res");
         foundit = strstr(line," ");
         foundit=strpbrk(foundit,"-.1234567890");
         /*---------------------------------------- read global node Id */
	 num = strtod(foundit,&end)-1;   
	 /*------------------------------------ determine local node Id */
	 num = globloc.a.iv[num];
         if (num<0) dserror("node number not valid!\n");
         for (k=0;k<numdf;k++)
         {
            foundit = strstr(foundit," ");
            foundit=strpbrk(foundit,"-.1234567890");
            val[num][k] = strtod(foundit,&end); 
         }
      }
      /*------------------------------------------ plausibility checkes */
      if (fgets(line,499,in)==NULL)
         dserror("An error occured reading a line from .flavia.res\n");
     /*------------------------------------- copy results to the nodes */
      for (i=0;i<numnp;i++)
      {
         actnode=&(actfield->dis[0].node[i]);      
	 num = actnode->Id_loc;
         for (k=0;k<numdf;k++)
            actnode->sol.a.da[counter][k+numdf]=val[num][k];
      }
      counter++;
      stepin+=incre;
      /*----------------------- check if there's another set of results */
      if (stepin<=laststep)
      for (i=0;i<incre;i++)
      {
         vis_frfind("# RESULT FSI LOADS on FIELD STRUCTURE");      
         if (fgets(line,499,in)==NULL)
            dserror("An error occured reading a line from .flavia.res\n");      
      }
   }
   if (counter!=ncols)
      dserror("an error occured reading .flavie.res!\n");
#endif
   amdel(&val_a);   
   amdel(&globloc);
   *FIRSTSTEP = firststep;
   *LASTSTEP  = laststep;
   *DSTEP     = incre;   
break;
case ale:   
#ifdef D_FSI
   numff = genprob.numff;
   numsf = genprob.numsf;
   if (numsf!=-1) structfield=&(field[numsf]);
   fluidfield =&(field[numff]);
   numnp = fluidfield->dis[0].numnp;
   if (numsf != -1) numnp_struct = structfield->dis[0].numnp;
   numnp_ale = actfield->dis[0].numnp;
   /*----------------------------------------------- create index array */
   for (i=0;i<numnp;i++)
   {
      actnode = &(fluidfield->dis[0].node[i]);  
      globloc.a.iv[actnode->Id] = actnode->Id_loc;
   }
   /* ---------------------------------------------------initialisation */
   fdyn = alldyn[genprob.numff].fdyn;
   numdf = fdyn->numdf-1;
   /*---------------------------------------- allocate temporary arrays */
   val = amdef("val",&val_a,numnp,numdf,"DA");
   /*----------------------------------------- allocate nodal sol field */
   for (i=0;i<numnp_ale;i++)
   {
      actnode=&(actfield->dis[0].node[i]);      
      amredef(&(actnode->sol),ncols,numdf,"DA");
   }
   /*------------------------------------ find firststep in .flavia.res */
   rewind(in);
   stepin=-1;
   eof=0;
   while (stepin<firststep)
   {
      vis_frfind("# RESULT DISPLACEMENTS on FIELD FSI");
      stepin++;
      if (eof==1)
         dserror("Cannot read from fluid_start.data: stepin non-existent!\n");
      if (fgets(line,499,in)==NULL)
         dserror("An error occured reading a line from fluid_start.data");
   }
   stepin = firststep;
   counter = 0;
   eof = 0;
   printf("\n     Reading ALE results ... \n");
   while (stepin<=laststep && eof==0)
   {
      printf("        STEP %d\n",stepin);
      /*------------------------------ find & read displacement results */
      vis_frfind("VALUES");
      /*------------ skip results of structural field if there are some */
      if (numsf!=-1)
      for (i=0;i<numnp_struct;i++)
      if (fgets(line,499,in)==NULL)
         dserror("An error occured reading a line from .flavia.res");
      /*----------------------------- read displacements of fluid nodes */
      for (i=0;i<numnp;i++)
      {
         if (fgets(line,499,in)==NULL)
            dserror("An error occured reading a line from .flavia.res");
         foundit = strstr(line," ");
         foundit=strpbrk(foundit,"-.1234567890");
         /*---------------------------------------- read global node Id */
	 num = strtod(foundit,&end)-1;   
	 /*------------------------------------ determine local node Id */
	 num = globloc.a.iv[num];
         if (num<0) dserror("node number not valid!\n");
         for (k=0;k<numdf;k++)
         {
            foundit = strstr(foundit," ");
            foundit=strpbrk(foundit,"-.1234567890");
            val[num][k] = strtod(foundit,&end); 
         }
      }
      /*------------------------------------------ plausibility checkes */
      if (fgets(line,499,in)==NULL)
         dserror("An error occured reading a line from .flavia.res\n");
      if (strstr(line,"END VALUES") == NULL)
         dserror("Number of Fluid nodes not correct in .flavia.res\n");
      /*------------------------------------- copy results to the nodes */
      for (i=0;i<numnp_ale;i++)
      {
         actnode=&(actfield->dis[0].node[i]);      
         actgnode = actnode->gnode;
	 actfnode = actgnode->mfcpnode[numff];
	 if (actfnode==NULL) continue;
	 num = actfnode->Id_loc;
         for (k=0;k<numdf;k++)
            actnode->sol.a.da[counter][k]=val[num][k];
      }
      counter++;
      stepin+=incre;
      /*----------------------- check if there's another set of results */
      if (stepin<=laststep)
      for (i=0;i<incre;i++)
      {
         vis_frfind("# RESULT DISPLACEMENTS on FIELD FSI");      
         if (fgets(line,499,in)==NULL)
            dserror("An error occured reading a line from .flavia.res\n");      
      }
   }
   if (counter!=ncols)
      dserror("an error occured reading .flavie.res!\n");
   *ntsteps=ncols;
   if (*ntsteps==0) dserror("no displ. results in .flavia.res\n");
   amdel(&val_a);   
   amdel(&globloc);
#else
dserror("FSI functions not compiled in!\n");
#endif
break;
default:
   dserror("fieldtyp not valid!\n");
}


#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of visual_readflaviares */
/*!---------------------------------------------------------------------
\brief find a character string in flavia.res                                            

<pre>                                                        genk 07/03 
searches for a given character string in flavia.res
</pre>
\param string   char[]   (i)   string to search for in .flavia.res  
\return void                                               

------------------------------------------------------------------------*/
void vis_frfind(char string[])
{

#ifdef DEBUG 
dstrc_enter("vis_frfind");
#endif

while ( strstr(line,string) == NULL )
{
   if (fgets(line,499,in)==NULL)
   {
      if (feof(in)!=0)
      {
         eof = 1;
         goto end;
      }
      else
      dserror("An error occured reading a line from .flavia.res");
   }
}
end:

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of vis_frfind */

#endif
