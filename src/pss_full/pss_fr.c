#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | ranks and communicators                                              |
 | This structure struct _PAR par; is defined in main_ccarat.c
 *----------------------------------------------------------------------*/
 extern struct _PAR   par;                      

/*----------------------------------------------------------------------*
 | init file reading system                               m.gee 8/00    |
 *----------------------------------------------------------------------*/
void frinit()
{
int     i=0;
int     linecount=0;
char   *remarkpointer;
allfiles.numcol=300;
if (par.myrank==0)
{
/*--------------------make a copy of the input file without commentars */
   rewind(allfiles.in_input);
/*----------------------------------------------------------read a line */
   while ( feof(allfiles.in_input) == 0 )
   {
      if ( (fgets(allfiles.line,499,allfiles.in_input) == NULL) &&
           (feof(allfiles.in_input) == 0) )
      {
         dserror("An error occured reading a line from input file");
      }
/* ----------------------------------check whether its a commentar line */
      if (strncmp(allfiles.line,"//",2)!=0)
      {
         linecount++; 
      }
   }
}
/*--------------------------------------------broadcast number of lines */
#ifdef PARALLEL 
if (par.nprocs>1)
{
   MPI_Bcast(&linecount,1,MPI_INT,0,MPI_COMM_WORLD);
}
#endif 
allfiles.numrows=linecount;
/*--------------------------------------allocate space for copy of file */
allfiles.input_file_hook=(char*)CALLOC(linecount*(allfiles.numcol),sizeof(char));
if (allfiles.input_file_hook==NULL) dserror("Allocation of memory failed");

allfiles.input_file=(char**)CALLOC(linecount,sizeof(char*));
if (allfiles.input_file==NULL) dserror("Allocation of memory failed");

for (i=0; i<linecount; i++)
{
   allfiles.input_file[i] = &(allfiles.input_file_hook[i*(allfiles.numcol)]);
}
/*-------------------know read the input file into allfiles.input_file */
if (par.myrank==0)
{
   rewind(allfiles.in_input);
   for (i=0; i<linecount; i++)
   {
      read:
      if (fgets(allfiles.input_file[i],499,allfiles.in_input)==NULL)
         dserror("An error occured reading a line from input file");
/*----------------------------------check whether its a commentar line */
      if ((strncmp(allfiles.input_file[i],"//",2)==0))
      {
         goto read;
      }
/*----------------------check whether there is a commentar in the line */
      remarkpointer=strchr(allfiles.input_file[i],'/');
      if (remarkpointer!=NULL)
      {
         strcpy(remarkpointer,"\n\0");
      }
   }
}
/*--------------------------------broadcast the copy of the input file */
#ifdef PARALLEL 
if (par.nprocs>1)
{
      MPI_Bcast(allfiles.input_file_hook,
                allfiles.numrows*allfiles.numcol,
                MPI_CHAR,
                0,
                MPI_COMM_WORLD);
}
#endif
/*---------------------give a copy of the "cleaned" input file on .err */
if (par.myrank==0)
{
   fprintf(allfiles.out_err,"===========================================\n");
   fprintf(allfiles.out_err,"broadcasted copy of input file:            \n");
   fprintf(allfiles.out_err,"===========================================\n");
   for (i=0; i<linecount; i++)
   {
      fprintf(allfiles.out_err,"%s",allfiles.input_file[i]);
   }
   fprintf(allfiles.out_err,"===========================================\n");
   fprintf(allfiles.out_err,"end of broadcasted copy of input file      \n");
   fprintf(allfiles.out_err,"===========================================\n");
   fflush(allfiles.out_err);
/*---------------------------close input file, 'cause no longer needed */
   fclose(allfiles.in_input);
}
/*--------------------------------set fr-system to begin of input_file */
allfiles.actrow=0;
allfiles.actplace=&(allfiles.input_file[0][0]);


return;
} /* end of frinit */




/*----------------------------------------------------------------------*
 | rewind the copy of the input file                      m.gee 8/00    |
 *----------------------------------------------------------------------*/
void frrewind()
{
allfiles.actrow=0;
allfiles.actplace=&(allfiles.input_file[0][0]);
return;
} /* end of frrewind */





/*----------------------------------------------------------------------*
 | find a character string                                m.gee 8/00    |
 | char string[] (input) character string to search for copy of inputfil|
 |                       terminates programm if not found               |
 *----------------------------------------------------------------------*/
void frfind(char string[])
{
char message[100];
int  i=0;

#ifdef DEBUG 
dstrc_enter("frfind");
#endif

frrewind();
while ( strstr(allfiles.input_file[i],string) == NULL )
{
   if ( i==allfiles.numrows )
   {
      
      sprintf(message,"frfind:  String %s is not in input file",string);
      dserror(message);
   }
   i++;
}
allfiles.actrow=i;
allfiles.actplace=strstr(allfiles.input_file[i],string);

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of frfind */





/*----------------------------------------------------------------------*
 |                                                        m.gee 8/00    |
 | sets a pointer to the next line in thecopy of input file on all procs|
 *----------------------------------------------------------------------*/
void frread()
{

#ifdef DEBUG 
dstrc_enter("frread");
#endif

allfiles.actrow+=1;
if (allfiles.actrow>=allfiles.numrows) dserror("Can't read line, end of input_file reached");
allfiles.actplace=&(allfiles.input_file[allfiles.actrow][0]);

#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of frread */


/*----------------------------------------------------------------------*
 | reads n integers from input_file                       m.gee 4/01    |
 | char string[] (input) keyword to search for in actual line           |
 | int *var      (output) adress of field to hold values read           |
 | int num       (input)  number of values to read                      |
 | int *ierr     (output) =0 keyword not found / =1 values read         |
 *----------------------------------------------------------------------*/
void frint_n(char string[],int *var,int num, int *ierr)
{
char message[100];
int  i=0;
char *foundit = NULL;
#ifdef DEBUG 
dstrc_enter("frint");
#endif
/*----------------------------------------------------------------------*/
foundit = strstr(allfiles.input_file[allfiles.actrow],string);
if (foundit != NULL)
{
   foundit += strlen(string);
   for (i=0; i<num; i++)
   {
      foundit=strpbrk(foundit,"-1234567890");
      var[i] = strtol(foundit,&foundit,10);
   }
   *ierr=1;
}
else
{
*ierr=0;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of frint */



/*----------------------------------------------------------------------*
 | reads an integer from input_file                       m.gee 8/00    |
 | starting in allfiles.actrow                                          |
 | char string[] (input) keyword to search for in actual line           |
 | int *var      (output)adress of variable to read to                  |
 | int *ierr     (output) flag to indicate success                      |
 | ierr=0 keyword not found                                             |
 | ierr=1 integer read                                                  |
 *----------------------------------------------------------------------*/
void frint(char string[],int *var, int *ierr)
{
char message[100];
int  i=0;
char *foundit = NULL;
#ifdef DEBUG 
dstrc_enter("frint");
#endif
/*----------------------------------------------------------------------*/
foundit = strstr(allfiles.input_file[allfiles.actrow],string);
if (foundit != NULL)
{
   foundit += strlen(string);
   foundit=strpbrk(foundit,"-1234567890");
   i=sscanf(foundit," %d ",var);
   if (i!=1 || i==EOF)
   {
      sprintf(message,"frint:  an error occured reading %s",string);
      dserror(message);
   }
   else
   {
   *ierr=1;
   }
}
else
{
*ierr=0;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of frint */


/*----------------------------------------------------------------------*
 | reads a double from input file line allfiles.line      m.gee 8/00    |
 | see frint_n                                                          |
 | ierr=0 kenner not found                                              |
 | ierr=1 double  read                                                  |
 *----------------------------------------------------------------------*/
void frdouble_n(char string[],double *var,int num, int *ierr)
{
int  i;
char *foundit = NULL;
char *end;
#ifdef DEBUG 
dstrc_enter("frdouble");
#endif
/*----------------------------------------------------------------------*/
foundit = strstr(allfiles.input_file[allfiles.actrow],string);
if (foundit != NULL) 
{
   foundit += strlen(string);
   for (i=0; i<num; i++)
   {
      foundit=strpbrk(foundit,"-.1234567890");
      var[i] = strtod(foundit,&foundit);
   }
   *ierr=1;
}
else
{
*ierr=0;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of frdouble */



/*----------------------------------------------------------------------*
 | reads a double from input file line allfiles.line      m.gee 8/00    |
 | see frint                                                            |
 | ierr=0 kenner not found                                              |
 | ierr=1 double  read                                                  |
 *----------------------------------------------------------------------*/
void frdouble(char string[],double *var, int *ierr)
{
char *foundit = NULL;
char *end;
#ifdef DEBUG 
dstrc_enter("frdouble");
#endif
/*----------------------------------------------------------------------*/
foundit = strstr(allfiles.input_file[allfiles.actrow],string);
if (foundit != NULL) 
{
   foundit += strlen(string);
   foundit=strpbrk(foundit,"-.1234567890");
   *var = strtod(foundit,&end);
   *ierr=1;
}
else
{
*ierr=0;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of frdouble */





/*----------------------------------------------------------------------*
 | reads a charstring from input file                                   |
 | user must assure, that the given charpointer space is long enough    |   
 | to hold the string                                                   |
 | ierr=0 keyword not found on line                                     |
 | ierr=1 char string read                                 m.gee 8/00   |
 *----------------------------------------------------------------------*/
void frchar(char string[],char *var, int *ierr)
{
int  i=0;
char *foundit = NULL;
char *end;

#ifdef DEBUG 
dstrc_enter("frchar");
#endif

foundit = strstr(allfiles.input_file[allfiles.actrow],string);
if (foundit != NULL) 
{
   foundit += strlen(string);
   sscanf(foundit," %s ",var);
   *ierr=1;
}
else
{
*ierr=0;
}

#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of frchar */



/*----------------------------------------------------------------------*
 | checks for a keyword in actual line                                  |
 | ierr=0 not found                                                     |
 | ierr=1 found                                            m.gee 8/00   |
 | char string[] (input) character string to check actual line for      |
 *----------------------------------------------------------------------*/
void frchk(char string[], int *ierr)
{
int  i=0;
char *foundit = NULL;
char *end;

foundit = strstr(allfiles.input_file[allfiles.actrow],string);
if (foundit != NULL) *ierr=1;
else *ierr=0;
return;
} /* end of frchk */



/*----------------------------------------------------------------------*
 | close and delete input file copy                        m.gee 4/01   |
 *----------------------------------------------------------------------*/
void frend()
{
/*----------------------------------------------------------------------*/
FREE(allfiles.input_file);
FREE(allfiles.input_file_hook);
allfiles.input_file=NULL;
allfiles.input_file_hook=NULL;
allfiles.actplace = NULL;
allfiles.numrows=0;
allfiles.numcol=0;
allfiles.actrow=0;
/*----------------------------------------------------------------------*/
return;
} /* end of frend */
