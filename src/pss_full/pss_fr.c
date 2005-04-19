/*!---------------------------------------------------------------------
\file
\brief file reading routines

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/

#include <ctype.h>
#include "../headers/standardtypes.h"
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;


/*!
\addtogroup FRSYSTEM
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>

*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


/* Very simple stack used to include files. */

#ifndef INCLUDE_STACK_DEPTH
#define INCLUDE_STACK_DEPTH 20
#endif

typedef struct
{
  FILE* file;
  CHAR* filename;
} include_stack_entry;


/*!---------------------------------------------------------------------
\brief init file reading system

<pre>                                                        m.gee 8/00
The routine open the inputfile allfiles.input on proc 0 and
reads it twice:
-Only counting lines, which are no comment (comment: //)
-Copying file to allfiles.input_file without comments
-Broadcasting allfiles.input_file to all procs
-Sets the fr-systems pointers used by the fr routines on all procs
</pre>

Additionally it's possible to #include outer files. The files are
merged in memory to one set of lines. Thus no other function is
concerned with this feature. Here, however, it adds a little
complexity.

------------------------------------------------------------------------*/
void frinit()
{
INT     i=0;
INT     linecount=0;
char   *remarkpointer;

static include_stack_entry include_stack[INCLUDE_STACK_DEPTH];
static INT include_stack_head = 0;
FILE* f;

allfiles.numcol=MAXNUMCOL;
if (par.myrank==0)
{
  /*-------------------- make a copy of the input file without comments */
  
  include_stack[0].file = allfiles.in_input;
  include_stack[0].filename = allfiles.inputfile_name;

  f = allfiles.in_input;
  rewind(f);
    
  while ( feof(f) == 0 )
  {
    INT p;

    for (;;)
    {
      /* see if we reached the end of a included file */
      if ((feof(f) != 0) && (include_stack_head>0))
      {
        fclose(include_stack[include_stack_head].file);
        CCAFREE(include_stack[include_stack_head].filename);
        include_stack_head--;
        f = include_stack[include_stack_head].file;
      }
      
      /* Read a line. Up to 499 characters. */
      if (fgets(allfiles.line,499,f) == NULL)
      {
        if ((feof(f) != 0) && (include_stack_head>0))
        {
          /* end of included file. handled above. */
          continue;
        }
        else if ((feof(f) != 0) && (include_stack_head==0))
        {
          /* end of main file. */
          break;
        }
        else
        {
          dserror("An error occured reading file '%s'",
                  include_stack[include_stack_head].filename);
        }
      }
      
      /* skip whitespaces */
      for (p=0; isspace(allfiles.line[p]); ++p)
      {}

      /* we accept this line only if it's not fully commented */
      if ((allfiles.line[p] != '/') ||
          (allfiles.line[p+1] != '/'))
      {

        /* look for includes */
        for (p=0; isspace(allfiles.line[p]); ++p)
        {}
        if (allfiles.line[p]=='#')
        {
          for (++p; isspace(allfiles.line[p]); ++p)
          {}
          if (strncmp("include", &(allfiles.line[p]), 7)==0)
          {
            CHAR* nameend;
            INT namelength;
        
            /* So we've got an include. Prepare the stack. Find the
             * filename. */

            include_stack_head++;
            if (include_stack_head >= INCLUDE_STACK_DEPTH)
              dserror("include stack overflow");
        
            for (p+=7; isspace(allfiles.line[p]); ++p)
            {}
            if (allfiles.line[p]!='"')
              dserror("'\"' expected in #include line");
            p += 1;
            nameend = strstr(&(allfiles.line[p]), "\"");
            if (nameend==NULL)
              dserror("terminating '\"' missing in #include line");

            namelength = nameend - &(allfiles.line[p]);
            include_stack[include_stack_head].filename =
              (CHAR*)CCACALLOC(namelength+1, 1);
            strncpy(include_stack[include_stack_head].filename,
                    &(allfiles.line[p]),
                    namelength);

            include_stack[include_stack_head].file =
              fopen(include_stack[include_stack_head].filename, "r");
            f = include_stack[include_stack_head].file;
            if (f==NULL)
              dserror("failed to open file '%s'",
                      include_stack[include_stack_head].filename);
            continue;
          }
        }
      }

      break;
    }

    /* count any non-comment lines */
    if (feof(f) == 0)
      linecount++;
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
allfiles.input_file_hook=(char*)CCACALLOC(linecount*(allfiles.numcol),sizeof(char));
if (allfiles.input_file_hook==NULL) dserror("Allocation of memory failed");

allfiles.input_file=(char**)CCACALLOC(linecount,sizeof(char*));
if (allfiles.input_file==NULL) dserror("Allocation of memory failed");

for (i=0; i<linecount; i++)
{
   allfiles.input_file[i] = &(allfiles.input_file_hook[i*(allfiles.numcol)]);
}
/*------------------- now read the input file into allfiles.input_file */
if (par.myrank==0)
{
  include_stack[0].file = allfiles.in_input;
  include_stack[0].filename = allfiles.inputfile_name;
  
  f = allfiles.in_input;
  rewind(f);

  for (i=0; i<linecount; i++)
  {
    INT p;

    for (;;)
    {
      /* see if we reached the end of a included file */
      if ((feof(f) != 0) && (include_stack_head>0))
      {
        fclose(include_stack[include_stack_head].file);
        CCAFREE(include_stack[include_stack_head].filename);
        include_stack_head--;
        f = include_stack[include_stack_head].file;
      }
      
      if (fgets(allfiles.input_file[i],499,f)==NULL)
      {
        if ((feof(f) != 0) && (include_stack_head>0))
        {
          /* end of included file. handled above. */
          continue;
        }
        else if ((feof(f) != 0) && (include_stack_head==0))
        {
          /* end of main file. */
          break;
        }
        else
        {
          dserror("An error occured reading file '%s'",
                  include_stack[include_stack_head].filename);
        }
      }

      /* skip whitespaces */
      for (p=0; isspace(allfiles.input_file[i][p]); ++p)
      {}

      /* we accept this line only if it's not fully commented */
      if ((allfiles.input_file[i][p] != '/') ||
          (allfiles.input_file[i][p+1] != '/'))
      {
        
        /* look for includes */
        for (p=0; isspace(allfiles.input_file[i][p]); ++p)
        {}
        if (allfiles.input_file[i][p]=='#')
        {
          for (++p; isspace(allfiles.input_file[i][p]); ++p)
          {}
          if (strncmp("include", &(allfiles.input_file[i][p]), 7)==0)
          {
            CHAR* nameend;
            INT namelength;
        
            /* So we've got an include. Prepare the stack. Find the
             * filename. */

            include_stack_head++;
            if (include_stack_head >= INCLUDE_STACK_DEPTH)
              dserror("include stack overflow");
        
            for (p+=7; isspace(allfiles.input_file[i][p]); ++p)
            {}
            if (allfiles.input_file[i][p]!='"')
              dserror("'\"' expected in #include line");
            p += 1;
            nameend = strstr(&(allfiles.input_file[i][p]), "\"");
            if (nameend==NULL)
              dserror("terminating '\"' missing in #include line");

            namelength = nameend - &(allfiles.input_file[i][p]);
            include_stack[include_stack_head].filename =
              (CHAR*)CCACALLOC(namelength+1, 1);
            strncpy(include_stack[include_stack_head].filename,
                    &(allfiles.input_file[i][p]),
                    namelength);

            include_stack[include_stack_head].file =
              fopen(include_stack[include_stack_head].filename, "r");
            f = include_stack[include_stack_head].file;
            if (f==NULL)
              dserror("failed to open file '%s'",
                      include_stack[include_stack_head].filename);
            continue;
          }
        }

        /* So this is not a comment line and on include here
         * either. We'll have to read this line the ordinary way. */
        break;
      }
    }
      
    /* cut comments that start inside the line */
    remarkpointer=strchr(allfiles.input_file[i],'/');
    if ((remarkpointer!=NULL) && (remarkpointer[1] == '/'))
    {
      remarkpointer[0] = '\n';
      remarkpointer[1] = '\0';
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
#ifdef DEBUG
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
}
#endif /* ifdef DEBUG */

/*---------------------------close input file, 'cause no longer needed */
if (par.myrank==0)
   fclose(allfiles.in_input);

/*--------------------------------set fr-system to begin of input_file */
allfiles.actrow=0;
allfiles.actplace=&(allfiles.input_file[0][0]);


return;
} /* end of frinit */




/*!---------------------------------------------------------------------
\brief rewind the copy of the input file

<pre>                                                        m.gee 8/00
sets pointer allfiles.actplace to allfiles.input_file[0][0]
</pre>
\return void

------------------------------------------------------------------------*/
void frrewind()
{
allfiles.actrow=0;
allfiles.actplace=&(allfiles.input_file[0][0]);
return;
} /* end of frrewind */





/*!---------------------------------------------------------------------
\brief find a character string

<pre>                                                        m.gee 8/00
searches for a given character string in allfiles.input_file
and sets allfiles.actplace to it
If string is not found it returns 0 otherwise 1.
</pre>
\param string   char[]   (i)   string to search for in input
\return INT

------------------------------------------------------------------------*/
INT frfind(char string[])
{
char message[100];
INT  i=0;

#ifdef DEBUG
dstrc_enter("frfind");
#endif

frrewind();
while ( strstr(allfiles.input_file[i],string) == NULL )
{
   if ( i==allfiles.numrows-1)
   {

      sprintf(message,"frfind:  String %s is not in input file",string);
      /*dserror(message);*/
      /*printf("%s\n",message);*/
      fprintf(allfiles.out_err,"%s\n",message);
      frrewind();
#ifdef DEBUG
      dstrc_exit();
#endif
      return(0);
   }
   i++;
}
allfiles.actrow=i;
allfiles.actplace=strstr(allfiles.input_file[i],string);

#ifdef DEBUG
dstrc_exit();
#endif
return(1);
} /* end of frfind */





/*!---------------------------------------------------------------------
\brief sets a pointer to the next line in the copy of input

<pre>                                                        m.gee 8/00
sets allfiles.actplace to the next line in allfiles.input_file
</pre>
\return void

------------------------------------------------------------------------*/
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


/*!---------------------------------------------------------------------
\brief reads n integers from input_file

<pre>                                                        m.gee 4/01
reads n integers from input_file
</pre>
\param string   char[]   (i)   string to search for in actual line
\param var      INT*     (o)   adress of field to hold values read
\param num      INT      (i)   number of values to read
\param ierr     INT*     (o)   =0  keyword not found / =1 values read
\return void
\sa frdouble_n() , frint()

------------------------------------------------------------------------*/
void frint_n(char string[],INT *var,INT num, INT *ierr)
{
INT  i=0;
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



/*!---------------------------------------------------------------------
\brief reads an integer from input_file

<pre>                                                        m.gee 8/00
reads an integer from input_file starting in allfiles.actrow
Searches for the keyword and reads the first thing consisting of -1234567890 behind it
</pre>
\param string   char[]   (i)   string to search for in actual line
\param var      INT*     (o)   adress of variable to hold value read
\param ierr     INT*     (o)   =0  keyword not found / =1 values read
\return void
\sa frdouble_n() , frint_n() , frdouble() , frchar()

------------------------------------------------------------------------*/
void frint(char string[],INT *var, INT *ierr)
{
char message[100];
INT  i=0;
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


/*!---------------------------------------------------------------------
\brief reads n doubles from input_file

<pre>                                                        m.gee 8/00
reads n doubles from input_file starting in allfiles.actrow
Searches for the keyword and reads the first thing consisting of -.1234567890 behind it.
All values to be read must then be continous
</pre>
\param string   char[]   (i)   string to search for in actual line
\param var      DOUBLE*  (o)   adress of field to hold values read
\param ierr     INT*     (o)   =0  keyword not found / =1 values read
\return void
\sa frdouble() , frint_n() , frint() , frchar()

------------------------------------------------------------------------*/
void frdouble_n(char string[],DOUBLE *var,INT num, INT *ierr)
{
INT  i;
char *foundit = NULL;
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



/*!---------------------------------------------------------------------
\brief reads a DOUBLE from input_file

<pre>                                                        m.gee 8/00
reads a DOUBLE from input_file starting in allfiles.actrow
Searches for the keyword and reads the first thing consisting of -.1234567890 behind it.
</pre>
\param string   char[]   (i)   string to search for in actual line
\param var      DOUBLE*  (o)   adress of variable to hold value read
\param ierr     INT*     (o)   =0  keyword not found / =1 value read
\return void
\sa frdouble_n() , frint_n() , frint() , frchar()

------------------------------------------------------------------------*/
void frdouble(char string[],DOUBLE *var, INT *ierr)
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





/*!---------------------------------------------------------------------
\brief reads a charstring from input file

<pre>                                                        m.gee 8/00
searches for keyword string and reads continous charstring behind it
user must assure, that the given charpointer space is long enough
to hold the string
</pre>
\param string   char[]   (i)   string to search for in actual line
\param var      char*    (o)   adress of variable to hold string read
\param ierr     INT*     (o)   =0  keyword not found / =1 value read
\return void
\sa frdouble_n() , frint_n() , frint()

------------------------------------------------------------------------*/
void frchar(char string[],char *var, INT *ierr)
{
char *foundit = NULL;

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



/*!---------------------------------------------------------------------
\brief checks for a keyword in actual line

<pre>                                                        m.gee 8/00
checks for a keyword in actual line of input_file
</pre>
\param string   char[]   (i)   string to search for in actual line
\param ierr     INT*     (o)   =0  keyword not found / =1 value read
\return void
\sa frdouble_n() , frint_n() , frint() , frchar()

------------------------------------------------------------------------*/
void frchk(char string[], INT *ierr)
{
char *foundit = NULL;

foundit = strstr(allfiles.input_file[allfiles.actrow],string);
if (foundit != NULL) *ierr=1;
else *ierr=0;
return;
} /* end of frchk */



/*!---------------------------------------------------------------------
\brief close and delete input file copy

<pre>                                                        m.gee 4/01
close and delete input file copy. All memory allocated to the
fr-system is freed
</pre>
\warning Nothing can be read after call to this routine
\return void
\sa frdouble_n() , frint_n() , frint() , frchar() , frchk()

------------------------------------------------------------------------*/
void frend()
{
/*----------------------------------------------------------------------*/
CCAFREE(allfiles.input_file);
CCAFREE(allfiles.input_file_hook);
allfiles.input_file=NULL;
allfiles.input_file_hook=NULL;
allfiles.actplace = NULL;
allfiles.numrows=0;
allfiles.numcol=0;
allfiles.actrow=0;
/*----------------------------------------------------------------------*/
return;
} /* end of frend */

/*! @} (documentation module close)*/
