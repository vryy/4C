#include "../headers/standardtypes.h"
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*!----------------------------------------------------------------------
\brief the tracing variable

<pre>                                                         m.gee 8/00
defined in pss_ds.c, declared in tracing.h                                                  
</pre>
*----------------------------------------------------------------------*/
#ifdef DEBUG
extern struct _TRACE         trace;
#endif

/*----------------------------------------------------------------------*
 | input of problem title         m.gee 8/00    |
 *----------------------------------------------------------------------*/
void inpctrhed()
{
int  linecount=0;
int  i=0;
#ifdef DEBUG 
dstrc_enter("inpctrhed");
#endif

frrewind();
frfind("-TITLE");
/*--------------------set input_file to next line */
allfiles.actrow+=1;
i=allfiles.actrow;
while ( strncmp(allfiles.input_file[i],"------",6) != 0)
{
   if (linecount<5)
   {
      strcpy(allfiles.title[linecount],allfiles.input_file[i]);
      linecount++;
      i++;
   }
   else dserror("Only 5 lines of title allowed");
}
frrewind();

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpctrhed */
/*----------------------------------------------------------------------*
 | input of tracing variable                              m.gee 8/00    |
 *----------------------------------------------------------------------*/
void inptrace()
{
int  linecount=0;
int  i=0;
char buffer[40];
int  ierr=0;
#ifdef DEBUG 

frrewind();
frfind("TRACE");
frchar("TRACE",buffer,&ierr);

if (
    strncmp(buffer,"secure",6)==NULL ||
    strncmp(buffer,"Secure",6)==NULL ||
    strncmp(buffer,"SECURE",6)==NULL 
   )
    trace.trace_on=1;
else
    trace.trace_on=0;
#endif
return;
} /* end of inptrace */
