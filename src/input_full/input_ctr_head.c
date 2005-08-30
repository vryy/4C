/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
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
extern struct _CCA_TRACE         trace;
#endif

/*----------------------------------------------------------------------*
 | input of problem title         m.gee 8/00    |
 *----------------------------------------------------------------------*/
void inpctrhed()
{
INT  linecount=0;
INT  i=0;
#ifdef DEBUG
dstrc_enter("inpctrhed");
#endif

frrewind();
if (frfind("-TITLE")==0) dserror("frfind: TITLE not in input file");
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


