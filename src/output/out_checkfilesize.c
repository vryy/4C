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
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par; 


/*!---------------------------------------------------------------------
\brief  check file size of .flavia.res

<pre>                                                         genk 01/03

It may happen, that the file .flavia.res get very big, even bigger than
by the system allowed. Then the calculation continues without writing
any results. This can be very annoying.
This function checks the size of the .flavia.res file by several system
commands, compares it with MAXFILESIZE of definitions.h. If the file
is bigger than allowed, .flavia.res is closed, renamed and a new file
is opened.
		     
</pre>

\param    opt         INT          (i)          evaluation option
\return void                                                                             
\warning this function is relatively slow

------------------------------------------------------------------------*/
void out_checkfilesize(INT opt)
{
#ifdef D_GENK
char  filename[100]=("                                                    ");
char  kenner[100]=("                                                      ");
char *charpointer;
char  command[120]=("du -sk                                               ");
char  remove[]=("rm resfilesize");
char  rename[120]=("mv                                                    ");
char  line [500]; 
char           *end;
long int   length;
INT   i,j;
long int   size;
FILE *resfilesize;

#ifdef DEBUG 
dstrc_enter("out_checkfilesize");
#endif

/*---------------------------------------------- get output file kenner */
strcpy(filename,allfiles.outputfile_kenner);

switch (opt)
{
case 1: /* flavia.res file */
   /*------------- determine name and length of actual .flavia.res file */
   if (allfiles.num_flaviaresfiles==1)
   {
      length = strlen(&filename[0]);
      charpointer=filename+length;
      strncpy(charpointer,"0.flavia.res",12);
      length+=12;
   }
   else
   {
      length = strlen(&filename[0]);
      charpointer=filename+length;
      strncpy(charpointer,"_file_",6);
      charpointer+=6;
      sprintf(charpointer,"%d",allfiles.num_flaviaresfiles);
      charpointer++;
      strncpy(charpointer,"_0.flavia.res",13);
      length+=20;
   }

   /*--------------- determine the size of the actual .flavia.res file 
     via a system command and pipe the result in the file "resfilesize" */
   for (i=0;i<length;i++)
   command[7+i] = filename[i];
   length+=7;
   charpointer = command+length;
   strncpy(charpointer," >resfilesize",13);
   system(&command[0]);

   /*------------ read the size of the actual .flavia.res from the file */
   if ( (resfilesize=fopen("resfilesize","r"))==NULL)
      dserror("Error opening file resfilesize\n");
   if (fgets(line,499,resfilesize)==NULL)
      dserror("An error occured reading a line from resfilesize");
   charpointer=line;
   size = strtod(charpointer,&end);

   /*------------------------------------ remove the file "resfilesize" */
   fclose(resfilesize);
   system(&remove[0]);

   if (size>MAXFILESIZE)
   {
      /*---------------------------------- close the actual .flavia.res */
      fclose(allfiles.gidres);
      /*--------------------------------- rename the actual .flavia.res */
      if (allfiles.num_flaviaresfiles==1)
      {
         strcpy(kenner,allfiles.outputfile_kenner);
         length = strlen(&kenner[0]);
         charpointer=kenner+length;
         strncpy(charpointer,"0.flavia.res",12);
         for (i=0;i<length+12;i++)
         rename[3+i] = kenner[i];
         strncpy(charpointer,"_file_1_0.flavia.res",20);
         for (j=0;j<length+20;j++)
         rename[4+i+j] = kenner[j];
         system(&rename[0]);
      }
      /*------------------------------ get name of the next .flavia.res */ 
      strcpy(filename,allfiles.outputfile_kenner);
      length = strlen(&filename[0]);
      charpointer=filename+length;
      strncpy(charpointer,"_file_",6);
      charpointer+=6;
      allfiles.num_flaviaresfiles++;
      sprintf(charpointer,"%d",allfiles.num_flaviaresfiles);
      charpointer++;
      strncpy(charpointer,"_0.flavia.res",13);
      if ( (allfiles.gidres=fopen(filename,"w"))==NULL)
         dserror("opening of .flavia.res file failes\n");

   }
   break;
   default:
      dserror("parameter opt out of range!\n");
}   
/*----------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_exit();
#endif
#endif
return;
} /* end of out_checkfilesize */
