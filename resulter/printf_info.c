#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "ctype.h"
#include "definitions.h"
#include "am.h"
extern int      dimension;
extern int      numnp;
extern int      numele;
extern int      nodeperele;
extern NODE    *node;
extern ELEMENT *element;
/*
resulttypes[0..number of resultfiles][0..number of DIFFERENT resultypes][resultname]
*/
extern char     resulttypes[150][50][100];
/*
lengthresult[0..number of resultfiles][0..number of DIFFERENT resultypes] how long is the string
*/
extern int      lengthresult[150][50];
/*
nresulttype[0..number of resultfiles] number of DIFFERENT resulttypes
*/
extern int      nresulttype[150];
/*
numresult[0..number of resultfiles][0..number of DIFFERENT resultypes] how often is the result
*/
extern int      numresult[150][50];
/*
firstresultnum[0..number of resultfiles][0..number of DIFFERENT resultypes] lowest number of the result
lastresultnum [0..number of resultfiles][0..number of DIFFERENT resultypes] highest number of the result
*/
extern int      firstresultnum[150][50];
extern int      lastresultnum[150][50];
/*----------------------------------------------------------------------*
 | printf info of a result file                           m.gee 5/03    |
 *----------------------------------------------------------------------*/
int print_info(FILE *info, FILE *result, char name[], int fnum)
{
int   i;
char  buffer[500],*bptr,*bptr2;
char  sign='"';
int   length;
int   number;
nresulttype[fnum] = 0;
rewind(result);
for (i=0; i<50; i++) 
{
   numresult[fnum][i]    = 0;
   lengthresult[fnum][i] = 0;
}

while(fgets(buffer,499,result)!=NULL)
{
   while( strncmp(buffer,"RESULT ",7)!=0 )
   {
      bptr = fgets(buffer,499,result);
      if (!bptr) goto endoffile;
   }
   /* set bptr to start of resulttype */
   bptr = buffer + 8;
   /* set bptr2 to end of resulttype */
   bptr2 = strpbrk(bptr,"\"");
   /* get length */
   length = (int)(bptr2 - bptr);
   /* get the number of the result */
   bptr2 = strpbrk(bptr2,"1234567890");
   number = strtol(bptr2,&bptr2,10);
   /* copy it to resulttypes */
   /* loop resultypes and check, whether it was there before */
   for (i=0; i<nresulttype[fnum]; i++)
   {
      if ( strncmp(resulttypes[fnum][i],bptr,length)==0 )
      {
         numresult[fnum][i]++;
         lastresultnum[fnum][i] = number;
         break;
      }
   }
   /* it was not there before */
   if (i==nresulttype[fnum])
   {
      strncpy(resulttypes[fnum][i],bptr,length);
      lengthresult[fnum][i] = length;
      firstresultnum[fnum][i] = number;
      numresult[fnum][i]++;
      nresulttype[fnum]++;
      if (nresulttype[fnum]==50) dserror("Overflow of nresulttype");
   }
}
endoffile:
/*---------------------------------- printf info to info file */
fprintf(info,"FILE %s\n",name);
for (i=0; i<nresulttype[fnum]; i++)
{
fprintf(info,"RESULT first %10d last %10d total %10d %s\n",
        firstresultnum[fnum][i],lastresultnum[fnum][i],numresult[fnum][i],resulttypes[fnum][i]);
}
fprintf(info,"\n");
fflush(info);







return(0);
} 
