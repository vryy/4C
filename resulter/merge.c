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
 |                                                        m.gee 5/03    |
 *----------------------------------------------------------------------*/
int merge_files(int nfiles, FILE *inputs[], FILE *outputs[], FILE *info, char *argv[])
{
int   i,j,fnum;
int   first,last;
int   head=0;
char  buffer[500],*bptr,*bptr2;
char  choice[100];
FILE *input,*output;
output = outputs[0];
rewind(output);

for (fnum=0; fnum<nfiles; fnum++)
{
   input = inputs[fnum];
   rewind(input);
   printf("======================FILE IS %s\n",argv[fnum+3]);
   for (j=0; j<nresulttype[fnum]; j++)
   {
      printf("RESULT first %10d last %10d total %10d %s\n",
             firstresultnum[fnum][j],lastresultnum[fnum][j],numresult[fnum][j],resulttypes[fnum][j]);
   
   }
   printf("\n");
   fflush(stdout);
   /*---------------------------------------- print head of result file */
   if (!head)
   {
      printf("PRINT HEAD OF FILE? (y/n)\n");
      fgets(choice,9,stdin);
      if (strncmp("y",choice,1)==0) 
      {
         fgets(buffer,499,input);
         while( strncmp(buffer,"# RESULT ",9)!=0 )
         {
            fputs(buffer,output);
            fgets(buffer,499,input);
         }
         head = 1;
      }
   
      fflush(output);
   }
   /*------------------------------------------------- loop resulttypes */
   for (j=0; j<nresulttype[fnum]; j++)
   {
      /* printf info about a certain type of result */
      printf("RESULTTYPE NOW IS first %10d last %10d total %10d %s\n",
             firstresultnum[fnum][j],lastresultnum[fnum][j],numresult[fnum][j],resulttypes[fnum][j]);
      /* skip this result or not */
      printf("SKIP? (y/n)\n");
      fgets(choice,9,stdin);
      if (strncmp("y",choice,1)==0) continue;
      /* check last step from previous */
      if (fnum>0)
         printf("LAST  STEP IN PREVIOUS FILE IS  %d\n",lastresultnum[fnum-1][j]);
         printf("FIRST STEP IN THIS     FILE IS  %d\n",firstresultnum[fnum][j]);
      /* check first step from next */
         printf("LAST  STEP IN THIS     FILE IS  %d\n",lastresultnum[fnum][j]);
      if (fnum<nfiles-1)
         printf("FIRST STEP IN NEXT     FILE IS  %d\n",firstresultnum[fnum+1][j]);
      /* get first step to print */
      printf("GIVE FIRST STEP YOU WANT:\n");
      fgets(choice,9,stdin);
      first = strtol(choice,&bptr,10);
      printf("GIVE LAST  STEP YOU WANT:\n");
      fgets(choice,9,stdin);
      last  = strtol(choice,&bptr,10);
      print_result(input,output,info,first,last,fnum,j);

   }
}
exit:
fflush(output);
return(0);
} 



/*----------------------------------------------------------------------*
 |                                                        m.gee 5/03    |
 *----------------------------------------------------------------------*/
int print_result(FILE *input, FILE *output, FILE *info, int first, int last, 
                 int fnum, int resnum)
{
int   i,number;
char *bptr;
char buffer[400];
char buffer1[400];
char buffer2[400];
char buffer3[400];
char buffer4[400];
char choice[10];
rewind(input);

/*----------------- search for next result and keep the last four lines */
goon:
bptr = fgets(buffer,499,input);
if (!bptr) goto exit;
/* find a result */
while(strncmp(buffer,"RESULT ",7)!=0 )
{
   strcpy(buffer4,buffer3);
   strcpy(buffer3,buffer2);
   strcpy(buffer2,buffer1);
   strcpy(buffer1,buffer);
   bptr = fgets(buffer,499,input);
   if (!bptr) goto exit;
}
/* check correct type of result */
bptr = strstr(buffer,resulttypes[fnum][resnum]);
if (!bptr) goto goon;
/* check correct number of result */
bptr   = strpbrk(buffer,"1234567890");
number = strtol(bptr,&bptr,10);
if (number < first || number > last) goto goon;
/* print result to output */
if (strncmp(buffer4,"#",1)==0)
fprintf(output,"%s",buffer4);
fprintf(output,"%s",buffer3);
fprintf(output,"%s",buffer2);
fprintf(output,"%s",buffer1);
fprintf(output,"%s",buffer);
while(strncmp(buffer,"END ",4)!=0 )
{
   bptr = fgets(buffer,499,input);
   fprintf(output,"%s",buffer);
}
fflush(output);
goto goon;



exit:
fflush(output);
return(0);
} 







