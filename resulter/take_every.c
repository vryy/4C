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
int take_every(int fnum, FILE *input, FILE *output, FILE *info, char name[])
{
int   i,j;
char  buffer[500],*bptr,*bptr2;
char  choice[10];
int   off,mod;
int   counter;
rewind(input);
rewind(output);
/*--------------------------------------------------- printf info again */
printf("\n======================FILE IS %s\n",name);
for (i=0; i<nresulttype[fnum]; i++)
{
printf("RESULT first %10d last %10d total %10d %s\n",
        firstresultnum[fnum][i],lastresultnum[fnum][i],numresult[fnum][i],resulttypes[fnum][i]);
}
printf("\n");
fflush(stdout);
/*-------------------------------------------------- print head of file */
fgets(buffer,499,input);
while( strncmp(buffer,"# RESULT ",9)!=0 )
{
   fputs(buffer,output);
   fgets(buffer,499,input);
}
/*------------------------------------------------ loop the resulttypes */
for (i=0; i<nresulttype[fnum]; i++)
{
   /* printf info about a certain type of result */
   printf("RESULTTYPE NOW IS first %10d last %10d total %10d %s\n",
          firstresultnum[fnum][i],lastresultnum[fnum][i],numresult[fnum][i],resulttypes[fnum][i]);
   /* skip this result or not */
   printf("SKIP? (y/n)\n");
   fgets(choice,9,stdin);
   if (strncmp("y",choice,1)==0) continue;
   /* ask for the offset */
   printf("GIVE THE OFFSET YOU WANT:\n");
   fgets(choice,9,stdin);
   off = strtol(choice,&bptr2,10);
   /* loop these results */
   rewind(input);
   /* ask and print these result to output */
   counter=0;
   for (j=0; j<numresult[fnum][i]; j++)
   {
      mod = counter % off;
      printresult_ask(input,output,info,resulttypes[fnum][i],mod);
      counter++;
   }
}

fflush(output);
return(0);
} 
/*----------------------------------------------------------------------*
 |                                                        m.gee 5/03    |
 *----------------------------------------------------------------------*/
int printresult_ask(FILE *input, FILE *output, FILE *info, char resulttype[], int mod)
{
int   i;
char *bptr;
char buffer[400];
char buffer1[400];
char buffer2[400];
char buffer3[400];
char buffer4[400];
char choice[10];
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
bptr = strstr(buffer,resulttype);
if (!bptr) goto goon;
/* check offset */
if (mod != 0) goto exit;
/* print result head to screen */
if (strncmp(buffer4,"#",1)==0)
printf("%s",buffer4);
printf("%s",buffer3);
printf("%s",buffer2);
printf("%s",buffer1);
printf("%s",buffer);
/*---------------------------------------------- print result to output */
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

exit:
fflush(output);
return(0);
} 




/*----------------------------------------------------------------------*
 |                                                        m.gee 5/03    |
 *----------------------------------------------------------------------*/
int take_every_s(int fnum, FILE *input, FILE *output, FILE *info, char name[])
{
int   i,j;
char  buffer[500],*bptr,*bptr2;
char  choice[10];
int   counter=0;
rewind(input);
rewind(output);
/*--------------------------------------------------- printf info again */
printf("\n======================FILE IS %s\n",name);
for (i=0; i<nresulttype[fnum]; i++)
{
printf("RESULT first %10d last %10d total %10d %s\n",
        firstresultnum[fnum][i],lastresultnum[fnum][i],numresult[fnum][i],resulttypes[fnum][i]);
}
printf("\n");
fflush(stdout);
/*-------------------------------------------------- print head of file */
fgets(buffer,499,input);
while( strncmp(buffer,"# RESULT ",9)!=0 )
{
   fputs(buffer,output);
   fgets(buffer,499,input);
}
/*------------------------------------------------ loop the resulttypes */
for (i=0; i<nresulttype[fnum]; i++)
{
   /* printf info about a certain type of result */
   printf("RESULTTYPE NOW IS first %10d last %10d total %10d %s\n",
          firstresultnum[fnum][i],lastresultnum[fnum][i],numresult[fnum][i],resulttypes[fnum][i]);
   /* skip this result or not */
   printf("SKIP? (y/n)\n");
   fgets(choice,9,stdin);
   if (strncmp("y",choice,1)==0) continue;
   /* loop these results */
   rewind(input);
   /* ask and print these result to output */
   for (j=0; j<numresult[fnum][i]; j++)
   {
      printresult_ask_s(input,output,info,resulttypes[fnum][i]);
   }
}

fflush(output);
return(0);
} 
/*----------------------------------------------------------------------*
 |                                                        m.gee 5/03    |
 *----------------------------------------------------------------------*/
int printresult_ask_s(FILE *input, FILE *output, FILE *info, char resulttype[])
{
int   i;
char *bptr;
char buffer[400];
char buffer1[400];
char buffer2[400];
char buffer3[400];
char buffer4[400];
char choice[10];
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
bptr = strstr(buffer,resulttype);
if (!bptr) goto goon;
/* print result head to screen */
if (strncmp(buffer4,"#",1)==0)
printf("%s",buffer4);
printf("%s",buffer3);
printf("%s",buffer2);
printf("%s",buffer1);
printf("%s",buffer);
printf("TAKE? (y/n)\n");
fgets(choice,9,stdin);
if (strncmp("n",choice,1)==0) goto goon;
/*---------------------------------------------- print result to output */
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

exit:
fflush(output);
return(0);
} 
