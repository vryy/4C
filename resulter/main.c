#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "ctype.h"
#include "definitions.h"
#include "am.h"
int      dimension; /* dimension of the problem */
int      numnp;     /* number of nodes */
int      numele;    /* number of elements */
int      nodeperele;/* nodes per element */
NODE    *node;      /* vector of nodes */
ELEMENT *element;   /* vecotr of elements */
/*
resulttypes[0..number of resultfiles][0..number of DIFFERENT resultypes][resultname]
*/
char     resulttypes[150][50][100];
/*
lengthresult[0..number of resultfiles][0..number of DIFFERENT resultypes] how long is the string
*/
int      lengthresult[150][50];
/*
nresulttype[0..number of resultfiles] number of DIFFERENT resulttypes
*/
int      nresulttype[150];
/*
numresult[0..number of resultfiles][0..number of DIFFERENT resultypes] how often is the result
*/
int      numresult[150][50];
/*
firstresultnum[0..number of resultfiles][0..number of DIFFERENT resultypes] lowest number of the result
lastresultnum [0..number of resultfiles][0..number of DIFFERENT resultypes] highest number of the result
*/
int      firstresultnum[150][50];
int      lastresultnum[150][50];
/*----------------------------------------------------------------------*
 | main                                                   m.gee 5/03    |
 | call:
 | resulter info .msh   1.res  2.res  3.res ....                        
 *----------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
int   i;
int   nresult;
char *bptr;
char  choice[10];
char  outputname[50];
FILE *inmesh;
FILE *input[200];
FILE *output[200];
FILE *info;
/*------------------------------------------------------ open info file */
bptr   = strstr(argv[1],".msh");
if (bptr) dserror("You gave a .msh as first argument instead of an info file");
bptr   = strstr(argv[1],".res");
if (bptr) dserror("You gave a .res as first argument instead of an info file");
info = fopen(argv[1],"w");
if (!info) dserror("Cannot open info file");
/*------------------------------------------------------ open mesh file */
bptr   = strstr(argv[2],".msh");
if (!bptr) dserror("First file must be a .msh mesh file");
inmesh = fopen(argv[2],"r");
if (!inmesh) dserror("Canot open mesh file");
/*--------------------------------------- result files start in argv[2] */
nresult = argc-3;
/*--------------------------------------------------- open result files */
for (i=0; i<nresult; i++)
{
   bptr = strstr(argv[3+i],".res");
   if (!bptr) dserror("files must be a .res result files");
   input[i] = fopen(argv[3+i],"r");
   if (!input[i]) dserror("Canot open result files");
}
/*--------------------------------------------------- open output files */
for (i=0; i<nresult; i++)
{
   sprintf(outputname,"resulter_out%d.flavia.res",i);
   output[i] = fopen(outputname,"w");
   if (!output[i]) dserror("Cannot open output file");
}
/*----------------------------------------------------------- read mesh */
/* not needed yet */
/*read_mesh(inmesh);*/
/*------------------------------------ scan result files and print info */
for (i=0; i<nresult; i++)
{
   print_info(info,input[i],argv[3+i],i);
}
/*==================================================== printf main menu */
mainmenu:
printf("=========================MAIN MENU============================\n");
printf("merge files                                          (m)\n");
printf("select steps from file by offset                     (o)\n");
printf("select single steps from file                        (s)\n");
printf("quit                                                 (q)\n");
printf("==============================================================\n");
/*choice = fgetc(stdin);*/
fgets(choice,9,stdin);
/*----------------------------------------------------------------- quit */
if (strncmp("q",choice,1)==0) goto exit;
/*---------------------------------------------------------- merge files */
if (strncmp("m",choice,1)==0)
{
   merge_files(nresult,input,output,info,argv);
}
/*------------------------ go through files and choose results by offset */
if (strncmp("o",choice,1)==0)
{
  for (i=0; i<nresult; i++)
  take_every(i,input[i],output[i],info,argv[3+i]);
}
/*------------------------ go through files and choose single results    */
if (strncmp("s",choice,1)==0)
{
  for (i=0; i<nresult; i++)
  take_every_s(i,input[i],output[i],info,argv[3+i]);
}
goto mainmenu;
/*=======================================================================*/










exit:
return;
} /* end of main */
