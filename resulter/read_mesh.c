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
/*----------------------------------------------------------------------*
 | read_mesh                                              m.gee 5/03    |
 *----------------------------------------------------------------------*/
int read_mesh(FILE *inmesh)
{
int   i;
char  buffer[500],*bptr;
/*------------------------------------------------- get the line "MESH" */
rewind(inmesh);
fgets(buffer,499,inmesh);
while( strncmp(buffer,"MESH",4)!=0 )
{
   fgets(buffer,499,inmesh);
}
/*------------------------------------------------------- get dimension */
bptr       = strstr(buffer,"DIMENSION");
bptr       = bptr + 10;
dimension  = strtol(bptr,&bptr,10);
/*------------------------------------------------------ get nodeperele */
bptr       = strstr(buffer,"NNODE");
bptr       = bptr + 6;
nodeperele = strtol(bptr,&bptr,10);
/*-------------------------------------- go where the coordinates start */
while( strncmp(buffer,"COORDINATES",11)!=0 )
{
   fgets(buffer,499,inmesh);
}
/*----------------------------------------------- count number of nodes */
numnp = 0;
while( strncmp(buffer,"END COORDINATES",15)!=0 )
{
   fgets(buffer,499,inmesh);
   numnp++;
}
/*--------------------------------------------- go where elements start */
while( strncmp(buffer,"ELEMENTS",8)!=0 )
{
   fgets(buffer,499,inmesh);
}
/*------------------------------------------------------ count elements */
numele = 0;
while( strncmp(buffer,"END ELEMENTS",12)!=0 )
{
   fgets(buffer,499,inmesh);
   numele++;
}
/*----------------------------------------- allocate nodes and elements */
/*
node = (NODE*)calloc(numnp,sizeof(NODE));
element = (ELEMENT*)calloc(numele,sizeof(ELEMENT));
*/
return(0);
} 
