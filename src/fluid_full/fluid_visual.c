/*!----------------------------------------------------------------------
\file
\brief writing and reading data for visual2/3 to/from pss-file

------------------------------------------------------------------------*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../headers/solution.h"
#include "../headers/visual.h"
#include "fluid_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;                     
 
/*!---------------------------------------------------------------------                                         
\brief output to pss-file

<pre>                                                         genk 07/02

in this routine the solution (node->sol) is written to the pss-file of
proc 0 used for visualisation with VISUAL2 / VISUAL3
all other data required by VISUAL2 / VISUAL3 are read from the input    
file
			     
</pre>   
\param *actfield	FIELD         (i)  actual field
\param  ntsteps         int           (i)  total num. of time steps
\param *time_a          ARRAY         (i)  time array 
\return void  

------------------------------------------------------------------------*/
void fluid_writevispss(FIELD  *actfield,  
                       int     ntsteps,  
		       ARRAY  *time_a	 
		      )
{
int                  i;             /* simply a counter                 */
int                  ierr;          /* error flag                       */
int                  numnp;         /* number of fluid nodes            */
long int            *node_handles;  /* handles for writing pss-file     */
long int             longdummy;     /* dummy                            */
FILE                *out;           /* pointer to pss-file              */
VISUAL_FLUID         vis;           /* visualisation data               */
NODE                *actnode;       /* actual node                      */

#ifdef DEBUG 
dstrc_enter("fluid_writevispss");
#endif

/*----------------------------------------------------------------------*/
out = allfiles.out_pss;
vis.step = ntsteps;

/*--------------------------------------------- write the ARRAY *time_a */
amredef(time_a,ntsteps,1,"DV");
pss_write_array(time_a,&(vis.time),out,&ierr);

/*----------------------------------------------------------------------*/
/* 
   now write the solution data, that are stored in the nodes: 
   
   actnode->sol
   
   so we need to store an array of numnp handles
*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
numnp = actfield->dis[0].numnp;
vis.node_handles = (long int*)malloc(numnp*sizeof(long int));
if (!vis.node_handles) dserror("Allocation of memory failed\n");
node_handles = vis.node_handles;

/*--------------- now we loop the nodes and each node writes his ARRAYs */
for (i=0; i<numnp; i++)
{
   actnode=&(actfield->dis[0].node[i]);
   pss_write_array(&(actnode->sol),&(node_handles[i]),out,&ierr);
   if (ierr != 1) dserror("Error writing visual data\n");
}

/*---------------------------------------------------------------------* 
   all nodal data is written, so we now write the node_handles and store
   the handle to this in vis.handle_of_node_handles
 *---------------------------------------------------------------------*/
pss_write("nod_hand",numnp,1,sizeof(long int),node_handles,&(vis.handle_of_node_handles),out,&ierr);
if (ierr != 1) dserror("Error writing visual data\n");

/*----------------- delete the vis.node_handles but keep the dimensions */
vis.node_fdim = numnp;
vis.node_sdim = 1;
free(vis.node_handles);

/*----------------------------------------------------------------------*/
/*
   the only thing to do is now write the VISUAL_FLUID itself with 
   its unique name "fluidvis"
   NOTE:
   names are limited to 9 characters.
*/   
pss_write("fluidvis",1,1,sizeof(VISUAL_FLUID),&vis,&longdummy,out,&ierr);
if (ierr != 1) dserror("Error writing visual data\n");

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_writevispss */ 


/*!---------------------------------------------------------------------                                         
\brief input from pss-file

<pre>                                                         genk 07/02

in this routine the solution (node->sol) is read from the pss-file of
proc 0 used for visualisation with VISUAL2 / VISUAL3
all other data required by VISUAL2 / VISUAL3 are read from the input    
file
			     
</pre>   
\param *actfield	FIELD         (i)  actual field
\param  ntsteps         int           (i)  total num. of time steps
\param *time_a          ARRAY         (i)  time array 
\return void  

------------------------------------------------------------------------*/
void fluid_readvispss(FIELD   *actfield, 
                      int     *ntsteps,  
		      ARRAY   *time_a	 
		     )
{
int                  i,j;              /* simply some counters          */
int                  ierr;             /* error flag                    */
int                  byte;             /* byte size                     */
int                  dims[3];          /* dims for reading an array     */
int                  numnp;            /* number of fluid nodes         */
long int             vishandle;        /* handle of "fluidvis"          */
long int            *node_handles;     /* handles for reading pss-file  */
FILE                *in;               /* pointer to pss-input file     */
VISUAL_FLUID         vis;              /* visualisation data            */
NODE                *actnode;          /* actual node                   */

#ifdef DEBUG 
dstrc_enter("fluid_readvispss");
#endif

/*--------------------------------------- set pointer to pss-input file */
in = allfiles.in_pss;

/*------------------------------ check if data exist that shall be read */
pss_chck("fluidvis",&vishandle,in,&ierr);
if (ierr != 1) dserror("Cannot visualise fluid, data don't exist in vis.pss-file\n");

/*------------------------ the structure fluidvis exists, so we read it */
pss_read_name_handle("fluidvis",&i,&i,&byte,&vis,&vishandle,in,&ierr);
if (ierr != 1) dserror("Restart structure exists, but cannot read it\n");
*ntsteps=vis.step;

/*------------------------------------------------- allocate time-array */
amdef("time",time_a,vis.step,1,"DV");

/*---------------------------------------------- read the ARRAY *time_a */
pss_read_array_name_handle(time_a->name,time_a,&(vis.time),in,&ierr);
if (ierr != 1) dserror("Cannot read visual data: time_a\n");

/*----------------------------------------------------------------------*/
/*
   now read the data that is stored in the nodes: 

   actnode->sol

   so we need to read an array of numnp handles first
*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
numnp = actfield->dis[0].numnp;
if (numnp != vis.node_fdim)
    dserror("Mismatch in number of nodes on reading visual data from file .pss\n");

/*----------------------------------------- define the array of handles */
vis.node_handles = (long int*)malloc(numnp*sizeof(long int));
if (!vis.node_handles) dserror("Allocation of memory failed\n");
node_handles = vis.node_handles;

/*------------------------------------------- read the array of handles */
pss_read_name_handle("nod_hand",&(vis.node_fdim),&(vis.node_sdim),&i,
                     node_handles,&vis.handle_of_node_handles,in,&ierr);
if (ierr != 1) dserror("Cannot read visual data\n");

/*---------------- now we loop the nodes and each node reads his ARRAYs */
for (i=0; i<numnp; i++)
{
   actnode = &(actfield->dis[0].node[i]);
   /*------------------------- check for the dimensions of actnode->sol */
   pss_getdims_name_handle(actnode->sol.name,&dims[0],&dims[1],&dims[2],&(node_handles[i]),in,&ierr);
   if (ierr != 1) dserror("Cannot read visual data\n");
   if (dims[2] != sizeof(double)) dserror("Cannot read visual data\n");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol.fdim ||
       dims[1] != actnode->sol.sdim)
   amredef(&(actnode->sol),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol.name,&(actnode->sol),&(node_handles[i]),in,&ierr);
   if (ierr != 1) dserror("Cannot read visual data\n");
   /*-------------------------------------------------------------------*/
} 

/*------------------------------- delete the handle array for the nodes */
free(vis.node_handles);

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_readvispss */ 

#endif
