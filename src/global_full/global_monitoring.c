/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127

General explanation:
   It might be interesting to have the results at a specific FE-nodes
   over the whole simulation time in one file. This file may be used to
   plot the results with another program like Excel or GNUplot.
   The three functions inp_monitor(), monitoring, out_monitor are used
   to create this information. out_monitor has to be called at the 
   beginning of the calculation to write the headers to the files
   (init=1)
   During the input phase the structure MONITOR is allocated for each
   field. There the information of the monitoring data is stored
   (which nodes, which DOFs, position of the node in the field).
   During the calculation at the end of a time step the function
   monitoring has to be called. It extracts the solution for the 
   specific monitoring nodes from the nodal sol_array and writes it
   to the <field>.mon file. For each field a separate output file 
   exists.
   
Input:
   The input is done manually (no input via GID) by deleting the
   // signs in the MONITORING block. The MONITORING
   block within the input file conains the following information:
   
--------------------------------------------------------MONITORING
FIELD globalId FLAG FLAG FLAG FLAG FLAG FLAG
FLUID 1 - 0 0 1 0 0 0 
STRUCTURE 2 - 0 2 0 0 0 0 
ALE 3 - 1 1 0 0 0 0 

   This means:
   From the fluid field the nodal solution (3rd dof) of the node 
   with the global (GID) ID 1 is written to *.fluid.mon
   From the strucutre field the nodal solution (2nd dof) of the node
   with the global (GID) ID 2 is written to *.structure.mon
   From the ALE field the nodal solution (1st and 2nd dof) of the node
   with the global (GID) ID 3 is written to *.ale.mon
   
   More Nodes can be added by copying a line!

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
extern struct _FILES        allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 | struct _FIELD         *field;                                        |
 *----------------------------------------------------------------------*/
extern struct _FIELD       *field;   

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*---------------------------------------------------------------------*
 | monotoring informations                                  genk 01/03 |
 *---------------------------------------------------------------------*/
extern struct _MONITOR *moni;

/*----------------------------------------------------------------------*
 | input of monitoring data                               genk 01/03    |
 *----------------------------------------------------------------------*/
void monitoring(
                  FIELD         *actfield,
                  INT            numf,
                  INT            actpos,
                  DOUBLE         time
               ) 
{
INT i,j;
INT numnp;
INT numr;
INT nodepos;
DOUBLE actval;
MONITOR *actmoni;
NODE *actnode;

#ifdef DEBUG 
dstrc_enter("monitoring");
#endif

actmoni = &(moni[numf]);
numnp   = actmoni->numnp;

if (numnp==0) goto end;

/*-------------------------------------------------------- store values */
for (i=0;i<numnp;i++)
{
   for (j=0;j<MAXDOFPERNODE;j++)
   {
      numr=actmoni->onoff.a.ia[i][j];
      if (numr==-1) continue;
      nodepos = actmoni->monnodes.a.ia[i][1];
      actnode = &(actfield->dis[0].node[nodepos]);
      dsassert(j<actnode->sol.sdim,"Monitoring fails!\n");
      actval = actnode->sol.a.da[actpos][j];
      actmoni->val.a.dv[numr] = actval;
   }
}

out_monitor(actfield,numf,time,0);

/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of monitoring */
