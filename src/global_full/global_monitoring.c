/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
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
void monitoring(FIELD *actfield,INT numf, INT actpos, INT actstep, DOUBLE time)
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

/*------------------------------------------------ enlarge if necessary */
if (actstep >= actmoni->val.sdim)
   amredef(&(actmoni->val),actmoni->val.fdim,actmoni->val.sdim+5,"DA");
if (actstep >= actmoni->time.fdim)   
   amredef(&(actmoni->time),actmoni->time.fdim+5,1,"DV");
/*-------------------------------------------------------- store values */
for (i=0;i<numnp;i++)
{
   for (j=0;j<MAXDOFPERNODE;j++)
   {
      numr=actmoni->onoff.a.ia[i][j];
      if (numr==-1) continue;
      nodepos = actmoni->monnodes.a.ia[i][1];
      actnode = &(actfield->dis[0].node[nodepos]);
      dsassert(j<actnode->numdf,"Monitoring fails!\n");
      actval = actnode->sol.a.da[actpos][j];
      actmoni->val.a.da[numr][actstep] = actval;
   }
}

actmoni->time.a.dv[actstep] = time;
actmoni->numstep+=1;

/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of monitoring */
