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
struct _FILES           allfiles;
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
struct _MONITOR *moni;

/*----------------------------------------------------------------------*
 | input of monitoring data                               genk 01/03    |
 *----------------------------------------------------------------------*/
void inp_monitor()
{
INT i,j;
INT ierr;
INT cs=0;
INT cf=0;
INT ca=0;
INT counter=0;
INT numsf,numff,numaf;
INT **snode=NULL, **fnode=NULL, **anode=NULL;
INT actnum;
INT **aonoff=NULL, **fonoff=NULL, **sonoff=NULL;
INT actId;
MONITOR *smoni=NULL, *fmoni=NULL, *amoni=NULL;
FIELD *actfield;
NODE *actnode;
char  *colptr;

#ifdef DEBUG 
dstrc_enter("inp_monitor");
#endif

/*----------------------------------------------------------------------*/
moni = (MONITOR*)CCACALLOC(genprob.numfld,sizeof(MONITOR));
if (!moni) dserror("Allocation of MONITOR failed");

for (i=0;i<genprob.numfld;i++)
{
   moni[i].numnp=0;
   moni[i].numval=0;
   moni[i].numstep=0;
}
if (ioflags.monitor==0) goto end;

numsf=genprob.numsf;
numff=genprob.numff;
numaf=genprob.numaf;

/*--------------------------------------------------- first count nodes */
if (frfind("--MONITORING")==1)
{
  frread();
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    frchk("STRUCTURE",&ierr);
    if (ierr==1) cs++;
    frchk("FLUID",&ierr);
    if (ierr==1) cf++;
    frchk("ALE",&ierr);
    if (ierr==1) ca++;
    frread();
  }
}


/*---------------------------------------------------------- check data */
/*----------------------------------------------------------- structure */
if (numsf>=0 && cs>0)
{
   smoni=&(moni[numsf]);
   smoni->numnp=cs;  
   snode = amdef("monnodes",&(moni[numsf].monnodes),cs,2,"IA");
   sonoff = amdef("onoff",&(moni[numsf].onoff),cs,MAXDOFPERNODE,"IA");
}
else if (numsf>=0 && cs==0)
{
   smoni=&(moni[numsf]);
   smoni->numnp=cs;
}
else if (numsf<0 && cs>0)
{
   dserror("No struct field - Monitoring fails!\n");
}

/*-------------------------------------------------------------- fluid */
if (numff>=0 && cf>0)
{
   fmoni=&(moni[numff]);
   fmoni->numnp=cf;  
   fnode = amdef("monnodes",&(moni[numff].monnodes),cf,2,"IA");
   fonoff = amdef("onoff",&(moni[numff].onoff),cf,MAXDOFPERNODE,"IA");
}
else if (numff>=0 && cf==0)
{
   fmoni=&(moni[numff]);
   fmoni->numnp=cf;
}
else if (numff<0 && cf>0)
{
   dserror("No fluid field - Monitoring fails!\n");
}

/*---------------------------------------------------------------- ale */
if (numaf>=0 && ca>0)
{
   amoni=&(moni[numaf]);
   amoni->numnp=ca;  
   anode = amdef("monnodes",&(moni[numaf].monnodes),ca,2,"IA");
   aonoff = amdef("onoff",&(moni[numaf].onoff),ca,MAXDOFPERNODE,"IA");
}
else if (numaf>=0 && ca==0)
{
   amoni=&(moni[numaf]);
   amoni->numnp=ca;  
}
else if (numaf<0 && ca>0)
{
   dserror("No ale field - Monitoring fails!\n");
}

cs=0;
cf=0;
ca=0;

/*---------------------------------------------- read in global node Ids*/
if (frfind("--MONITORING")==1)
{
  frread();
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    frint("STRUCTURE",&actnum,&ierr);
    if (ierr==1)
    {
      actnum--;
      snode[cs][0] = actnum;
      /*------------------------------ move pointer behind the "-" sign */
      colptr = strstr(allfiles.actplace,"-");
      dsassert(colptr!=NULL,"Cannot read monitoring ");
      colptr++;      
      for (i=0; i<MAXDOFPERNODE; i++)     
        sonoff[cs][i] = strtol(colptr,&colptr,10);
      cs++;
    }
    frint("FLUID",&actnum,&ierr);
    if (ierr==1)
    {
      actnum--;
      fnode[cf][0] = actnum;
      /*------------------------------ move pointer behind the "-" sign */
      colptr = strstr(allfiles.actplace,"-");
      dsassert(colptr!=NULL,"Cannot read monitoring ");
      colptr++;      
      for (i=0; i<MAXDOFPERNODE; i++)     
        fonoff[cf][i] = strtol(colptr,&colptr,10);
      cf++;
    }
    frint("ALE",&actnum,&ierr);
    if (ierr==1)
    {
      actnum--;
      anode[ca][0] = actnum;
      /*------------------------------ move pointer behind the "-" sign */
      colptr = strstr(allfiles.actplace,"-");
      dsassert(colptr!=NULL,"Cannot read monitoring ");
      colptr++;      
      for (i=0; i<MAXDOFPERNODE; i++)     
        aonoff[ca][i] = strtol(colptr,&colptr,10);
      ca++;
    }
    frread();
  }
}

/*-------------------------------------------------------- count values */
if (numsf>=0)
{
   counter = 0;
   for (i=0;i<smoni->numnp; i++) 
   {
      for (j=0;j<MAXDOFPERNODE; j++)
      {
         if (sonoff[i][j]!=0) 
         {         
	    sonoff[i][j]=counter;
	    counter++;      	
         }
         else sonoff[i][j]=-1;
      }
   }
   smoni->numval = counter;
}

if (numff>=0)
{
   counter = 0;
   for (i=0;i<fmoni->numnp; i++) 
   {
      for (j=0;j<MAXDOFPERNODE; j++)
      {
         if (fonoff[i][j]!=0) 
         {         
	    fonoff[i][j]=counter;
   	    counter++;      	
         }
         else fonoff[i][j]=-1;
      }
   }
   fmoni->numval = counter;
}

if (numaf>=0)
{
   counter = 0;
   for (i=0;i<amoni->numnp; i++)
   {
      for (j=0;j<MAXDOFPERNODE; j++)
      {
         if (aonoff[i][j]!=0) 
         {         
	    aonoff[i][j]=counter;
	    counter++;      	
         }
         else aonoff[i][j]=-1;
      }
   }
   amoni->numval = counter;
}

/*-------------------------------------------- determine local node Ids */
if (numsf>=0 && smoni->numnp>0)
{
   amdef("val",&(moni[numsf].val),smoni->numval,5,"DA");
   amdef("time",&(moni[numsf].time),5,1,"DV");   
   actfield = &(field[numsf]);
   for (i=0;i<actfield->dis[0].numnp;i++)
   {
      actnode = &(actfield->dis[0].node[i]);
      actId = actnode->Id;
      for (j=0;j<smoni->numnp;j++)
      {
         if (snode[j][0] == actId)
	 {
	    snode[j][1] = i;
	    break;
	 }
      }
   }
}
if (numff>=0 && fmoni->numnp>0)
{
   amdef("val",&(moni[numff].val),fmoni->numval,5,"DA");
   amdef("time",&(moni[numff].time),5,1,"DV"); 
   actfield = &(field[numff]);
   for (i=0;i<actfield->dis[0].numnp;i++)
   {
      actnode = &(actfield->dis[0].node[i]);
      actId = actnode->Id;
      for (j=0;j<fmoni->numnp;j++)
      {
         if (fnode[j][0] == actId)
	 {
	    fnode[j][1] = i;
	    break;
	 }
      }
   }
}
if (numaf>=0 && amoni->numnp>0)
{
   amdef("val",&(moni[numaf].val),amoni->numval,5,"DA");
   amdef("time",&(moni[numaf].time),5,1,"DV"); 
   actfield = &(field[numaf]);
   for (i=0;i<actfield->dis[0].numnp;i++)
   {
      actnode = &(actfield->dis[0].node[i]);
      actId = actnode->Id;
      for (j=0;j<amoni->numnp;j++)
      {
         if (anode[j][0] == actId)
	 {
	    anode[j][1] = i;
	    break;
	 }
      }
   }
}

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_monitor */
