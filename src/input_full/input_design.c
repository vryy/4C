/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
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
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design                                           |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;

/*----------------------------------------------------------------------*
 | input of design                                        m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inpdesign()
{
#ifdef DEBUG 
dstrc_enter("inpctrdesign");
#endif

design = (DESIGN*)CCACALLOC(1,sizeof(DESIGN));
if (design==NULL) dserror("Unable to allocate DESIGN");
/*----------------------- input of design size (number of dnodes etc... */
inp_designsize();
/*------------------------------------------------input of design nodes */
inp_dnode();
/*------------------------------------------------input of design lines */
inp_dline();
/*---------------------------------------------input of design surfaces */
inp_dsurface();
/*----------------------------------------------input of design volumes */
inp_dvolume();
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpctrdesign */




/*----------------------------------------------------------------------*
 | input of design size                                   m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_designsize()
{
INT  i,ierr;
INT  numdnode;
INT  numdline;
INT  numdsurf;
INT  numdvol;
#ifdef DEBUG 
dstrc_enter("inp_designsize");
#endif

/*-------------------------------------- count number of design objects */
if (frfind("--DESIGN DESCRIPTION")==0) dserror("frfind: DESIGN DESCRIPTION is not in input file");
frread();
frint("NDPOINT",&numdnode,&ierr);
if (!ierr) dserror("Cannot read design");
frread();
frint("NDLINE",&numdline,&ierr);
if (!ierr) dserror("Cannot read design");
frread();
frint("NDSURF",&numdsurf,&ierr);
if (!ierr) dserror("Cannot read design");
frread();
frint("NDVOL",&numdvol,&ierr);
if (!ierr) dserror("Cannot read design");
/*-----------------------------------------------put sizes to structure */
design->ndnode = numdnode;
design->ndline = numdline;
design->ndsurf = numdsurf;
design->ndvol  = numdvol;
/*----------------------------------------------allocate design vectors */
design->dnode = (DNODE*)CCACALLOC(design->ndnode,sizeof(DNODE));
if (!design->dnode) dserror("Allocation of design nodes failed");
for (i=0; i<design->ndnode; i++) design->dnode[i].Id = i;

design->dline = (DLINE*)CCACALLOC(design->ndline,sizeof(DLINE));
if (!design->dline) dserror("Allocation of design lines failed");
for (i=0; i<design->ndline; i++) design->dline[i].Id = i;

design->dsurf = (DSURF*)CCACALLOC(design->ndsurf,sizeof(DSURF));
if (!design->dsurf) dserror("Allocation of design surfaces failed");
for (i=0; i<design->ndsurf; i++) design->dsurf[i].Id = i;

design->dvol = (DVOL*)CCACALLOC(design->ndvol,sizeof(DVOL));
if (!design->dvol) dserror("Allocation of design volumes failed");
for (i=0; i<design->ndvol; i++) design->dvol[i].Id = i;
/*-------------------------------------------------------------- rewind */
frrewind();
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_designsize */






/*----------------------------------------------------------------------*
 | input of design nodes                                  m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_dnode()
{
INT    i;
INT    readID;
DNODE *actdnode;
#ifdef DEBUG 
dstrc_enter("inp_dnode");
#endif
/*-------------------------------------------------------------- rewind */
frrewind();
/*------------------------ now read the description of the design nodes */
if (frfind("--DESIGN POINTS")==0) dserror("frfind: DESIGN POINTS is not in input file");
frread();
for (i=0; i<design->ndnode; i++)
{
   actdnode = &(design->dnode[i]);
   /*------------------------------------- GID writes fortran style Ids */
   readID   = actdnode->Id+1;
   /*---------------------------------------------- find the word POINT */
   read_1_dnode(actdnode,readID);
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_dnode */

/*----------------------------------------------------------------------*
 | input of one design nodes                              m.gee 1/02    |
 *----------------------------------------------------------------------*/
void read_1_dnode(DNODE *dnode, INT readId)
{
INT    i,ierr;
#ifdef DEBUG 
dstrc_enter("read_1_dnode");
#endif
/*----------------------------------------------------------------------*/
ierr=0;
frchk("POINT",&ierr);
while (!ierr) {frread(); frchk("POINT",&ierr);}

frchk("Num:",&ierr);
while (!ierr) {frread(); frchk("Num:",&ierr);}
frint("Num:",&i,&ierr);
if (!ierr) dserror("Cannot read DNODE");

if (i != readId) dserror("DNODEs got mixed up");

frint("conditions:",&(dnode->ncond),&ierr);
if (!ierr) dserror("Cannot read DNODE");

frchk("Coord:",&ierr);
while (!ierr) {frread(); frchk("Coord:",&ierr);}
frdouble_n("Coord:",dnode->x,3,&ierr);
if (!ierr) dserror("Cannot read DNODE");

frchk("END POINT",&ierr);
while (!ierr) {frread(); frchk("END POINT",&ierr);}
frread();
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of read_1_dnode */








/*----------------------------------------------------------------------*
 | input of design lines                                  m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_dline()
{
INT    i;
DLINE *actdline;
INT    readID;
#ifdef DEBUG 
dstrc_enter("inp_dline");
#endif

/*------------------------ now read the description of the design lines */
if (frfind("--DESIGN LINES")==0) dserror("frfind: DESIGN LINES is not in input file");
frread();

for (i=0; i<design->ndline; i++)
{
   actdline = &(design->dline[i]);
   /*------------------------------------- GID writes fortran style Ids */
   readID   = actdline->Id+1;
   /*---------------------------------------------- find the word POINT */
   read_1_dline(actdline,readID);
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_dline */


/*----------------------------------------------------------------------*
 | input of one design line                               m.gee 1/02    |
 *----------------------------------------------------------------------*/
void read_1_dline(DLINE *dline, INT readId)
{
INT    i,ierr;
INT    isnurb;

#ifdef DEBUG 
dstrc_enter("read_1_dline");
#endif
/*----------------------------------------------------------------------*/
ierr=0;
/* ------------------------------------------------------ read line typ */
frchk("LINE",&ierr);
while (!ierr) {frread(); frchk("LINE",&ierr);}
frchk("STLINE",&ierr);
if (ierr==1) dline->typ = stline;
frchk("NURBLINE",&ierr);
if (ierr==1) dline->typ = nurbline;
frchk("ARCLINE",&ierr);
if (ierr==1) dline->typ = arcline;

/* --------------------------------------------------- read line number */
frchk("Num:",&ierr);
while (!ierr) {frread(); frchk("Num:",&ierr);}
frint("Num:",&i,&ierr);
if (!ierr) dserror("Cannot read DLINE");
if (i != readId) dserror("DLINEs got mixed up");

/* ----------------------------- read number of conditions to this line */
frint("conditions:",&(dline->ncond),&ierr);
if (!ierr) dserror("Cannot read DLINE");

/* -------------------------------- read dpoints connected to this line */
frchk("Points:",&ierr);
while (!ierr) {frread(); frchk("Points:",&ierr);}
frint_n("Points:",dline->my_dnodeId,2,&ierr);
if (!ierr) dserror("Cannot read DLINE");
dline->my_dnodeId[0]--;
dline->my_dnodeId[1]--;
dline->ndnode=2;

/* -------------------------------------------- read arcline properties */
if(dline->typ == arcline)
{
  dline->props.arcline = (ARCLINE*)CCACALLOC(1,sizeof(ARCLINE));
  if (dline->props.arcline==NULL) dserror("Allocation of dline failed");
  /* ------------------------------------------------- read line number */
  frchk("2D center",&ierr);
  while (!ierr) {frread(); frchk("2D center",&ierr);}
  /* -------------------------------------- read radius of this arcline */
  frdouble("radius=",&(dline->props.arcline->radius),&ierr);
  if (!ierr) dserror("Cannot read DLINE");
  /* ------------------------------------- read initang of this arcline */
  frdouble("initang=",&(dline->props.arcline->initang),&ierr);
  if (!ierr) dserror("Cannot read DLINE");
  /* -------------------------------------- read endang of this arcline */
  frdouble("endang=",&(dline->props.arcline->endang),&ierr);
  if (!ierr) dserror("Cannot read DLINE");
  /* ------------------------------------------------ find total length */
  frchk("TotalLength",&ierr);
  while (!ierr) {frread(); frchk("TotalLength",&ierr);}
  /* ------------------------------------------------ read total length */
  frdouble("TotalLength=",&(dline->props.arcline->total_length),&ierr);
  if (!ierr) dserror("Cannot read DLINE");
}
/* -------------------------------------------- read arcline properties */
if(dline->typ == stline)
{
  dline->props.stline = (STLINE*)CCACALLOC(1,sizeof(STLINE));
  if (dline->props.stline==NULL) dserror("Allocation of dline failed");
}

frchk("END",&ierr);
while (!ierr) {frread(); frchk("END",&ierr);}
frread();
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of read_1_dline */




/*----------------------------------------------------------------------*
 | input of design surfaces                               m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_dsurface()
{
INT    i;
DSURF *actdsurf;
INT    readId;
#ifdef DEBUG 
dstrc_enter("inp_dsurface");
#endif
/*------------------------ now read the description of the design lines */
if (frfind("--DESIGN SURFACES")==0) dserror("frfind: DESIGN SURFACES is not in input file");
frread();
for (i=0; i<design->ndsurf; i++)
{
   actdsurf = &(design->dsurf[i]);
   /*------------------------------------- GID writes fortran style Ids */
   readId   = actdsurf->Id+1;
   /*---------------------------------------------- find the word POINT */
   read_1_dsurf(actdsurf,readId);
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_dsurface */



/*----------------------------------------------------------------------*
 | input of one design surface                            m.gee 1/02    |
 *----------------------------------------------------------------------*/
void read_1_dsurf(DSURF *dsurf, INT readId)
{
INT    i,ierr;
char   buffer[100];
#ifdef DEBUG 
dstrc_enter("read_1_dsurf");
#endif
/*----------------------------------------------------------------------*/
ierr=0;
frchk("NURBSURFACE",&ierr);
while (!ierr) {frread(); frchk("NURBSURFACE",&ierr);}

frchk("Num:",&ierr);
while (!ierr) {frread(); frchk("Num:",&ierr);}
frint("Num:",&i,&ierr);
if (!ierr) dserror("Cannot read DSURF");

if (i != readId) dserror("DSURF got mixed up");

frint("conditions:",&(dsurf->ncond),&ierr);
if (!ierr) dserror("Cannot read DSURF");

frchk("NumLines:",&ierr);
while (!ierr) {frread(); frchk("NumLines:",&ierr);}
frint("NumLines:",&(dsurf->ndline),&ierr);
if (!ierr) dserror("Cannot read DSURF");
amdef("my_dlineIDs",&(dsurf->my_dlineId),dsurf->ndline,2,"IA");

frchk("Line:",&ierr);
while (!ierr) {frread(); frchk("Line:",&ierr);}
for (i=0; i<dsurf->ndline; i++)
{
   frint("Line:",&(dsurf->my_dlineId.a.ia[i][0]),&ierr);
   if (!ierr) dserror("Cannot read DSURF");
   dsurf->my_dlineId.a.ia[i][0]--;
   frchar("Orientation:",buffer,&ierr);
   if (!ierr) dserror("Cannot read DSURF");
   if (strncmp("SAME1ST",buffer,7)==0) dsurf->my_dlineId.a.ia[i][1]=0;
   else                                dsurf->my_dlineId.a.ia[i][1]=1;
   frread();
}

frchk("END NURBSURFACE",&ierr);
while (!ierr) {frread(); frchk("END NURBSURFACE",&ierr);}
frread();
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of read_1_dsurf */





/*----------------------------------------------------------------------*
 | input of design volumes                                m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_dvolume()
{
INT   i;
DVOL *actdvol;
INT   readId;
#ifdef DEBUG 
dstrc_enter("inp_dvolume");
#endif

/*---------------------- now read the description of the design volumes */
if (frfind("--DESIGN VOLUMES")==0) dserror("frfind: DESIGN VOLUMES is not in input file");
frread();
for (i=0; i<design->ndvol; i++)
{
   actdvol = &(design->dvol[i]);
   /*------------------------------------- GID writes fortran style Ids */
   readId   = actdvol->Id+1;
   /*---------------------------------------------- find the word POINT */
   read_1_dvol(actdvol,readId);
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_dvolume */

/*----------------------------------------------------------------------*
 | input of one design volume                             m.gee 1/02    |
 *----------------------------------------------------------------------*/
void read_1_dvol(DVOL *dvol, INT readId)
{
INT    i,ierr;
char   buffer[100];
#ifdef DEBUG 
dstrc_enter("read_1_dvol");
#endif
/*----------------------------------------------------------------------*/
ierr=0;
frchk("VOLUME",&ierr);
while (!ierr) {frread(); frchk("VOLUME",&ierr);}

frchk("Num:",&ierr);
while (!ierr) {frread(); frchk("Num:",&ierr);}
frint("Num:",&i,&ierr);
if (!ierr) dserror("Cannot read DVOL");

if (i != readId) dserror("DVOLs got mixed up");

frint("conditions:",&(dvol->ncond),&ierr);
if (!ierr) dserror("Cannot read DVOL");

frchk("NumSurfaces:",&ierr);
while (!ierr) {frread(); frchk("NumSurfaces:",&ierr);}
frint("NumSurfaces:",&(dvol->ndsurf),&ierr);
if (!ierr) dserror("Cannot read DVOL");
amdef("my_dsurfIDs",&(dvol->my_dsurfId),dvol->ndsurf,2,"IA");

frchk("Surface:",&ierr);
while (!ierr) {frread(); frchk("Surface:",&ierr);}
for (i=0; i<dvol->ndsurf; i++)
{
   frint("Surface:",&(dvol->my_dsurfId.a.ia[i][0]),&ierr);
   if (!ierr) dserror("Cannot read DVOL");
   dvol->my_dsurfId.a.ia[i][0]--;
   frchar("Orientation:",buffer,&ierr);
   if (!ierr) dserror("Cannot read DVOL");
   if (strncmp("SAME1ST",buffer,7)==0) dvol->my_dsurfId.a.ia[i][1]=0;
   else                                dvol->my_dsurfId.a.ia[i][1]=1;
   frread();
}
frchk("END VOLUME",&ierr);
while (!ierr) {frread(); frchk("END VOLUME",&ierr);}
frread();
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of read_1_dvol */
