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
int  ierr;
#ifdef DEBUG 
dstrc_enter("inpctrdesign");
#endif

design = (DESIGN*)CALLOC(1,sizeof(DESIGN));
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
int  i,ierr;
int  numdnode;
int  numdline;
int  numdsurf;
int  numdvol;
#ifdef DEBUG 
dstrc_enter("inp_designsize");
#endif

/*-------------------------------------- count number of design objects */
frfind("--DESIGN DESCRIPTION");
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
design->dnode = (DNODE*)CALLOC(design->ndnode,sizeof(DNODE));
if (!design->dnode) dserror("Allocation of design nodes failed");
for (i=0; i<design->ndnode; i++) design->dnode[i].Id = i;

design->dline = (DLINE*)CALLOC(design->ndline,sizeof(DLINE));
if (!design->dline) dserror("Allocation of design lines failed");
for (i=0; i<design->ndline; i++) design->dline[i].Id = i;

design->dsurf = (DSURF*)CALLOC(design->ndsurf,sizeof(DSURF));
if (!design->dsurf) dserror("Allocation of design surfaces failed");
for (i=0; i<design->ndsurf; i++) design->dsurf[i].Id = i;

design->dvol = (DVOL*)CALLOC(design->ndvol,sizeof(DVOL));
if (!design->dvol) dserror("Allocation of design volumes failed");
for (i=0; i<design->ndvol; i++) design->dvol[i].Id = i;
/*-------------------------------------------------------------- rewind */
frrewind();
/*----------------------------------------------------------------------*/
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
int    i,ierr;
int    maxdnode=0;
int    dnode=0;  
int    readID;
DNODE *actdnode;
#ifdef DEBUG 
dstrc_enter("inp_dnode");
#endif
/*-------------------------------------------------------------- rewind */
frrewind();
/*------------------------ now read the description of the design nodes */
frfind("--DESIGN POINTS");
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
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_dnode */

/*----------------------------------------------------------------------*
 | input of one design nodes                              m.gee 1/02    |
 *----------------------------------------------------------------------*/
void read_1_dnode(DNODE *dnode, int readId)
{
int    i,ierr;
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
int    i,ierr,counter;
int    maxdline=0;
int    dline=0;  
DLINE *actdline;
int    readID;
#ifdef DEBUG 
dstrc_enter("inp_dline");
#endif

/*------------------------ now read the description of the design lines */
frfind("--DESIGN LINES");
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
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_dline */
/*----------------------------------------------------------------------*
 | input of one design line                               m.gee 1/02    |
 *----------------------------------------------------------------------*/
void read_1_dline(DLINE *dline, int readId)
{
int    i,ierr;
int    isnurb;

#ifdef DEBUG 
dstrc_enter("read_1_dline");
#endif
/*----------------------------------------------------------------------*/
ierr=0;
frchk("LINE",&ierr);
while (!ierr) {frread(); frchk("LINE",&ierr);}
frchk("STLINE",&isnurb);
frchk("NURBLINE",&isnurb);

frchk("Num:",&ierr);
while (!ierr) {frread(); frchk("Num:",&ierr);}
frint("Num:",&i,&ierr);
if (!ierr) dserror("Cannot read DLINE");

if (i != readId) dserror("DLINEs got mixed up");

frint("conditions:",&(dline->ncond),&ierr);
if (!ierr) dserror("Cannot read DLINE");

frchk("Points:",&ierr);
while (!ierr) {frread(); frchk("Points:",&ierr);}
frint_n("Points:",dline->my_dnodeId,2,&ierr);
if (!ierr) dserror("Cannot read DLINE");
dline->my_dnodeId[0]--;
dline->my_dnodeId[1]--;
dline->ndnode=2;

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
int    i,ierr,counter;
int    maxdsurf=0;
int    dsurf=0;  
DSURF *actdsurf;
int    readId;
#ifdef DEBUG 
dstrc_enter("inp_dsurface");
#endif
/*------------------------ now read the description of the design lines */
frfind("--DESIGN SURFACES");
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
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_dsurface */

/*----------------------------------------------------------------------*
 | input of one design surface                            m.gee 1/02    |
 *----------------------------------------------------------------------*/
void read_1_dsurf(DSURF *dsurf, int readId)
{
int    i,ierr;
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
int   i,ierr,counter;
int   maxdvol=0;
int   dvol=0;  
DVOL *actdvol;
int   readId;
#ifdef DEBUG 
dstrc_enter("inp_dvolume");
#endif

/*---------------------- now read the description of the design volumes */
frfind("--DESIGN VOLUMES");
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
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_dvolume */

/*----------------------------------------------------------------------*
 | input of one design volume                             m.gee 1/02    |
 *----------------------------------------------------------------------*/
void read_1_dvol(DVOL *dvol, int readId)
{
int    i,ierr;
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
