#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
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
int inp_designsize()
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

design->dline = (DLINE*)CALLOC(design->ndline,sizeof(DLINE));
if (!design->dline) dserror("Allocation of design lines failed");

design->dsurf = (DSURF*)CALLOC(design->ndsurf,sizeof(DSURF));
if (!design->dsurf) dserror("Allocation of design surfaces failed");

design->dvol = (DVOL*)CALLOC(design->ndvol,sizeof(DVOL));
if (!design->dvol) dserror("Allocation of design volumes failed");
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
/*------------------------------- find fe-nodes belonging to this dnode */
for (i=0; i<design->ndnode; i++)
{
   design->dnode[i].Id=i;
   frfind("-DESIGN-FE TOPOLOGY");
   frread();
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("DNODE",&dnode,&ierr);
      if (ierr==1)
      {
         if (dnode==i+1)
         {
            frint("NODE",&(design->dnode[i].mynode),&ierr);
            (design->dnode[i].mynode)--;
            goto nextdnode;
         }
      }
      frread();
   }
   nextdnode:
   frrewind();
}
/*------------------------ now read the description of the design nodes */
frfind("--DESIGN DESCRIPTION");
for (i=0; i<4; i++) frread();

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
frchk("{POINT",&ierr);
while (!ierr) {frread(); frchk("{POINT",&ierr);}

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

/*-------------------------------------------------------------- rewind */
frrewind();
/*---------------------------------- count number of nodes on this line */
for (i=0; i<design->ndline; i++)
{
   design->dline[i].Id = i;
   counter=0;
   frfind("-DESIGN-FE TOPOLOGY");
   frread();
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("DLINE",&dline,&ierr);
      if (ierr==1)
      {
         if (dline==i+1) counter++;
      }
      frread();
   }
   amdef("dline_node",&(design->dline[i].mynode),counter,1,"IV");
   frrewind();
}
/*------------------------------- find fe-nodes belonging to this dline */
for (i=0; i<design->ndline; i++)
{
   frfind("-DESIGN-FE TOPOLOGY");
   frread();
   counter=0;
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("DLINE",&dline,&ierr);
      if (ierr==1)
      {
         if (dline==i+1)
         {
            frint("NODE",&(design->dline[i].mynode.a.iv[counter]),&ierr);
            (design->dline[i].mynode.a.iv[counter])--;
            counter++;
         } 
      }
      frread();
   }
   frrewind();
}
/*------------------------ now read the description of the design lines */
frfind("--DESIGN DESCRIPTION");
for (i=0; i<4; i++) frread();

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
frchk("{STLINE",&isnurb);
frchk("{NURBLINE",&isnurb);

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

/*-------------------------------------------------------------- rewind */
frrewind();
/*---------------------------------- count number of nodes on this surf */
for (i=0; i<design->ndsurf; i++)
{
   design->dsurf[i].Id = i;
   counter=0;
   frfind("-DESIGN-FE TOPOLOGY");
   frread();
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("DSURFACE",&dsurf,&ierr);
      if (ierr==1)
      {
         if (dsurf==i+1) counter++;
      }
      frread();
   }
   amdef("dsurf_node",&(design->dsurf[i].mynode),counter,1,"IV");
   frrewind();
}
/*------------------------------- find fe-nodes belonging to this dsurf */
for (i=0; i<design->ndsurf; i++)
{
   frfind("-DESIGN-FE TOPOLOGY");
   frread();
   counter=0;
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("DSURFACE",&dsurf,&ierr);
      if (ierr==1)
      {
         if (dsurf==i+1)
         {
            frint("NODE",&(design->dsurf[i].mynode.a.iv[counter]),&ierr);
            (design->dsurf[i].mynode.a.iv[counter])--;
            counter++;
         } 
      }
      frread();
   }
   frrewind();
}
/*------------------------ now read the description of the design lines */
frfind("--DESIGN DESCRIPTION");
for (i=0; i<4; i++) frread();

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
int    nline;
char   buffer[100];
#ifdef DEBUG 
dstrc_enter("read_1_dsurf");
#endif
/*----------------------------------------------------------------------*/
ierr=0;
frchk("{NURBSURFACE",&ierr);
while (!ierr) {frread(); frchk("{NURBSURFACE",&ierr);}

frchk("Num:",&ierr);
while (!ierr) {frread(); frchk("Num:",&ierr);}
frint("Num:",&i,&ierr);
if (!ierr) dserror("Cannot read DSURF");

if (i != readId) dserror("DSURF got mixed up");

frint("conditions:",&(dsurf->ncond),&ierr);
if (!ierr) dserror("Cannot read DSURF");

frchk("NumLines:",&ierr);
while (!ierr) {frread(); frchk("NumLines:",&ierr);}
frint("NumLines:",&nline,&ierr);
if (!ierr) dserror("Cannot read DSURF");
amdef("my_dlineIDs",&(dsurf->my_dlineId),nline,2,"IA");

frchk("Line:",&ierr);
while (!ierr) {frread(); frchk("Line:",&ierr);}
for (i=0; i<nline; i++)
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

/*-------------------------------------------------------------- rewind */
frrewind();
/*---------------------------------- count number of nodes on this vol */
for (i=0; i<design->ndvol; i++)
{
   design->dvol[i].Id = i;
   counter=0;
   frfind("-DESIGN-FE TOPOLOGY");
   frread();
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("DVOLUME",&dvol,&ierr);
      if (ierr==1)
      {
         if (dvol==i+1) counter++;
      }
      frread();
   }
   amdef("dvol_node",&(design->dvol[i].mynode),counter,1,"IV");
   frrewind();
}
/*------------------------------- find fe-nodes belonging to this vol */
for (i=0; i<design->ndvol; i++)
{
   frfind("-DESIGN-FE TOPOLOGY");
   frread();
   counter=0;
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("DVOLUME",&dvol,&ierr);
      if (ierr==1)
      {
         if (dvol==i+1)
         {
            frint("NODE",&(design->dvol[i].mynode.a.iv[counter]),&ierr);
            (design->dvol[i].mynode.a.iv[counter])--;
            counter++;
         } 
      }
      frread();
   }
   frrewind();
}
/*---------------------- now read the description of the design volumes */
frfind("--DESIGN DESCRIPTION");
for (i=0; i<4; i++) frread();

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
int    nsurf;
char   buffer[100];
#ifdef DEBUG 
dstrc_enter("read_1_dvol");
#endif
/*----------------------------------------------------------------------*/
ierr=0;
frchk("{VOLUME",&ierr);
while (!ierr) {frread(); frchk("{VOLUME",&ierr);}

frchk("Num:",&ierr);
while (!ierr) {frread(); frchk("Num:",&ierr);}
frint("Num:",&i,&ierr);
if (!ierr) dserror("Cannot read DVOL");

if (i != readId) dserror("DVOLs got mixed up");

frint("conditions:",&(dvol->ncond),&ierr);
if (!ierr) dserror("Cannot read DVOL");

frchk("NumSurfaces:",&ierr);
while (!ierr) {frread(); frchk("NumSurfaces:",&ierr);}
frint("NumSurfaces:",&nsurf,&ierr);
if (!ierr) dserror("Cannot read DVOL");
amdef("my_dsurfIDs",&(dvol->my_dsurfId),nsurf,2,"IA");

frchk("Surface:",&ierr);
while (!ierr) {frread(); frchk("Surface:",&ierr);}
for (i=0; i<nsurf; i++)
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
