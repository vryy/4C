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

design = (DESIGN*)calloc(1,sizeof(DESIGN));
if (design==NULL) dserror("Unable to allocate DESIGN");
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
 | input of design nodes                                  m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_dnode()
{
int  i,ierr;
int  maxdnode=0;
int  dnode=0;  
#ifdef DEBUG 
dstrc_enter("inp_dnode");
#endif

/*---------------------------------------- count number of design nodes */
frfind("-DESIGN-FE TOPOLOGY");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frint("DNODE",&dnode,&ierr);
   if (ierr==1)
   {
      if (dnode>maxdnode) maxdnode=dnode;
   }
   frread();
}
/*--------------------------------------------- lowest dnodenumber is 0 */
design->ndnode = maxdnode+1;
/*------------------------------------------------allocate design nodes */
design->dnode = (DNODE*)calloc(design->ndnode,sizeof(DNODE));
if (design->dnode==NULL) dserror("Allocation of design nodes failed");
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
         if (dnode==i)
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
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_dnode */





/*----------------------------------------------------------------------*
 | input of design lines                                  m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_dline()
{
int  i,ierr,counter;
int  maxdline=0;
int  dline=0;  
#ifdef DEBUG 
dstrc_enter("inp_dline");
#endif

/*---------------------------------------- count number of design lines */
frfind("-DESIGN-FE TOPOLOGY");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frint("DLINE",&dline,&ierr);
   if (ierr==1)
   {
      if (dline>maxdline) maxdline=dline;
   }
   frread();
}
/*--------------------------------------------- lowest dlinenumber is 0 */
design->ndline = maxdline+1;
/*------------------------------------------------allocate design nodes */
design->dline = (DLINE*)calloc(design->ndline,sizeof(DLINE));
if (design->dline==NULL) dserror("Allocation of design lines failed");
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
         if (dline==i) counter++;
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
         if (dline==i)
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
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_dline */






/*----------------------------------------------------------------------*
 | input of design surfaces                               m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_dsurface()
{
int  i,ierr,counter;
int  maxdsurf=0;
int  dsurf=0;  
#ifdef DEBUG 
dstrc_enter("inp_dsurface");
#endif

/*---------------------------------------- count number of design lines */
frfind("-DESIGN-FE TOPOLOGY");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frint("DSURFACE",&dsurf,&ierr);
   if (ierr==1)
   {
      if (dsurf>maxdsurf) maxdsurf=dsurf;
   }
   frread();
}
/*--------------------------------------------- lowest dsurfnumber is 0 */
design->ndsurf = maxdsurf+1;
/*------------------------------------------------allocate design nodes */
design->dsurf  = (DSURF*)calloc(design->ndsurf,sizeof(DSURF));
if (design->dsurf==NULL) dserror("Allocation of design surfaces failed");
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
         if (dsurf==i) counter++;
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
         if (dsurf==i)
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
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_dsurface */





/*----------------------------------------------------------------------*
 | input of design volumes                                m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_dvolume()
{
int  i,ierr,counter;
int  maxdvol=0;
int  dvol=0;  
#ifdef DEBUG 
dstrc_enter("inp_dvolume");
#endif

/*---------------------------------------- count number of design lines */
frfind("-DESIGN-FE TOPOLOGY");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frint("DVOLUME",&dvol,&ierr);
   if (ierr==1)
   {
      if (dvol>maxdvol) maxdvol=dvol;
   }
   frread();
}
/*--------------------------------------------- lowest dsurfnumber is 0 */
design->ndvol = maxdvol+1;
/*------------------------------------------------allocate design nodes */
design->dvol  = (DVOL*)calloc(design->ndvol,sizeof(DVOL));
if (design->dvol==NULL) dserror("Allocation of design volumes failed");
/*-------------------------------------------------------------- rewind */
frrewind();
/*---------------------------------- count number of nodes on this surf */
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
         if (dvol==i) counter++;
      }
      frread();
   }
   amdef("dvol_node",&(design->dvol[i].mynode),counter,1,"IV");
   frrewind();
}
/*------------------------------- find fe-nodes belonging to this dsurf */
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
         if (dvol==i)
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
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_dvolume */

