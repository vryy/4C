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
/*----------------------------------------------------------------------*
 | prototypes of functions callable only in this file                   |
 *----------------------------------------------------------------------*/
static void inp_read_locsys(char *string, INT actnum);
static void inpdesign_nodal_locsys(void);
static void inpdesign_line_locsys(void);
static void inpdesign_surf_locsys(void);
static void inpdesign_vol_locsys(void);

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
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;
/*----------------------------------------------------------------------*
 |                                                        genk 04/04    |
 | number of local co-ordinate systems                                  |
 | vector of structures of local co-ordinate systems                    |
 | defined in input_locsys.c                                            |
 *----------------------------------------------------------------------*/
INT            numlocsys;
struct _LOCSYS *locsys;
/*!--------------------------------------------------------------------- 
\brief input local co-ordinate system

<pre>                                                         genk 04/04
			     
</pre>
\return void                                               
                                 
------------------------------------------------------------------------*/
void inp_cond_locsys()
{
INT  ierr;
INT  i;
#ifdef DEBUG 
dstrc_enter("inp_cond_locsys");
#endif

/*------------------------ count the number of different locsys (max=5) */
numlocsys=0;
if (frfind("--LOCSYS1")==1)
{
  frread();
  frchk("---",&ierr);
  if (ierr==0) (numlocsys)++;
}

if (frfind("--LOCSYS2")==1)
{
  frread();
  frchk("---",&ierr);
  if (ierr==0) (numlocsys)++;
}

if (frfind("--LOCSYS3")==1)
{
  frread();
  frchk("---",&ierr);
  if (ierr==0) (numlocsys)++;
}

if (frfind("--LOCSYS4")==1)
{
  frread();
  frchk("---",&ierr);
  if (ierr==0) (numlocsys)++;
}

if (frfind("--LOCSYS5")==1)
{
  frread();
  frchk("---",&ierr);
  if (ierr==0) (numlocsys)++;
}

/*------------------------------------------------- allocate the cosys */
locsys = (LOCSYS*)CCACALLOC(numlocsys,sizeof(LOCSYS));
for (i=0; i<numlocsys; i++) locsys[i].Id=i;

/*----------------------------------------------------- read the cosys */
inp_read_locsys("--LOCSYS1",1);
inp_read_locsys("--LOCSYS2",2);
inp_read_locsys("--LOCSYS3",3);
inp_read_locsys("--LOCSYS4",4);
inp_read_locsys("--LOCSYS5",5);

/*--------------------------------------- no read the design conditions */
inpdesign_nodal_locsys();
inpdesign_line_locsys();
inpdesign_surf_locsys();
inpdesign_vol_locsys();

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_cond_locsys */

/*!--------------------------------------------------------------------- 
\brief input local co-ordinate system

<pre>                                                         genk 04/04

read the local co-system from the input file. The so defined system has
 to be cartesian up to now!
			     
</pre>
\param   *string     char           string to find
\param    actnum     INT            
\return void                                               
                                 
------------------------------------------------------------------------*/
static void inp_read_locsys(char *string, INT actnum)
{
INT      ierr;
INT      i;
INT      dlineId;
DOUBLE   norm;
DOUBLE   *xloc, *yloc, *zloc;
LOCSYS   *actlocsys;
DLINE    *actdline1, *actdline2;

#ifdef DEBUG 
dstrc_enter("inp_read_locsys");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------- check whether there is info on this locsys */
if (frfind(string)==0) goto end;
frread();
frchk("---",&ierr);
if (ierr==1) goto end;
/*-------------------------------------------------- define actual curve */
dsassert(actnum<=numlocsys,"cannot read local co-ordinate system!\n");
actlocsys = &(locsys[actnum-1]);

/*-------------------------------------------------------- start reading */
if (frfind(string)==0) goto end;
frread();

/*------------------------------------------------ allocate base vectors */
xloc=amdef("xloc",&actlocsys->xloc,3,1,"DV");
yloc=amdef("yloc",&actlocsys->yloc,3,1,"DV");
zloc=amdef("zloc",&actlocsys->zloc,3,1,"DV");

/*--------------------------------------------------- read typ of locsys */   
frchk("None",&ierr);
if (ierr==1) 
{
   actlocsys->locsystyp = locsys_none;
}
frchk("BASEVEC",&ierr);
if (ierr==1) 
{
   actlocsys->locsystyp = locsys_basevec;

   /* in this case the three base vector are explicitely specified in the
      input file
   */
   /* read the base vectors */
   frdouble_n("XLOC",xloc,3,&ierr);
   frdouble_n("YLOC",yloc,3,&ierr);
   frdouble_n("ZLOC",zloc,3,&ierr);
}   
frchk("LINEPLANE",&ierr);
if (ierr==1) 
{
   actlocsys->locsystyp = locsys_line_plane;
   /* in this case the system is created through a given DLine and a plane:
      e.g. LX_PXY means:
         the given line lies in the xy-plane and represents the x-direction
         the y-direction is chosen orthogonal in such a way that the
         resulting z-direction is positive
   */

   /* read the dlineId */
   frint("LINEID",&dlineId,&ierr);
   dlineId--;
   /* find dline */
   actdline1=NULL;
   for (i=0;i<design->ndline;i++)
   {
      if(design->dline[i].Id == dlineId)
      {
         actdline1= &(design->dline[i]);
         break;
      }
   }
   if (actdline1==NULL)
      dserror("error reading locsys\n");   
   frchk("LX_PXY",&ierr);
   if (ierr==1)
   {
      for (i=0;i<3;i++)
         xloc[i] = actdline1->dnode[1]->x[i] - actdline1->dnode[0]->x[i];
      if (FABS(xloc[2]) > EPS8)
         dserror("LX not in PXY\n");
      /* define yloc orthogonal in xy-plane to xloc */
      yloc[0] = -xloc[1];
      yloc[1] =  xloc[0];
      yloc[2] = ZERO;
      /* vector product to check zloc */
      zloc[0] = xloc[1]*yloc[2] - xloc[2]*yloc[1];
      zloc[1] = xloc[2]*yloc[0] - xloc[0]*yloc[2];
      zloc[2] = xloc[0]*yloc[1] - xloc[1]*yloc[0];
      if (zloc[2]<ZERO) /* change sign of yloc */
      {
         yloc[0] *=-ONE;
         yloc[1] *=-ONE;
         zloc[2] *=-ONE;
      }
   }
   else
      dserror("error reading locsys: given plane not implemented yet!\n");
}
frchk("LINELINE",&ierr);
if (ierr==1) 
{
   actlocsys->locsystyp = locsys_line_line;
   /* in this case the system is created through two given DLines which
      have to be orthogonal!
   */

   /* read the dlineIdx */
   frint("LINEIDX",&dlineId,&ierr);
   dlineId--;
   /* find dline */
   actdline1=NULL;
   for (i=0;i<design->ndline;i++)
   {
      if(design->dline[i].Id == dlineId)
      {
         actdline1= &(design->dline[i]);
         break;
      }
   }
   /* read the dlineIdy */
   frint("LINEIDY",&dlineId,&ierr);
   dlineId--;
   /* find dline */
   actdline2=NULL;
   for (i=0;i<design->ndline;i++)
   {
      if(design->dline[i].Id == dlineId)
      {
         actdline2= &(design->dline[i]);
         break;
      }
   }
   /* read the dlineIdz*/
   frint("LINEIDZ",&dlineId,&ierr);
   if (dlineId<=0)
      dserror("locsys via LINEDZ not implemented yet\n");
   if (actdline1==NULL || actdline2==NULL)
      dserror("error reading locsys: none or only one LineID given\n");
      
   /*--------------------------------------------------- x-base vector */
   for (i=0;i<3;i++)
      xloc[i] = actdline1->dnode[1]->x[i] - actdline1->dnode[0]->x[i];
   /*--------------------------------------------------- y-base vector */
   for (i=0;i<3;i++)
      yloc[i] = actdline2->dnode[1]->x[i] - actdline2->dnode[0]->x[i];

   /* vector product to get zloc */
   zloc[0] = xloc[1]*yloc[2] - xloc[2]*yloc[1];
   zloc[1] = xloc[2]*yloc[0] - xloc[0]*yloc[2];
   zloc[2] = xloc[0]*yloc[1] - xloc[1]*yloc[0];
}
frchk("FLUIDMASSCONS",&ierr);
if (ierr==1)
{
   /* this co-sys type allows to define normal and tangential BCs along
      a curved line (see Gresho & Sani "Incompressible Flow and the
      Finite Element Method (chapter 3.13.1 e) 
      The co-sys is defined by the normal and tangent at the points - 
      this has to be done later                                        */
   actlocsys->locsystyp = locsys_fmc;
   goto end;
}
/*----------------------------------------- normalise the base vectors */
norm = sqrt(DSQR(xloc[0])+DSQR(xloc[1])+DSQR(xloc[2]));
for (i=0;i<3;i++) xloc[i] /=norm;
norm = sqrt(DSQR(yloc[0])+DSQR(yloc[1])+DSQR(yloc[2]));
for (i=0;i<3;i++) yloc[i] /=norm;
norm = sqrt(DSQR(zloc[0])+DSQR(zloc[1])+DSQR(zloc[2]));
for (i=0;i<3;i++) zloc[i] /=norm;

/*------------------------------------------------ plausibility checks */
/* vectors have to be orthogonal: check scalar products */
if (FABS(xloc[0]*yloc[0]+xloc[1]*yloc[1]+xloc[2]*yloc[2])>EPS8)
   dserror("locsys base vectors are not orthogonal!\n");
if (FABS(xloc[0]*zloc[0]+xloc[1]*zloc[1]+xloc[2]*zloc[2])>EPS8)
   dserror("locsys base vectors are not orthogonal!\n");
if (FABS(zloc[0]*yloc[0]+zloc[1]*yloc[1]+zloc[2]*yloc[2])>EPS8)
   dserror("locsys base vectors are not orthogonal!\n");

/* we need a right hand system: check vector product */
if (FABS(zloc[0]-(xloc[1]*yloc[2] - xloc[2]*yloc[1]))
   +FABS(zloc[1]-(xloc[2]*yloc[0] - xloc[0]*yloc[2]))
   +FABS(zloc[2]-(xloc[0]*yloc[1] - xloc[1]*yloc[0]))>EPS8)
    dserror("chosen local cosys no right hand system!\n");

/*------------------------------------------ store the direction cosini */
actlocsys->lXx=xloc[0];
actlocsys->lXy=yloc[0];  
actlocsys->lXz=zloc[0];
actlocsys->lYx=xloc[1];
actlocsys->lYy=yloc[1];  
actlocsys->lYz=zloc[1];
actlocsys->lZx=xloc[2];
actlocsys->lZy=yloc[2]; 
actlocsys->lZz=zloc[2];
   
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_read_locsys */

/*!--------------------------------------------------------------------- 
\brief input local co-ordinate system

<pre>                                                         genk 04/04
			     
</pre>
\return void                                               
                                 
------------------------------------------------------------------------*/
static void inpdesign_nodal_locsys()
{
INT    i;
INT    ierr;
INT    ndnode;
INT    dnodeId;
char  *colptr;
DNODE *actdnode;

#ifdef DEBUG 
dstrc_enter("inpdesign_nodal_locsys");
#endif

/*---------------------------------------------------------- initialise */
for (i=0;i<design->ndnode;i++) design->dnode[i].locsysId=0;

/*-------------------- find the beginning of nodal dirichlet conditions */
if (frfind("--DESIGN POINT LOCSYS CONDITIONS")==0) goto end;
frread();
/*------------------------ read number of design points with conditions */
frint("DPOINT",&ndnode,&ierr);
dsassert(ierr==1,"Cannot read design-nodal locsys conditions");
frread();

/*-------------------------------------- start reading the design nodes */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design node Id */
   frint("E",&dnodeId,&ierr);
   dsassert(ierr==1,"Cannot read design-nodal locsys conditions");
   dnodeId--;
   /*--------------------------------------------------- find the dnode */
   actdnode=NULL;
   for (i=0; i<design->ndnode; i++)
   {
      if (design->dnode[i].Id ==  dnodeId) 
      {
         actdnode = &(design->dnode[i]);
         break;
      }
   }
   dsassert(actdnode!=NULL,"Cannot read design-nodal locsys conditions");
  /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-nodal locsys conditions");
   colptr++;
   actdnode->locsysId = strtod(colptr,&colptr);
   frread();
}

/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_nodal_locsys */


/*!--------------------------------------------------------------------- 
\brief input local co-ordinate system

<pre>                                                         genk 04/04
			     
</pre>
\return void                                               
                                 
------------------------------------------------------------------------*/
static void inpdesign_line_locsys()
{
INT    i;
INT    ierr;
INT    ndline;
INT    dlineId;
char  *colptr;
DLINE *actdline;

#ifdef DEBUG 
dstrc_enter("inpdesign_line_locsys");
#endif

/*---------------------------------------------------------- initialise */
for (i=0;i<design->ndline;i++) design->dline[i].locsysId=0;

/*-------------------- find the beginning of nodal dirichlet conditions */
if (frfind("--DESIGN LINE LOCSYS CONDITIONS")==0) goto end;
frread();
/*------------------------ read number of design points with conditions */
frint("DLINE",&ndline,&ierr);
dsassert(ierr==1,"Cannot read design-line locsys conditions");
frread();

/*-------------------------------------- start reading the design lines */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design line Id */
   frint("E",&dlineId,&ierr);
   dsassert(ierr==1,"Cannot read design-line locsys conditions");
   dlineId--;
   /*--------------------------------------------------- find the dline */
   actdline=NULL;
   for (i=0; i<design->ndline; i++)
   {
      if (design->dline[i].Id ==  dlineId) 
      {
         actdline = &(design->dline[i]);
         break;
      }
   }
   dsassert(actdline!=NULL,"Cannot read design-line locsys conditions");
  /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-line locsys conditions");
   colptr++;
   actdline->locsysId = strtod(colptr,&colptr);
   frread();
}

/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_line_locsys */


/*!--------------------------------------------------------------------- 
\brief input local co-ordinate system

<pre>                                                         genk 04/04
			     
</pre>
\return void                                               
                                 
------------------------------------------------------------------------*/
static void inpdesign_surf_locsys()
{
INT    i;
INT    ierr;
INT    ndsurf;
INT    dsurfId;
char  *colptr;
DSURF *actdsurf;

#ifdef DEBUG 
dstrc_enter("inpdesign_surf_locsys");
#endif

/*---------------------------------------------------------- initialise */
for (i=0;i<design->ndsurf;i++) design->dsurf[i].locsysId=0;

/*-------------------- find the beginning of nodal dirichlet conditions */
if (frfind("--DESIGN SURF LOCSYS CONDITIONS")==0) goto end;
frread();
/*------------------------ read number of design points with conditions */
frint("DSURF",&ndsurf,&ierr);
dsassert(ierr==1,"Cannot read design-surf locsys conditions");
frread();

/*-------------------------------------- start reading the design surfs */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design surf Id */
   frint("E",&dsurfId,&ierr);
   dsassert(ierr==1,"Cannot read design-surf locsys conditions");
   dsurfId--;
   /*--------------------------------------------------- find the dsurf */
   actdsurf=NULL;
   for (i=0; i<design->ndsurf; i++)
   {
      if (design->dsurf[i].Id ==  dsurfId) 
      {
         actdsurf = &(design->dsurf[i]);
         break;
      }
   }
   dsassert(actdsurf!=NULL,"Cannot read design-surf locsys conditions");
  /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-surf locsys conditions");
   colptr++;
   actdsurf->locsysId = strtod(colptr,&colptr);
   frread();
}

/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_surf_locsys */

/*!--------------------------------------------------------------------- 
\brief input local co-ordinate system

<pre>                                                         genk 04/04
			     
</pre>
\return void                                               
                                 
------------------------------------------------------------------------*/
static void inpdesign_vol_locsys()
{
INT    i;
INT    ierr;
INT    ndvol;
INT    dvolId;
char  *colptr;
DVOL *actdvol;

#ifdef DEBUG 
dstrc_enter("inpdesign_vol_locsys");
#endif

/*---------------------------------------------------------- initialise */
for (i=0;i<design->ndvol;i++) design->dvol[i].locsysId=0;

/*-------------------- find the beginning of nodal dirichlet conditions */
if (frfind("--DESIGN VOL LOCSYS CONDITIONS")==0) goto end;
frread();
/*------------------------ read number of design points with conditions */
frint("DVOL",&ndvol,&ierr);
dsassert(ierr==1,"Cannot read design-vol locsys conditions");
frread();

/*-------------------------------------- start reading the design vols */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design vol Id */
   frint("E",&dvolId,&ierr);
   dsassert(ierr==1,"Cannot read design-vol locsys conditions");
   dvolId--;
   /*--------------------------------------------------- find the dvol */
   actdvol=NULL;
   for (i=0; i<design->ndvol; i++)
   {
      if (design->dvol[i].Id ==  dvolId) 
      {
         actdvol = &(design->dvol[i]);
         break;
      }
   }
   dsassert(actdvol!=NULL,"Cannot read design-vol locsys conditions");
  /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-vol locsys conditions");
   colptr++;
   actdvol->locsysId = strtod(colptr,&colptr);
   frread();
}

/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_vol_locsys */
