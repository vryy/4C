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
/*! 
\addtogroup INPUT 
*//*! @{ (documentation module open)*/
#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;



/*!----------------------------------------------------------------------
\brief 

<pre>                                                         oezdem 8/03
</pre>
*----------------------------------------------------------------------*/
#ifdef WALLCONTACT
extern struct _WALL_CONTACT contact;
#endif





/*----------------------------------------------------------------------*
 | prototypes of functions callable only in this file                   |
 *----------------------------------------------------------------------*/
static void inpdesign_nodal_dirich(void);
static void inpdesign_line_dirich(void);
static void inpdesign_surf_dirich(void);
static void inpdesign_vol_dirich(void);

static void inpdesign_nodal_neum(void);
static void inpdesign_line_neum(void);
static void inpdesign_surf_neum(void);
static void inpdesign_vol_neum(void);

static void inpdesign_nodal_couple(void);
static void inpdesign_line_couple(void);
static void inpdesign_surf_couple(void);
static void inpdesign_vol_couple(void);

#ifdef D_FSI
static void inpdesign_line_fsicouple(void);
#endif
static void inpdesign_nodal_freesurf(void);
static void inpdesign_line_freesurf(void);
static void inpdesign_line_liftdrag(void);
static void inpdesign_surf_liftdrag(void);
static void inpdesign_surf_stability(void);
static void inpdesign_vol_stability(void);

#ifdef D_AXISHELL
static void inpdesign_line_thickness(void);
static void inpdesign_line_axishellload(void);
static void inpdesign_point_axishellcos(void);
#endif

#ifdef WALLCONTACT
static void inpdesign_line_contact(void);
#endif

/*----------------------------------------------------------------------*
 | input of conditions                                    m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_conditions()
{
#ifdef DEBUG 
dstrc_enter("inp_conditions");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------ input of time curves */
inp_cond_curve();
/* input of spatial functions */
inp_cond_funct();
/*--------------------------------- input of nodal dirichlet conditions */
inpdesign_nodal_dirich();
/*---------------------------------- input of line dirichlet conditions */
inpdesign_line_dirich();
/*------------------------------- input of surface dirichlet conditions */
inpdesign_surf_dirich();
/*-------------------------------- input of volume dirichlet conditions */
inpdesign_vol_dirich();
/*----------------------------------- input of nodal neumann conditions */
inpdesign_nodal_neum();
/*------------------------------------ input of line neumann conditions */
inpdesign_line_neum();
/*--------------------------------- input of surface neumann conditions */
inpdesign_surf_neum();
/*---------------------------------- input of volume neumann conditions */
inpdesign_vol_neum();
/*---------------------------------- input of nodal coupling conditions */
inpdesign_nodal_couple();
/*----------------------------------- input of line coupling conditions */
inpdesign_line_couple();
/*----------------------------------- input of surf coupling conditions */
inpdesign_surf_couple();
/*------------------------------------ input of vol coupling conditions */
inpdesign_vol_couple();
/*----------------------------------------- input of wall contact stuff */
#ifdef WALLCONTACT
inpdesign_line_contact();
#endif

/*------------------------------- input of line FSI coupling conditions */
#ifdef D_FSI
if (genprob.probtyp==prb_fsi)
{
   inpdesign_line_fsicouple();
}
#endif
#ifdef D_FLUID
if (genprob.probtyp==prb_fluid || genprob.probtyp==prb_fsi)
{
/*------------------------------ input of nodal free surface conditions */
   inpdesign_nodal_freesurf();
/*------------------------------- input of line free surface conditions */   
   inpdesign_line_freesurf();
/*---------------------------------- input of line lift&drag definition */   
   inpdesign_line_liftdrag();
/*---------------------------------- input of surf lift&drag definition */   
   inpdesign_surf_liftdrag();
/*------------------------------- input of surface stability definition */   
   inpdesign_surf_stability();
/*-------------------------------- input of volume stability definition */   
   inpdesign_vol_stability();
}
#endif
#ifdef D_AXISHELL
   /* -----------------------input of line axishell thickness condition */
   inpdesign_line_thickness();
   inpdesign_line_axishellload();
   inpdesign_point_axishellcos();
#endif
#ifdef D_LS
if (genprob.probtyp==prb_twophase)
{
/*------------------------------- input of surface stability definition */   
   inpdesign_surf_stability();
}
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_conditions */








/*----------------------------------------------------------------------*
 | input of nodal dirichlet conditions on design          m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_nodal_dirich()
{
INT    i;
INT    ierr;
INT    ndnode;
INT    dnodeId;
char  *colptr;
char   buffer[200];
DNODE *actdnode;
#ifdef DEBUG 
dstrc_enter("inpdesign_nodal_dirich");
#endif
/*----------------------------------------------------------------------*/
/*-------------------- find the beginning of nodal dirichlet conditions */
if (frfind("--DESIGN POINT DIRICH CONDITIONS")==0) goto end;
frread();
/*------------------------ read number of design points with conditions */
frint("DPOINT",&ndnode,&ierr);
dsassert(ierr==1,"Cannot read design-nodal dirichlet conditions");
frread();
/*-------------------------------------- start reading the design nodes */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design node Id */
   frint("E",&dnodeId,&ierr);
   dsassert(ierr==1,"Cannot read design-nodal dirichlet conditions");
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
   dsassert(actdnode!=NULL,"Cannot read design-nodal dirichlet conditions");
   /*----------- allocate space for a dirichlet condition in this dnode */
   actdnode->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
   if (!actdnode->dirich) dserror("Allocation of memory failed");
   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-nodal dirichlet conditions");
   colptr++;
   /*--- now read the 6 flags for the dirichlet conditions */ 
   /*---------------------------- and the 6 values of them */

   amdef("onoff",&(actdnode->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
   amdef("val",&(actdnode->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
   amdef("curve",&(actdnode->dirich->curve),MAXDOFPERNODE,1,"IV");
   amdef("function",&(actdnode->dirich->funct),MAXDOFPERNODE,1,"IV");

   /* Initialize some arrays to ZERO if MAXDOFPERNODE != 6 */
   amzero(&(actdnode->dirich->dirich_onoff));
   amzero(&(actdnode->dirich->dirich_val));
   amzero(&(actdnode->dirich->curve));
   amzero(&(actdnode->dirich->funct));

   /* NOTE: number of read values = 6  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdnode->dirich->dirich_onoff.a.iv[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdnode->dirich->dirich_val.a.dv[i] = strtod(colptr,&colptr);
     else
       strtod(colptr,&colptr);
   /*---------------------------------- read the curvenumbers or "none" */
   for (i=0; i<6; i++)
   {
     ierr=sscanf(colptr," %s ",buffer);   
     dsassert(ierr==1,"Cannot read design-nodal dirichlet conditions");
     if (strncmp(buffer,"none",4)==0) 
     {
       if (i < MAXDOFPERNODE)
         actdnode->dirich->curve.a.iv[i] = 0;
       colptr = strstr(colptr,"none");
       dsassert(colptr!=NULL,"Cannot read design-nodal dirichlet conditions");
       colptr += 4;
     }
     else
     {
       if (i < MAXDOFPERNODE)
         ierr=sscanf(colptr," %d ",&(actdnode->dirich->curve.a.iv[i]));
       dsassert(ierr==1,"Cannot read design-nodal dirichlet conditions");
       colptr = strpbrk(colptr,"1234567890");
       colptr++;
     }   
   } 
   /* read function number */
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdnode->dirich->funct.a.iv[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_nodal_dirich */


/*----------------------------------------------------------------------*
 | input of line dirichlet conditions on design           m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_line_dirich()
{
INT    i;
INT    ierr;
INT    ndline;
INT    dlineId;
char  *colptr;
char   buffer[200];
DLINE *actdline;
#ifdef DEBUG 
dstrc_enter("inpdesign_line_dirich");
#endif
/*----------------------------------------------------------------------*/
/*-------------------- find the beginning of line dirichlet conditions */
if (frfind("--DESIGN LINE DIRICH CONDITIONS")==0) goto end;
frread();
/*------------------------- read number of design lines with conditions */
frint("DLINE",&ndline,&ierr);
dsassert(ierr==1,"Cannot read design-line dirichlet conditions");
frread();
/*-------------------------------------- start reading the design lines */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design line Id */
   frint("E",&dlineId,&ierr);
   dsassert(ierr==1,"Cannot read design-line dirichlet conditions");
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
   dsassert(actdline!=NULL,"Cannot read design-line dirichlet conditions");
   /*----------- allocate space for a dirichlet condition in this dline */
   actdline->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
   if (!actdline->dirich) dserror("Allocation of memory failed");
   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-line dirichlet conditions");
   colptr++;
   /*----------------------------------- read the curvenumber or "none" */
   /*--- now read the 6 flags for the dirichlet conditions */ 
   /*---------------------------- and the 6 values of them */

   amdef("onoff",&(actdline->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
   amdef("val",&(actdline->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
   amdef("curve",&(actdline->dirich->curve),MAXDOFPERNODE,1,"IV");
   amdef("function",&(actdline->dirich->funct),MAXDOFPERNODE,1,"IV");

   /* Initialize some arrays to ZERO if MAXDOFPERNODE != 6 */
   amzero(&(actdline->dirich->dirich_onoff));
   amzero(&(actdline->dirich->dirich_val));
   amzero(&(actdline->dirich->curve));
   amzero(&(actdline->dirich->funct));

   /* NOTE: number of read values = 6  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdline->dirich->dirich_onoff.a.iv[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdline->dirich->dirich_val.a.dv[i] = strtod(colptr,&colptr);
     else
       strtod(colptr,&colptr);
   /*---------------------------------- read the curvenumbers or "none" */
   for (i=0; i<6; i++)
   {
     ierr=sscanf(colptr," %s ",buffer);
     dsassert(ierr==1,"Cannot read design-line dirichlet conditions");
     if (strncmp(buffer,"none",4)==0) 
     {
       if (i < MAXDOFPERNODE)
         actdline->dirich->curve.a.iv[i] = 0;
       colptr = strstr(colptr,"none");
       dsassert(colptr!=NULL,"Cannot read design-line dirichlet conditions");
       colptr += 4;
     }
     else
     {
       if (i < MAXDOFPERNODE)
         ierr=sscanf(colptr," %d ",&(actdline->dirich->curve.a.iv[i]));
       dsassert(ierr==1,"Cannot read design-line dirichlet conditions");
       colptr = strpbrk(colptr,"1234567890");
       colptr++;
     }
   }   
   /* read function number */
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdline->dirich->funct.a.iv[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_line_dirich */


/*----------------------------------------------------------------------*
 | input of surf dirichlet conditions on design           m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_surf_dirich()
{
INT    i;
INT    ierr;
INT    ndsurf;
INT    dsurfId;
char  *colptr;
char   buffer[200];
DSURF *actdsurf;
#ifdef DEBUG 
dstrc_enter("inpdesign_surf_dirich");
#endif
/*----------------------------------------------------------------------*/
/*------------------ find the beginning of surface dirichlet conditions */
if (frfind("--DESIGN SURF DIRICH CONDITIONS")==0) goto end;
frread();
/*------------------------- read number of design surfs with conditions */
frint("DSURF",&ndsurf,&ierr);
dsassert(ierr==1,"Cannot read design-surf dirichlet conditions");
frread();
/*-------------------------------------- start reading the design surfs */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design surf Id */
   frint("E",&dsurfId,&ierr);
   dsassert(ierr==1,"Cannot read design-surf dirichlet conditions");
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
   dsassert(actdsurf!=NULL,"Cannot read design-surf dirichlet conditions");
   /*----------- allocate space for a dirichlet condition in this dsurf */
   actdsurf->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
   if (!actdsurf->dirich) dserror("Allocation of memory failed");
   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-surf dirichlet conditions");
   colptr++;
   /*--- now read the 6 flags for the dirichlet conditions */ 
   /*---------------------------- and the 6 values of them */

   amdef("onoff",&(actdsurf->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
   amdef("val",&(actdsurf->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
   amdef("curve",&(actdsurf->dirich->curve),MAXDOFPERNODE,1,"IV");
   amdef("function",&(actdsurf->dirich->funct),MAXDOFPERNODE,1,"IV");

   /* Initialize some arrays to ZERO if MAXDOFPERNODE != 6 */
   amzero(&(actdsurf->dirich->dirich_onoff));
   amzero(&(actdsurf->dirich->dirich_val));
   amzero(&(actdsurf->dirich->curve));
   amzero(&(actdsurf->dirich->funct));

   /* NOTE: number of read values = 6  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdsurf->dirich->dirich_onoff.a.iv[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdsurf->dirich->dirich_val.a.dv[i] = strtod(colptr,&colptr);   
     else
       strtod(colptr,&colptr);   
   /*----------------------------------- read the curvenumber or "none" */
   for (i=0; i<6; i++)
   {
     ierr=sscanf(colptr," %s ",buffer);
     dsassert(ierr==1,"Cannot read design-surf dirichlet conditions");
     if (strncmp(buffer,"none",4)==0) 
     {
       if (i < MAXDOFPERNODE)
         actdsurf->dirich->curve.a.iv[i] = 0;
       colptr = strstr(colptr,"none");
       dsassert(colptr!=NULL,"Cannot read design-surf dirichlet conditions");
       colptr += 4;
     }
     else
     {
       if (i < MAXDOFPERNODE)
         ierr=sscanf(colptr," %d ",&(actdsurf->dirich->curve.a.iv[i]));
       dsassert(ierr==1,"Cannot read design-surf dirichlet conditions");
       colptr = strpbrk(colptr,"1234567890");
       colptr++;
     }
   } 
   /* read function number */
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdsurf->dirich->funct.a.iv[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_surf_dirich */


/*----------------------------------------------------------------------*
 | input of vol  dirichlet conditions on design           m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_vol_dirich()
{
INT    i;
INT    ierr;
INT    ndvol;
INT    dvolId;
char  *colptr;
char   buffer[200];
DVOL  *actdvol;
#ifdef DEBUG 
dstrc_enter("inpdesign_vol_dirich");
#endif
/*----------------------------------------------------------------------*/
/*------------------ find the beginning of surface dirichlet conditions */
if (frfind("--DESIGN VOL DIRICH CONDITIONS")==0) goto end;
frread();
/*----------------------- read number of design volumes with conditions */
frint("DVOL",&ndvol,&ierr);
dsassert(ierr==1,"Cannot read design-vol dirichlet conditions");
frread();
/*--------------------------------------- start reading the design vols */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------- read the design vol Id */
   frint("E",&dvolId,&ierr);
   dsassert(ierr==1,"Cannot read design-vol dirichlet conditions");
   dvolId--;
   /*---------------------------------------------------- find the dvol */
   actdvol=NULL;
   for (i=0; i<design->ndvol; i++)
   {
      if (design->dvol[i].Id ==  dvolId) 
      {
         actdvol = &(design->dvol[i]);
         break;
      }
   }
   dsassert(actdvol!=NULL,"Cannot read design-vol dirichlet conditions");
   /*------------ allocate space for a dirichlet condition in this dvol */
   actdvol->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
   if (!actdvol->dirich) dserror("Allocation of memory failed");
   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-vol dirichlet conditions");
   colptr++;
   /*--- now read the 6 flags for the dirichlet conditions */ 
   /*---------------------------- and the 6 values of them */

   amdef("onoff",&(actdvol->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
   amdef("val",&(actdvol->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
   amdef("curve",&(actdvol->dirich->curve),MAXDOFPERNODE,1,"IV"); 
   amdef("function",&(actdvol->dirich->funct),MAXDOFPERNODE,1,"IV");

   /* Initialize some arrays to ZERO if MAXDOFPERNODE != 6 */
   amzero(&(actdvol->dirich->dirich_onoff));
   amzero(&(actdvol->dirich->dirich_val));
   amzero(&(actdvol->dirich->curve));
   amzero(&(actdvol->dirich->funct));

   /* NOTE: number of read values = 6  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdvol->dirich->dirich_onoff.a.iv[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdvol->dirich->dirich_val.a.dv[i] = strtod(colptr,&colptr);   
     else
       strtod(colptr,&colptr);   
   /*----------------------------------- read the curvenumber or "none" */
   for (i=0; i<6; i++)
   {
     ierr=sscanf(colptr," %s ",buffer);
     dsassert(ierr==1,"Cannot read design-vol dirichlet conditions");
     if (strncmp(buffer,"none",4)==0) 
     {
       if (i < MAXDOFPERNODE)
         actdvol->dirich->curve.a.iv[i] = 0;
       colptr = strstr(colptr,"none");
       dsassert(colptr!=NULL,"Cannot read design-vol dirichlet conditions");
       colptr += 4;
     }
     else
     {
       if (i < MAXDOFPERNODE)
         ierr=sscanf(colptr," %d ",&(actdvol->dirich->curve.a.iv[i]));
       dsassert(ierr==1,"Cannot read design-vol dirichlet conditions");
       colptr = strpbrk(colptr,"1234567890");
       colptr++;
     }  
   }   
   /* read function number */
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdvol->dirich->funct.a.iv[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_vol_dirich */




/*----------------------------------------------------------------------*
 | input of nodal neumann   conditions on design          m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_nodal_neum()
{
INT    i;
INT    ierr;
INT    ndnode;
INT    dnodeId;
char  *colptr;
char   buffer[200];
DNODE *actdnode;
#ifdef DEBUG 
dstrc_enter("inpdesign_nodal_neum");
#endif
/*----------------------------------------------------------------------*/
/*-------------------- find the beginning of nodal neumann conditions */
if (frfind("--DESIGN POINT NEUMANN CONDITIONS")==0) goto end;
frread();
/*------------------------ read number of design points with conditions */
frint("DPOINT",&ndnode,&ierr);
dsassert(ierr==1,"Cannot read design-nodal neumann conditions");
frread();
/*-------------------------------------- start reading the design nodes */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design node Id */
   frint("E",&dnodeId,&ierr);
   dsassert(ierr==1,"Cannot read design-nodal neumann conditions");
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
   dsassert(actdnode!=NULL,"Cannot read design-nodal neumann conditions");
   /*----------- allocate space for a neumann condition in this dnode */
   actdnode->neum = (NEUM_CONDITION*)CCACALLOC(1,sizeof(NEUM_CONDITION));
   if (!actdnode->neum) dserror("Allocation of memory failed");
   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-nodal neumann conditions");
   colptr++;
   /*----------------------------------- read the curvenumber or "none" */
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot read design-nodal neumann conditions");
   if (strncmp(buffer,"none",4)==0) 
   {
      actdnode->neum->curve = 0;
      colptr = strstr(allfiles.actplace,"none");
      dsassert(colptr!=NULL,"Cannot read design-nodal neumann conditions");
      colptr += 4;
   }
   else
   {
      ierr=sscanf(colptr," %d ",&(actdnode->neum->curve));
      dsassert(ierr==1,"Cannot read design-nodal neumann conditions");
      colptr = strpbrk(colptr,"1234567890");
      colptr++;
   }
   /*--- now read the 6 flags for the neumann conditions */ 
   /*-------------------------- and the 6 values of them */
   amdef("onoff",&(actdnode->neum->neum_onoff),MAXDOFPERNODE,1,"IV");
   amdef("val",&(actdnode->neum->neum_val),MAXDOFPERNODE,1,"DV");

   /* Initialize some arrays to ZERO if MAXDOFPERNODE != 6 */
   amzero(&(actdnode->neum->neum_onoff));
   amzero(&(actdnode->neum->neum_val));

   /* NOTE: number of read values = 6  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdnode->neum->neum_onoff.a.iv[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdnode->neum->neum_val.a.dv[i] = strtod(colptr,&colptr);
     else
       strtod(colptr,&colptr);
   /*----------- read if load is applied on surface -> shell elements */
   frchk("Mid",&ierr);
   if (ierr) actdnode->neum->neum_surf = mid;
   frchk("Top",&ierr);
   if (ierr) actdnode->neum->neum_surf = top;
   frchk("Bot",&ierr);   
   if (ierr) actdnode->neum->neum_surf = bot;
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_nodal_neum */


/*----------------------------------------------------------------------*
 | input of line neumann conditions on design           m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_line_neum()
{
INT    i;
INT    ierr;
INT    ndline;
INT    dlineId;
char  *colptr;
char   buffer[200];
DLINE *actdline;
#ifdef DEBUG 
dstrc_enter("inpdesign_line_neum");
#endif
/*----------------------------------------------------------------------*/
/*-------------------- find the beginning of line neumann conditions */
if (frfind("--DESIGN LINE NEUMANN CONDITIONS")==0) goto end;
frread();
/*------------------------- read number of design lines with conditions */
frint("DLINE",&ndline,&ierr);
dsassert(ierr==1,"Cannot read design-line neumann conditions");
frread();
/*-------------------------------------- start reading the design lines */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design line Id */
   frint("E",&dlineId,&ierr);
   dsassert(ierr==1,"Cannot read design-line neumann conditions");
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
   dsassert(actdline!=NULL,"Cannot read design-line neumann conditions");
   /*----------- allocate space for a neumann condition in this dline */
   actdline->neum = (NEUM_CONDITION*)CCACALLOC(1,sizeof(NEUM_CONDITION));
   if (!actdline->neum) dserror("Allocation of memory failed");
   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-line neumann conditions");
   colptr++;
   /*----------------------------------- read the curvenumber or "none" */
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot read design-line neumann conditions");
   if (strncmp(buffer,"none",4)==0) 
   {
      actdline->neum->curve = 0;
      colptr = strstr(allfiles.actplace,"none");
      dsassert(colptr!=NULL,"Cannot read design-line neumann conditions");
      colptr += 4;
   }
   else
   {
      ierr=sscanf(colptr," %d ",&(actdline->neum->curve));
      dsassert(ierr==1,"Cannot read design-line neumann conditions");
      colptr = strpbrk(colptr,"1234567890");
      colptr++;
   }
   /*--- now read the 6 flags for the neumann conditions */ 
   /*-------------------------- and the 6 values of them */
   amdef("onoff",&(actdline->neum->neum_onoff),MAXDOFPERNODE,1,"IV");
   amdef("val",&(actdline->neum->neum_val),MAXDOFPERNODE,1,"DV");

   /* Initialize some arrays to ZERO if MAXDOFPERNODE != 6  */
   amzero(&(actdline->neum->neum_onoff));
   amzero(&(actdline->neum->neum_val));

   /* NOTE: number of read values = 6  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdline->neum->neum_onoff.a.iv[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdline->neum->neum_val.a.dv[i] = strtod(colptr,&colptr);
     else
       strtod(colptr,&colptr);
   /*----------- read if load is applied on surface -> shell elements */
   frchk("Mid",&ierr);
   if (ierr) actdline->neum->neum_surf = mid;
   frchk("Top",&ierr);
   if (ierr) actdline->neum->neum_surf = top;
   frchk("Bot",&ierr);   
   if (ierr) actdline->neum->neum_surf = bot;
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_line_neum */


/*----------------------------------------------------------------------*
 | input of surf neumann conditions on design           m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_surf_neum()
{
INT    i;
INT    ierr;
INT    ndsurf;
INT    dsurfId;
char  *colptr;
char   buffer[200];
DSURF *actdsurf;
#ifdef DEBUG 
dstrc_enter("inpdesign_surf_neum");
#endif
/*----------------------------------------------------------------------*/
/*------------------ find the beginning of surface neumann conditions */
if (frfind("--DESIGN SURF NEUMANN CONDITIONS")==0) goto end;
frread();
/*------------------------- read number of design surfs with conditions */
frint("DSURF",&ndsurf,&ierr);
dsassert(ierr==1,"Cannot read design-surf neumann conditions");
frread();
/*-------------------------------------- start reading the design surfs */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design surf Id */
   frint("E",&dsurfId,&ierr);
   dsassert(ierr==1,"Cannot read design-surf neumann conditions");
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
   dsassert(actdsurf!=NULL,"Cannot read design-surf neumann conditions");
   /*----------- allocate space for a neumann condition in this dsurf */
   actdsurf->neum = (NEUM_CONDITION*)CCACALLOC(1,sizeof(NEUM_CONDITION));
   if (!actdsurf->neum) dserror("Allocation of memory failed");
   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-surf neumann conditions");
   colptr++;
   /*----------------------------------- read the curvenumber or "none" */
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot read design-surf neumann conditions");
   if (strncmp(buffer,"none",4)==0) 
   {
      actdsurf->neum->curve = 0;
      colptr = strstr(allfiles.actplace,"none");
      dsassert(colptr!=NULL,"Cannot read design-surf neumann conditions");
      colptr += 4;
   }
   else
   {
      ierr=sscanf(colptr," %d ",&(actdsurf->neum->curve));
      dsassert(ierr==1,"Cannot read design-surf neumann conditions");
      colptr = strpbrk(colptr,"1234567890");
      colptr++;
   }
   /*--- now read the 6 flags for the neumann conditions */ 
   /*-------------------------- and the 6 values of them */
   amdef("onoff",&(actdsurf->neum->neum_onoff),MAXDOFPERNODE,1,"IV");
   amdef("val",&(actdsurf->neum->neum_val),MAXDOFPERNODE,1,"DV");

   /* Initialize some arrays to ZERO if MAXDOFPERNODE != 6  */
   amzero(&(actdsurf->neum->neum_onoff));
   amzero(&(actdsurf->neum->neum_val));

   /* NOTE: number of read values = 6  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdsurf->neum->neum_onoff.a.iv[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdsurf->neum->neum_val.a.dv[i] = strtod(colptr,&colptr);
     else
       strtod(colptr,&colptr);
   /*----------- read if load is applied on surface -> shell elements */
   frchk("Mid",&ierr);
   if (ierr) actdsurf->neum->neum_surf = mid;
   frchk("Top",&ierr);
   if (ierr) actdsurf->neum->neum_surf = top;
   frchk("Bot",&ierr);   
   if (ierr) actdsurf->neum->neum_surf = bot;
   /*----------------------------------- read type of neumann condition */
   frchk("Live",&ierr);
   if (ierr) actdsurf->neum->neum_type = neum_live;
   frchk("Dead",&ierr);
   if (ierr) actdsurf->neum->neum_type = neum_dead;
   frchk("PrescribedDomainLoad",&ierr);   
   if (ierr) actdsurf->neum->neum_type = pres_domain_load;
   frchk("constHydro_z",&ierr);
   if (ierr) actdsurf->neum->neum_type = neum_consthydro_z;
   frchk("increaseHydro_z",&ierr);
   if (ierr) actdsurf->neum->neum_type = neum_increhydro_z;
   frchk("orthopressure",&ierr);
   if (ierr) actdsurf->neum->neum_type = neum_orthopressure;
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_surf_neum */


/*----------------------------------------------------------------------*
 | input of vol  neumann conditions on design           m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_vol_neum()
{
INT    i;
INT    ierr;
INT    ndvol;
INT    dvolId;
char  *colptr;
char   buffer[200];
DVOL  *actdvol;
#ifdef DEBUG 
dstrc_enter("inpdesign_vol_neum");
#endif
/*----------------------------------------------------------------------*/
/*------------------ find the beginning of surface neumann conditions */
if (frfind("--DESIGN VOL NEUMANN CONDITIONS")==0) goto end;
frread();
/*----------------------- read number of design volumes with conditions */
frint("DVOL",&ndvol,&ierr);
dsassert(ierr==1,"Cannot read design-vol neumann conditions");
frread();
/*--------------------------------------- start reading the design vols */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------- read the design vol Id */
   frint("E",&dvolId,&ierr);
   dsassert(ierr==1,"Cannot read design-vol neumann conditions");
   dvolId--;
   /*---------------------------------------------------- find the dvol */
   actdvol=NULL;
   for (i=0; i<design->ndvol; i++)
   {
      if (design->dvol[i].Id ==  dvolId) 
      {
         actdvol = &(design->dvol[i]);
         break;
      }
   }
   dsassert(actdvol!=NULL,"Cannot read design-vol neumann conditions");
   /*------------ allocate space for a neumann condition in this dvol */
   actdvol->neum = (NEUM_CONDITION*)CCACALLOC(1,sizeof(NEUM_CONDITION));
   if (!actdvol->neum) dserror("Allocation of memory failed");
   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-vol neumann conditions");
   colptr++;
   /*----------------------------------- read the curvenumber or "none" */
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot read design-vol neumann conditions");
   if (strncmp(buffer,"none",4)==0) 
   {
      actdvol->neum->curve = 0;
      colptr = strstr(allfiles.actplace,"none");
      dsassert(colptr!=NULL,"Cannot read design-vol neumann conditions");
      colptr += 4;
   }
   else
   {
      ierr=sscanf(colptr," %d ",&(actdvol->neum->curve));
      dsassert(ierr==1,"Cannot read design-vol neumann conditions");
      colptr = strpbrk(colptr,"1234567890");
      colptr++;
   }
   /*--- now read the 6 flags for the neumann conditions */ 
   /*-------------------------- and the 6 values of them */
   amdef("onoff",&(actdvol->neum->neum_onoff),MAXDOFPERNODE,1,"IV");
   amdef("val",&(actdvol->neum->neum_val),MAXDOFPERNODE,1,"DV");

   /* Initialize some arrays to ZERO if MAXDOFPERNODE != 6 */
   amzero(&(actdvol->neum->neum_onoff));
   amzero(&(actdvol->neum->neum_val));

   /* NOTE: number of read values = 6  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdvol->neum->neum_onoff.a.iv[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       actdvol->neum->neum_val.a.dv[i] = strtod(colptr,&colptr);
     else
       strtod(colptr,&colptr);
   /*----------------------------------- read type of neumann condition */
   frchk("Dead",&ierr);
   if (ierr) actdvol->neum->neum_type = neum_dead;   
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_vol_neum */













/*----------------------------------------------------------------------*
 | input of nodal coupling conditions on design           m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_nodal_couple()
{
INT    i;
INT    ierr;
INT    ndnode;
INT    dnodeId;
char  *colptr;
char   buffer[200];
DNODE *actdnode;
INT    coupleId;
INT    dofcouple;
INT    dofflags[MAXDOFPERNODE];
INT    geocouple;
INT    geoflags[MAXDOFPERNODE];
#ifdef DEBUG 
dstrc_enter("inpdesign_nodal_couple");
#endif
/*----------------------------------------------------------------------*/
/*-------------------- find the beginning of nodal coupling conditions */
if (frfind("--DESIGN COUPLING POINT CONDITIONS")==0) goto end;
frread();
/*------------------------ read number of design points with conditions */
frint("DPOINT",&ndnode,&ierr);
dsassert(ierr==1,"Cannot read design-nodal coupling conditions");
frread();
/*-------------------------------------- start reading the design nodes */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design node Id */
   frint("E",&dnodeId,&ierr);
   dsassert(ierr==1,"Cannot read design-nodal coupling conditions");
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
   dsassert(actdnode!=NULL,"Cannot read design-nodal coupling conditions");
   /*----------- allocate space for a coupling condition in this dnode */
   actdnode->couple = (COUPLE_CONDITION*)CCACALLOC(1,sizeof(COUPLE_CONDITION));
   if (!actdnode->couple) dserror("Allocation of memory failed");
   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-nodal coupling conditions");
   colptr++;
   /*----------------------------------------- read the coupling set Id */
   coupleId=-1;
   coupleId = strtol(colptr,&colptr,10);
   dsassert(coupleId>0,"Parameter out of range: CoupleId<=0\n");
   /*------------------------------------------------ read the fieldtyp */
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot reading coupling conditions");
   if (strncmp(buffer,"structure",9)==0) 
   {
       actdnode->couple->fieldtyp=structure;
       colptr = strstr(colptr,"structure");
       dsassert(colptr!=NULL,"Cannot reading coupling conditions");
       colptr+=9;
   }
   if (strncmp(buffer,"fluid",5)==0) 
   {
      actdnode->couple->fieldtyp=fluid;
       colptr = strstr(colptr,"fluid");
       dsassert(colptr!=NULL,"Cannot reading coupling conditions");
       colptr+=5;
   }
   if (strncmp(buffer,"ale",3)==0) 
   {
      actdnode->couple->fieldtyp=ale;
       colptr = strstr(colptr,"ale");
       dsassert(colptr!=NULL,"Cannot reading coupling conditions");
       colptr+=3;
   }
   /*-------------------------------------------------- allocate arrays */
   amdef("couple",&(actdnode->couple->couple),MAXDOFPERNODE,4,"IA");
   amzero(&(actdnode->couple->couple));
   /*------------------------------------------ read the dofcouple flag */
   /* dofcouple Id */
   /* NOTE: number of read values = 6  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   dofcouple = strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       dofflags[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
   {
     if (i < MAXDOFPERNODE)
       if (dofflags[i]!=0) 
         actdnode->couple->couple.a.ia[i][1]=dofcouple;
   }
   /*------------------------------------------ read the geocouple flag */
   geocouple = strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       geoflags[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
   {
     if (i < MAXDOFPERNODE)
       if (geoflags[i]!=0) 
         actdnode->couple->couple.a.ia[i][0]=geocouple;
   }
   /*------------------------------------------------ read the fsi flag */   
/*   actdnode->couple->fsi_iscoupled = strtol(colptr,&colptr,10);*/
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_nodal_couple */


/*----------------------------------------------------------------------*
 | input of line coupling conditions on design            m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_line_couple()
{
INT    i;
INT    ierr;
INT    ndline;
INT    dlineId;
char  *colptr;
char   buffer[200];
DLINE *actdline;
INT    coupleId;
INT    dofcouple;
INT    dofflags[MAXDOFPERNODE];
INT    geocouple;
INT    geoflags[MAXDOFPERNODE];
#ifdef DEBUG 
dstrc_enter("inpdesign_line_couple");
#endif
/*----------------------------------------------------------------------*/
/*--------------------- find the beginning of line coupling conditions */
if (frfind("--DESIGN COUPLING LINE CONDITIONS")==0) goto end;
frread();
/*------------------------ read number of design lines with conditions */
frint("DLINE",&ndline,&ierr);
dsassert(ierr==1,"Cannot read design-line coupling conditions");
frread();
/*-------------------------------------- start reading the design lines */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design line Id */
   frint("E",&dlineId,&ierr);
   dsassert(ierr==1,"Cannot read design-line coupling conditions");
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
   dsassert(actdline!=NULL,"Cannot read design-line coupling conditions");
   /*----------- allocate space for a coupling condition in this dline */
   actdline->couple = (COUPLE_CONDITION*)CCACALLOC(1,sizeof(COUPLE_CONDITION));
   if (!actdline->couple) dserror("Allocation of memory failed");
   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-line coupling conditions");
   colptr++;
   /*----------------------------------------- read the coupling set Id */
   coupleId=-1;
   coupleId = strtol(colptr,&colptr,10);
   dsassert(coupleId>0,"Cannot reading coupling conditions");
   /*------------------------------------------------ read the fieldtyp */
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot reading coupling conditions");
   if (strncmp(buffer,"structure",9)==0) 
   {
       actdline->couple->fieldtyp=structure;
       colptr = strstr(colptr,"structure");
       dsassert(colptr!=NULL,"Cannot reading coupling conditions");
       colptr+=9;
   }
   if (strncmp(buffer,"fluid",5)==0) 
   {
      actdline->couple->fieldtyp=fluid;
       colptr = strstr(colptr,"fluid");
       dsassert(colptr!=NULL,"Cannot reading coupling conditions");
       colptr+=5;
   }
   if (strncmp(buffer,"ale",3)==0) 
   {
      actdline->couple->fieldtyp=ale;
       colptr = strstr(colptr,"ale");
       dsassert(colptr!=NULL,"Cannot reading coupling conditions");
       colptr+=3;
   }
   /*-------------------------------------------------- allocate arrays */
   amdef("couple",&(actdline->couple->couple),MAXDOFPERNODE,4,"IA");
   amzero(&(actdline->couple->couple));
   /*------------------------------------------ read the dofcouple flag */
   /* dofcouple Id */
   /* NOTE: number of read values = 6  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   dofcouple = strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       dofflags[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
   {
     if (i < MAXDOFPERNODE)
       if (dofflags[i]!=0) 
         actdline->couple->couple.a.ia[i][1]=dofcouple;
   }
   /*------------------------------------------ read the geocouple flag */
   geocouple = strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       geoflags[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
   {
     if (i < MAXDOFPERNODE)
       if (geoflags[i]!=0) 
         actdline->couple->couple.a.ia[i][0]=geocouple;
   }
   /*------------------------------------------------ read the fsi flag */   
/*   actdline->couple->fsi_iscoupled = strtol(colptr,&colptr,10);*/
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_line_couple */



/*----------------------------------------------------------------------*
 | input of surface coupling conditions on design           m.gee 3/02  |
 *----------------------------------------------------------------------*/
static void inpdesign_surf_couple()
{
INT    i;
INT    ierr;
INT    ndsurf;
INT    dsurfId;
char  *colptr;
char   buffer[200];
DSURF *actdsurf;
INT    coupleId;
INT    dofcouple;
INT    dofflags[MAXDOFPERNODE];
INT    geocouple;
INT    geoflags[MAXDOFPERNODE];
#ifdef DEBUG 
dstrc_enter("inpdesign_surf_couple");
#endif
/*----------------------------------------------------------------------*/
/*-------------------- find the beginning of surf coupling conditions */
if (frfind("--DESIGN COUPLING SURF CONDITIONS")==0) goto end;
frread();
/*------------------------ read number of design points with conditions */
frint("DSURF",&ndsurf,&ierr);
dsassert(ierr==1,"Cannot read design-nodal coupling conditions");
frread();
/*-------------------------------------- start reading the design surfs */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design surf Id */
   frint("E",&dsurfId,&ierr);
   dsassert(ierr==1,"Cannot read design-nodal coupling conditions");
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
   dsassert(actdsurf!=NULL,"Cannot read design-nodal coupling conditions");
   /*----------- allocate space for a coupling condition in this dsurf */
   actdsurf->couple = (COUPLE_CONDITION*)CCACALLOC(1,sizeof(COUPLE_CONDITION));
   if (!actdsurf->couple) dserror("Allocation of memory failed");
   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-nodal coupling conditions");
   colptr++;
   /*----------------------------------------- read the coupling set Id */
   coupleId=-1;
   coupleId = strtol(colptr,&colptr,10);
   dsassert(coupleId>0,"Cannot reading coupling conditions");
   /*------------------------------------------------ read the fieldtyp */
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot reading coupling conditions");
   if (strncmp(buffer,"structure",9)==0) 
   {
       actdsurf->couple->fieldtyp=structure;
       colptr = strstr(colptr,"structure");
       dsassert(colptr!=NULL,"Cannot reading coupling conditions");
       colptr+=9;
   }
   if (strncmp(buffer,"fluid",5)==0) 
   {
      actdsurf->couple->fieldtyp=fluid;
       colptr = strstr(colptr,"fluid");
       dsassert(colptr!=NULL,"Cannot reading coupling conditions");
       colptr+=5;
   }
   if (strncmp(buffer,"ale",3)==0) 
   {
      actdsurf->couple->fieldtyp=ale;
       colptr = strstr(colptr,"ale");
       dsassert(colptr!=NULL,"Cannot reading coupling conditions");
       colptr+=3;
   }
   /*-------------------------------------------------- allocate arrays */
   amdef("couple",&(actdsurf->couple->couple),MAXDOFPERNODE,4,"IA");
   amzero(&(actdsurf->couple->couple));
   /*------------------------------------------ read the dofcouple flag */
   /* dofcouple Id */
   /* NOTE: number of read values = 6  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   dofcouple = strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       dofflags[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
   {
     if (i < MAXDOFPERNODE)
       if (dofflags[i]!=0) 
         actdsurf->couple->couple.a.ia[i][1]=dofcouple;
   }
   /*------------------------------------------ read the geocouple flag */
   geocouple = strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       geoflags[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
   {
     if (i < MAXDOFPERNODE)
       if (geoflags[i]!=0) 
         actdsurf->couple->couple.a.ia[i][0]=geocouple;
   }
   /*------------------------------------------------ read the fsi flag */   
/*   actdsurf->couple->fsi_iscoupled = strtol(colptr,&colptr,10);*/
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_surf_couple */



/*----------------------------------------------------------------------*
 | input of volume coupling conditions on design           m.gee 3/02   |
 *----------------------------------------------------------------------*/
static void inpdesign_vol_couple()
{
INT    i;
INT    ierr;
INT    ndvol;
INT    dvolId;
char  *colptr;
char   buffer[200];
DVOL *actdvol;
INT    coupleId;
INT    dofcouple;
INT    dofflags[MAXDOFPERNODE];
INT    geocouple;
INT    geoflags[MAXDOFPERNODE];
#ifdef DEBUG 
dstrc_enter("inpdesign_vol_couple");
#endif
/*----------------------------------------------------------------------*/
/*-------------------- find the beginning of volume coupling conditions */
if (frfind("--DESIGN COUPLING VOL CONDITIONS")==0) goto end;
frread();
/*------------------------ read number of design points with conditions */
frint("DVOL",&ndvol,&ierr);
dsassert(ierr==1,"Cannot read design-nodal coupling conditions");
frread();
/*-------------------------------------- start reading the design vols */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design vol Id */
   frint("E",&dvolId,&ierr);
   dsassert(ierr==1,"Cannot read design-nodal coupling conditions");
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
   dsassert(actdvol!=NULL,"Cannot read design-nodal coupling conditions");
   /*----------- allocate space for a coupling condition in this dvol */
   actdvol->couple = (COUPLE_CONDITION*)CCACALLOC(1,sizeof(COUPLE_CONDITION));
   if (!actdvol->couple) dserror("Allocation of memory failed");
   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-nodal coupling conditions");
   colptr++;
   /*----------------------------------------- read the coupling set Id */
   coupleId=-1;
   coupleId = strtol(colptr,&colptr,10);
   dsassert(coupleId>0,"Cannot reading coupling conditions");
   /*------------------------------------------------ read the fieldtyp */
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot reading coupling conditions");
   if (strncmp(buffer,"structure",9)==0) 
   {
       actdvol->couple->fieldtyp=structure;
       colptr = strstr(colptr,"structure");
       dsassert(colptr!=NULL,"Cannot reading coupling conditions");
       colptr+=9;
   }
   if (strncmp(buffer,"fluid",5)==0) 
   {
      actdvol->couple->fieldtyp=fluid;
       colptr = strstr(colptr,"fluid");
       dsassert(colptr!=NULL,"Cannot reading coupling conditions");
       colptr+=5;
   }
   if (strncmp(buffer,"ale",3)==0) 
   {
      actdvol->couple->fieldtyp=ale;
       colptr = strstr(colptr,"ale");
       dsassert(colptr!=NULL,"Cannot reading coupling conditions");
       colptr+=3;
   }
   /*-------------------------------------------------- allocate arrays */
   amdef("couple",&(actdvol->couple->couple),MAXDOFPERNODE,4,"IA");
   amzero(&(actdvol->couple->couple));
   /*------------------------------------------ read the dofcouple flag */
   /* dofcouple Id */
   /* NOTE: number of read values = 6  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   dofcouple = strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       dofflags[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
   {
     if (i < MAXDOFPERNODE)
       if (dofflags[i]!=0) 
         actdvol->couple->couple.a.ia[i][1]=dofcouple;
   }
   /*------------------------------------------ read the geocouple flag */
   geocouple = strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
     if (i < MAXDOFPERNODE)
       geoflags[i] = strtol(colptr,&colptr,10);
     else
       strtol(colptr,&colptr,10);
   for (i=0; i<6; i++)
   {
     if (i < MAXDOFPERNODE)
       if (geoflags[i]!=0) 
         actdvol->couple->couple.a.ia[i][0]=geocouple;
   }
   /*------------------------------------------------ read the fsi flag */   
/*   actdvol->couple->fsi_iscoupled = strtol(colptr,&colptr,10);*/
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_vol_couple */

#ifdef D_FSI
/*----------------------------------------------------------------------*
 | input of line  fsi coupling conditions on design          genk 10/02 |
 *----------------------------------------------------------------------*/
static void inpdesign_line_fsicouple()
{
INT    i;
INT    ierr;
INT    ndline;
INT    dlineId;
char  *colptr;
char   buffer[200];
DLINE *actdline;
INT    coupleId;

#ifdef DEBUG 
dstrc_enter("inpdesign_line_fsicouple");
#endif

/*----------------------------------------------------------------------*/
/*------------------ find the beginning of line fsi coupling conditions */
if (frfind("--DESIGN FSI COUPLING LINE CONDITIONS")==0) goto end;
frread();
/*------------------------ read number of design lines with conditions */
frint("DLINE",&ndline,&ierr);
dsassert(ierr==1,"Cannot read design-line fsi coupling conditions");
frread();
/*-------------------------------------- start reading the design lines */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design line Id */
   frint("E",&dlineId,&ierr);
   dsassert(ierr==1,"Cannot read design-line fsi coupling conditions");
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
   dsassert(actdline!=NULL,"Cannot read design-line fsi coupling conditions");
   /*----------- allocate space for a coupling condition in this dline */
   actdline->fsicouple = (FSI_COUPLE_CONDITION*)CCACALLOC(1,sizeof(FSI_COUPLE_CONDITION));
   if (!actdline->fsicouple) dserror("Allocation of memory failed");
   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-line fsi coupling conditions");
   colptr++;
   /*----------------------------------------- read the  fsi couplingId */
   coupleId=-1;
   coupleId = strtol(colptr,&colptr,10);
   dsassert(coupleId>0,"Cannot reading fsi coupling conditions");
   actdline->fsicouple->fsi_coupleId=coupleId;
   /*------------------------------------------------ read the fieldtyp */
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot readingfsi  coupling conditions");
   if (strncmp(buffer,"structure",9)==0) 
   {
       actdline->fsicouple->fieldtyp=structure;
       colptr = strstr(colptr,"structure");
       dsassert(colptr!=NULL,"Cannot reading fsi coupling conditions");
       colptr+=9;
   }
   if (strncmp(buffer,"fluid",5)==0) 
   {
      actdline->fsicouple->fieldtyp=fluid;
       colptr = strstr(colptr,"fluid");
       dsassert(colptr!=NULL,"Cannot reading fsi coupling conditions");
       colptr+=5;
   }
   if (strncmp(buffer,"ale",3)==0) 
   {
      actdline->fsicouple->fieldtyp=ale;
       colptr = strstr(colptr,"ale");
       dsassert(colptr!=NULL,"Cannot reading fsi coupling conditions");
       colptr+=3;
   }
   /*-------------------------------------------------- read the mesh */
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot reading fsi coupling conditions");
   if (strncmp(buffer,"conforming",10)==0) 
   {
       actdline->fsicouple->fsi_mesh=conforming;
       colptr = strstr(colptr,"conforming");
       dsassert(colptr!=NULL,"Cannot reading coupling conditions");
       colptr+=10;
   }
   if (strncmp(buffer,"non-conforming",14)==0) 
   {
       actdline->fsicouple->fsi_mesh=non_conforming;
       colptr = strstr(colptr,"non-conforming");
       dsassert(colptr!=NULL,"Cannot reading fsi coupling conditions");
       colptr+=14;
   }
   /*----------------------------------------------------- read the typ */
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot reading fsi coupling conditions");
   if (strncmp(buffer,"fsi_real",10)==0) 
   {
       actdline->fsicouple->fsi_typ=fsi_real;
       colptr = strstr(colptr,"fsi_real");
       dsassert(colptr!=NULL,"Cannot reading coupling conditions");
       colptr+=8;
   }
   if (strncmp(buffer,"fsi_pseudo",14)==0) 
   {
       actdline->fsicouple->fsi_typ=fsi_pseudo;
       colptr = strstr(colptr,"fsi_pseudo");
       dsassert(colptr!=NULL,"Cannot reading fsi coupling conditions");
       colptr+=10;
   }

   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_line_couple */
#endif


#ifdef D_FLUID
/*----------------------------------------------------------------------*
 | input of nodal fluid freesurface conditions on design     genk 01/03 |
 *----------------------------------------------------------------------*/
static void inpdesign_nodal_freesurf()
{
INT    i;
INT    ierr;
INT    ndnode;
INT    dnodeId;
char  *colptr;
char   buffer[200];
DNODE *actdnode;

#ifdef DEBUG 
dstrc_enter("inpdesign_nodal_freesurf");
#endif

/*----------------------------------------------------------------------*/
/*---- find the beginning of line fluid freesurface coupling conditions */
if (frfind("--DESIGN FLUID FREE SURFACE POINT CONDITIONS")==0) goto end;
frread();
/*------------------------ read number of design lines with conditions */
frint("DPOINT",&ndnode,&ierr);
dsassert(ierr==1,"Cannot read design-line fsi coupling conditions");
frread();
/*-------------------------------------- start reading the design lines */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design line Id */
   frint("E",&dnodeId,&ierr);
   dsassert(ierr==1,"Cannot read design-point fsi coupling conditions");
   dnodeId--;
   /*--------------------------------------------------- find the dline */
   actdnode=NULL;
   for (i=0; i<design->ndnode; i++)
   {
      if (design->dnode[i].Id ==  dnodeId) 
      {
         actdnode = &(design->dnode[i]);
         break;
      }
   }
   dsassert(actdnode!=NULL,"Cannot read design-point freesurface conditions");
   /*----------- allocate space for a coupling condition in this dline */
   actdnode->freesurf = (FLUID_FREESURF_CONDITION*)CCACALLOC(1,sizeof(FLUID_FREESURF_CONDITION));
   if (!actdnode->freesurf) dserror("Allocation of memory failed");

   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-line freesurf conditions");
   colptr++;
   /*------------------------------------------------ read the fieldtyp */
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot reading freesurface  coupling conditions");
   if (strncmp(buffer,"fluid",5)==0) 
   {
       actdnode->freesurf->fieldtyp=fluid;
       colptr = strstr(colptr,"fluid");
       dsassert(colptr!=NULL,"Cannot reading freesurface coupling conditions");
       colptr+=5;
   }
   if (strncmp(buffer,"ale",3)==0) 
   {
       actdnode->freesurf->fieldtyp=ale;
       colptr = strstr(colptr,"ale");
       dsassert(colptr!=NULL,"Cannot reading freesurface coupling conditions");
       colptr+=3;
   }
   /*---- now read the MAXDOFPERNODE flags for the local slippage conditions */ 
   amdef("fixed",&(actdnode->freesurf->fixed_onoff),MAXDOFPERNODE,1,"IV");
   for (i=0; i<MAXDOFPERNODE; i++)
   actdnode->freesurf->fixed_onoff.a.iv[i] = strtol(colptr,&colptr,10);

   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_line_freesurf */



/*----------------------------------------------------------------------*
 | input of line  fluid freesurface conditions on design     genk 01/03 |
 *----------------------------------------------------------------------*/
static void inpdesign_line_freesurf()
{
INT    i;
INT    ierr;
INT    ndline;
INT    dlineId;
char  *colptr;
char   buffer[200];
DLINE *actdline;

#ifdef DEBUG 
dstrc_enter("inpdesign_line_freesurf");
#endif

/*----------------------------------------------------------------------*/
/*---- find the beginning of line fluid freesurface coupling conditions */
if (frfind("--DESIGN FLUID FREE SURFACE LINE CONDITIONS")==0) goto end;
frread();
/*------------------------ read number of design lines with conditions */
frint("DLINE",&ndline,&ierr);
dsassert(ierr==1,"Cannot read design-line freesurface conditions");
frread();
/*-------------------------------------- start reading the design lines */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design line Id */
   frint("E",&dlineId,&ierr);
   dsassert(ierr==1,"Cannot read design-line freesurface conditions");
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
   dsassert(actdline!=NULL,"Cannot read design-line freesurface conditions");
   /*----------- allocate space for a coupling condition in this dline */
   actdline->freesurf = (FLUID_FREESURF_CONDITION*)CCACALLOC(1,sizeof(FLUID_FREESURF_CONDITION));
   if (!actdline->freesurf) dserror("Allocation of memory failed");

   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-line freesurf conditions");
   colptr++;
   /*------------------------------------------------ read the fieldtyp */
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot read freesurface conditions");
   if (strncmp(buffer,"fluid",5)==0) 
   {
       actdline->freesurf->fieldtyp=fluid;
       colptr = strstr(colptr,"fluid");
       dsassert(colptr!=NULL,"Cannot reading free surface conditions");
       colptr+=5;
   }
   if (strncmp(buffer,"ale",3)==0) 
   {
       actdline->freesurf->fieldtyp=ale;
       colptr = strstr(colptr,"ale");
       dsassert(colptr!=NULL,"Cannot reading free surface conditions");
       colptr+=3;
   }
  
   /*---- now read the MAXDOFPERNODE flags for the local slippage conditions */ 
   amdef("fixed",&(actdline->freesurf->fixed_onoff),MAXDOFPERNODE,1,"IV");
   for (i=0; i<MAXDOFPERNODE; i++)
   actdline->freesurf->fixed_onoff.a.iv[i] = strtol(colptr,&colptr,10);
   
   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_line_freesurf */



/*!----------------------------------------------------------------------
\brief input of line lift&drag definition on design

<pre>                                                            mn 03/04
This routine reads an line lift&drag definition on design.
</pre>

\warning There is nothing special to this routine
\return void                                               
\sa

*----------------------------------------------------------------------*/
static void inpdesign_line_liftdrag()
{
  INT    i;
  INT    ierr;
  INT    ndline;
  INT    dlineId;
  char  *colptr;
  DLINE *actdline;

#ifdef DEBUG 
  dstrc_enter("inpdesign_line_liftdrag");
#endif

  /* find the beginning of line fluid liftdrag conditions */
  if (frfind("--DESIGN FLUID LINE LIFT&DRAG") == 0) goto end;
  frread();

  /* read number of design lines with conditions */
  frint("DLINE",&ndline,&ierr);
  dsassert(ierr==1,"Cannot read design-line lift&drag definition");
  frread();


  /* start reading the design lines */
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /* read the design line Id */
    frint("E",&dlineId,&ierr);
    dsassert(ierr==1,"Cannot read design-line left&drag definition");
    dlineId--;

    /* find the dline */
    actdline=NULL;
    for (i=0; i<design->ndline; i++)
    {
      if (design->dline[i].Id ==  dlineId) 
      {
        actdline = &(design->dline[i]);
        break;
      }
    }
    dsassert(actdline!=NULL,"Cannot read designline-lift&drag definition");

    /* allocate space for a liftdrag condition in this dline */
    actdline->liftdrag = (FLUID_LIFTDRAG_CONDITION*)CCACALLOC(1,sizeof(FLUID_LIFTDRAG_CONDITION));
    if (!actdline->liftdrag) dserror("Allocation of memory failed");

    /* move pointer behind the "-" sign */
    colptr = strstr(allfiles.actplace,"-");
    dsassert(colptr!=NULL,"Cannot read design-line left&drag definition");
    colptr++;

    /* read the number */
    actdline->liftdrag->liftdrag = strtol(colptr,&colptr,10);
    if (actdline->liftdrag->liftdrag >= FLUID_NUM_LD)
      dserror("Lift&Drag id larger than FLUID_NUM_LD");

    /* read center coordinates */
    actdline->liftdrag->ld_center[0] = strtod(colptr,&colptr);
    actdline->liftdrag->ld_center[1] = strtod(colptr,&colptr);
    actdline->liftdrag->ld_center[2] = 0.0;
    /* read corresponding ALE-DLINE */
    actdline->liftdrag->aledline = strtol(colptr,&colptr,10);
    actdline->liftdrag->aledline--;

    frread();
  }
end:
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of inpdesign_line_liftdrag */



/*!----------------------------------------------------------------------
\brief input of surface lift&drag definition on design

<pre>                                                            mn 03/04
This routine reads an surface lift&drag definition on design.
</pre>

\warning There is nothing special to this routine
\return void                                               
\sa

*----------------------------------------------------------------------*/
static void inpdesign_surf_liftdrag()
{
  INT    i;
  INT    ierr;
  INT    ndsurf;
  INT    dsurfId;
  char  *colptr;
  DSURF *actdsurf;

#ifdef DEBUG 
  dstrc_enter("inpdesign_surf_liftdrag");
#endif

  /* find the beginning of surf fluid liftdrag conditions */
  if (frfind("--DESIGN FLUID SURF LIFT&DRAG") == 0) goto end;
  frread();

  /* read number of design surfs with conditions */
  frint("DSURF",&ndsurf,&ierr);
  dsassert(ierr==1,"Cannot read designsurf lift&drag definition");
  frread();


  /* start reading the design surfs */
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /* read the design surf Id */
    frint("E",&dsurfId,&ierr);
    dsassert(ierr==1,"Cannot read design-surf lift&drag definition");
    dsurfId--;

    /* find the dsurf */
    actdsurf=NULL;
    for (i=0; i<design->ndsurf; i++)
    {
      if (design->dsurf[i].Id ==  dsurfId) 
      {
        actdsurf = &(design->dsurf[i]);
        break;
      }
    }
    dsassert(actdsurf!=NULL,"Cannot read designsurf-lift&drag definition");

    /* allocate space for a liftdrag condition in this dsurf */
    actdsurf->liftdrag = (FLUID_LIFTDRAG_CONDITION*)CCACALLOC(1,sizeof(FLUID_LIFTDRAG_CONDITION));
    if (!actdsurf->liftdrag) dserror("Allocation of memory failed");

    /* move pointer behind the "-" sign */
    colptr = strstr(allfiles.actplace,"-");
    dsassert(colptr!=NULL,"Cannot read design-surf left&drag definition");
    colptr++;

    /* read the number */
    actdsurf->liftdrag->liftdrag = strtol(colptr,&colptr,10);
    if (actdsurf->liftdrag->liftdrag >= FLUID_NUM_LD)
      dserror("Lift&Drag id larger than FLUID_NUM_LD");

    /* read center coordinates */
    actdsurf->liftdrag->ld_center[0] = strtod(colptr,&colptr);
    actdsurf->liftdrag->ld_center[1] = strtod(colptr,&colptr);
    actdsurf->liftdrag->ld_center[2] = strtod(colptr,&colptr);
    /* read corresponding ALE-DLINE */
    actdsurf->liftdrag->aledline = strtol(colptr,&colptr,10);
    actdsurf->liftdrag->aledline--;

    frread();
  }
end:
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of inpdesign_surf_liftdrag */


/*!----------------------------------------------------------------------
\brief input of surface stability parameters

<pre>                                                       chfoe 01/04
This routine reads the stability parameters as a condition assigned to
a design surface for 2D fluid elements.
The single parameters within one condition line have to be in prescribed
order while additional free spaces do not matter.

The parameters accociated with gls stabilisation as well as the way they
are stored is due to genk (as it was converted from elemental input 
'f2_inp') and previous historical developments.
</pre>

\warning There is nothing special to this routine
\return void                                               
\sa

*-----------------------------------------------------------------------*/
static void inpdesign_surf_stability()
{
INT      i;
INT      ierr;
INT      ndsurf;
INT      dsurfId;
INT      ndum;          /* dummy value                                  */
INT      ihelem;        /* temporary variable for ihelem                */
INT      itaumu;        /*                                              */ 
INT      itaump;        /*                                              */
INT      itauc;         /* element flags                                */
char    *colptr;
char     buffer[200];
DSURF   *actdsurf;
STAB_PAR_GLS *gls;      /* pointer to GLS stabilisation parameters      */

#ifdef DEBUG 
dstrc_enter("inpdesign_surf_stability");
#endif

/*----------------------------------------------------------------------*/
/*---- find the beginning of line fluid freesurface coupling conditions */
if (frfind("--DESIGN FLUID SURFACE STABILISING CONDITIONS") == 0) goto end;
frread();
/*-------------------------- read number of surfaces with conditions ---*/
frint("DSURF",&ndsurf,&ierr);
dsassert(ierr==1,"Cannot read design surface for stabilisation");
frread();
/*-------------------------------------- start reading the design surfs */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
  /*------------------------------------------ read the design surf Id */
  frint("E",&dsurfId,&ierr);
  dsassert(ierr==1,"Cannot read design surface within stabilistation condition");
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
  dsassert(actdsurf!=NULL,"Cannot read design surface stabilisation definition");

  /*--------------------------------- move pointer behind the "-" sign */
  colptr = strstr(allfiles.actplace,"-");
  dsassert(colptr!=NULL,"Cannot read design surface stabilisation definition");
  colptr++;
  
  /*-------------------------------- read the type of stabilisation ---*/
  ierr=sscanf(colptr," %s ",buffer);
  dsassert(ierr==1,"Cannot read design surface stabilisation conditions");
  if (strncmp(buffer,"GLS",3)==0) 
  {
    actdsurf->stab_type = stab_gls;
    /* move pointer after "GLS" */
    while (colptr[0] == ' ')
      colptr++;
    colptr += 3;
  }
  else if (strncmp(buffer,"PresPro",7)==0)
  {
    actdsurf->stab_type = stab_prespro;
    /* move pointer after "PresPro" */
    while (colptr[0] == ' ')
      colptr++;
    colptr += 7;
  }
  else dserror("Unknown stabilisation type!");  

  /*---- All the following is done for eigther stabilisation type! ----*/
  switch (actdsurf->stab_type)
  {
  case stab_gls:
    /*-- allocate space for stabilisation parameters at this dsurf ---*/
    actdsurf->stabi.gls = (STAB_PAR_GLS*)CCACALLOC(1,sizeof(STAB_PAR_GLS));
    gls = actdsurf->stabi.gls;
    
    /*------------------------------------------------ read ISTABI ---*/
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design surface stabilisation ISTABI");
    if (strncmp(buffer,"yes",3)==0)
    {
      gls->istabi=1;
      /* move pointer after "yes" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 3;
    }
    else if (strncmp(buffer,"no",2)==0)
    {
      gls->istabi=0;
      /* move pointer after "no" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 2;
    }
    else
      dserror("Cannot read design surface stabilisation ISTABI");     
    
    /*------------------------------------------------ read IADVEC ---*/
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design surface stabilisation IADVEC");
    if (strncmp(buffer,"yes",3)==0)
    {
      gls->iadvec=1;
      /* move pointer after "yes" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 3;
    }
    else if (strncmp(buffer,"no",2)==0)
    {
      gls->iadvec=0;
      /* move pointer after "no" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 2;
    }
    else
      dserror("Cannot read design surface stabilisation IADVEC");     
    
    /*------------------------------------------------- read IPRES ---*/
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design surface stabilisation IPRES");
    if (strncmp(buffer,"yes",3)==0)
    {
      gls->ipres=1;
      /* move pointer after "yes" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 3;
    }
    else if (strncmp(buffer,"no",2)==0)
    {
      gls->ipres=0;
      /* move pointer after "no" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 2;
    }
    else
      dserror("Cannot read design surface stabilisation IPRES");
       
    /*------------------------------------------------- read IVISC ---*/
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design surface stabilisation IVISC");
    if (strncmp(buffer,"GLS-",4)==0)
    {
      gls->ivisc=1;
      /* move pointer after "GLS-" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 4;
    }
    else if (strncmp(buffer,"GLS+",4)==0)
    {
      gls->ivisc=2;
      /* move pointer after "GLS+" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 4;
    }
    else if (strncmp(buffer,"USFEM",5)==0)
    {
      gls->ivisc=3;
      /* move pointer after "USFEM" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 5;
    }
    else if (strncmp(buffer,"no_GLS",6)==0)
    {
      gls->ivisc=0;
      /* move pointer after "no_GLS" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 6;
    }
    else if (strncmp(buffer,"no_USFEM",8)==0)
    {
      gls->ivisc=-1;
      /* move pointer after "no_USFEM" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 8;
    }
    else
      dserror("Cannot read design surface stabilisation IVISC");     
    
    /*------------------------------------------------- read ICONT ---*/
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design surface stabilisation ICONT");
    if (strncmp(buffer,"yes",3)==0)
    {
      gls->icont=1;
      /* move pointer after "yes" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 3;
    }
    else if (strncmp(buffer,"no",2)==0)
    {
      gls->icont=0;
      /* move pointer after "no" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 2;
    }
    else
      dserror("Cannot read design surface stabilisation ICONT");
       
    /*------------------------------------------------ read ISTAPA ---*/
    gls->istapa = strtol(colptr,&colptr,10);
    
    /*------------------------------------------------ read NORM_P ---*/
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design surface stabilisation NORM_P");
    if (strncmp(buffer,"L_2",3)==0)
    {
      gls->norm_p=2;
      /* move pointer after "L_2" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 3;
    }
    else if (strncmp(buffer,"L_1",2)==0)
    {
      gls->norm_p=1;
      /* move pointer after "L_1" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 3;
    }
    else
      dserror("Cannot read design surface stabilisation NORM_P");
       
    /*---------------------------------------------------- read MK ---*/
    gls->mk = strtol(colptr,&colptr,10);
     
    /*------------------------------------------------ read IHELEM ---*/
    ihelem = strtol(colptr,&colptr,10);
    
    /*------------------------------------------------ read NINTHS ---*/
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design surface stabilisation NINTHS");
    if (strncmp(buffer,"at_center",9)==0)
    {
      gls->ninths=1;
      /* move pointer after "at_center" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 9;
    }
    else if (strncmp(buffer,"every_intpt",11)==0)
    {
      gls->ninths=2;
      /* move pointer after "every_intpt" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 11;
    }
    else
      dserror("Cannot read design surface stabilisation NINTHS");
       
    /*------------------------------------------------ read ISTAPC ---*/
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design surface stabilisation ISTAPC");
    if (strncmp(buffer,"at_center",9)==0)
    {
      gls->istapc=1;
      /* move pointer after "at_center" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 9;
    }
    else if (strncmp(buffer,"every_intpt",11)==0)
    {
      gls->istapc=2;
      /* move pointer after "every_intpt" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 11;
    }
    else
      dserror("Cannot read design surface stabilisation ISTAPC");

    /*------------------------------------- read stab const C_LAMB ---*/
    gls->clamb = strtod(colptr,&colptr);
    
    /*------------------------ initialise some stabilisation flags ---*/
    gls->istrle   = 0;
    gls->iareavol = 0;
    gls->iduring  = 0;
    gls->idiaxy   = 0;
    itaumu        = 0;
    itaump        = 0;
    itauc         = 0;

    math_intextract(ihelem,&ndum, 
                    &(gls->ihele[0]), 
		    &(gls->ihele[1]),
	            &(gls->ihele[2]));

    for(i=0;i<3;i++)
    {
      if (gls->ihele[i]==5)
      {
        gls->istrle  = 1;
        if (gls->iadvec!=0 && gls->ninths==1)
          itaumu = -1;
        if (gls->iadvec!=0 && gls->ninths!=1)
          itaumu = 1;
        if (gls->ipres!=0 && gls->ninths==1)
          itaump = -1;
        if (gls->ipres!=0 && gls->ninths!=1)
          itaump = 1;
        if (gls->icont!=0 && gls->ninths==1)
          itauc = -1;
        if (gls->icont!=0 && gls->ninths!=1)
          itauc = 1;     
      }
      else if (gls->ihele[i]!=0)
      {
        gls->iareavol = 1;
        if (gls->iadvec!=0 && gls->istapc==1)
          itaumu = -1;
        if (gls->iadvec!=0 && gls->istapc!=1)
          itaumu = 1;
        if (gls->ipres!=0 && gls->istapc==1)
          itaump = -1;
        if (gls->ipres!=0 && gls->istapc!=1)
          itaump = 1;
        if (gls->icont!=0 && gls->istapc==1)
          itauc = -1;
        if (gls->icont!=0 && gls->istapc!=1)
          itauc = 1;
	 
        if (gls->ihele[i]==4)
          gls->idiaxy = 1;         
      }
    }

    if (gls->istrle==1 && gls->ninths!=1)
      gls->iduring = 1;
    if (gls->iareavol==1 && gls->istapc!=1)
      gls->iduring = 1;

    /*------------------------------ store data within the condition ---*/
    gls->itau[0] = itaumu;
    gls->itau[1] = itaump;
    gls->itau[2] = itauc;

  break;
  case stab_prespro:
    dserror("Pressure Projection type of stabilisation has not yet been implemented!");
  break;
  default:
    dserror("Unknown stabilisation type!");
  }
  frread();
}
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_surf_stability */


/*!----------------------------------------------------------------------
\brief input of surface stability parameters

<pre>                                                       chfoe 01/04
This routine reads the stability parameters as a condition assigned to
a design volume for 3D fluid elements.
The single parameters within one condition line have to be in prescribed
order while additional free spaces do not matter.

The parameters accociated with gls stabilisation as well as the way they
are stored is due to genk (as it was converted from elemental input 
'f3_inp') and previous historical developments.
</pre>

\warning There is nothing special to this routine
\return void                                               
\sa

*-----------------------------------------------------------------------*/
static void inpdesign_vol_stability()
{

INT      i;
INT      ierr;
INT      ndvol;
INT      dvolId;
INT      ndum;          /* dummy value                                  */
INT      ihelem;        /* temporary variable for ihelem                */
INT      itaumu;        /*                                              */ 
INT      itaump;        /*                                              */
INT      itauc;         /* element flags                                */
char    *colptr;
char     buffer[200];
DVOL   *actdvol;
STAB_PAR_GLS *gls;      /* pointer to GLS stabilisation parameters      */

#ifdef DEBUG 
dstrc_enter("inpdesign_vol_stability");
#endif
/*----------------------------------------------------------------------*/
/*---- find the beginning of line fluid freesurface coupling conditions */
if (frfind("--DESIGN FLUID VOLUME STABILISING CONDITIONS") == 0) goto end;
frread();
/*-------------------------- read number of surfaces with conditions ---*/
frint("DVOL",&ndvol,&ierr);
dsassert(ierr==1,"Cannot read design volume for stabilisation");
frread();
/*-------------------------------------- start reading the design surfs */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
  /*------------------------------------------ read the design surf Id */
  frint("E",&dvolId,&ierr);
  dsassert(ierr==1,"Cannot read design volume within stabilistation condition");
  dvolId--;
  
  /*--------------------------------------------------- find the dsurf */
  actdvol=NULL; 
  for (i=0; i<design->ndvol; i++)
  {
    if (design->dvol[i].Id == dvolId) 
    {
      actdvol = &(design->dvol[i]);
      break;
    }
  }
  dsassert(actdvol!=NULL,"Cannot read design volume stabilisation definition");

  /*--------------------------------- move pointer behind the "-" sign */
  colptr = strstr(allfiles.actplace,"-");
  dsassert(colptr!=NULL,"Cannot read design volume stabilisation definition");
  colptr++;
  
  /*-------------------------------- read the type of stabilisation ---*/
  ierr=sscanf(colptr," %s ",buffer);
  dsassert(ierr==1,"Cannot read design volume stabilisation conditions");
  if (strncmp(buffer,"GLS",3)==0) 
  {
    actdvol->stab_type = stab_gls;
    /* move pointer after "GLS" */
    while (colptr[0] == ' ')
      colptr++;
    colptr += 3;
  }
  else if (strncmp(buffer,"PresPro",7)==0)
  {
    actdvol->stab_type = stab_prespro;
    /* move pointer after "PresPro" */
    while (colptr[0] == ' ')
      colptr++;
    colptr += 7;
  }
  else dserror("Unknown stabilisation type!");  
  
  /*---- All the following is done for eigther stabilisation type! ----*/
  switch (actdvol->stab_type)
  {
  case stab_gls:
    /*-- allocate space for stabilisation parameters at this dsurf ---*/
    actdvol->stabi.gls = (STAB_PAR_GLS*)CCACALLOC(1,sizeof(STAB_PAR_GLS));
    gls = actdvol->stabi.gls;
    
    /*------------------------------------------------ read ISTABI ---*/
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design volume stabilisation ISTABI");
    if (strncmp(buffer,"yes",3)==0)
    {
      gls->istabi=1;
      /* move pointer after "yes" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 3;
    }
    else if (strncmp(buffer,"no",2)==0)
    {
      gls->istabi=0;
      /* move pointer after "no" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 2;
    }
    else
      dserror("Cannot read design volume stabilisation ISTABI");     
    
    /*------------------------------------------------ read IADVEC ---*/
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design volume stabilisation IADVEC");
    if (strncmp(buffer,"yes",3)==0)
    {
      gls->iadvec=1;
      /* move pointer after "yes" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 3;
    }
    else if (strncmp(buffer,"no",2)==0)
    {
      gls->iadvec=0;
      /* move pointer after "no" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 2;
    }
    else
      dserror("Cannot read design volume stabilisation IADVEC");     
    
    /*------------------------------------------------- read IPRES ---*/
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design volume stabilisation IPRES");
    if (strncmp(buffer,"yes",3)==0)
    {
      gls->ipres=1;
      /* move pointer after "yes" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 3;
    }
    else if (strncmp(buffer,"no",2)==0)
    {
      gls->ipres=0;
      /* move pointer after "no" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 2;
    }
    else
      dserror("Cannot read design volume stabilisation IPRES");
       
    /*------------------------------------------------- read IVISC ---*/
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design volume stabilisation IVISC");
    if (strncmp(buffer,"GLS-",4)==0)
    {
      gls->ivisc=1;
      /* move pointer after "GLS-" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 4;
    }
    else if (strncmp(buffer,"GLS+",4)==0)
    {
      gls->ivisc=2;
      /* move pointer after "GLS+" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 4;
    }
    else if (strncmp(buffer,"USFEM",5)==0)
    {
      gls->ivisc=3;
      /* move pointer after "USFEM" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 5;
    }
    else if (strncmp(buffer,"no_GLS",6)==0)
    {
      gls->ivisc=0;
      /* move pointer after "no_GLS" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 6;
    }
    else if (strncmp(buffer,"no_USFEM",8)==0)
    {
      gls->ivisc=-1;
      /* move pointer after "no_USFEM" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 8;
    }
    else
      dserror("Cannot read design volume stabilisation IVISC");     
    
    /*------------------------------------------------- read ICONT ---*/
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design volume stabilisation ICONT");
    if (strncmp(buffer,"yes",3)==0)
    {
      gls->icont=1;
      /* move pointer after "yes" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 3;
    }
    else if (strncmp(buffer,"no",2)==0)
    {
      gls->icont=0;
      /* move pointer after "no" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 2;
    }
    else
      dserror("Cannot read design volume stabilisation ICONT");
       
    /*------------------------------------------------ read ISTAPA ---*/
    gls->istapa = strtol(colptr,&colptr,10);
    
    /*------------------------------------------------ read NORM_P ---*/
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design volume stabilisation NORM_P");
    if (strncmp(buffer,"L_2",3)==0)
    {
      gls->norm_p=2;
      /* move pointer after "L_2" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 3;
    }
    else if (strncmp(buffer,"L_1",2)==0)
    {
      gls->norm_p=1;
      /* move pointer after "L_1" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 3;
    }
    else
      dserror("Cannot read design volume stabilisation NORM_P");
       
    /*---------------------------------------------------- read MK ---*/
    gls->mk = strtol(colptr,&colptr,10);
     
    /*------------------------------------------------ read IHELEM ---*/
    ihelem = strtol(colptr,&colptr,10);
    
    /*------------------------------------------------ read NINTHS ---*/
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design volume stabilisation NINTHS");
    if (strncmp(buffer,"at_center",9)==0)
    {
      gls->ninths=1;
      /* move pointer after "at_center" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 9;
    }
    else if (strncmp(buffer,"every_intpt",11)==0)
    {
      gls->ninths=2;
      /* move pointer after "every_intpt" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 11;
    }
    else
      dserror("Cannot read design volume stabilisation NINTHS");
       
    /*------------------------------------------------ read ISTAPC ---*/
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design volume stabilisation ISTAPC");
    if (strncmp(buffer,"at_center",9)==0)
    {
      gls->istapc=1;
      /* move pointer after "at_center" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 9;
    }
    else if (strncmp(buffer,"every_intpt",11)==0)
    {
      gls->istapc=2;
      /* move pointer after "every_intpt" */
      while (colptr[0] == ' ')
        colptr++;
      colptr += 11;
    }
    else
      dserror("Cannot read design volume stabilisation ISTAPC");

    /*------------------------------------- read stab const C_LAMB ---*/
    gls->clamb = strtod(colptr,&colptr);
    
    /*------------------------ initialise some stabilisation flags ---*/
    gls->istrle   = 0;
    gls->iareavol = 0;
    gls->iduring  = 0;
    gls->idiaxy   = 0;
    itaumu        = 0;
    itaump        = 0;
    itauc         = 0;

    math_intextract(ihelem,&ndum, 
                    &(gls->ihele[0]), 
		    &(gls->ihele[1]),
	            &(gls->ihele[2]));

    for(i=0;i<3;i++)
    {
      if (gls->ihele[i]==5)
      {
        gls->istrle  = 1;
        if (gls->iadvec!=0 && gls->ninths==1)
          itaumu = -1;
        if (gls->iadvec!=0 && gls->ninths!=1)
          itaumu = 1;
        if (gls->ipres!=0 && gls->ninths==1)
          itaump = -1;
        if (gls->ipres!=0 && gls->ninths!=1)
          itaump = 1;
        if (gls->icont!=0 && gls->ninths==1)
          itauc = -1;
        if (gls->icont!=0 && gls->ninths!=1)
          itauc = 1;     
      }
      else if (gls->ihele[i]!=0)
      {
        gls->iareavol = 1;
        if (gls->iadvec!=0 && gls->istapc==1)
          itaumu = -1;
        if (gls->iadvec!=0 && gls->istapc!=1)
          itaumu = 1;
        if (gls->ipres!=0 && gls->istapc==1)
          itaump = -1;
        if (gls->ipres!=0 && gls->istapc!=1)
          itaump = 1;
        if (gls->icont!=0 && gls->istapc==1)
          itauc = -1;
        if (gls->icont!=0 && gls->istapc!=1)
          itauc = 1;
	 
        if (gls->ihele[i]==4)
          gls->idiaxy = 1;         
      }
    }

    if (gls->istrle==1 && gls->ninths!=1)
      gls->iduring = 1;
    if (gls->iareavol==1 && gls->istapc!=1)
      gls->iduring = 1;

    /*------------------------------ store data within the condition ---*/
    gls->itau[0] = itaumu;
    gls->itau[1] = itaump;
    gls->itau[2] = itauc;

  break;
  case stab_prespro:
    dserror("Pressure Projection type of stabilisation has not yet been implemented!");
  break;
  default:
    dserror("Unknown stabilisation type!");
  }
  frread();
}
 
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_vol_stability */
#endif



#ifdef D_AXISHELL
/*!----------------------------------------------------------------------
\brief input of line axishell thickness conditions on design

<pre>                                                              mn 05/03
This routine reads an axishell thickness conditions.
</pre>

\warning There is nothing special to this routine
\return void                                               
\sa

*----------------------------------------------------------------------*/
static void inpdesign_line_thickness()
{
INT    i;
INT    ierr;
INT    ndline;
INT    dlineId;
char  *colptr;
char   buffer[200];
DLINE *actdline;

#ifdef DEBUG 
dstrc_enter("inpdesign_line_thickness");
#endif

/*----------------------------------------------------------------------*/
/*---- find the beginning of line fluid freesurface coupling conditions */
if (frfind("--DESIGN AXISHELL THICKNESS LINE CONDITIONS") == 0) goto end;
frread();
/*------------------------ read number of design lines with conditions */
frint("DLINE",&ndline,&ierr);
dsassert(ierr==1,"Cannot read design-line thickness conditions");
frread();
/*-------------------------------------- start reading the design lines */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design line Id */
   frint("E",&dlineId,&ierr);
   dsassert(ierr==1,"Cannot read design-line thickness conditions");
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
   dsassert(actdline!=NULL,"Cannot read design-line thickness conditions");
   /*----------- allocate space for a coupling condition in this dline */
   actdline->thickness = (SAXI_THICK_CONDITION*)CCACALLOC(1,sizeof(SAXI_THICK_CONDITION));
   if (!actdline->thickness) dserror("Allocation of memory failed");

   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-line thickness conditions");
   colptr++;
   /*------------------------------------------------ read the fieldtyp */
  
   /*---- now read the MAXDOFPERNODE flags for the local slippage conditions */ 
   for (i=0; i<2; i++)
   actdline->thickness->value[i] = strtod(colptr,&colptr);
   
   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_line_thickness */



/*!----------------------------------------------------------------------
\brief input of line axishell load conditions on design

<pre>                                                              mn 05/03
This routine reads an axishell load conditions.
</pre>

\warning There is nothing special to this routine
\return void                                               
\sa

*----------------------------------------------------------------------*/
static void inpdesign_line_axishellload()
{
INT    i;
INT    ierr;
INT    ndline;
INT    dlineId;
char  *colptr;
char   buffer[200];
DLINE *actdline;

#ifdef DEBUG 
dstrc_enter("inpdesign_line_axishellload");
#endif

/*----------------------------------------------------------------------*/
/*---- find the beginning of line fluid freesurface coupling conditions */
if (frfind("--DESIGN AXISHELL LOAD LINE CONDITIONS") == 0) goto end;
frread();
/*------------------------ read number of design lines with conditions */
frint("DLINE",&ndline,&ierr);
dsassert(ierr==1,"Cannot read design-line thickness conditions");
frread();
/*-------------------------------------- start reading the design lines */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design line Id */
   frint("E",&dlineId,&ierr);
   dsassert(ierr==1,"Cannot read design-line thickness conditions");
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
   dsassert(actdline!=NULL,"Cannot read design-line thickness conditions");
   /*----------- allocate space for a coupling condition in this dline */
   actdline->axishellload = (SAXI_LOAD_CONDITION*)CCACALLOC(1,sizeof(SAXI_LOAD_CONDITION));
   if (!actdline->axishellload) dserror("Allocation of memory failed");

   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-line axishell-load conditions");
   colptr++;
   /*------------------------------------------------ read the fieldtyp */
  
   /*---- now read the MAXDOFPERNODE flags for the local slippage conditions */ 
   for (i=0; i<2; i++)
     actdline->axishellload->pv[i]     = strtod(colptr,&colptr);
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot read axishell interpolation type");
   if (strncmp(buffer,"arclength",9)==0) 
   {
     actdline->axishellload->interpol_pv = 0;
     colptr = strstr(colptr,"arclength");
     dsassert(colptr!=NULL,"Cannot read axishell interpolation type");
     colptr+=9;
   }
   if (strncmp(buffer,"vert._axis",10)==0) 
   {
     actdline->axishellload->interpol_pv = 1;
     colptr = strstr(colptr,"vert._axis");
     dsassert(colptr!=NULL,"Cannot read axishell interpolation type");
     colptr+=10;
   }


   for (i=0; i<2; i++)
     actdline->axishellload->ph[i]     = strtod(colptr,&colptr);
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot read axishell interpolation type");
   if (strncmp(buffer,"arclength",9)==0) 
   {
     actdline->axishellload->interpol_ph = 0;
     colptr = strstr(colptr,"arclength");
     dsassert(colptr!=NULL,"Cannot read axishell interpolation type");
     colptr+=9;
   }
   if (strncmp(buffer,"vert._axis",10)==0) 
   {
     actdline->axishellload->interpol_ph = 1;
     colptr = strstr(colptr,"vert._axis");
     dsassert(colptr!=NULL,"Cannot read axishell interpolation type");
     colptr+=10;
   }


   for (i=0; i<2; i++)
     actdline->axishellload->px[i]     = strtod(colptr,&colptr);
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot read axishell interpolation type");
   if (strncmp(buffer,"arclength",9)==0) 
   {
     actdline->axishellload->interpol_px = 0;
     colptr = strstr(colptr,"arclength");
     dsassert(colptr!=NULL,"Cannot read axishell interpolation type");
     colptr+=9;
   }
   if (strncmp(buffer,"vert._axis",10)==0) 
   {
     actdline->axishellload->interpol_px = 1;
     colptr = strstr(colptr,"vert._axis");
     dsassert(colptr!=NULL,"Cannot read axishell interpolation type");
     colptr+=10;
   }


   for (i=0; i<2; i++)
     actdline->axishellload->pw[i]     = strtod(colptr,&colptr);
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot read axishell interpolation type");
   if (strncmp(buffer,"arclength",9)==0) 
   {
     actdline->axishellload->interpol_pw = 0;
     colptr = strstr(colptr,"arclength");
     dsassert(colptr!=NULL,"Cannot read axishell interpolation type");
     colptr+=9;
   }
   if (strncmp(buffer,"vert._axis",10)==0) 
   {
     actdline->axishellload->interpol_pw = 1;
     colptr = strstr(colptr,"vert._axis");
     dsassert(colptr!=NULL,"Cannot read axishell interpolation type");
     colptr+=10;
   }

   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_line_axishellload */


/*!----------------------------------------------------------------------
\brief input of line axishell coordinate system conditions on design

<pre>                                                              mn 05/03
This routine reads an axishell coordinate system conditions.
</pre>

\warning There is nothing special to this routine
\return void                                               
\sa

*----------------------------------------------------------------------*/
static void inpdesign_point_axishellcos()
{
INT     i;
INT     ierr;
INT     ndnode;
INT     dnodeId;
char   *colptr;
char    buffer[200];
DNODE *actdnode;
#ifdef DEBUG 
dstrc_enter("inpdesign_point_axishellcos");
#endif

/*----------------------------------------------------------------------*/
/*---- find the beginning of line fluid freesurface coupling conditions */
if (frfind("--DESIGN AXISHELL COS POINT CONDITIONS") == 0) goto end;
frread();
/*------------------------ read number of design points with conditions */
frint("DPOINT",&ndnode,&ierr);
dsassert(ierr==1,"Cannot read design-point cos conditions");



frread();
/*-------------------------------------- start reading the design lines */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design line Id */
   frint("E",&dnodeId,&ierr);
   dsassert(ierr==1,"Cannot read design-point cos conditions");
   dnodeId--;
   /*--------------------------------------------------- find the dline */
   actdnode=NULL;
   for (i=0; i<design->ndnode; i++)
   {
      if (design->dnode[i].Id ==  dnodeId) 
      {
         actdnode = &(design->dnode[i]);
         break;
      }
   }
   dsassert(actdnode!=NULL,"Cannot read design-point cos conditions");
   /*--------------------------------- move pointer behind the "-" sign */
   colptr = strstr(allfiles.actplace,"-");
   dsassert(colptr!=NULL,"Cannot read design-point cos conditions");
   colptr++;
   /*------------------------------------------------ read the fieldtyp */
   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot read axishell cos type");
   if (strncmp(buffer,"global",6)==0) 
   {
     actdnode->cos_type = 0;
     colptr = strstr(colptr,"global");
     dsassert(colptr!=NULL,"Cannot read axishell cos type");
     colptr+=6;
   }
   if (strncmp(buffer,"local",5)==0) 
   {
     actdnode->cos_type = 1;
     colptr = strstr(colptr,"local");
     dsassert(colptr!=NULL,"Cannot read axishell cos type");
     colptr+=5;
   }
  


   frread();
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_line_axishellload */
#endif






#ifdef WALLCONTACT
/*----------------------------------------------------------------------*
 |                                                         oezdem 08/03 |
 *----------------------------------------------------------------------*/
static void inpdesign_line_contact()
{
INT    i,counter=0;
INT    ierr;
INT    ndline;
INT    dlineId;
char  *colptr;
char   buffer[200];
DLINE *actdline;
#ifdef DEBUG 
dstrc_enter("inpdesign_line_contact");
#endif
/*----------------------------------------------------------------------*/
/*----------------------- find the beginning of line contact conditions */
frfind("----CONTACT CONDITIONS");
frread();
/*------------------------- read number of design lines with conditions */
frint("DLINE",&ndline,&ierr);
dsassert(ierr==1,"Cannot read design-line contact conditions");
contact.ndline = ndline;
contact.dline = (DLINE**)CCACALLOC(ndline,sizeof(DLINE*));
if (!(contact.dline)) dserror("Allocation of memory failed");
frread();
/*-------------------------------------- start reading the design lines */
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   /*------------------------------------------ read the design line Id */
   frint("E",&dlineId,&ierr);
   dsassert(ierr==1,"Cannot read design-line contact conditions");
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
   dsassert(actdline!=NULL,"Cannot read design-line contact conditions");
   contact.dline[counter] = actdline;
   counter++;
   actdline->contype = contact_none;
   frchk("Master",&ierr);
   if (ierr) actdline->contype = contact_master;
   frchk("Slave",&ierr);
   if (ierr) actdline->contype = contact_slave;
   frchk("Self",&ierr);
   if (ierr) actdline->contype = contact_self;


   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_line_contact */
#endif
/*! @} (documentation module close)*/
