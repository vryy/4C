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
/*----------------------------------------------------------------------*
 | input of conditions                                    m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_conditions()
{
int  ierr;
int  i;
#ifdef DEBUG 
dstrc_enter("inp_conditions");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------ input of time curves */
inp_cond_curve();
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
/*---------------------------------- input of line coupling conditions */
inpdesign_line_couple();
/*---------------------------------- input of surf coupling conditions */
inpdesign_surf_couple();
/*---------------------------------- input of vol coupling conditions */
inpdesign_vol_couple();
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
int    i,j;
int    ierr;
int    ndnode;
int    dnodeId;
int    foundit;
char  *colptr;
char   buffer[200];
DNODE *actdnode;
#ifdef DEBUG 
dstrc_enter("inpdesign_nodal_dirich");
#endif
/*----------------------------------------------------------------------*/
/*-------------------- find the beginning of nodal dirichlet conditions */
frfind("--DESIGN POINT DIRICH CONDITIONS");
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
   /*----------------------------------- read the curvenumber or "none" */
/*   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot read design-nodal dirichlet conditions");
   if (strncmp(buffer,"none",4)==0) 
   {
      actdnode->dirich->curve = 0;
      colptr = strstr(allfiles.actplace,"none");
      dsassert(colptr!=NULL,"Cannot read design-nodal dirichlet conditions");
      colptr += 4;
   }
   else
   {
      ierr=sscanf(colptr," %d ",&(actdnode->dirich->curve));
      dsassert(ierr==1,"Cannot read design-nodal dirichlet conditions");
      colptr = strpbrk(colptr,"1234567890");
      colptr++;
   }  */
   /*--- now read the MAXCONDPERNODE flags for the dirichlet conditions */ 
   /*---------------------------- and the MAXCONDPERNODE values of them */
   amdef("onoff",&(actdnode->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
   amdef("val",&(actdnode->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
   amdef("curve",&(actdnode->dirich->curve),MAXDOFPERNODE,1,"IV");

   /* Initialize some arrays to ZERO if MAXDOFPERNODE != MAXCONDPERNODE -> shell9 */
   #ifdef D_SHELL9
      amzero(&(actdnode->dirich->dirich_onoff));
      amzero(&(actdnode->dirich->dirich_val));
      amzero(&(actdnode->dirich->curve));
   #endif /*D_SHELL9*/

   /* NOTE: number of read values = MAXCONDPERNODE  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   for (i=0; i<MAXCONDPERNODE; i++)
   actdnode->dirich->dirich_onoff.a.iv[i] = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   actdnode->dirich->dirich_val.a.dv[i] = strtod(colptr,&colptr);
   /*---------------------------------- read the curvenumbers or "none" */
   for (i=0; i<MAXCONDPERNODE; i++)
   {
      ierr=sscanf(colptr," %s ",buffer);   
      dsassert(ierr==1,"Cannot read design-nodal dirichlet conditions");
      if (strncmp(buffer,"none",4)==0) 
      {
         actdnode->dirich->curve.a.iv[i] = 0;
         colptr = strstr(colptr,"none");
         dsassert(colptr!=NULL,"Cannot read design-nodal dirichlet conditions");
         colptr += 4;
      }
      else
      {
         ierr=sscanf(colptr," %d ",&(actdnode->dirich->curve.a.iv[i]));
         dsassert(ierr==1,"Cannot read design-nodal dirichlet conditions");
         colptr = strpbrk(colptr,"1234567890");
         colptr++;
      }   
   }
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/
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
int    i,j;
int    ierr;
int    ndline;
int    dlineId;
int    foundit;
char  *colptr;
char   buffer[200];
DLINE *actdline;
#ifdef DEBUG 
dstrc_enter("inpdesign_line_dirich");
#endif
/*----------------------------------------------------------------------*/
/*-------------------- find the beginning of line dirichlet conditions */
frfind("--DESIGN LINE DIRICH CONDITIONS");
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
/*   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot read design-line dirichlet conditions");
   if (strncmp(buffer,"none",4)==0) 
   {
      actdline->dirich->curve = 0;
      colptr = strstr(allfiles.actplace,"none");
      dsassert(colptr!=NULL,"Cannot read design-line dirichlet conditions");
      colptr += 4;
   }
   else
   {
      ierr=sscanf(colptr," %d ",&(actdline->dirich->curve));
      dsassert(ierr==1,"Cannot read design-line dirichlet conditions");
      colptr = strpbrk(colptr,"1234567890");
      colptr++;
   }         */
   /*--- now read the MAXCONDPERNODE flags for the dirichlet conditions */ 
   /*---------------------------- and the MAXCONDPERNODE values of them */
   amdef("onoff",&(actdline->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
   amdef("val",&(actdline->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
   amdef("curve",&(actdline->dirich->curve),MAXDOFPERNODE,1,"IV");

   /* Initialize some arrays to ZERO if MAXDOFPERNODE != MAXCONDPERNODE -> shell9 */
   #ifdef D_SHELL9
      amzero(&(actdline->dirich->dirich_onoff));
      amzero(&(actdline->dirich->dirich_val));
      amzero(&(actdline->dirich->curve));
   #endif /*D_SHELL9*/

   /* NOTE: number of read values = MAXCONDPERNODE  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   for (i=0; i<MAXCONDPERNODE; i++)
   actdline->dirich->dirich_onoff.a.iv[i] = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   actdline->dirich->dirich_val.a.dv[i] = strtod(colptr,&colptr);
   /*---------------------------------- read the curvenumbers or "none" */
   for (i=0; i<MAXCONDPERNODE; i++)
   {
      ierr=sscanf(colptr," %s ",buffer);
      dsassert(ierr==1,"Cannot read design-line dirichlet conditions");
      if (strncmp(buffer,"none",4)==0) 
      {
         actdline->dirich->curve.a.iv[i] = 0;
         colptr = strstr(colptr,"none");
         dsassert(colptr!=NULL,"Cannot read design-line dirichlet conditions");
         colptr += 4;
      }
      else
      {
         ierr=sscanf(colptr," %d ",&(actdline->dirich->curve.a.iv[i]));
         dsassert(ierr==1,"Cannot read design-line dirichlet conditions");
         colptr = strpbrk(colptr,"1234567890");
         colptr++;
      }
   }           
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/
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
int    i,j;
int    ierr;
int    ndsurf;
int    dsurfId;
int    foundit;
char  *colptr;
char   buffer[200];
DSURF *actdsurf;
#ifdef DEBUG 
dstrc_enter("inpdesign_surf_dirich");
#endif
/*----------------------------------------------------------------------*/
/*------------------ find the beginning of surface dirichlet conditions */
frfind("--DESIGN SURF DIRICH CONDITIONS");
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
   /*----------------------------------- read the curvenumber or "none" */
/*   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot read design-surf dirichlet conditions");
   if (strncmp(buffer,"none",4)==0) 
   {
      actdsurf->dirich->curve = 0;
      colptr = strstr(allfiles.actplace,"none");
      dsassert(colptr!=NULL,"Cannot read design-surf dirichlet conditions");
      colptr += 4;
   }
   else
   {
      ierr=sscanf(colptr," %d ",&(actdsurf->dirich->curve));
      dsassert(ierr==1,"Cannot read design-surf dirichlet conditions");
      colptr = strpbrk(colptr,"1234567890");
      colptr++;
   }   */
   /*--- now read the MAXCONDPERNODE flags for the dirichlet conditions */ 
   /*---------------------------- and the MAXCONDPERNODE values of them */
   amdef("onoff",&(actdsurf->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
   amdef("val",&(actdsurf->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
   amdef("curve",&(actdsurf->dirich->curve),MAXDOFPERNODE,1,"IV");

   /* Initialize some arrays to ZERO if MAXDOFPERNODE != MAXCONDPERNODE -> shell9 */
   #ifdef D_SHELL9
      amzero(&(actdsurf->dirich->dirich_onoff));
      amzero(&(actdsurf->dirich->dirich_val));
      amzero(&(actdsurf->dirich->curve));
   #endif /*D_SHELL9*/

   /* NOTE: number of read values = MAXCONDPERNODE  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   for (i=0; i<MAXCONDPERNODE; i++)
   actdsurf->dirich->dirich_onoff.a.iv[i] = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   actdsurf->dirich->dirich_val.a.dv[i] = strtod(colptr,&colptr);   
   /*----------------------------------- read the curvenumber or "none" */
   for (i=0; i<MAXCONDPERNODE; i++)
   {
      ierr=sscanf(colptr," %s ",buffer);
      dsassert(ierr==1,"Cannot read design-surf dirichlet conditions");
      if (strncmp(buffer,"none",4)==0) 
      {
         actdsurf->dirich->curve.a.iv[i] = 0;
         colptr = strstr(colptr,"none");
         dsassert(colptr!=NULL,"Cannot read design-surf dirichlet conditions");
         colptr += 4;
      }
      else
      {
         ierr=sscanf(colptr," %d ",&(actdsurf->dirich->curve.a.iv[i]));
         dsassert(ierr==1,"Cannot read design-surf dirichlet conditions");
         colptr = strpbrk(colptr,"1234567890");
         colptr++;
      }
   }
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/
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
int    i,j;
int    ierr;
int    ndvol;
int    dvolId;
int    foundit;
char  *colptr;
char   buffer[200];
DVOL  *actdvol;
#ifdef DEBUG 
dstrc_enter("inpdesign_vol_dirich");
#endif
/*----------------------------------------------------------------------*/
/*------------------ find the beginning of surface dirichlet conditions */
frfind("--DESIGN VOL DIRICH CONDITIONS");
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
   /*----------------------------------- read the curvenumber or "none" */
/*   ierr=sscanf(colptr," %s ",buffer);
   dsassert(ierr==1,"Cannot read design-vol dirichlet conditions");
   if (strncmp(buffer,"none",4)==0) 
   {
      actdvol->dirich->curve = 0;
      colptr = strstr(allfiles.actplace,"none");
      dsassert(colptr!=NULL,"Cannot read design-vol dirichlet conditions");
      colptr += 4;
   }
   else
   {
      ierr=sscanf(colptr," %d ",&(actdvol->dirich->curve));
      dsassert(ierr==1,"Cannot read design-vol dirichlet conditions");
      colptr = strpbrk(colptr,"1234567890");
      colptr++;
   }  */
   /*--- now read the MAXCONDPERNODE flags for the dirichlet conditions */ 
   /*---------------------------- and the MAXCONDPERNODE values of them */
   amdef("onoff",&(actdvol->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
   amdef("val",&(actdvol->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
   amdef("curve",&(actdvol->dirich->curve),MAXDOFPERNODE,1,"IV"); 

   /* Initialize some arrays to ZERO if MAXDOFPERNODE != MAXCONDPERNODE -> shell9 */
   #ifdef D_SHELL9
      amzero(&(actdvol->dirich->dirich_onoff));
      amzero(&(actdvol->dirich->dirich_val));
      amzero(&(actdvol->dirich->curve));
   #endif /*D_SHELL9*/

   /* NOTE: number of read values = MAXCONDPERNODE  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   for (i=0; i<MAXCONDPERNODE; i++)
   actdvol->dirich->dirich_onoff.a.iv[i] = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   actdvol->dirich->dirich_val.a.dv[i] = strtod(colptr,&colptr);   
   /*----------------------------------- read the curvenumber or "none" */
   for (i=0; i<MAXCONDPERNODE; i++)
   {
        ierr=sscanf(colptr," %s ",buffer);
        dsassert(ierr==1,"Cannot read design-vol dirichlet conditions");
        if (strncmp(buffer,"none",4)==0) 
        {
           actdvol->dirich->curve.a.iv[i] = 0;
           colptr = strstr(colptr,"none");
           dsassert(colptr!=NULL,"Cannot read design-vol dirichlet conditions");
           colptr += 4;
        }
        else
        {
           ierr=sscanf(colptr," %d ",&(actdvol->dirich->curve.a.iv[i]));
           dsassert(ierr==1,"Cannot read design-vol dirichlet conditions");
           colptr = strpbrk(colptr,"1234567890");
           colptr++;
        }  
    }
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/
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
int    i,j;
int    ierr;
int    ndnode;
int    dnodeId;
int    foundit;
char  *colptr;
char   buffer[200];
DNODE *actdnode;
#ifdef DEBUG 
dstrc_enter("inpdesign_nodal_neum");
#endif
/*----------------------------------------------------------------------*/
/*-------------------- find the beginning of nodal neumann conditions */
frfind("--DESIGN POINT NEUMANN CONDITIONS");
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
   /*--- now read the MAXCONDPERNODE flags for the neumann conditions */ 
   /*-------------------------- and the MAXCONDPERNODE values of them */
   amdef("onoff",&(actdnode->neum->neum_onoff),MAXDOFPERNODE,1,"IV");
   amdef("val",&(actdnode->neum->neum_val),MAXDOFPERNODE,1,"DV");

   /* Initialize some arrays to ZERO if MAXDOFPERNODE != MAXCONDPERNODE -> shell9 */
   #ifdef D_SHELL9
      amzero(&(actdnode->neum->neum_onoff));
      amzero(&(actdnode->neum->neum_val));
   #endif /*D_SHELL9*/

   /* NOTE: number of read values = MAXCONDPERNODE  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   for (i=0; i<MAXCONDPERNODE; i++)
   actdnode->neum->neum_onoff.a.iv[i] = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   actdnode->neum->neum_val.a.dv[i] = strtod(colptr,&colptr);
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
int    i,j;
int    ierr;
int    ndline;
int    dlineId;
int    foundit;
char  *colptr;
char   buffer[200];
DLINE *actdline;
#ifdef DEBUG 
dstrc_enter("inpdesign_line_neum");
#endif
/*----------------------------------------------------------------------*/
/*-------------------- find the beginning of line neumann conditions */
frfind("--DESIGN LINE NEUMANN CONDITIONS");
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
   /*--- now read the MAXCONDPERNODE flags for the neumann conditions */ 
   /*-------------------------- and the MAXCONDPERNODE values of them */
   amdef("onoff",&(actdline->neum->neum_onoff),MAXDOFPERNODE,1,"IV");
   amdef("val",&(actdline->neum->neum_val),MAXDOFPERNODE,1,"DV");

   /* Initialize some arrays to ZERO if MAXDOFPERNODE != MAXCONDPERNODE -> shell9 */
   #ifdef D_SHELL9
      amzero(&(actdline->neum->neum_onoff));
      amzero(&(actdline->neum->neum_val));
   #endif /*D_SHELL9*/

   /* NOTE: number of read values = MAXCONDPERNODE  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   for (i=0; i<MAXCONDPERNODE; i++)
   actdline->neum->neum_onoff.a.iv[i] = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   actdline->neum->neum_val.a.dv[i] = strtod(colptr,&colptr);
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
int    i,j;
int    ierr;
int    ndsurf;
int    dsurfId;
int    foundit;
char  *colptr;
char   buffer[200];
DSURF *actdsurf;
#ifdef DEBUG 
dstrc_enter("inpdesign_surf_neum");
#endif
/*----------------------------------------------------------------------*/
/*------------------ find the beginning of surface neumann conditions */
frfind("--DESIGN SURF NEUMANN CONDITIONS");
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
   /*--- now read the MAXCONDPERNODE flags for the neumann conditions */ 
   /*-------------------------- and the MAXCONDPERNODE values of them */
   amdef("onoff",&(actdsurf->neum->neum_onoff),MAXDOFPERNODE,1,"IV");
   amdef("val",&(actdsurf->neum->neum_val),MAXDOFPERNODE,1,"DV");

   /* Initialize some arrays to ZERO if MAXDOFPERNODE != MAXCONDPERNODE -> shell9 */
   #ifdef D_SHELL9
      amzero(&(actdsurf->neum->neum_onoff));
      amzero(&(actdsurf->neum->neum_val));
   #endif /*D_SHELL9*/

   /* NOTE: number of read values = MAXCONDPERNODE  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   for (i=0; i<MAXCONDPERNODE; i++)
   actdsurf->neum->neum_onoff.a.iv[i] = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   actdsurf->neum->neum_val.a.dv[i] = strtod(colptr,&colptr);
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
   frchk("PredescrDomainLoad",&ierr);   
   if (ierr) actdsurf->neum->neum_type = pres_domain_load;
   frchk("constHydro_z",&ierr);
   if (ierr) actdsurf->neum->neum_type = neum_consthydro_z;
   frchk("increaseHydro_z",&ierr);
   if (ierr) actdsurf->neum->neum_type = neum_increhydro_z;
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/
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
int    i,j;
int    ierr;
int    ndvol;
int    dvolId;
int    foundit;
char  *colptr;
char   buffer[200];
DVOL  *actdvol;
#ifdef DEBUG 
dstrc_enter("inpdesign_vol_neum");
#endif
/*----------------------------------------------------------------------*/
/*------------------ find the beginning of surface neumann conditions */
frfind("--DESIGN VOL NEUMANN CONDITIONS");
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
   /*--- now read the MAXCONDPERNODE flags for the neumann conditions */ 
   /*-------------------------- and the MAXCONDPERNODE values of them */
   amdef("onoff",&(actdvol->neum->neum_onoff),MAXDOFPERNODE,1,"IV");
   amdef("val",&(actdvol->neum->neum_val),MAXDOFPERNODE,1,"DV");

   /* Initialize some arrays to ZERO if MAXDOFPERNODE != MAXCONDPERNODE -> shell9 */
   #ifdef D_SHELL9
      amzero(&(actdvol->neum->neum_onoff));
      amzero(&(actdvol->neum->neum_val));
   #endif /*D_SHELL9*/

   /* NOTE: number of read values = MAXCONDPERNODE  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   for (i=0; i<MAXCONDPERNODE; i++)
   actdvol->neum->neum_onoff.a.iv[i] = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   actdvol->neum->neum_val.a.dv[i] = strtod(colptr,&colptr);
   /*----------------------------------- read type of neumann condition */
   frchk("Dead",&ierr);
   if (ierr) actdvol->neum->neum_type = neum_dead;   
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/
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
int    i,j;
int    ierr;
int    ndnode;
int    dnodeId;
int    foundit;
char  *colptr;
char   buffer[200];
DNODE *actdnode;
int    coupleId;
int    dofcouple;
int    dofflags[MAXDOFPERNODE];
int    geocouple;
int    geoflags[MAXDOFPERNODE];
#ifdef DEBUG 
dstrc_enter("inpdesign_nodal_couple");
#endif
/*----------------------------------------------------------------------*/
/*-------------------- find the beginning of nodal coupling conditions */
frfind("--DESIGN COUPLING POINT CONDITIONS");
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
   dsassert(coupleId>0,"Cannot reading coupling conditions");
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
   /* NOTE: number of read values = MAXCONDPERNODE  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   dofcouple = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   dofflags[i] = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   {
      if (dofflags[i]!=0) 
      actdnode->couple->couple.a.ia[i][1]=dofcouple;
   }
   /*------------------------------------------ read the geocouple flag */
   geocouple = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   geoflags[i] = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   {
      if (geoflags[i]!=0) 
      actdnode->couple->couple.a.ia[i][0]=geocouple;
   }
   /*------------------------------------------------ read the fsi flag */   
   actdnode->couple->fsi_iscoupled = strtol(colptr,&colptr,10);
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/
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
int    i,j;
int    ierr;
int    ndline;
int    dlineId;
int    foundit;
char  *colptr;
char   buffer[200];
DLINE *actdline;
int    coupleId;
int    dofcouple;
int    dofflags[MAXDOFPERNODE];
int    geocouple;
int    geoflags[MAXDOFPERNODE];
#ifdef DEBUG 
dstrc_enter("inpdesign_line_couple");
#endif
/*----------------------------------------------------------------------*/
/*--------------------- find the beginning of line coupling conditions */
frfind("--DESIGN COUPLING LINE CONDITIONS");
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
   /* NOTE: number of read values = MAXCONDPERNODE  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   dofcouple = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   dofflags[i] = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   {
      if (dofflags[i]!=0) 
      actdline->couple->couple.a.ia[i][1]=dofcouple;
   }
   /*------------------------------------------ read the geocouple flag */
   geocouple = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   geoflags[i] = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   {
      if (geoflags[i]!=0) 
      actdline->couple->couple.a.ia[i][0]=geocouple;
   }
   /*------------------------------------------------ read the fsi flag */   
   actdline->couple->fsi_iscoupled = strtol(colptr,&colptr,10);
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/
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
int    i,j;
int    ierr;
int    ndsurf;
int    dsurfId;
int    foundit;
char  *colptr;
char   buffer[200];
DSURF *actdsurf;
int    coupleId;
int    dofcouple;
int    dofflags[MAXDOFPERNODE];
int    geocouple;
int    geoflags[MAXDOFPERNODE];
#ifdef DEBUG 
dstrc_enter("inpdesign_surf_couple");
#endif
/*----------------------------------------------------------------------*/
/*-------------------- find the beginning of surf coupling conditions */
frfind("--DESIGN COUPLING SURF CONDITIONS");
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
   /* NOTE: number of read values = MAXCONDPERNODE  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   dofcouple = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   dofflags[i] = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   {
      if (dofflags[i]!=0) 
      actdsurf->couple->couple.a.ia[i][1]=dofcouple;
   }
   /*------------------------------------------ read the geocouple flag */
   geocouple = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   geoflags[i] = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   {
      if (geoflags[i]!=0) 
      actdsurf->couple->couple.a.ia[i][0]=geocouple;
   }
   /*------------------------------------------------ read the fsi flag */   
   actdsurf->couple->fsi_iscoupled = strtol(colptr,&colptr,10);
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/
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
int    i,j;
int    ierr;
int    ndvol;
int    dvolId;
int    foundit;
char  *colptr;
char   buffer[200];
DVOL *actdvol;
int    coupleId;
int    dofcouple;
int    dofflags[MAXDOFPERNODE];
int    geocouple;
int    geoflags[MAXDOFPERNODE];
#ifdef DEBUG 
dstrc_enter("inpdesign_vol_couple");
#endif
/*----------------------------------------------------------------------*/
/*-------------------- find the beginning of volume coupling conditions */
frfind("--DESIGN COUPLING VOL CONDITIONS");
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
   /* NOTE: number of read values = MAXCONDPERNODE  does not need to be */
   /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
   dofcouple = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   dofflags[i] = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   {
      if (dofflags[i]!=0) 
      actdvol->couple->couple.a.ia[i][1]=dofcouple;
   }
   /*------------------------------------------ read the geocouple flag */
   geocouple = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   geoflags[i] = strtol(colptr,&colptr,10);
   for (i=0; i<MAXCONDPERNODE; i++)
   {
      if (geoflags[i]!=0) 
      actdvol->couple->couple.a.ia[i][0]=geocouple;
   }
   /*------------------------------------------------ read the fsi flag */   
   actdvol->couple->fsi_iscoupled = strtol(colptr,&colptr,10);
   /*--------------------------------------------------- read next line */
   frread();
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_vol_couple */
