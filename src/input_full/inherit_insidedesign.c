#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;
/*----------------------------------------------------------------------*
 | functions only accessible in this file                    m.gee 3/02 |
 *----------------------------------------------------------------------*/
static void inherit_dvoldsurf_dirich(void);
static void inherit_dsurfdline_dirich(void);
static void inherit_dlinednode_dirich(void);
static void inherit_dvoldsurf_couple(void);
static void inherit_dsurfdline_couple(void);
static void inherit_dlinednode_couple(void);
/*----------------------------------------------------------------------*
 | inherit boundary conditions inside design                 m.gee 3/02 |
 *----------------------------------------------------------------------*/
void inherit_dirich_coup_indesign()
{
#ifdef DEBUG 
dstrc_enter("inherit_dirich_coup_indesign");
#endif
/*----------------------------------------------------------------------*/
/* 
dirichlet conditions are inherited as follows:
DVOL inherits to its DSURFS if the DSURF does not have its own 
DSURF inherits to its DLINEs if the DLINE does not have its own 
DLINE inherits to its DNODEs if the DNODE does not have its own 
*/
inherit_dvoldsurf_dirich();
inherit_dsurfdline_dirich();
inherit_dlinednode_dirich();
/* 
coupling conditions are inherited as follows:
DVOL inherits to its DSURFS if the DSURF does not have its own 
DSURF inherits to its DLINEs if the DLINE does not have its own 
DLINE inherits to its DNODEs if the DNODE does not have its own 
*/
inherit_dvoldsurf_couple();
inherit_dsurfdline_couple();
inherit_dlinednode_couple();
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inherit_dirich_coup_indesign */



/*----------------------------------------------------------------------*
 | inherit boundary conditions DVOL to DSURF                 m.gee 3/02 |
 *----------------------------------------------------------------------*/
static void inherit_dvoldsurf_dirich()
{
int             i,j;
DVOL           *actdvol;
DSURF          *actdsurf;
#ifdef DEBUG 
dstrc_enter("inherit_dvoldsurf_dirich");
#endif
/* 
dirichlet conditions are inherited as follows:
DVOL inherits to its DSURFS if the DSURF does not have its own 
*/
/*----------------------------------------------------------------------*/
/*---------------------------------------------------------- loop DVOLs */
for (i=0; i<design->ndvol; i++)
{
   actdvol = &(design->dvol[i]);
   /*-------------------- do nothing if DVOL has no dirichlet condition */
   if (actdvol->dirich==NULL) continue;
   /*------------------------------- loop the DSURFs related to actdvol */
   for (j=0; j<actdvol->ndsurf; j++)
   {
      actdsurf = actdvol->dsurf[j];
      /*-------- if actdsurf has its own dirichlet condition do nothing */
      if (actdsurf->dirich != NULL) continue;
      /*------ inherit the dirichlet condition from actdvol to actdsurf */
      actdsurf->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
      if (!actdsurf->dirich) dserror("Allocation of memory failed");
/*      actdsurf->dirich->curve = actdvol->dirich->curve;*/
      am_alloc_copy(&(actdvol->dirich->dirich_onoff),&(actdsurf->dirich->dirich_onoff));  
      am_alloc_copy(&(actdvol->dirich->dirich_val),&(actdsurf->dirich->dirich_val));  
      am_alloc_copy(&(actdvol->dirich->curve),&(actdsurf->dirich->curve));
   }/* loop j over dsurfs */
}/* loop i over dvols */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inherit_dvoldsurf_dirich */
/*----------------------------------------------------------------------*
 | inherit boundary conditions DSURF to DLINE                m.gee 3/02 |
 *----------------------------------------------------------------------*/
static void inherit_dsurfdline_dirich()
{
int             i,j;
DSURF          *actdsurf;
DLINE          *actdline;
#ifdef DEBUG 
dstrc_enter("inherit_dsurfdline_dirich");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------------------------------------- loop DSURFs */
for (i=0; i<design->ndsurf; i++)
{
   actdsurf = &(design->dsurf[i]);
   /*-------------------- do nothing if DSURF has no dirichlet condition */
   if (actdsurf->dirich==NULL) continue;
   /*------------------------------- loop the DLINEs related to actsurf */
   for (j=0; j<actdsurf->ndline; j++)
   {
      actdline = actdsurf->dline[j];
      /*-------- if actdline has its own dirichlet condition do nothing */
      if (actdline->dirich != NULL) continue;
      /*------ inherit the dirichlet condition from actdsurf to actdline */
      actdline->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
      if (!actdline->dirich) dserror("Allocation of memory failed");
/*      actdline->dirich->curve = actdsurf->dirich->curve;*/
      am_alloc_copy(&(actdsurf->dirich->dirich_onoff),&(actdline->dirich->dirich_onoff));  
      am_alloc_copy(&(actdsurf->dirich->dirich_val),&(actdline->dirich->dirich_val));  
      am_alloc_copy(&(actdsurf->dirich->curve),&(actdline->dirich->curve));
   }/* loop j over dlines */
}/* loop i over dsurfs */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inherit_dsurfdline_dirich */
/*----------------------------------------------------------------------*
 | inherit boundary conditions DLINE to DNODE                m.gee 3/02 |
 *----------------------------------------------------------------------*/
static void inherit_dlinednode_dirich()
{
int             i,j;
DNODE          *actdnode;
DLINE          *actdline;
#ifdef DEBUG 
dstrc_enter("inherit_dlinednode_dirich");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------------------------------------- loop DLINE */
for (i=0; i<design->ndline; i++)
{
   actdline = &(design->dline[i]);
   /*-------------------- do nothing if DLINE has no dirichlet condition */
   if (actdline->dirich==NULL) continue;
   /*------------------------------- loop the DNODEs related to actdline */
   for (j=0; j<actdline->ndnode; j++)
   {
      actdnode = actdline->dnode[j];
      /*-------- if actdnode has its own dirichlet condition do nothing */
      if (actdnode->dirich != NULL) continue;
      /*------ inherit the dirichlet condition from actdline to actdnode */
      actdnode->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
      if (!actdnode->dirich) dserror("Allocation of memory failed");
/*      actdnode->dirich->curve = actdline->dirich->curve;*/
      am_alloc_copy(&(actdline->dirich->dirich_onoff),&(actdnode->dirich->dirich_onoff));  
      am_alloc_copy(&(actdline->dirich->dirich_val),&(actdnode->dirich->dirich_val));  
      am_alloc_copy(&(actdline->dirich->curve),&(actdnode->dirich->curve));
   }/* loop j over dnodes */
}/* loop i over dlines */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inherit_dlinednode_dirich */



/*----------------------------------------------------------------------*
 | inherit coupling conditions DVOL to DSURF                 m.gee 3/02 |
 *----------------------------------------------------------------------*/
static void inherit_dvoldsurf_couple()
{
int             i,j;
DVOL           *actdvol;
DSURF          *actdsurf;
#ifdef DEBUG 
dstrc_enter("inherit_dvoldsurf_couple");
#endif
/* 
couple conditions are inherited as follows:
DVOL inherits to its DSURFS if the DSURF does not have its own 
*/
/*----------------------------------------------------------------------*/
/*---------------------------------------------------------- loop DVOLs */
for (i=0; i<design->ndvol; i++)
{
   actdvol = &(design->dvol[i]);
   /*-------------------- do nothing if DVOL has no couple condition */
   if (actdvol->couple==NULL) continue;
   /*------------------------------- loop the DSURFs related to actdvol */
   for (j=0; j<actdvol->ndsurf; j++)
   {
      actdsurf = actdvol->dsurf[j];
      /*-------- if actdsurf has its own couple condition do nothing */
      if (actdsurf->couple != NULL) continue;
      /*------ inherit the couple condition from actdvol to actdsurf */
      actdsurf->couple = (COUPLE_CONDITION*)CCACALLOC(1,sizeof(COUPLE_CONDITION));
      if (!actdsurf->couple) dserror("Allocation of memory failed");
      actdsurf->couple->fieldtyp      = actdvol->couple->fieldtyp;
      actdsurf->couple->fsi_iscoupled = actdvol->couple->fsi_iscoupled;
      am_alloc_copy(&(actdvol->couple->couple),&(actdsurf->couple->couple));  
   }/* loop j over dsurfs */
}/* loop i over dvols */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inherit_dvoldsurf_couple */
/*----------------------------------------------------------------------*
 | inherit coupling conditions DSURF to DLINE                m.gee 3/02 |
 *----------------------------------------------------------------------*/
static void inherit_dsurfdline_couple()
{
int             i,j;
DSURF          *actdsurf;
DLINE          *actdline;
#ifdef DEBUG 
dstrc_enter("inherit_dsurfdline_couple");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------------------------------------- loop DSURFs */
for (i=0; i<design->ndsurf; i++)
{
   actdsurf = &(design->dsurf[i]);
   /*-------------------- do nothing if DSURF has no couple condition */
   if (actdsurf->couple==NULL) continue;
   /*------------------------------- loop the DLINEs related to actsurf */
   for (j=0; j<actdsurf->ndline; j++)
   {
      actdline = actdsurf->dline[j];
      /*-------- if actdline has its own couple condition do nothing */
      if (actdline->couple != NULL) continue;
      /*------ inherit the couple condition from actdsurf to actdline */
      actdline->couple = (COUPLE_CONDITION*)CCACALLOC(1,sizeof(COUPLE_CONDITION));
      if (!actdline->couple) dserror("Allocation of memory failed");
      actdline->couple->fieldtyp      = actdsurf->couple->fieldtyp;
      actdline->couple->fsi_iscoupled = actdsurf->couple->fsi_iscoupled;
      am_alloc_copy(&(actdsurf->couple->couple),&(actdline->couple->couple));  
   }/* loop j over dlines */
}/* loop i over dsurfs */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inherit_dsurfdline_couple */
/*----------------------------------------------------------------------*
 | inherit coupling conditions DLINE to DNODE                m.gee 3/02 |
 *----------------------------------------------------------------------*/
static void inherit_dlinednode_couple()
{
int             i,j;
DNODE          *actdnode;
DLINE          *actdline;
#ifdef DEBUG 
dstrc_enter("inherit_dlinednode_couple");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------------------------------------- loop DLINE */
for (i=0; i<design->ndline; i++)
{
   actdline = &(design->dline[i]);
   /*-------------------- do nothing if DLINE has no couple condition */
   if (actdline->couple==NULL) continue;
   /*------------------------------- loop the DNODEs related to actdline */
   for (j=0; j<actdline->ndnode; j++)
   {
      actdnode = actdline->dnode[j];
      /*-------- if actdnode has its own couple condition do nothing */
      if (actdnode->couple != NULL) continue;
      /*------ inherit the couple condition from actdline to actdnode */
      actdnode->couple = (COUPLE_CONDITION*)CCACALLOC(1,sizeof(COUPLE_CONDITION));
      if (!actdnode->couple) dserror("Allocation of memory failed");
      actdnode->couple->fieldtyp       = actdline->couple->fieldtyp;
       actdnode->couple->fsi_iscoupled = actdline->couple->fsi_iscoupled;
     am_alloc_copy(&(actdline->couple->couple),&(actdnode->couple->couple));  
   }/* loop j over dnodes */
}/* loop i over dlines */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inherit_dlinednode_couple */
