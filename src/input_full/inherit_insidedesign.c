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
#ifdef D_FSI
static void inherit_dlinednode_fsicouple(void);
#endif
static void inherit_dlinednode_freesurf(void);
#ifdef D_AXISHELL
static void inherit_dlinednode_thickness(void);
static void inherit_dlinednode_axishellload(void);
#endif
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
/* 
fsi coupling conditions are inherited as follows:
DLINE inherits to its DNODEs if the DNODE does not have its own 
*/
#ifdef D_FSI
inherit_dlinednode_fsicouple();
#endif
/* 
freesruface conditions are inherited as follows:
DLINE inherits to its DNODEs if the DNODE does not have its own 
*/
#ifdef D_FLUID
inherit_dlinednode_freesurf(); 
#endif
/*
 * thickness and axishellload conditions are inherited as follows:
 * DLINE inherits to its DNODEs if the DNODE does not have its own
 */
#ifdef D_AXISHELL
inherit_dlinednode_thickness(); 
inherit_dlinednode_axishellload(); 
#endif
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
INT             i,j;
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
      actdsurf->dirich->dirich_type = actdvol->dirich->dirich_type;
      am_alloc_copy(&(actdvol->dirich->dirich_onoff),&(actdsurf->dirich->dirich_onoff));  
      am_alloc_copy(&(actdvol->dirich->dirich_val),&(actdsurf->dirich->dirich_val));  
      am_alloc_copy(&(actdvol->dirich->curve),&(actdsurf->dirich->curve));
      am_alloc_copy(&(actdvol->dirich->funct),&(actdsurf->dirich->funct));
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
INT             i,j;
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
      am_alloc_copy(&(actdsurf->dirich->funct),&(actdline->dirich->funct));
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
INT             i,j;
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
/*      actdnode->dirich->curve = actdline->dirich->curve;*/
      am_alloc_copy(&(actdline->dirich->dirich_onoff),&(actdnode->dirich->dirich_onoff));  
      am_alloc_copy(&(actdline->dirich->dirich_val),&(actdnode->dirich->dirich_val));  
      am_alloc_copy(&(actdline->dirich->curve),&(actdnode->dirich->curve));
      am_alloc_copy(&(actdline->dirich->funct),&(actdnode->dirich->funct));
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
INT             i,j;
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
INT             i,j;
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
INT             i,j;
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
     am_alloc_copy(&(actdline->couple->couple),&(actdnode->couple->couple));  
   }/* loop j over dnodes */
}/* loop i over dlines */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inherit_dlinednode_couple */



#ifdef D_FSI
/*----------------------------------------------------------------------*
 | inherit fsicoupling conditions DLINE to DNODE             genk 10/02 |
 *----------------------------------------------------------------------*/
static void inherit_dlinednode_fsicouple()
{
INT             i,j;
DNODE          *actdnode;
DLINE          *actdline;
#ifdef DEBUG 
dstrc_enter("inherit_dlinednode_fsicouple");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------------------------------------- loop DLINE */
for (i=0; i<design->ndline; i++)
{
   actdline = &(design->dline[i]);
   /*-------------------- do nothing if DLINE has no dirichlet condition */
   if (actdline->fsicouple==NULL) continue;
   /*------------------------------- loop the DNODEs related to actdline */
   for (j=0; j<actdline->ndnode; j++)
   {
      actdnode = actdline->dnode[j];
      /*--------------- if actdnode has its own fsi condition do nothing */
      if (actdnode->fsicouple != NULL) continue;
      /*------ inherit the dirichlet condition from actdline to actdnode */
      actdnode->fsicouple = (FSI_COUPLE_CONDITION*)CCACALLOC(1,sizeof(FSI_COUPLE_CONDITION));
      if (!actdnode->fsicouple) dserror("Allocation of memory failed");
      actdnode->fsicouple->fieldtyp = actdline->fsicouple->fieldtyp;
      actdnode->fsicouple->fsi_coupleId = actdline->fsicouple->fsi_coupleId;
      actdnode->fsicouple->fsi_mesh = actdline->fsicouple->fsi_mesh;
   }/* loop j over dnodes */
}/* loop i over dlines */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inherit_dlinednode_fsicouple */
#endif


/*----------------------------------------------------------------------*
 | inherit freesurface conditions DLINE to DNODE             genk 01/03 |
 *----------------------------------------------------------------------*/
static void inherit_dlinednode_freesurf()
{
#ifdef D_FLUID
INT             i,j;
DNODE          *actdnode;
DLINE          *actdline;
#ifdef DEBUG 
dstrc_enter("inherit_dlinednode_freesurf");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------------------------------------- loop DLINE */
for (i=0; i<design->ndline; i++)
{
   actdline = &(design->dline[i]);
   /*-------------------- do nothing if DLINE has no dirichlet condition */
   if (actdline->freesurf==NULL) continue;
   /*------------------------------- loop the DNODEs related to actdline */
   for (j=0; j<actdline->ndnode; j++)
   {
      actdnode = actdline->dnode[j];
      /*--------- if actdnode has its own freesurf condition do nothing */
      if (actdnode->freesurf != NULL) continue;
      /*------ inherit the freesurf condition from actdline to actdnode */
      actdnode->freesurf = (FLUID_FREESURF_CONDITION*)CCACALLOC(1,sizeof(FLUID_FREESURF_CONDITION));
      if (!actdnode->freesurf) dserror("Allocation of memory failed");
      actdnode->freesurf->fieldtyp = actdline->freesurf->fieldtyp;
      am_alloc_copy(&(actdline->freesurf->fixed_onoff),&(actdnode->freesurf->fixed_onoff)); 
   }/* loop j over dnodes */
}/* loop i over dlines */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif
return;
} /* end of inherit_dlinednode_freesurf */




#ifdef D_AXISHELL
/*!----------------------------------------------------------------------
\brief inherit thickness conditions DLINE to DNODE

<pre>                                                              mn 05/03
inherit thickness conditions DLINE to DNODE
</pre>

\warning There is nothing special to this routine
\return void                                               
\sa

*----------------------------------------------------------------------*/
static void inherit_dlinednode_thickness()
{
INT             i,j;
DNODE          *actdnode;
DLINE          *actdline;
#ifdef DEBUG 
dstrc_enter("inherit_dlinednode_thickness");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------------------------------------- loop DLINE */
for (i=0; i<design->ndline; i++)
{
   actdline = &(design->dline[i]);
   /*-------------------- do nothing if DLINE has no thickness condition */
   if (actdline->thickness==NULL) continue;
   /*------------------------------- loop the DNODEs related to actdline */
   for (j=0; j<actdline->ndnode; j++)
   {
      actdnode = actdline->dnode[j];
      /*-------- if actdnode has its own thickness condition check value */
      if (actdnode->thickness != NULL) 
      {
        if (actdnode->thickness->value[0] !=  actdline->thickness->value[j])
          dserror("Thickness ust be equal at both sides of a dnode!");
      }
      /*------ inherit the dirichlet condition from actdline to actdnode */
      actdnode->thickness = (SAXI_THICK_CONDITION*)CCACALLOC(1,sizeof(SAXI_THICK_CONDITION));
      if (!actdnode->thickness) dserror("Allocation of memory failed");
      actdnode->thickness->value[0] = actdline->thickness->value[j];
   }/* loop j over dnodes */
}/* loop i over dlines */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inherit_dlinednode_thickness */
#endif




#ifdef D_AXISHELL 
/*!----------------------------------------------------------------------
\brief inherit axishell_load conditions DLINE to DNODE

<pre>                                                              mn 05/03
inherit axishell_load conditions DLINE to DNODE
</pre>

\warning There is nothing special to this routine
\return void                                               
\sa

*----------------------------------------------------------------------*/
static void inherit_dlinednode_axishellload()
{
INT             i,j;
DNODE          *actdnode;
DLINE          *actdline;
#ifdef DEBUG 
dstrc_enter("inherit_dlinednode_axishellload");
#endif
/*---------------------------------------------------------------------*/
/*--------------------------------------------------------- loop DLINE */
for (i=0; i<design->ndline; i++)
{
   actdline = &(design->dline[i]);
   /*-------------------- do nothing if DLINE has no axishell_load condition */
   if (actdline->axishellload==NULL) continue;
   /*----------------------------- loop the DNODEs related to actdline */
   for (j=0; j<actdline->ndnode; j++)
   {
      actdnode = actdline->dnode[j];
      /*-- if actdnode has its own axishell_load condition check value */
      if (actdnode->axishellload != NULL) 
      {
        if (actdnode->axishellload->pv[0] !=  actdline->axishellload->pv[j] ||
            actdnode->axishellload->ph[0] !=  actdline->axishellload->ph[j] ||
            actdnode->axishellload->px[0] !=  actdline->axishellload->px[j] ||
            actdnode->axishellload->pw[0] !=  actdline->axishellload->pw[j] )
          dserror("Loadvalue must be equal at both sides of a dnode!");
      }
      /*------ inherit the dirichlet condition from actdline to actdnode */
      actdnode->axishellload = (SAXI_LOAD_CONDITION*)CCACALLOC(1,sizeof(SAXI_LOAD_CONDITION));
      if (!actdnode->axishellload) dserror("Allocation of memory failed");
      actdnode->axishellload->pv[0] = actdline->axishellload->pv[j];
      actdnode->axishellload->ph[0] = actdline->axishellload->ph[j];
      actdnode->axishellload->px[0] = actdline->axishellload->px[j];
      actdnode->axishellload->pw[0] = actdline->axishellload->pw[j];
   }/* loop j over dnodes */
}/* loop i over dlines */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inherit_dlinednode_axishellload */
#endif



