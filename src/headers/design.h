/*!---------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | variables needed by design                             m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _DESIGN
{
   INT ndnode;                /* number of design nodes */
   struct _DNODE *dnode;      /* vector of design nodes */

   INT ndline;                /* number of design lines */
   struct _DLINE *dline;      /* vector of design lines */

   INT ndsurf;                /* number of design surfaces */
   struct _DSURF *dsurf;      /* vector of design surfaces */

   INT ndvol;                 /* number of design volumes */
   struct _DVOL *dvol;        /* vector of design volumes */
} DESIGN;


/*----------------------------------------------------------------------*
 | design node                                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _DNODE
{
   INT                       Id;           /* Id of this design node */
   INT                       ncond;        /* number of conditions associated with me (not used at the moment) */
   DOUBLE                    x[3];         /* the coordinates */
   /*------- design topology section */
   INT                       ndline;       /* number of lines connected to me */
   struct _DLINE           **dline;        /* vector of pointers to these line */

   /* boundary and coupling conditions */
   struct _NEUM_CONDITION   *neum;         /* neumann conditions to this DNODE, else NULL */
   struct _DIRICH_CONDITION *dirich;       /* dirichlet conditions to this DNODE, else NULL */
   struct _COUPLE_CONDITION *couple;       /* coupling conditions to this DNODE, else NULL */
#ifdef D_FSI
   struct _FSI_COUPLE_CONDITION *fsicouple;
#endif
#ifdef D_FLUID
   struct _FLUID_FREESURF_CONDITION *freesurf;
#endif
#ifdef D_AXISHELL
   struct _SAXI_THICK_CONDITION   *thickness;
   struct _SAXI_LOAD_CONDITION    *axishellload;
   INT                             cos_type;  /* type of cos used for this node:
                                                 0: global (default)
                                                 1: local */
#endif
} DNODE;


/*----------------------------------------------------------------------*
 | design line                                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _DLINE
{
   INT                       Id;            /* Id of this design line */
   enum _DLINE_TYP           typ   ;        /* typ of the dline (st,nurb,arc) */
   INT                       ncond;         /* number of conditions associated with me (not used at the moment) */

   union                                    /* union pointer to arcline properties */
   {
     struct _ARCLINE        *arcline;       /* arcline */
     struct _STLINE         *stline;        /* stline */
   }                          props;        /* name of union */ 

   /*------- design topology section */
   INT                       my_dnodeId[2]; /* IDs of design nodes to me */

   INT                       ndnode;        /* number of design nodes to this line (=2)*/
   struct _DNODE           **dnode;         /* vector of pointers to these design nodes */

   INT                       ndsurf;        /* number of surfaces connected to me */
   struct _DSURF           **dsurf;         /* vector of pointers to these surfaces */

   /*------------fe topology section */
   /* boundary and coupling conditions */
   struct _NEUM_CONDITION   *neum;          /* neumann conditions to this DLINE, else NULL */
   struct _DIRICH_CONDITION *dirich;        /* dirichlet conditions to this DLINE, else NULL */
   struct _COUPLE_CONDITION *couple;        /* coupling conditions to this DLINE, else NULL */
#ifdef D_FSI
   struct _FSI_COUPLE_CONDITION *fsicouple;
#endif
#ifdef D_FLUID
   struct _FLUID_FREESURF_CONDITION *freesurf;
   INT                       liftdrag;
#endif
#ifdef D_AXISHELL
   struct _SAXI_THICK_CONDITION   *thickness;
   struct _SAXI_LOAD_CONDITION    *axishellload;
#endif
#ifdef WALLCONTACT
   enum   _CONTACTTYPE       contype;
#endif
} DLINE;


/*!----------------------------------------------------------------------
\brief structure ARCLINE                                            

<pre>                                                              mn 05/03
This structure contains all geometric properties of an arcline. 
</pre>

*----------------------------------------------------------------------*/
typedef struct _ARCLINE
{
   DOUBLE                     radius;            /* radius of the arc */
   DOUBLE                     initang;           /* angle in radians */
   DOUBLE                     endang;            /* angle in radians */
   DOUBLE                     center[2];         /* coords of the center point */
   DOUBLE                     total_length;      /* total arc length */

} ARCLINE;

/*!----------------------------------------------------------------------
\brief structure STLINE                                            

<pre>                                                              mn 05/03
This structure contains all geometric properties of a straight line. 
</pre>

*----------------------------------------------------------------------*/
typedef struct _STLINE
{
   DOUBLE                     total_length;      /* total length of the line */

} STLINE;


/*----------------------------------------------------------------------*
 | design surf                                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _DSURF
{
   INT                       Id;          /* Id of this design surface */
   INT                       ncond;       /* number of conditions associated with me (not used at the moment) */
   /*------- design topology section */
   ARRAY                     my_dlineId;  /* Id's and orientation of the design lines to me */
                                          /* first row: Ids of the design lines of this surface */
                                          /* scnd row:  orientation: 0=same orientation as first dline, 1=opposite */
   INT                       ndline;      /* number of DLINEs to me */
   struct _DLINE           **dline;       /* vector of pointers to my dlines */

   INT                       ndvol;       /* number of volumes connected to me */
   struct _DVOL            **dvol;        /* vector of pointers to these volumes */

   /*------------fe topology section */
   /* boundary and coupling conditions */
   struct _NEUM_CONDITION   *neum;        /* neumann conditions to this DSURF, else NULL */
   struct _DIRICH_CONDITION *dirich;      /* dirichlet conditions to this DSURF, else NULL */
   struct _COUPLE_CONDITION *couple;      /* coupling conditions to this DSURF, else NULL */
} DSURF;


/*----------------------------------------------------------------------*
 | design vol                                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _DVOL
{
   INT                       Id;          /* Id of this design volume */
   INT                       ncond;       /* number of conditions associated with me (not used at the moment) */
   /*------- design topology section */
   ARRAY                     my_dsurfId;  /* Id's and orientation of the design surfaces to me */
                                          /* first row: Ids of the design surfs of this volume */
                                          /* scnd row:  orientation: 0=same orientation as first dsurf, 1=opposite */
   INT                       ndsurf;      /* number of DSURFs to me */
   struct _DSURF           **dsurf;       /* vector of pointers to these surfaces */
   /*------------fe topology section */
   /* boundary and coupling conditions */
   struct _NEUM_CONDITION   *neum;        /* neumann conditions to this DVOL, else NULL */
   struct _DIRICH_CONDITION *dirich;      /* dirichlet conditions to this DVOL, else NULL */
   struct _COUPLE_CONDITION *couple;      /* coupling conditions to this DVOL, else NULL */
} DVOL;
