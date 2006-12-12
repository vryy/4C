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



   INT *ndnode_fenode;
   INT **dnode_fenode;
   INT **dnode_fenode2;

   INT *ndline_fenode;
   INT **dline_fenode;
   INT **dline_fenode2;

   INT *ndsurf_fenode;
   INT **dsurf_fenode;
   INT **dsurf_fenode2;

   INT *ndvol_fenode;
   INT **dvol_fenode;
   INT **dvol_fenode2;

#ifdef CCADISCRET    /* new discretization management module */
   void* ccadesign;  /* is of type RefCountPtr<DRT::Design>* */
#endif

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
   struct _DIRICH_CONDITION     *ale_dirich;   /* dirichlet conditions for ale dis */
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
#ifdef D_SSI
   /*enum _SSI_COUPTYP         ssi_couptyp;*/
   struct _SSI_COUPLE_CONDITION *ssicouple;
#endif
   INT                     locsysId;
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
     struct _NURBLINE       *nurbline;      /* nurbline */
   }                          props;        /* name of union */

   /*------- design topology section */
   INT                       my_dnodeId[2]; /* IDs of design nodes to me */

   INT                       ndnode;        /* number of design nodes to this line (=2)*/
  /*struct _DNODE           **dnode;*/
   struct _DNODE            *dnode[2];      /* vector of pointers to these design nodes */

   INT                       ndsurf;        /* number of surfaces connected to me */
   struct _DSURF           **dsurf;         /* vector of pointers to these surfaces */

   /*------------fe topology section */
   /* boundary and coupling conditions */
   struct _NEUM_CONDITION   *neum;          /* neumann conditions to this DLINE, else NULL */
   struct _DIRICH_CONDITION *dirich;        /* dirichlet conditions to this DLINE, else NULL */
   struct _COUPLE_CONDITION *couple;        /* coupling conditions to this DLINE, else NULL */

#ifdef D_FSI
   struct _FSI_COUPLE_CONDITION *fsicouple;
   struct _SLIPDIRICH_CONDITION *slipdirich;
   struct _DIRICH_CONDITION     *ale_dirich;   /* dirichlet conditions for ale dis */
#endif

#ifdef D_FLUID
   struct _FLUID_FREESURF_CONDITION *freesurf;
   struct _FLUID_LIFTDRAG_CONDITION *liftdrag;
#endif

#ifdef D_AXISHELL
   struct _SAXI_THICK_CONDITION   *thickness;
   struct _SAXI_LOAD_CONDITION    *axishellload;
#endif

#ifdef WALLCONTACT
   enum   _CONTACTTYPE       contype;
#endif

#ifdef D_SSI
   struct _SSI_COUPLE_CONDITION *ssicouple;
#endif

   INT                     locsysId;
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
   DOUBLE                     trans_matrix[4][4];/* transformation matrix */

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


/*!----------------------------------------------------------------------
\brief structure NURBLINE

<pre>                                                              mn 05/03
This structure contains all geometric properties of a straight line.
</pre>

*----------------------------------------------------------------------*/
typedef struct _NURBLINE
{
  INT                        degree;
  INT                        num_cp;
  DOUBLE                   **cp;
  INT                        num_knots;
  DOUBLE                    *knots;
  INT                        rational;
  DOUBLE                    *weights;

  DOUBLE                     total_length;      /* total length of the line */

} NURBLINE;


/*----------------------------------------------------------------------*
 | design surf                                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _DSURF
{
   INT                       Id;          /* Id of this design surface */
   enum _DSURF_TYP           typ;         /* typ of the dline (st,nurb,arc) */
   INT                       ncond;       /* number of conditions associated with me (not used at the moment) */

   union                                    /* union pointer to surf properties */
   {
     struct _NURBSURF       *nurbsurf;      /* nurbsurf */
   }                         props;         /* name of union */


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
#ifdef D_FSI
   struct _FSI_COUPLE_CONDITION *fsicouple;
   struct _DIRICH_CONDITION     *ale_dirich;   /* dirichlet conditions for ale dis */
#endif
#ifdef D_FLUID
   struct _FLUID_FREESURF_CONDITION *freesurf;
   struct _FLUID_LIFTDRAG_CONDITION *liftdrag;
   /* fluid2 stabilisation is condition assigned to dsurface */
   enum _STABILISATION_TYP	stab_type;	/* enum of stabilisation	*/
   union
   {
     struct _STAB_PAR_GLS  *gls;/*! pointer to stabilisation parameters	*/
   /*  struct _STAB_PRES_PRO *pp; */
   } stabi;
#endif
#ifdef SURFACE_ENERGY
  INT     interface_flag;          /* =1 if dsurf=surface */
  INT     surface_energy_flag;     /* =0 for surfactant,
                             * =1 for const. surface tension */
  DOUBLE  const_gamma;      /* for constant surface tension only */
  DOUBLE  k1;               /* adsorption coefficient (surfactant) */
  DOUBLE  k2;               /* desorption coefficient (surfactant) */
  DOUBLE  C;                /* bulk concentration of surfactant */
  DOUBLE  m1;               /* isotherm slope regime 1 (surfactant) */
  DOUBLE  m2;               /* isotherm slope regime 2 (surfactant) */
  DOUBLE  gamma_0;
  DOUBLE  gamma_min;        /* minimal surface stress (surfactant) */
  DOUBLE  gamma_min_eq ;    /* minimum equilibrium surface stress
                             * (surfactant) */
#endif
   INT                     locsysId;
} DSURF;




/*!----------------------------------------------------------------------
\brief structure NURBSURF

<pre>                                                              mn 05/03
This structure contains all geometric properties of a nurb surface.
</pre>

*----------------------------------------------------------------------*/
typedef struct _NURBSURF
{
  INT                        degree[2];
  INT                        num_cp[2];
  DOUBLE                  ***cp;
  INT                        num_knots[2];
  DOUBLE                    *knots[2];
  INT                        rational;
  DOUBLE                   **weights;
  INT                        istrimmed;
  DOUBLE                     center[3];
  DOUBLE                     normal[3];

} NURBSURF;



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

#ifdef D_FLUID /* fluid3 stabilisation is condition assigned to dsurface */
   enum _STABILISATION_TYP	stab_type;	/* enum of stabilisation	*/
   union
   {
     struct _STAB_PAR_GLS  *gls;/*! pointer to stabilisation parameters	*/
   /*  struct _STAB_PRES_PRO *pp; */
   } stabi;
#endif

#ifdef D_FSI
   struct _DIRICH_CONDITION     *ale_dirich;   /* dirichlet conditions for ale dis */
#endif

   INT                     locsysId;
} DVOL;
