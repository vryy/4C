/*----------------------------------------------------------------------*
 | variables needed by design                             m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _DESIGN
{
   int ndnode;                /* number of design nodes */
   struct _DNODE *dnode;      /* vector of design nodes */

   int ndline;                /* number of design lines */
   struct _DLINE *dline;      /* vector of design lines */

   int ndsurf;                /* number of design surfaces */
   struct _DSURF *dsurf;      /* vector of design surfaces */

   int ndvol;                 /* number of design volumes */
   struct _DVOL *dvol;        /* vector of design volumes */
} DESIGN;


/*----------------------------------------------------------------------*
 | design node                                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _DNODE
{
   int                Id;           /* Id of this design node */
   int                ncond;        /* number of conditions associated with me (not used at the moment) */
   double             x[3];         /* the coordinates */
   /* design topology section */
   int                ndline;       /* number of lines connected to me */
   struct _DLINE    **dline;        /* vector of pointers to these line */
   /*     fe topology section */
   int                mynode;       /* global Id of my FE-node */
   struct _NODE      *node;         /* ptr to my FE-node */
   struct _FIELD     *field;        /* ptr to FIELD this design belongs to */
} DNODE;


/*----------------------------------------------------------------------*
 | design line                                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _DLINE
{
   int                 Id;            /* Id of this design line */
   int                 isnurb;        /* flag to show whether it's a standard or a nurb line */
   int                 ncond;         /* number of conditions associated with me (not used at the moment) */
   /* design topology section */
   int                 my_dnodeId[2]; /* IDs of design nodes to me */
   struct _DNODE     **dnode;         /* vector of pointers to these design nodes */
   int                 ndsurf;        /* number of surfaces connected to me */
   struct _DSURF     **dsurf;         /* vector of pointers to these surfaces */
   /*     fe topology section */
   ARRAY               mynode;        /* vector of global Ids of FE-nodes */
   struct _NODE      **node;          /* vector of ptrs to FE-nodes */
   struct _FIELD      *field;         /* ptr to FIELD this design belongs to */
} DLINE;


/*----------------------------------------------------------------------*
 | design surf                                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _DSURF
{
   int              Id;          /* Id of this design surface */
   int              ncond;       /* number of conditions associated with me (not used at the moment) */
   /* design topology section */
   ARRAY            my_dlineId;  /* Id's and orientation of the design lines to me */
                                 /* first row: Ids of the design lines of this surface */
                                 /* scnd row:  orientation: 0=same orientation as first dline, 1=opposite */
   struct _DLINE  **dline;       /* vector of pointers to my dlines */
   int              ndvol;       /* number of volumes connected to me */
   struct _DVOL   **dvol;        /* vector of pointers to these volumes */
   /*     fe topology section */
   ARRAY            mynode;      /* vector of global Ids of FE-nodes */
   struct _NODE   **node;        /* vector of ptrs to FE-nodes */
   struct _FIELD   *field;       /* ptr to FIELD this design belongs to */
} DSURF;


/*----------------------------------------------------------------------*
 | design vol                                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _DVOL
{
   int              Id;          /* Id of this design volume */
   int              ncond;       /* number of conditions associated with me (not used at the moment) */
   /* design topology section */
   ARRAY            my_dsurfId;  /* Id's and orientation of the design surfaces to me */
                                 /* first row: Ids of the design surfs of this volume */
                                 /* scnd row:  orientation: 0=same orientation as first dsurf, 1=opposite */
   struct _DSURF  **dsurf;       /* vector of pointers to these surfaces */
   /*     fe topology section */
   ARRAY            mynode;      /* vector of global Ids of FE-nodes */
   struct _NODE   **node;        /* vector of ptrs to FE-nodes */
   struct _FIELD   *field;       /* ptr to FIELD this design belongs to */
} DVOL;
