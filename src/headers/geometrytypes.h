/*----------------------------------------------------------------------*
 | type definitions of geometry based information         m.gee 8/00    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | 1 NODE                                                 m.gee 8/00    |
 *----------------------------------------------------------------------*/
typedef struct _NODE
{
     int                        Id;            /* global Id (numbering starts with 0)*/              
     int                        Id_loc;        /* field-local Id  (numbering starts with 0)*/
     int                        proc;          /* my owner intra-proc */

     double                     x[3];          /* my coordinates */

     struct _ARRAY              sol;           /* my solution history */
     struct _ARRAY              sol_increment; /* my incremental solution */
     struct _ARRAY              sol_residual ; /* my residual solution */

     int                        numdf;         /* my number of degrees of freedom */
     int                       *dof;           /* my dof-numbers  */

     int                        numele;        /* number of elements to me */
     struct _ELEMENT          **element;       /* ptrs to elements to me */

     struct _GNODE             *gnode;         /* ptr to my gnode */
} NODE;


/*----------------------------------------------------------------------*
 | 1 ELEMENT                                              m.gee 8/00    |
 *----------------------------------------------------------------------*/
typedef struct _ELEMENT
{  
     int                        Id;             /* global Id (numbering starts with 0)*/                          
     int                        Id_loc;         /* field-local Id (numbering starts with 0)*/
     int                        proc;           /* my owner intra-proc */

     int                        numnp;          /* number of nodes to me */
     int                       *lm;             /* only used for reading from input (this will be eliminated)*/
     struct _NODE             **node;           /* ptrs to my nodes */

     int                        mat;            /* number of material law associated with me */

     enum _ELEMENT_TYP          eltyp;          /* my element type */
     enum _DIS_TYP              distyp;         /* my actual discretization type */

     union                                      /* union pointer to elementformulation */
     {
     struct _SHELL8   *s8;                      /* shell8 element */
     struct _BRICK1   *b1;                      /* structural volume element */
     struct _WALL1    *w1;                      /* 2D plane stress - plane strain element */
     struct _FLUID2   *f2;                      /* 2D fluid element */
     struct _FLUID3   *f3;                      /* 3D fluid element */
     struct _ALE      *ale;                     /* pseudo structural 2D or 3D ale element */
     }                          e;              /* name of union */ 

     union
     {
     struct _GSURF    *gsurf;                   /* my gsurf, if I am a 2D element */
     struct _GVOL     *gvol;                    /* my gvol,  if I am a 3D element */
     }                          g;              /* name of union */           
} ELEMENT;



/*----------------------------------------------------------------------*
 | 1 GNODE                                                 m.gee 3/02   |
 *----------------------------------------------------------------------*/
typedef struct _GNODE
{
#ifdef DEBUG 
     int                           Id;          /* for debugging only, do not use in code !*/
#endif
   /*------------fe topology section */
     struct _NODE                 *node;        /* pointer to my node */
     int                           ngline;      /* number of GLINEs to me */
     struct _GLINE               **gline;       /* pointers to the GLINEs to me */
   /*------- design topology section */
     enum
     {
        ondnothing,                             /* I am NOT positioned on any design object (error!) */
        ondnode,                                /* I am on a DNODE */
        ondline,                                /* I am on a DLINE */
        ondsurf,                                /* I am on a DSURF */
        ondvol                                  /* I am inside a DVOL */
     }                             ondesigntyp; 
     union
     {
        struct _DNODE      *dnode;              
        struct _DLINE      *dline;
        struct _DSURF      *dsurf;
        struct _DVOL       *dvol;
     }                             d;           /* ptr to the design object I am positioned on */
   /* boundary and coupling conditions */
     struct _DIRICH_CONDITION     *dirich;      /* a dirichlet condition on this gnode, else NULL */
     struct _COUPLE_CONDITION     *couple;      /* a coupling conditions on this gnode, else NULL */
     struct _NEUM_CONDITION       *neum;        /* a neumann condition on this gnode, else NULL */
} GNODE;


/*----------------------------------------------------------------------*
 | 1 GLINE                                                 m.gee 3/02   |
 *----------------------------------------------------------------------*/
typedef struct _GLINE
{
#ifdef DEBUG 
     int                        Id;             /* for debugging only, do not use in code !*/
#endif
   /*------------fe topology section */
     int                        ngnode;         /* number of gnodes on me */
     struct _GNODE            **gnode;          /* vector of ptrs to these gnodes */

     int                        ngsurf;         /* number of gsurfs to me */
     struct _GSURF            **gsurf;          /* vector of ptrs to these gsurfs */

   /*------- design topology section */
     struct _DLINE             *dline;          /* The DLINE I am on, else NULL */

   /*----------- boundary conditions */
     struct _NEUM_CONDITION       *neum;        /* neumann conditions to this GLINE, else NULL */
} GLINE;
/*----------------------------------------------------------------------*
 | 1 GSURF                                                 m.gee 3/02   |
 *----------------------------------------------------------------------*/
typedef struct _GSURF
{
#ifdef DEBUG 
     int                        Id;             /* for debugging only, do not use in code !*/
#endif
   /*------------fe topology section */
     struct _ELEMENT           *element;        /* ptr to my ELEMENT, if I am a 2D element, else NULL */

     int                        ngnode;         /* number of GNODEs to me */
     struct _GNODE            **gnode;          /* ptrs to these GNODEs */

     int                        ngline;         /* number of GLINEs to me */
     struct _GLINE            **gline;          /* ptrs to these GLINEs */

     int                        ngvol;          /* number of GVOLs to me */
     struct _GVOL             **gvol;           /* ptrs to these GVOLs, else NULL */

   /*------- design topology section */
     struct _DSURF             *dsurf;          /* DSURF I am on, else NULL */

   /*----------- boundary conditions */
     struct _NEUM_CONDITION       *neum;        /* neumann conditions to this GSURF, else NULL */
} GSURF;
/*----------------------------------------------------------------------*
 | 1 GVOL                                                  m.gee 3/02   |
 *----------------------------------------------------------------------*/
typedef struct _GVOL
{
#ifdef DEBUG 
     int                        Id;             /* for debugging only, do not use in code !*/
#endif
   /*------------fe topology section */
     struct _ELEMENT           *element;        /* ptr to my ELEMENT */

     int                        ngline;         /* number of GLINEs to me */
     struct _GLINE            **gline;          /* ptrs to these GLINEs */

     int                        ngsurf;         /* number of GSURFs to me */
     struct _GSURF            **gsurf;          /* ptrs to these GSURFs */

   /*------- design topology section */
     struct _DVOL              *dvol;           /* the DVOL I am placed in */

   /*----------- boundary conditions */
     struct _NEUM_CONDITION       *neum;        /* neumann conditions to this GVOL, else NULL */
} GVOL;




