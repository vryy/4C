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

     struct _COND_NODE         *c;             /* my conditions, if any, else NULL */

     enum                                      /* type of design object that owns me */
       {
          not_owned
          dnode_owned,
          dline_owned,
          dsurf_owned,
          dvol_owned
       }                        downertyp;
     union                                     /* union pointer to design object that owns me */
       {
       struct _DNODE  *dnode;
       struct _DLINE  *dline;
       struct _DSURF  *dsurf;
       struct _DVOL   *dvol;
       }                        d;
       
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
     struct _WALL1    *w1;                      /* 2D plain stress - plain strain element */
     struct _FLUID3   *f3;                      /* 3D fluid element */
     struct _ALE      *ale;                     /* pseudo structural 2D or 3D ale element */
     }                          e;              /* name of union */ 
     struct _COND_ELEMENT      *c;              /* my conditions, if any, else NULL */
} ELEMENT;







