/*----------------------------------------------------------------------*
 | type definitions of geometry based information         m.gee 8/00    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | 1 NODE                                                 m.gee 8/00    |
 *----------------------------------------------------------------------*/
typedef struct _NODE
{
     int                        Id;            /* global Id */              
     int                        Id_loc;        /* field-local Id */
     int                        proc;          /* my master intra-proc */
     double                     x[3];          /* my coordinates */
     struct _ARRAY              sol;           /* my solution history */
     struct _ARRAY              sol_increment; /* my incremental solution */
     struct _ARRAY              sol_residual ; /* my residual solution */
     int                        numdf;         /* my numdf */
     int                       *dof;           /* my dof-numbers */

     int                        numele;        /* number of elements to me */
     struct _ELEMENT          **element;       /* ptrs to elements to me */

     struct _COND_NODE         *c;             /* my conditions, if any */

     struct _ARRAY              db_access;     /* my bunker access Ids (not used) */     
     
} NODE;


/*----------------------------------------------------------------------*
 | 1 ELEMENT                                              m.gee 8/00    |
 *----------------------------------------------------------------------*/
typedef struct _ELEMENT
{  
     int                        Id;             /* global Id */                          
     int                        Id_loc;         /* field-local Id */
     int                        proc;           /* my master intra-proc */
     int                        numnp;          /* number of nodes to me */
     int                       *lm;             /* only used for reading from input*/
     struct _NODE             **node;           /* ptrs to my nodes */
     int                        mat;
     
     enum _ELEMENT_TYP          eltyp;          /* my element type */
     enum _DIS_TYP              distyp;         /* my actual discretization type */

     union 
     {
     struct _SHELL8   *s8;
     struct _BRICK1   *b1;
     struct _FLUID3   *f3;
     struct _ALE      *ale;
     }                          e;              /* element specific data */ 
     struct _COND_ELEMENT      *c;



     struct _ARRAY              db_access;     /* my bunker access Ids (not used) */  

} ELEMENT;






