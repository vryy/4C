/*----------------------------------------------------------------------*
 | type definitions conditions                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |  nodal conditions                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _COND_NODE
{
     int                       curve;            /* number of the curve to a time dependent condition */

     int                       isneum;           /* flag for neumann condition */
     struct _ARRAY             neum_onoff;       /* array holding flags for each dof of node */
     struct _ARRAY             neum_val;         /* array holding values for each dof of node */

     int                       isdirich;         /* flag for dirichlet condition on this node */
     struct _ARRAY             dirich_onoff;     /* array holding flags for each dof of node */
     struct _ARRAY             dirich_val;       /* aray holding values for each dof of node */

     int                       iscoupled;        /* flag whether this node is coupled to another node */                  
     struct _ARRAY             couple;           /* array holding flags foreach dof of node */

     int                       fsi_iscoupled;    /* flag for an fsi coupled node */                      

} COND_NODE;

/*----------------------------------------------------------------------*
 | element neumann/dirichlet condition types             m.gee 10/01    |
 *----------------------------------------------------------------------*/
typedef enum _COND_ELE_TYP
{
                       ne_live,   /* given neumann or dirichlet values */
                       ne_dead    /* neuman value is computed from the specific weight of element */
} COND_ELE_TYP;
/*----------------------------------------------------------------------*
 |  element conditions                                    m.gee 5/01    |
 *----------------------------------------------------------------------*/
typedef struct _COND_ELEMENT
{
     enum _COND_ELE_TYP        condtyp;      /* type of element condition */
     int                       curve;        /* number of time curve for time dependent condition */

     int                       isneum;       /* flag for neumann condition */
     struct _ARRAY             neum_onoff;   /* array holding flags for each dof of node */
     struct _ARRAY             neum_val;     /* array holding values for each dof of node */

} COND_ELEMENT;



