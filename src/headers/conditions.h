/*----------------------------------------------------------------------*
 | neumann condition                                      m.gee 3/02    |
 |                                                                      |
 | this structure holds a neumann condition                             |
 | depend on number of dofs to a fe-node and type of element it is      |
 | is connected to, the arrays can be defined in several styles         |
 *----------------------------------------------------------------------*/
typedef struct _NEUM_CONDITION
{
     int                       curve;        /* number of load curve associated with this conditions */        

     struct _ARRAY             neum_onoff;   /* array of on-off flags */
     struct _ARRAY             neum_val;     /* values of this condition */    

} NEUM_CONDITION;
/*----------------------------------------------------------------------*
 | dirichlet condition                                    m.gee 3/02    |
 |                                                                      |
 | this structure holds a dirichlet condition                           |
 | depend on number of dofs to a fe-node and type of element it is      |
 | is connected to, the arrays can be defined in several styles         |
 *----------------------------------------------------------------------*/
typedef struct _DIRICH_CONDITION
{
     int                       curve;        /* number of load curve associated with this conditions */                 

     struct _ARRAY             dirich_onoff; /* array of on-off flags */
     struct _ARRAY             dirich_val;   /* values of this condition */    


} DIRICH_CONDITION;
/*----------------------------------------------------------------------*
 | coupling condition                                     m.gee 3/02    |
 |                                                                      |
 | this structure is assigned to nodes, which are coupled in some or    |
 | all of their dofs                                                    |
 *----------------------------------------------------------------------*/
typedef struct _COUPLE_CONDITION
{
     enum _FIELDTYP            fieldtyp;        /* type of field this structure is in */

     struct _ARRAY             couple;          /* array of the coupling conditions */
     int                       fsi_iscoupled;   /* flag for fsi coupling */

} COUPLE_CONDITION;
