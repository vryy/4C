/*----------------------------------------------------------------------*
 | variables needed for parallel comp.                    m.gee 8/00    |
 *----------------------------------------------------------------------*/
typedef struct _PAR
{
int               myrank;                /* the individual processor number */
int               nprocs;                /* total number of processors */
#ifdef PARALLEL 
int               numfld;                /* number of intra-communicators */
struct _INTRA    *intra;                 /* vector of intra-communicator-structures */
#endif
} PAR;

/*----------------------------------------------------------------------*
 | one intra-communicator            .                    m.gee 9/01    |
 *----------------------------------------------------------------------*/
typedef struct _INTRA
{
enum   _FIELDTYP    intra_fieldtyp;      /* type of field */
int                 intra_rank;          /* proc's intra-rank */
int                 intra_nprocs;        /* number of procs in this intracomm. */
#ifdef PARALLEL 
MPI_Comm            MPI_INTRA_COMM;      /* the intra-communicator itself */
MPI_Group           MPI_INTRA_GROUP;     /* ? */
#endif
} INTRA;

/*----------------------------------------------------------------------*
 | one proc's info about his partition                    m.gee 8/00    |
 *----------------------------------------------------------------------*/
typedef struct _PARTITION
{
enum   _FIELDTYP    fieldtyp;            /* type of field */

int                 numnp;               /* total number of nodes on this parition */
int                 numele;              /* total number of elements on this partition */
struct _NODE      **node;                /* ptrs to nodes on this part. */
struct _ELEMENT   **element;             /* ptrs to elements on this partition */

int                 inner_numnp;         /* number of pure inner-nodes */
struct _NODE      **inner_node;          /* ptrs to pure inner-nodes */
int                 bou_numnp;           /* number of boundary nodes */
struct _NODE      **bou_node;            /* ptrs to boundary nodes */

int                 inner_numele;        /* number of pure inner elements */
struct _ELEMENT   **inner_element;       /* pts to pure inner elements */
int                 bou_numele;          /* number of boundary elements */
struct _ELEMENT   **bou_element;         /* ptrs to boundary elements */

struct _ARRAY       coupledofs;          /* number of coupled dofs, which dofs, */
                                         /* who is master owner and who is slave owner */ 
                                         /* for details, see global_mask_matrices.c */ 
struct _ARRAY       db_access;           /* bunker access handles (not used) */
} PARTITION;



