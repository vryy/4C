/*!---------------------------------------------------------------------
\file
\brief domain decomposition and metis structures

---------------------------------------------------------------------*/

/*! 
\addtogroup PARALLEL 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
typedef struct _PAR
{
int               myrank;                /*!< the individual processor number */
int               nprocs;                /*!< total number of processors */
#ifdef PARALLEL 
int               numfld;                /*!< number of intra-communicators == number of fields */
struct _INTRA    *intra;                 /*!< vector of intra-communicator-structures correspondent to vector of FIELDs */
#endif
} PAR;

/*!----------------------------------------------------------------------
\brief one intra-communicator

<pre>                                                         m.gee 9/01
-holds data associated with one intra-communicator
-the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
typedef struct _INTRA
{
enum   _FIELDTYP    intra_fieldtyp;      /*!< type of field */
int                 intra_rank;          /*!< proc's intra-rank */
int                 intra_nprocs;        /*!< number of procs in this intracomm. */
#ifdef PARALLEL 
MPI_Comm            MPI_INTRA_COMM;      /*!< the intra-communicator itself */
MPI_Group           MPI_INTRA_GROUP;     /*!< not needed, but for some reason you cannot have an intra-communicator without group */
#endif
} INTRA;

/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
typedef struct _PARTITION
{
enum   _FIELDTYP     fieldtyp;            /*!< type of field */
int                  ndis;                /*!< number of discretization in this field */
struct _PARTDISCRET *pdis;                /*!< vector of partitions of discretizations of this field */
} PARTITION;

/*!----------------------------------------------------------------------
\brief one proc's info about his partition of one discret.

<pre>                                                         m.gee 2/02
-the partition of one proc of one discretization
-the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
typedef struct _PARTDISCRET
{
int                 numnp;               /*!< total number of nodes on this parition including boundary nodes */
int                 numele;              /*!< total number of elements on this partition */
struct _NODE      **node;                /*!< ptrs to nodes on this part. */
struct _ELEMENT   **element;             /*!< ptrs to elements on this partition */

int                 inner_numnp;         /*!< number of pure inner-nodes */
struct _NODE      **inner_node;          /*!< ptrs to pure inner-nodes */
int                 bou_numnp;           /*!< number of boundary nodes */
struct _NODE      **bou_node;            /*!< ptrs to boundary nodes */

int                 inner_numele;        /*!< number of pure inner elements */
struct _ELEMENT   **inner_element;       /*!< pts to pure inner elements */
int                 bou_numele;          /*!< number of boundary elements */
struct _ELEMENT   **bou_element;         /*!< ptrs to boundary elements */

struct _ARRAY       coupledofs;          /*!< number of coupled dofs, which dofs, */

} PARTDISCRET;

/*! @} (documentation module close)*/


