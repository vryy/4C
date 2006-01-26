/*!---------------------------------------------------------------------
\file
\brief domain decomposition and metis structures

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

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
INT               myrank;                /*!< the individual processor number */
INT               nprocs;                /*!< total number of processors */
#ifdef PARALLEL
INT               numfld;                /*!< number of intra-communicators
                                           == number of fields */
struct _INTRA    *intra;                 /*!< vector of intra-communicator-structures
                                           correspondent to vector of FIELDs */
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
INT                 intra_rank;          /*!< proc's intra-rank */
INT                 intra_nprocs;        /*!< number of procs in this intracomm. */
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
INT                  ndis;                /*!< number of discretization in this field */
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
INT                 numnp;               /*!< total number of nodes on this parition including boundary nodes */
INT                 numele;              /*!< total number of elements on this partition */
struct _NODE      **node;                /*!< ptrs to nodes on this part. */
struct _ELEMENT   **element;             /*!< ptrs to elements on this partition */

  INT numlele;                  /* number of local elements on this partition */

INT                 inner_numnp;         /*!< number of pure inner-nodes */
struct _NODE      **inner_node;          /*!< ptrs to pure inner-nodes */
INT                 bou_numnp;           /*!< number of boundary nodes */
struct _NODE      **bou_node;            /*!< ptrs to boundary nodes */

INT                 inner_numele;        /*!< number of pure inner elements */
struct _ELEMENT   **inner_element;       /*!< pts to pure inner elements */
INT                 bou_numele;          /*!< number of boundary elements */
struct _ELEMENT   **bou_element;         /*!< ptrs to boundary elements */

struct _ARRAY       coupledofs;          /*!< number of coupled dofs, which dofs, */


#ifdef D_FLUID3_F
  INT                  num_fele;         /* number of sets of fast elements */
  struct _FAST_ELES   *fast_eles;        /* pointers to the sets of fast elements */
#endif


} PARTDISCRET;



#ifdef D_FLUID3_F
/*----------------------------------------------------------------------*/
/*!
  \brief structure for a set of fast elements

  This structure represents a set of fast elements. These elements are
  collected as pointers in a vector.

  \author mn
  \date 10/04

*/
/*----------------------------------------------------------------------*/
typedef struct _FAST_ELES
{
  FAST_ELE_TYP       fast_ele_typ;   /* the typ of fast elements collected
                                        in this set */
  INT                aloopl;         /* the number of elements in this sets */
  struct _ELEMENT  **ele_vec;        /* a vector containing pointers to the
                                        elements in this set */
} FAST_ELES;
#endif

/*! @} (documentation module close)*/


