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
int             Id;

int             mynode;       /* global Id of my FE-node */
struct _NODE    *node;        /* ptr to my FE-node */
struct _FIELD   *field;       /* ptr to FIELD this design belongs to */
} DNODE;

/*----------------------------------------------------------------------*
 | design line                                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _DLINE
{
int              Id;

ARRAY            mynode;      /* vector of global Ids of FE-nodes */
struct _NODE   **node;        /* vector of ptrs to FE-nodes */
struct _FIELD   *field;       /* ptr to FIELD this design belongs to */
} DLINE;

/*----------------------------------------------------------------------*
 | design surf                                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _DSURF
{
int              Id;

ARRAY            mynode;      /* vector of global Ids of FE-nodes */
struct _NODE   **node;        /* vector of ptrs to FE-nodes */
struct _FIELD   *field;       /* ptr to FIELD this design belongs to */
} DSURF;

/*----------------------------------------------------------------------*
 | design vol                                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _DVOL
{
int              Id;

ARRAY            mynode;      /* vector of global Ids of FE-nodes */
struct _NODE   **node;        /* vector of ptrs to FE-nodes */
struct _FIELD   *field;       /* ptr to FIELD this design belongs to */
} DVOL;
